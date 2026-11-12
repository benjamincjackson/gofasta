package pipeline

import (
	"bytes"
	"context"
	"errors"
	"fmt"
	"io"
	"math/rand"
	"reflect"
	"runtime"
	"strings"
	"sync/atomic"
	"testing"
	"time"
)

// counter returns a source of the integers 0..n-1.
func counter(n int) func() (int, error) {
	i := 0
	return func() (int, error) {
		if i == n {
			return 0, io.EOF
		}
		i++
		return i - 1, nil
	}
}

// jitter makes workers finish out of order.
func jitter(x int) (int, error) {
	time.Sleep(time.Duration(rand.Intn(200)) * time.Microsecond)
	return x * 10, nil
}

func TestRunDeliversInInputOrder(t *testing.T) {
	for _, threads := range []int{0, 1, 2, 7} {
		var got []int
		err := Run(context.Background(), threads, counter(500), jitter, func(r Result[int]) error {
			if r.Err != nil {
				t.Fatalf("unexpected error %v", r.Err)
			}
			if r.Val != r.Idx*10 {
				t.Fatalf("record %d has value %d", r.Idx, r.Val)
			}
			got = append(got, r.Idx)
			return nil
		})
		if err != nil {
			t.Fatal(err)
		}
		if len(got) != 500 {
			t.Fatalf("threads %d: got %d results", threads, len(got))
		}
		for i, idx := range got {
			if idx != i {
				t.Fatalf("threads %d: position %d has record %d", threads, i, idx)
			}
		}
	}
}

func TestRunEmptyInput(t *testing.T) {
	called := false
	err := Run(context.Background(), 4, counter(0), jitter, func(Result[int]) error { called = true; return nil })
	if err != nil || called {
		t.Errorf("got err %v, sink called %v", err, called)
	}
}

// Every record reaches the sink, failed ones included, in order.
func TestRunWorkErrorsReachTheSink(t *testing.T) {
	bad := errors.New("bad record")
	work := func(x int) (int, error) {
		if x%3 == 0 {
			return 0, bad
		}
		return jitter(x)
	}
	var seen, failed int
	err := Run(context.Background(), 4, counter(100), work, func(r Result[int]) error {
		if r.Idx != seen {
			t.Fatalf("got record %d, want %d", r.Idx, seen)
		}
		seen++
		if (r.Err != nil) != (r.Idx%3 == 0) {
			t.Fatalf("record %d: err %v", r.Idx, r.Err)
		}
		if r.Err != nil {
			failed++
		}
		return nil
	})
	if err != nil || seen != 100 || failed != 34 {
		t.Errorf("err %v, seen %d, failed %d", err, seen, failed)
	}
}

// A read error arrives after every record before it, and ends the input.
func TestRunReadErrorArrivesInOrder(t *testing.T) {
	broken := errors.New("unreadable")
	i := 0
	next := func() (int, error) {
		if i == 20 {
			return 0, broken
		}
		i++
		return i - 1, nil
	}
	var seen []Result[int]
	err := Run(context.Background(), 4, next, jitter, func(r Result[int]) error {
		seen = append(seen, r)
		return nil
	})
	if err != nil {
		t.Fatal(err)
	}
	var readErr *ReadError
	if len(seen) != 21 || !errors.Is(seen[20].Err, broken) || !errors.As(seen[20].Err, &readErr) || seen[20].Idx != 20 {
		t.Fatalf("got %d results, last %+v", len(seen), seen[len(seen)-1])
	}
	for _, r := range seen[:20] {
		if r.Err != nil {
			t.Fatalf("record %d: %v", r.Idx, r.Err)
		}
	}
}

// With FirstError, the error of the first failed record is returned, and later records do not reach f.
func TestFirstErrorStopsAtTheFirstFailure(t *testing.T) {
	work := func(x int) (int, error) {
		if x == 40 || x == 41 || x == 90 {
			return 0, fmt.Errorf("record %d", x)
		}
		return jitter(x)
	}
	var written []int
	err := Run(context.Background(), 8, counter(1000), work, FirstError(func(v int) error {
		written = append(written, v)
		return nil
	}))
	if err == nil || err.Error() != "record 40" {
		t.Fatalf("got %v, want the error of record 40", err)
	}
	if len(written) != 40 {
		t.Fatalf("wrote %d records, want 40", len(written))
	}
}

// After the sink fails, Run returns without leaking goroutines, even though the input is not exhausted.
func TestRunSinkErrorStopsEverything(t *testing.T) {
	before := runtime.NumGoroutine()
	var read atomic.Int64
	next := func() (int, error) {
		read.Add(1)
		return int(read.Load()), nil // endless input
	}
	stop := errors.New("disk full")
	err := Run(context.Background(), 4, next, jitter, func(r Result[int]) error {
		if r.Idx == 10 {
			return stop
		}
		return nil
	})
	if err != stop {
		t.Fatalf("got %v", err)
	}
	deadline := time.Now().Add(2 * time.Second)
	for runtime.NumGoroutine() > before && time.Now().Before(deadline) {
		time.Sleep(time.Millisecond)
	}
	if n := runtime.NumGoroutine(); n > before {
		t.Errorf("%d goroutines still running (had %d)", n, before)
	}
	if n := read.Load(); n > 10+4*inFlightPerThread+1 {
		t.Errorf("read %d records after the sink failed at record 10", n)
	}
}

// A record that takes a long time does not let the others pile up without limit.
func TestRunBoundsRecordsInFlight(t *testing.T) {
	const threads = 2
	var read atomic.Int64
	next := func() (int, error) {
		n := read.Add(1)
		if n > 10000 {
			return 0, io.EOF
		}
		return int(n - 1), nil
	}
	release := make(chan struct{})
	work := func(x int) (int, error) {
		if x == 0 {
			<-release
		}
		return x, nil
	}
	done := make(chan error)
	go func() { done <- Run(context.Background(), threads, next, work, func(Result[int]) error { return nil }) }()
	time.Sleep(50 * time.Millisecond)
	if n := read.Load(); n > threads*inFlightPerThread+1 {
		t.Errorf("read %d records while the first was still being processed", n)
	}
	close(release)
	if err := <-done; err != nil {
		t.Fatal(err)
	}
}

// Cancelling the caller's context stops the pipeline: Run returns ctx.Err() without leaking goroutines.
func TestRunStopsWhenTheContextIsCancelled(t *testing.T) {
	before := runtime.NumGoroutine()
	ctx, cancel := context.WithCancel(context.Background())
	var read atomic.Int64
	next := func() (int, error) {
		return int(read.Add(1) - 1), nil // endless input
	}
	delivered := 0
	err := Run(ctx, 4, next, jitter, func(r Result[int]) error {
		delivered++
		if r.Idx == 50 {
			cancel()
		}
		return nil
	})
	if err != context.Canceled {
		t.Fatalf("got %v, want context.Canceled", err)
	}
	if delivered != 51 {
		t.Errorf("sink called %d times after cancelling at record 50", delivered)
	}
	waitForGoroutines(t, before)
}

// A context cancelled after the last record has reached the sink does not make a finished run fail.
func TestRunCompleteDespiteLateCancel(t *testing.T) {
	ctx, cancel := context.WithCancel(context.Background())
	err := Run(ctx, 4, counter(10), jitter, func(r Result[int]) error {
		if r.Idx == 9 {
			cancel()
		}
		return nil
	})
	if err != nil {
		t.Errorf("got %v", err)
	}
}

// A panic in work becomes that record's error, as a *RecordPanicError; the other records are processed
// and every record reaches the sink in order.
func TestRunPanicInWorkBecomesTheRecordsError(t *testing.T) {
	before := runtime.NumGoroutine()
	cause := errors.New("record 30 is cursed")
	work := func(x int) (int, error) {
		if x == 30 {
			panic(cause)
		}
		return jitter(x)
	}
	seen := 0
	err := Run(context.Background(), 4, counter(1000), work, func(r Result[int]) error {
		if r.Idx != seen {
			t.Fatalf("got record %d, want %d", r.Idx, seen)
		}
		seen++
		if r.Idx == 30 {
			var pe *RecordPanicError
			if !errors.As(r.Err, &pe) || !errors.Is(r.Err, cause) || len(pe.Stack) == 0 {
				t.Errorf("record 30: got %#v", r.Err)
			}
			if r.Err.Error() != "panic: record 30 is cursed" {
				t.Errorf("record 30: got message %q", r.Err.Error())
			}
			return nil
		}
		if r.Err != nil || r.Val != r.Idx*10 {
			t.Errorf("record %d: got %v, %v", r.Idx, r.Val, r.Err)
		}
		return nil
	})
	if err != nil || seen != 1000 {
		t.Errorf("got %v after %d records", err, seen)
	}
	waitForGoroutines(t, before)
}

// A Handler warns about and skips a record whose work panicked, like any other failed record, and a
// strict Handler stops with its error.
func TestHandlerSkipsPanickedRecords(t *testing.T) {
	work := func(x int) (int, error) {
		if x == 2 {
			panic("boom")
		}
		return x, nil
	}
	for _, strict := range []bool{false, true} {
		var warn bytes.Buffer
		var got []int
		h := NewHandler(strict, &warn, Noun{One: "record", Many: "records"}, func(x int) error {
			got = append(got, x)
			return nil
		})
		err := Run(context.Background(), 2, counter(5), work, h.Sink)
		if strict {
			if err == nil || err.Error() != "(record 3): panic: boom" {
				t.Errorf("strict: got error %v", err)
			}
			continue
		}
		if err != nil {
			t.Fatal(err)
		}
		if err := h.Finish(); err != nil {
			t.Fatal(err)
		}
		if !reflect.DeepEqual(got, []int{0, 1, 3, 4}) {
			t.Errorf("got records %v", got)
		}
		want := "warning: skipping (record 3): panic: boom\nwarning: skipped 1 of 5 records\n"
		if warn.String() != want {
			t.Errorf("got warnings %q, want %q", warn.String(), want)
		}
	}
}

// namedInt is a record that knows its name
type namedInt int

func (n namedInt) RecordName() string { return fmt.Sprintf("rec%d", int(n)) }

// Errors from work on a record with a RecordName, including panics, are named with it unless they already
// have a name, so warnings and strict errors say which record failed.
func TestRecordErrorsAreNamed(t *testing.T) {
	newNext := func() func() (namedInt, error) {
		i := 0
		return func() (namedInt, error) {
			if i == 5 {
				return 0, io.EOF
			}
			i++
			return namedInt(i - 1), nil
		}
	}
	work := func(n namedInt) (int, error) {
		switch n {
		case 1:
			return 0, errors.New("bad")
		case 2:
			panic("boom")
		case 3:
			return 0, Named("custom", errors.New("named already"))
		}
		return int(n), nil
	}

	var warn bytes.Buffer
	h := NewHandler(false, &warn, Noun{One: "record", Many: "records"}, func(int) error { return nil })
	if err := Run(context.Background(), 2, newNext(), work, h.Sink); err != nil {
		t.Fatal(err)
	}
	if err := h.Finish(); err != nil {
		t.Fatal(err)
	}
	want := "warning: skipping rec1 (record 2): bad\n" +
		"warning: skipping rec2 (record 3): panic: boom\n" +
		"warning: skipping custom (record 4): named already\n" +
		"warning: skipped 3 of 5 records\n"
	if warn.String() != want {
		t.Errorf("got warnings %q, want %q", warn.String(), want)
	}

	h = NewHandler(true, io.Discard, Noun{One: "record", Many: "records"}, func(int) error { return nil })
	err := Run(context.Background(), 2, newNext(), work, h.Sink)
	if err == nil || err.Error() != "rec1 (record 2): bad" {
		t.Errorf("strict: got %v", err)
	}
}

// A panic in next is re-raised in the caller too, and an error value stays reachable with errors.As.
func TestRunReraisesPanicsInNext(t *testing.T) {
	cause := errors.New("reader fell over")
	i := 0
	next := func() (int, error) {
		if i == 5 {
			panic(cause)
		}
		i++
		return i - 1, nil
	}
	p := catchPanic(func() {
		_ = Run(context.Background(), 2, next, jitter, func(Result[int]) error { return nil })
	})
	pe, ok := p.(*PanicError)
	if !ok || !errors.Is(pe, cause) {
		t.Fatalf("recovered %#v", p)
	}
}

func catchPanic(f func()) (p any) {
	defer func() { p = recover() }()
	f()
	return nil
}

func waitForGoroutines(t *testing.T, want int) {
	t.Helper()
	deadline := time.Now().Add(2 * time.Second)
	for runtime.NumGoroutine() > want && time.Now().Before(deadline) {
		time.Sleep(time.Millisecond)
	}
	if n := runtime.NumGoroutine(); n > want {
		t.Errorf("%d goroutines still running (had %d)", n, want)
	}
}

// handlerRun runs 10 records through a Handler; records 3 and 7 fail in work (7 with a name), and if
// readFailAt >= 0 reading fails at that record.
func handlerRun(t *testing.T, strict bool, readFailAt int) (got []int, warnings string, runErr, finishErr error) {
	t.Helper()
	i := 0
	next := func() (int, error) {
		if i == 10 {
			return 0, io.EOF
		}
		if i == readFailAt {
			return 0, errors.New("unreadable")
		}
		i++
		return i - 1, nil
	}
	work := func(x int) (int, error) {
		switch x {
		case 3:
			return 0, errors.New("bad three")
		case 7:
			return 0, Named("seven", errors.New("bad seven"))
		}
		return x, nil
	}
	var warn strings.Builder
	h := NewHandler(strict, &warn, Noun{"record", "records"}, func(v int) error {
		got = append(got, v)
		return nil
	})
	runErr = Run(context.Background(), 3, next, work, h.Sink)
	if runErr == nil {
		finishErr = h.Finish()
	}
	return got, warn.String(), runErr, finishErr
}

func TestHandlerSkipsFailedRecords(t *testing.T) {
	got, warnings, runErr, finishErr := handlerRun(t, false, -1)
	if runErr != nil || finishErr != nil {
		t.Fatal(runErr, finishErr)
	}
	if fmt.Sprint(got) != "[0 1 2 4 5 6 8 9]" {
		t.Errorf("got %v", got)
	}
	want := "warning: skipping (record 4): bad three\n" +
		"warning: skipping seven (record 8): bad seven\n" +
		"warning: skipped 2 of 10 records\n"
	if warnings != want {
		t.Errorf("warnings:\n%s\nwant:\n%s", warnings, want)
	}
}

func TestHandlerStrictStopsAtTheFirstFailure(t *testing.T) {
	got, warnings, runErr, _ := handlerRun(t, true, -1)
	if runErr == nil || runErr.Error() != "(record 4): bad three" || fmt.Sprint(got) != "[0 1 2]" || warnings != "" {
		t.Errorf("got %v, err %v, warnings %q", got, runErr, warnings)
	}
}

// A read error stops the pipeline even when failed records are being skipped.
func TestHandlerReadErrorsAreFatal(t *testing.T) {
	got, _, runErr, _ := handlerRun(t, false, 5)
	if runErr == nil || runErr.Error() != "unreadable" || fmt.Sprint(got) != "[0 1 2 4]" {
		t.Errorf("got %v, err %v", got, runErr)
	}
}

func TestHandlerEveryRecordFailed(t *testing.T) {
	var warn strings.Builder
	h := NewHandler(false, &warn, Noun{"query", "queries"}, func(int) error { return nil })
	work := func(int) (int, error) { return 0, errors.New("bad") }
	if err := Run(context.Background(), 2, counter(3), work, h.Sink); err != nil {
		t.Fatal(err)
	}
	err := h.Finish()
	if err == nil || err.Error() != "none of the 3 queries could be processed" {
		t.Errorf("got %v", err)
	}
	if !strings.HasSuffix(warn.String(), "warning: skipped 3 of 3 queries\n") || h.Total() != 3 || h.Skipped() != 3 {
		t.Errorf("warnings %q, total %d, skipped %d", warn.String(), h.Total(), h.Skipped())
	}
	// no records at all is not a failure
	h = NewHandler(false, &warn, Noun{"query", "queries"}, func(int) error { return nil })
	if err := Run(context.Background(), 2, counter(0), work, h.Sink); err != nil || h.Finish() != nil {
		t.Errorf("empty input: %v", err)
	}
}
