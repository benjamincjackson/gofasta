package updown

import (
	"bytes"
	"fmt"
	"os"
	"runtime/debug"
	"strconv"
	"strings"
	"testing"
	"time"
)

// fuzzMode is one set of topranking's selection flags, taken from the bits of mode
type fuzzMode struct {
	sizeTotal, sizeUp, sizeSide, distAll, distUp, distPush int
	noFill, table                                          bool
	thresholdPair                                          float32
}

func modeFromBits(mode uint8) fuzzMode {
	var m fuzzMode
	switch mode % 4 {
	case 0:
		m.sizeTotal = 4
	case 1:
		m.distAll = 3
	case 2:
		m.distPush, m.distUp = 1, 1
	case 3:
		m.sizeUp, m.sizeSide, m.noFill = 2, 1, true
	}
	m.table = mode&4 != 0
	m.thresholdPair = 0.1
	if mode&8 != 0 {
		m.thresholdPair = 1
	}
	return m
}

// topRanking runs TopRanking with m on query and target of the given types ("csv" or "fasta")
func topRanking(query, target []byte, qtype, ttype string, root []byte, m fuzzMode, threads int) (string, error) {
	var out bytes.Buffer
	err := TopRanking(bytes.NewReader(query), bytes.NewReader(target), bytes.NewReader(root), &out, m.table,
		qtype, ttype, nil,
		m.sizeTotal, m.sizeUp, 0, m.sizeSide, 0,
		m.distAll, m.distUp, 0, 0,
		m.thresholdPair, 10000, m.noFill, m.distPush, threads, false)
	return out.String(), err
}

// FuzzTopRankingCSV: arbitrary bytes as the --query and --target csv files (updown list's format, header
// included) must give an error or output, never a panic or a hang, and the output must not depend on the
// number of threads.
func FuzzTopRankingCSV(f *testing.F) {
	h := listHeader + "\n"
	f.Add([]byte(h+"q1,A3C|G10T,20-25|40,2,7\n"), []byte(h+"t1,A3C,,1,0\nt2,A3C|G10T|C12A,1-5,3,5\nt3,,30,0,1\n"), uint8(0))
	f.Add([]byte(h+"q1,G10T|A3C,40|20-25,2,7\n"), []byte(h+"t1,A3G|G10T,10,2,1\nt2,C12A|A3C,5-1,2,5\n"), uint8(5))
	f.Add([]byte(h+"q1,,,0,0\nq2,A1C,1-100,1,100\n"), []byte(h+"t1,A1C|A1G,,2,0\n"), uint8(10))
	f.Add([]byte("not,a,header\n"), []byte(h), uint8(3))
	f.Fuzz(func(t *testing.T, query, target []byte, mode uint8) {
		quietStderr(t)
		m := modeFromBits(mode)
		collect(append(append(append([]byte{}, query...), 0), target...), func() {
			runWithTimeout(t, func() {
				_, _, _ = InputType(bytes.NewReader(query))
				one, err1 := topRanking(query, target, "csv", "csv", nil, m, 1)
				four, err4 := topRanking(query, target, "csv", "csv", nil, m, 4)
				if (err1 == nil) != (err4 == nil) || one != four {
					t.Fatalf("1 thread: %v\n%s\n4 threads: %v\n%s", err1, one, err4, four)
				}
			})
		})
	})
}

// fuzzRoot is the root of FuzzListMatchesTopRanking: 47 bases, so the last seven sites of a sequence are
// not compared eight at a time
const fuzzRoot = "ATGAAACCCGGGTTTAAACCCGGGTTTAAACCCGGGTTTAAATAACG"

// fuzzAlignment builds eight sequences of the root's width from data: each the root changed at the
// (position, character) pairs of the bytes, taken in turn
func fuzzAlignment(data []byte) []string {
	const alphabet = "ACGTACGTACGTNRYK-?acgtn"
	rows := make([][]byte, 8)
	for i := range rows {
		rows[i] = []byte(fuzzRoot)
	}
	for i := 0; i+1 < len(data); i += 2 {
		row := rows[(i/2)%len(rows)]
		row[int(data[i])%len(row)] = alphabet[int(data[i+1])%len(alphabet)]
	}
	seqs := make([]string, len(rows))
	for i, r := range rows {
		seqs[i] = string(r)
	}
	return seqs
}

func fastaOf(names []string, seqs []string) []byte {
	var b bytes.Buffer
	for i := range names {
		b.WriteString(">" + names[i] + "\n" + seqs[i] + "\n")
	}
	return b.Bytes()
}

// FuzzListMatchesTopRanking: for alignments built from the fuzz input (with, from the bits of mode, an
// ambiguous root), topranking on the fasta files must give exactly what it gives on updown list's csv of
// the same files, whatever the number of threads.
func FuzzListMatchesTopRanking(f *testing.F) {
	f.Add([]byte{}, uint8(0))
	f.Add([]byte{3, 1, 10, 2, 20, 3, 30, 12, 40, 16, 3, 1, 44, 13, 7, 0}, uint8(1))
	f.Add([]byte{0, 1, 1, 2, 2, 3, 3, 0, 46, 12, 45, 12, 44, 12, 43, 12, 5, 14, 5, 13}, uint8(18))
	f.Add([]byte{9, 2, 9, 2, 9, 2, 9, 3, 9, 3, 9, 3, 9, 1, 9, 1, 12, 16, 12, 17}, uint8(6))
	f.Fuzz(func(t *testing.T, data []byte, mode uint8) {
		quietStderr(t)
		root := []byte(fuzzRoot)
		if mode&16 != 0 {
			root[5], root[20], root[33] = 'N', 'R', '-'
		}
		rootFasta := []byte(">root\n" + string(root) + "\n")
		seqs := fuzzAlignment(data)
		names := []string{"s0", "s1", "s2", "s3", "s4", "s5", "s6", "s7"}
		query := fastaOf(names[:3], seqs[:3])
		target := fastaOf(append(names[3:], names[0]), append(seqs[3:], seqs[0]))
		m := modeFromBits(mode)
		collect(data, func() {
			runWithTimeout(t, func() {
				var qcsv, tcsv bytes.Buffer
				if err := List(bytes.NewReader(rootFasta), bytes.NewReader(query), &qcsv, 1, true); err != nil {
					t.Fatalf("list --query: %v", err)
				}
				if err := List(bytes.NewReader(rootFasta), bytes.NewReader(target), &tcsv, 1, true); err != nil {
					t.Fatalf("list --target: %v", err)
				}
				fromFasta, err := topRanking(query, target, "fasta", "fasta", rootFasta, m, 1)
				if err != nil {
					t.Fatalf("topranking on fasta: %v", err)
				}
				fromCSV, err := topRanking(qcsv.Bytes(), tcsv.Bytes(), "csv", "csv", nil, m, 4)
				if err != nil {
					t.Fatalf("topranking on csv: %v", err)
				}
				mixed, err := topRanking(query, tcsv.Bytes(), "fasta", "csv", rootFasta, m, 4)
				if err != nil {
					t.Fatalf("topranking on fasta and csv: %v", err)
				}
				if fromFasta != fromCSV || fromFasta != mixed {
					t.Fatalf("fasta:\n%s\ncsv:\n%s\nfasta query, csv target:\n%s\nquery:\n%s\ntarget:\n%s\nroot %s",
						fromFasta, fromCSV, mixed, query, target, root)
				}
			})
		})
	})
}

// collect runs fn on data. Normally a panic fails the fuzz test. With GOFASTA_FUZZ_COLLECT=DIR it is
// recovered and the first input seen for each distinct panic site is written to DIR instead, so one run
// can survey every distinct crash rather than stopping at the first.
func collect(data []byte, fn func()) {
	dir := os.Getenv("GOFASTA_FUZZ_COLLECT")
	if dir == "" {
		fn()
		return
	}
	defer func() {
		if r := recover(); r != nil {
			site := "?"
			for _, l := range strings.Split(string(debug.Stack()), "\n") {
				if i := strings.Index(l, "/gofasta/pkg/"); i >= 0 && strings.HasPrefix(l, "\t") && !strings.Contains(l, "_test.go") {
					site = strings.Fields(l[i+len("/gofasta/"):])[0]
					break
				}
			}
			name := strings.NewReplacer("/", "_", ":", "_", ".", "_").Replace(fmt.Sprint(r, " ", site))
			if len(name) > 120 {
				name = name[:120]
			}
			_ = os.MkdirAll(dir, 0o755)
			_ = os.WriteFile(dir+"/"+name, data, 0o644)
		}
	}()
	fn()
}

// runWithTimeout fails the test if fn does not return soon (a hang). The limit is 10 s, or
// GOFASTA_FUZZ_HANG_TIMEOUT seconds if that is set.
func runWithTimeout(t *testing.T, fn func()) {
	limit := 10 * time.Second
	if s, err := strconv.Atoi(os.Getenv("GOFASTA_FUZZ_HANG_TIMEOUT")); err == nil && s > 0 {
		limit = time.Duration(s) * time.Second
	}
	done := make(chan struct{})
	go func() { fn(); close(done) }()
	select {
	case <-done:
	case <-time.After(limit):
		t.Fatalf("hang: did not return within %v", limit)
	}
}

// quietStderr discards what is written to stderr (warnings about skipped records) until the test ends
func quietStderr(t *testing.T) {
	devnull, err := os.OpenFile(os.DevNull, os.O_WRONLY, 0)
	if err != nil {
		t.Fatal(err)
	}
	old := os.Stderr
	os.Stderr = devnull
	t.Cleanup(func() { os.Stderr = old; devnull.Close() })
}
