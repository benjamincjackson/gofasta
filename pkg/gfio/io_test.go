package gfio

import (
	"errors"
	"os"
	"testing"
)

func TestOpenIn(t *testing.T) {
	_, err := OpenIn("not/a/file.whatever", "-r / --reference")
	if err == nil || err.Error() != "open -r / --reference not/a/file.whatever: no such file or directory" {
		t.Error(err)
	}
	if f, err := OpenIn("stdin", "x"); err != nil || f != os.Stdin {
		t.Errorf("stdin: got %v, %v", f, err)
	}
}

func TestOpenOutStdout(t *testing.T) {
	if f, err := OpenOut("stdout"); err != nil || f != os.Stdout {
		t.Errorf("stdout: got %v, %v", f, err)
	}
}

func TestCloseOutReportsCloseError(t *testing.T) {
	f, err := os.CreateTemp(t.TempDir(), "x")
	if err != nil {
		t.Fatal(err)
	}
	f.Close() // already closed, so the next Close fails

	var got error
	CloseOut(f, &got)
	if got == nil {
		t.Errorf("expected the close error")
	}

	first := errors.New("first")
	got = first
	CloseOut(f, &got)
	if got != first {
		t.Errorf("an existing error must not be overwritten")
	}
}
