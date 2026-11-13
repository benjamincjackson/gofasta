package gfio

import (
	"bytes"
	"compress/gzip"
	"errors"
	"io"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

func TestOpenIn(t *testing.T) {
	_, err := OpenIn("not/a/file.whatever", "-r / --reference")
	if err == nil || err.Error() != "open -r / --reference not/a/file.whatever: no such file or directory" {
		t.Error(err)
	}
	if f, err := OpenIn("stdin", "x"); err != nil || f.(readCloser).Closer != os.Stdin {
		t.Errorf("stdin: got %v, %v", f, err)
	}
}

// gzipped writes text as gzip, in members of at most n bytes (several members, as bgzip writes)
func gzipped(t *testing.T, text string, n int) []byte {
	t.Helper()
	var b bytes.Buffer
	for i := 0; i < len(text); i += n {
		zw := gzip.NewWriter(&b)
		zw.Write([]byte(text[i:min(i+n, len(text))]))
		if err := zw.Close(); err != nil {
			t.Fatal(err)
		}
	}
	return b.Bytes()
}

// A gzipped file is decompressed whatever it is called, can be read again from its start, and an
// uncompressed one is returned as its *os.File
func TestOpenInGzip(t *testing.T) {
	dir := t.TempDir()
	text := ">a\nACGT\n>b\nTTTT\n"
	for name, content := range map[string][]byte{
		"plain.fasta":   []byte(text),
		"one.fasta.gz":  gzipped(t, text, len(text)),
		"many.fasta.gz": gzipped(t, text, 5),
		"no_suffix.fa":  gzipped(t, text, len(text)),
		"short.fasta":   []byte(">"),
		"empty.fasta":   nil,
	} {
		path := filepath.Join(dir, name)
		if err := os.WriteFile(path, content, 0o644); err != nil {
			t.Fatal(err)
		}
		f, err := OpenIn(path, "--x")
		if err != nil {
			t.Fatalf("%s: %v", name, err)
		}
		want := text
		if !bytes.HasPrefix(content, []byte{0x1f, 0x8b}) {
			want = string(content)
			if _, ok := f.(*os.File); !ok {
				t.Errorf("%s: uncompressed, but not returned as an *os.File", name)
			}
		}
		for pass := 1; pass <= 2; pass++ {
			got, err := io.ReadAll(f)
			if err != nil || string(got) != want {
				t.Errorf("%s, read %d: got %q, %v, want %q", name, pass, got, err, want)
			}
			if _, err := f.(io.Seeker).Seek(0, io.SeekStart); err != nil {
				t.Errorf("%s: seeking to the start: %v", name, err)
			}
		}
		f.Close()
	}
}

// A corrupt gzip header is an error naming the flag
func TestOpenInBadGzip(t *testing.T) {
	path := filepath.Join(t.TempDir(), "bad.gz")
	if err := os.WriteFile(path, []byte{0x1f, 0x8b, 0, 0}, 0o644); err != nil {
		t.Fatal(err)
	}
	if _, err := OpenIn(path, "-q / --query"); err == nil || !strings.Contains(err.Error(), "-q / --query") {
		t.Errorf("got %v, want an error naming -q / --query", err)
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
