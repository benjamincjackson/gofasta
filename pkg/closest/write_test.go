package closest

import (
	"errors"
	"testing"
)

var errFull = errors.New("disk full")

// failingWriter accepts the first ok writes and fails every one after
type failingWriter struct{ ok int }

func (w *failingWriter) Write(p []byte) (int, error) {
	if w.ok == 0 {
		return 0, errFull
	}
	w.ok--
	return len(p), nil
}

// A write that fails after the header is reported, by every way of writing the results
func TestWriteErrorsAreReported(t *testing.T) {
	results := []resultsStruct{{qname: "q", tname: "t", distance: 1}}
	catchments := []catchmentStruct{{qname: "q", catchment: results}}
	for name, write := range map[string]func(w *failingWriter) error{
		"closest raw":        func(w *failingWriter) error { return writeClosest(results, "raw", w) },
		"closest snp":        func(w *failingWriter) error { return writeClosest(results, "snp", w) },
		"closest tn93":       func(w *failingWriter) error { return writeClosest(results, "tn93", w) },
		"closest -n":         func(w *failingWriter) error { return writeClosestN(catchments, w) },
		"closest -n --table": func(w *failingWriter) error { return writeClosestNTable(catchments, w, "raw") },
		"closest -n --table snp": func(w *failingWriter) error {
			return writeClosestNTable(catchments, w, "snp")
		},
	} {
		if err := write(&failingWriter{ok: 1}); !errors.Is(err, errFull) {
			t.Errorf("%s: got %v, want the write error", name, err)
		}
	}
}
