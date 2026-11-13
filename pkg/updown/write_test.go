package updown

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

// A write that fails is reported, the header's included (with no rows after it to fail as well)
func TestWriteUpdownTableReportsWriteErrors(t *testing.T) {
	rows := []updownCatchmentStruct{{qname: "q", up: updownCatchmentSubStruct{catchment: []resultsStruct{{tname: "t", distance: 1}}}}}
	for _, tc := range []struct {
		name    string
		results []updownCatchmentStruct
		ok      int
	}{
		{"header, no rows", nil, 0},
		{"header", rows, 0},
		{"row", rows, 1},
	} {
		if err := writeUpdownTable(&failingWriter{ok: tc.ok}, tc.results); !errors.Is(err, errFull) {
			t.Errorf("%s: got %v, want the write error", tc.name, err)
		}
	}
}
