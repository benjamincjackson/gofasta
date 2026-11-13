package updown

import (
	"io"
	"strings"
	"testing"
)

// The type comes from the content, and the returned reader gives all of the input
func TestInputType(t *testing.T) {
	for _, tc := range []struct{ in, want string }{
		{">a\nACGT\n", "fasta"},
		{"\n\n  >a\nACGT\n", "fasta"},
		{listHeader + "\nq1,A1T,,1,0\n", "csv"},
		{listHeader, "csv"},
		{"query,something,else\n", ""},
		{"ACGT\n", ""},
		{"", ""},
	} {
		r, got, err := InputType(strings.NewReader(tc.in))
		if got != tc.want || (tc.want == "") != (err != nil) {
			t.Errorf("%q: got %q, %v, want %q", tc.in, got, err, tc.want)
		}
		if all, _ := io.ReadAll(r); string(all) != tc.in {
			t.Errorf("%q: reader gave %q", tc.in, all)
		}
	}
}
