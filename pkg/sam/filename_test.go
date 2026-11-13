package sam

import (
	"strings"
	"testing"
	"unicode/utf8"
)

// File names replace '/' and are cut, at a character boundary, so that they fit in 255 bytes
func TestPairFileName(t *testing.T) {
	for _, tc := range []struct {
		query, want string
		cut         bool
	}{
		{"a/b", "a_b.fasta", false},
		{"hCoV-19/England/X/2020", "hCoV-19_England_X_2020.fasta", false},
		{strings.Repeat("x", 249), strings.Repeat("x", 249) + ".fasta", false},
		{strings.Repeat("x", 250), strings.Repeat("x", 249) + ".fasta", true},
		{strings.Repeat("é", 150), strings.Repeat("é", 124) + ".fasta", true},            // 2 bytes each: byte 249 is mid-character
		{"x" + strings.Repeat("€", 100), "x" + strings.Repeat("€", 82) + ".fasta", true}, // 3 bytes each
	} {
		got, cut := pairFileName(tc.query)
		if got != tc.want || cut != tc.cut {
			t.Errorf("%.20q...: got %.30q... (cut %v), want %.30q... (cut %v)", tc.query, got, cut, tc.want, tc.cut)
		}
		if !utf8.ValidString(got) || len(got) > 255 {
			t.Errorf("%.20q...: %d bytes, valid UTF-8 %v", tc.query, len(got), utf8.ValidString(got))
		}
	}
}
