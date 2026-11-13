package variants

import (
	"math/rand"
	"reflect"
	"strings"
	"testing"

	"github.com/virus-evolution/gofasta/pkg/encoding"
)

// randomAlignment returns an aligned reference and query of width n (as IUPAC text): runs of matching
// columns, substitutions, ambiguity codes and Ns, deletions (gaps in the query), insertions (gaps in the
// reference) and columns gapped in both, and gaps at either end of the query
func randomAlignment(rng *rand.Rand, n int) (string, string) {
	bases, ambig := "ACGT", "RYKMSWN"
	ref, que := make([]byte, 0, n), make([]byte, 0, n)
	for len(ref) < n {
		run := 1 + rng.Intn(4)
		kind := rng.Intn(20)
		for k := 0; k < run && len(ref) < n; k++ {
			b := bases[rng.Intn(4)]
			switch kind {
			case 0: // substitution
				ref, que = append(ref, b), append(que, bases[rng.Intn(4)])
			case 1: // ambiguity code or N
				ref, que = append(ref, b), append(que, ambig[rng.Intn(len(ambig))])
			case 2: // deletion
				ref, que = append(ref, b), append(que, '-')
			case 3: // insertion
				ref, que = append(ref, '-'), append(que, b)
			case 4: // gapped in both
				ref, que = append(ref, '-'), append(que, '-')
			default:
				ref, que = append(ref, b), append(que, b)
			}
		}
	}
	for i := 0; i < rng.Intn(10) && i < n; i++ {
		que[i] = '-'
	}
	for i := n - min(n, rng.Intn(10)); i < n; i++ {
		que[i] = '-'
	}
	return string(ref), string(que)
}

func encodeString(s string) []byte {
	EA := encoding.MakeEncodingArray()
	b := []byte(s)
	for i := range b {
		b[i] = EA[b[i]]
	}
	return b
}

// The one-pass calls are those of GetVariantsPair with no CDSs, whatever the alignment
func TestNoAnnotationMatchesGetVariantsPair(t *testing.T) {
	rng := rand.New(rand.NewSource(1))
	for rep := 0; rep < 3000; rep++ {
		n := 1 + rng.Intn(200)
		refRow, queRow := randomAlignment(rng, n)
		ref, que := encodeString(refRow), encodeString(queRow)
		degapped := strings.Count(refRow, "-")
		inter := make([]int, n-degapped)
		for i := range inter {
			inter[i] = i + 1
		}
		off, _ := GetMSAOffsets(ref)
		AS, err := GetVariantsPair(ref, que, "ref", "q", 0, nil, inter, off)
		if err != nil {
			t.Fatal(err)
		}
		if got := getVariantsNoAnnotation(ref, que); !reflect.DeepEqual(got, AS.Vs) {
			t.Fatalf("ref %s\nque %s\ngot  %v\nwant %v", refRow, queRow, got, AS.Vs)
		}
	}
}

// Some cases by hand: an insertion and a substitution at the same position, an insertion at the very
// start, a deletion across a gap column, and gaps at the ends
func TestNoAnnotationByHand(t *testing.T) {
	for _, tc := range []struct{ ref, que, want string }{
		{"ACGTACGTACGT", "ACGTACGTACGT", ""},
		{"ACGTACGTACGT", "ACGTTCGTACGT", "nuc:A5T"},
		{"ACG--TACGTACGT", "ACCTTTACGTACGT", "ins:3:2|nuc:G3C"}, // sorted as GetVariantsPair: ins before nuc
		{"--ACGTACGT", "GGACGTACGT", "ins:0:2"},
		{"ACGTA-CGTACGT", "ACG---CGTACGT", "del:4:2"},
		{"ACGTACGTACGT", "---TACGTA---", ""},
		{"ACGTACGTACGTACGTACGT", "ACGTACGTACGTNCGTACGT", ""},
	} {
		got := make([]string, 0)
		for _, v := range getVariantsNoAnnotation(encodeString(tc.ref), encodeString(tc.que)) {
			s, _ := FormatVariant(v, false)
			got = append(got, s)
		}
		if g := strings.Join(got, "|"); g != tc.want {
			t.Errorf("ref %s que %s: got %q, want %q", tc.ref, tc.que, g, tc.want)
		}
	}
}
