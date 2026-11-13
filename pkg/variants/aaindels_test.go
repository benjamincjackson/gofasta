package variants

import (
	"strings"
	"testing"

	"github.com/virus-evolution/gofasta/pkg/encoding"
)

// aaIndelRef has a CDS at 1..27 (ATG GAG TTC AGA CAT CAT CAT GGT TAA: M E F R H H H G *) and three
// non-coding bases after it.
const aaIndelRef = "ATGGAGTTCAGACATCATCATGGTTAAAAA"

func plusCDS() Region {
	pos := make([]int, 27)
	for i := range pos {
		pos[i] = i + 1
	}
	return Region{Name: "g", Start: 1, Stop: 27, Strand: 1, Positions: pos, Translation: "MEFRHHHG*"}
}

// minusCDS is the same CDS on the reverse complement of aaIndelRef: plus-strand position p is 31-p
func minusCDS() Region {
	pos := make([]int, 27)
	for i := range pos {
		pos[i] = 30 - i
	}
	return Region{Name: "g", Start: 4, Stop: 30, Strand: -1, Positions: pos, Translation: "MEFRHHHG*"}
}

func revComp(s string) string {
	comp := map[byte]byte{'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', '-': '-'}
	b := make([]byte, len(s))
	for i := range s {
		b[len(s)-1-i] = comp[s[i]]
	}
	return string(b)
}

// callPair runs GetVariantsPair on one pairwise alignment and returns the formatted variants
func callPair(t *testing.T, refRow, queRow string, r Region, appendSNP bool) string {
	t.Helper()
	EA := encoding.MakeEncodingArray()
	ref, que := []byte(refRow), []byte(queRow)
	refLen := 0
	for i := range ref {
		if ref[i] != '-' {
			refLen++
		}
		ref[i] = EA[ref[i]]
		que[i] = EA[que[i]]
	}
	off, _ := GetMSAOffsets(ref)
	AS, err := GetVariantsPair(ref, que, "ref", "q", 0, []Region{r}, codes([]Region{r}, refLen), off)
	if err != nil {
		t.Fatal(err)
	}
	out := make([]string, 0)
	for _, v := range AS.Vs {
		s, err := FormatVariant(v, appendSNP)
		if err != nil {
			t.Fatal(err)
		}
		out = append(out, s)
	}
	return strings.Join(out, "|")
}

// del replaces the 1-based positions a..b of s with gaps; ins puts x after 1-based position p and returns
// the reference and query rows
func del(s string, a, b int) string { return s[:a-1] + strings.Repeat("-", b-a+1) + s[b:] }
func ins(s string, p int, x string) (string, string) {
	return s[:p] + strings.Repeat("-", len(x)) + s[p:], s[:p] + x + s[p:]
}
func snp(s string, p int, b string) string { return s[:p-1] + b + s[p:] }

func TestInFrameIndelAAsPlusStrand(t *testing.T) {
	r := plusCDS()
	ref := aaIndelRef
	insRef, insQue := ins(ref, 12, "GCT")
	repRef, repQue := ins(ref, 12, "CAT")
	stopRef, stopQue := ins(ref, 11, "CTG")
	stopCodonRef, stopCodonQue := ins(ref, 12, "TAA")
	for _, tc := range []struct {
		name, refRow, queRow, want string
	}{
		{"codon deletion", ref, del(ref, 7, 9), "aa:g:F3-|del:7:3"},
		{"deletion across codons changes a residue", ref, del(ref, 5, 10), "aa:g:E2G|del:5:6|aa:g:F3-|aa:g:R4-"},
		{"deletion in a repeat is placed 3'", ref, del(ref, 13, 15), "del:13:3|aa:g:H7-"},
		{"deletion in a repeat, other copy", ref, del(ref, 16, 18), "del:16:3|aa:g:H7-"},
		{"codon insertion", insRef, insQue, "aa:g:4A|ins:12:3"},
		{"insertion in a repeat is placed 3'", repRef, repQue, "ins:12:3|aa:g:7H"},
		{"insertion inside a codon with a stop", stopRef, stopQue, "aa:g:R4S|aa:g:4*|ins:11:3"},
		{"inserted stop codon: nothing after it", stopCodonRef, stopCodonQue, "ins:12:3|aa:g:H5*"},
		{"frameshift: no amino acid call", ref, del(ref, 7, 7), "del:7:1"},
		{"after a frameshift: no amino acid call", ref, del(del(ref, 7, 7), 13, 15), "del:7:1|del:13:3"},
		{"touches the stop codon", ref, del(ref, 22, 27), "del:22:6"},
		{"touches the start codon", ref, del(ref, 2, 4), "del:2:3"},
		{"SNP in the window is part of the call", ref, snp(del(ref, 5, 10), 4, "T"), "aa:g:E2*|del:5:6|aa:g:F3-|aa:g:R4-"},
		{"synonymous SNP in the window stays a nuc", insRef, snp(insQue, 18, "C"), "aa:g:4A|ins:12:3|nuc:T15C"},
		{"uncovered start is not a deletion", ref, del(ref, 1, 9), ""},
		{"deletion next to an uncovered start", ref, del(del(ref, 1, 6), 13, 15), "del:13:3|aa:g:H7-"},
		{"N in the window: no amino acid call", ref, snp(del(ref, 5, 10), 11, "N"), "del:5:6"},
		{"deletion is not shifted across an N", ref, snp(del(ref, 13, 15), 20, "N"), "del:13:3|aa:g:H6-"},
	} {
		if got := callPair(t, tc.refRow, tc.queRow, r, false); got != tc.want {
			t.Errorf("%s: got %q, want %q", tc.name, got, tc.want)
		}
	}
}

// --append-snps gives the nucleotide indels (and SNPs) behind an amino acid change caused by an indel
func TestInFrameIndelAAsAppendSNPs(t *testing.T) {
	r := plusCDS()
	ref := aaIndelRef
	if got, want := callPair(t, ref, del(ref, 5, 10), r, true), "aa:g:E2G(del:5:6)|del:5:6|aa:g:F3-(del:5:6)|aa:g:R4-(del:5:6)"; got != want {
		t.Errorf("got %q, want %q", got, want)
	}
	if got, want := callPair(t, ref, snp(del(ref, 5, 10), 4, "T"), r, true), "aa:g:E2*(del:5:6;nuc:G4T)|del:5:6|aa:g:F3-(del:5:6;nuc:G4T)|aa:g:R4-(del:5:6;nuc:G4T)"; got != want {
		t.Errorf("got %q, want %q", got, want)
	}
}

// The same changes on the minus strand give the same amino acid calls
func TestInFrameIndelAAsMinusStrand(t *testing.T) {
	r := minusCDS()
	ref := revComp(aaIndelRef)
	// plus-strand insertion of GCT after 12 is AGC between minus-strand 18 and 19
	insRef, insQue := ins(ref, 18, "AGC")
	for _, tc := range []struct {
		name, refRow, queRow, want string
	}{
		{"codon deletion", ref, del(ref, 22, 24), "del:22:3|aa:g:F3-"},
		{"deletion across codons changes a residue", ref, del(ref, 21, 26), "aa:g:R4-|del:21:6|aa:g:F3-|aa:g:E2G"},
		{"deletion in a repeat is placed 3'", ref, del(ref, 16, 18), "aa:g:H7-|del:16:3"},
		{"codon insertion", insRef, insQue, "ins:18:3|aa:g:4A"},
	} {
		if got := callPair(t, tc.refRow, tc.queRow, r, false); got != tc.want {
			t.Errorf("%s: got %q, want %q", tc.name, got, tc.want)
		}
	}
}
