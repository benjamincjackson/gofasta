package variants

import (
	"bytes"
	"strings"
	"testing"

	"github.com/virus-evolution/gofasta/pkg/encoding"
)

// annotatePair runs GetVariantsPair on one pairwise alignment against the CDS r
func annotatePair(t *testing.T, refRow, queRow string, r Region) AnnoStructs {
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
	return AS
}

// writeGrouped writes one query's variants with a LineWriter (groupAA, within start..end) and returns the
// mutations column
func writeGrouped(t *testing.T, AS AnnoStructs, appendSNP bool, start, end int) string {
	t.Helper()
	var out bytes.Buffer
	lw, err := NewLineWriter(&out, start, end, appendSNP, true, "ref")
	if err != nil {
		t.Fatal(err)
	}
	if err := lw.Write(AS); err != nil {
		t.Fatal(err)
	}
	line := strings.TrimSuffix(strings.Split(out.String(), "\n")[1], "\n")
	return strings.TrimPrefix(line, "q,")
}

// With --group-aa the amino acid changes from one indel window are one entry; single changes and pure
// insertions are written as without it, and nucleotide entries are unchanged
func TestGroupAA(t *testing.T) {
	r := plusCDS()
	ref := aaIndelRef
	insRef, insQue := ins(ref, 12, "GCT")
	stopRef, stopQue := ins(ref, 11, "CTG")
	for _, tc := range []struct {
		name, refRow, queRow, want string
	}{
		{"deletion across codons", ref, del(ref, 5, 10), "aa:g:EFR2G--|del:5:6"},
		{"two residues deleted in a repeat, placed 3'", ref, del(ref, 13, 18), "del:13:6|aa:g:HH6--"},
		{"insertion next to a changed residue", stopRef, stopQue, "aa:g:R4S*|ins:11:3"},
		{"frame restored 3 codons on", ref, del(del(ref, 7, 7), 17, 18), "aa:g:FRHH3SDI-|del:7:1|del:17:2"},
		{"single deletion: as without the flag", ref, del(ref, 7, 9), "aa:g:F3-|del:7:3"},
		{"pure insertion: as without the flag", insRef, insQue, "aa:g:4A|ins:12:3"},
		{"SNP alone: as without the flag", ref, snp(ref, 4, "T"), "aa:g:E2*"},
	} {
		if got := writeGrouped(t, annotatePair(t, tc.refRow, tc.queRow, r), false, -1, -1); got != tc.want {
			t.Errorf("%s: got %q, want %q", tc.name, got, tc.want)
		}
	}
}

// The grouped entry carries the window's SNPs with --append-snps, and is written if any of its residues'
// codons is within --start/--end, once
func TestGroupAAAppendSNPsAndRange(t *testing.T) {
	r := plusCDS()
	AS := annotatePair(t, aaIndelRef, snp(del(aaIndelRef, 5, 10), 4, "T"), r)
	if got, want := writeGrouped(t, AS, true, -1, -1), "aa:g:EFR2*--(del:5:6;nuc:G4T)|del:5:6"; got != want {
		t.Errorf("--append-snps: got %q, want %q", got, want)
	}
	// only the codon of R4 (10-12) is in range: the group is still written, the deletion (5-10) too
	if got, want := writeGrouped(t, AS, false, 10, 12), "aa:g:EFR2*--|del:5:6"; got != want {
		t.Errorf("--start 10 --end 12: got %q, want %q", got, want)
	}
	if got, want := writeGrouped(t, AS, false, 13, -1), ""; got != want {
		t.Errorf("--start 13: got %q, want %q", got, want)
	}
}

// An unchanged residue between two changed ones is part of the grouped entry, with the query's residue
// repeated (A10X, B unchanged, C12 deleted: ABC10XB-)
func TestGroupWindowUnchangedResidue(t *testing.T) {
	pos := make([]int, 21)
	for i := range pos {
		pos[i] = i + 1
	}
	region := Region{Name: "g", Strand: 1, Positions: pos, Translation: "MABCDE*"}
	entries := []Variant{
		{Changetype: "aa", Feature: "g", RefAl: "A", QueAl: "X", Residue: 2},
		{Changetype: "aa", Feature: "g", RefAl: "C", QueAl: "-", Residue: 4},
	}
	groupWindow(entries, region, region.Translation, 1, "XB")
	if entries[0].Group == nil || entries[0].Group != entries[1].Group {
		t.Fatalf("the entries should share one group: %+v", entries)
	}
	got, err := FormatVariant(*entries[0].Group, false)
	if err != nil {
		t.Fatal(err)
	}
	if want := "aa:g:ABC2XB-"; got != want {
		t.Errorf("got %q, want %q", got, want)
	}
}

// --aggregate with --group-aa counts the grouped entries
func TestGroupAAAggregate(t *testing.T) {
	r := plusCDS()
	var out bytes.Buffer
	agg, err := NewAggregator(&out, -1, -1, false, true, "ref")
	if err != nil {
		t.Fatal(err)
	}
	for i, q := range []string{del(aaIndelRef, 5, 10), del(aaIndelRef, 5, 10), del(aaIndelRef, 7, 9)} {
		AS := annotatePair(t, aaIndelRef, q, r)
		AS.Queryname, AS.Idx = "q"+string(rune('1'+i)), i
		if err := agg.Add(AS); err != nil {
			t.Fatal(err)
		}
	}
	if err := agg.Write(0); err != nil {
		t.Fatal(err)
	}
	want := "mutation,frequency\n" +
		"aa:g:EFR2G--,0.666666667\n" +
		"del:5:6,0.666666667\n" +
		"aa:g:F3-,0.333333333\n" +
		"del:7:3,0.333333333\n"
	if out.String() != want {
		t.Errorf("got\n%s\nwant\n%s", out.String(), want)
	}
}
