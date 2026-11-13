package variants

import (
	"bytes"
	"errors"
	"reflect"
	"strings"
	"testing"

	"github.com/virus-evolution/gofasta/pkg/alphabet"
	"github.com/virus-evolution/gofasta/pkg/genbank"
	"github.com/virus-evolution/gofasta/pkg/gff"
)

// phaseTestRef is the reference used by the tests below (20 bases)
const phaseTestRef = "GGGGTTAAACCCGTTTTTTT"

// regionsFromGFFRows parses the given CDS rows (tab-separated, on seqid ref1) and returns the CDS regions
func regionsFromGFFRows(t *testing.T, rows ...string) ([]Region, []int, error) {
	t.Helper()
	text := "##gff-version 3\n##sequence-region ref1 1 20\n" + strings.Join(rows, "\n") + "\n"
	anno, err := gff.ReadGFF(strings.NewReader(text))
	if err != nil {
		t.Fatal(err)
	}
	return RegionsFromGFF(anno, phaseTestRef)
}

// A CDS with phase 2 is translated from its third base, and the partial codon at its 3' end is left out.
// Changes in the bases left out are reported as nucleotide changes.
func TestGFFPhaseIsApplied(t *testing.T) {
	// CDS 5..13, phase 2: bases 5-6 come before the first complete codon, so the codons are 7-9 (AAA, K)
	// and 10-12 (CCC, P), and base 13 is a partial codon.
	gffText := []byte("##gff-version 3\n##sequence-region ref1 1 20\n" +
		"ref1\ttest\tCDS\t5\t13\t.\t+\t2\tID=cds1;Name=gene1\n" +
		"##FASTA\n>ref1\n" + phaseTestRef + "\n")
	// the query has T5C (a phase base), A7C (AAA -> CAA, K1Q) and G13A (the partial codon)
	msa := []byte(">q\nGGGGCTCAACCCATTTTTTT\n")

	out := new(bytes.Buffer)
	err := Variants(bytes.NewReader(msa), false, "", bytes.NewReader(gffText), "gff", out, -1, -1, false, 0.0, false, 1, true)
	if err != nil {
		t.Fatal(err)
	}
	for _, want := range []string{"aa:gene1:K1Q", "nuc:T5C", "nuc:G13A"} {
		if !strings.Contains(out.String(), want) {
			t.Errorf("want %s in the output, got:\n%s", want, out.String())
		}
	}
}

// On the minus strand the phase is that of the 5' segment, the one with the highest coordinates; the
// phases of the other segments are not used.
func TestGFFPhaseMinusStrandUsesThe5PrimeSegment(t *testing.T) {
	cds, inter, err := regionsFromGFFRows(t,
		"ref1\ttest\tCDS\t3\t8\t.\t-\t2\tID=c1;Name=g1",
		"ref1\ttest\tCDS\t12\t20\t.\t-\t1\tID=c1;Name=g1",
	)
	if err != nil {
		t.Fatal(err)
	}
	if len(cds) != 1 {
		t.Fatalf("got %d CDS, want 1", len(cds))
	}
	// 5' to 3': 20..12 then 8..3 (15 bases); phase 1 drops 20, and the partial codon 4..3 is left out
	wantPos := []int{19, 18, 17, 16, 15, 14, 13, 12, 8, 7, 6, 5}
	if !reflect.DeepEqual(cds[0].Positions, wantPos) {
		t.Errorf("got positions %v, want %v", cds[0].Positions, wantPos)
	}
	if cds[0].Start != 5 || cds[0].Stop != 19 || len(cds[0].Translation) != 4 {
		t.Errorf("got start %d, stop %d, translation %q", cds[0].Start, cds[0].Stop, cds[0].Translation)
	}
	// the bases left out are not coding
	for _, p := range []int{3, 4, 20} {
		found := false
		for _, q := range inter {
			found = found || q == p
		}
		if !found {
			t.Errorf("position %d is not in the non-coding positions", p)
		}
	}
}

// A CDS with phase 0 must be a whole number of codons, and a CDS must have a complete codon after its phase.
func TestGFFPhaseErrors(t *testing.T) {
	if _, _, err := regionsFromGFFRows(t, "ref1\ttest\tCDS\t1\t10\t.\t+\t0\tID=c1;Name=g1"); !errors.Is(err, alphabet.ErrorCDSNotModThree) {
		t.Errorf("phase 0, 10 bases: want ErrorCDSNotModThree, got %v", err)
	}
	if _, _, err := regionsFromGFFRows(t, "ref1\ttest\tCDS\t1\t4\t.\t+\t2\tID=c1;Name=g1"); err == nil {
		t.Errorf("phase 2, 4 bases: want an error, got none")
	}
}

// GFF phase p and GenBank codon_start p+1 give the same coding positions, on either strand.
func TestGFFPhaseMatchesGenbankCodonStart(t *testing.T) {
	for _, tc := range []struct {
		strand, location string
	}{
		{"+", "4..14"},
		{"-", "complement(4..14)"},
	} {
		for phase := 0; phase <= 2; phase++ {
			end := []string{"15", "14", "15"}[phase] // a whole number of codons when the phase is 0
			loc := strings.Replace(tc.location, "14", end, 1)
			cds, _, err := regionsFromGFFRows(t, "ref1\ttest\tCDS\t4\t"+end+"\t.\t"+tc.strand+"\t"+string(rune('0'+phase))+"\tID=c1;Name=g1")
			if err != nil {
				t.Fatalf("%s phase %d: %v", tc.strand, phase, err)
			}
			gb, err := CDSRegionfromGenbank(genbank.GenbankFeature{
				Feature:  "CDS",
				Location: genbank.Location{Representation: loc},
				Info:     map[string][]string{"gene": {"g1"}, "codon_start": {string(rune('1' + phase))}, "translation": {"XXXX"}},
			})
			if err != nil {
				t.Fatalf("%s codon_start %d: %v", tc.strand, phase+1, err)
			}
			if !reflect.DeepEqual(cds[0].Positions, gb.Positions) || cds[0].Strand != gb.Strand {
				t.Errorf("%s phase %d: GFF positions %v (strand %d), GenBank %v (strand %d)", tc.strand, phase, cds[0].Positions, cds[0].Strand, gb.Positions, gb.Strand)
			}
		}
	}
}
