package variants

import (
	"os"
	"strings"
	"testing"

	"github.com/virus-evolution/gofasta/pkg/fasta"
)

const e2eInputs = "../../cmd/testdata/e2e/inputs/"

func readInput(t *testing.T, name string) string {
	t.Helper()
	b, err := os.ReadFile(e2eInputs + name)
	if err != nil {
		t.Fatal(err)
	}
	return string(b)
}

func loadRef(t *testing.T, name string) fasta.EncodedRecord {
	t.Helper()
	refs, err := fasta.LoadEncodeAlignment(strings.NewReader(readInput(t, name)), false, false, false)
	if err != nil {
		t.Fatal(err)
	}
	return refs[0]
}

func TestLoadAnnotation(t *testing.T) {
	gb := readInput(t, "ref.gb")
	gffNoFasta := readInput(t, "ref.gff")
	gffText := gffNoFasta + "##FASTA\n" + readInput(t, "ref.fasta")
	gffTwoFasta := gffText + strings.Replace(readInput(t, "ref.fasta"), ">ref", ">second", 1)

	tests := []struct {
		name    string
		anno    string
		format  string
		ref     fasta.EncodedRecord
		wantErr string // "" for success
	}{
		{"genbank, its own sequence", gb, "gb", fasta.EncodedRecord{}, ""},
		{"genbank, given reference", gb, "gb", loadRef(t, "ref.fasta"), ""},
		{"genbank, letter that is not IUPAC", strings.Replace(gb, "gacagaagct", "gacagaxgct", 1), "gb", fasta.EncodedRecord{}, `the genbank annotation's sequence: invalid nucleotide in fasta record ("x")`},
		{"genbank, short reference", gb, "gb", loadRef(t, "ref_short.fasta"), "is not the same length as the genbank annotation's sequence (120 bases)"},
		{"gff, its own sequence", gffText, "gff", fasta.EncodedRecord{}, ""},
		{"gff, given reference", gffText, "gff", loadRef(t, "ref.fasta"), ""},
		{"gff, short reference", gffText, "gff", loadRef(t, "ref_short.fasta"), "is not the same length as the gff annotation's sequence (120 bases)"},
		{"gff without sequence, given reference", gffNoFasta, "gff", loadRef(t, "ref.fasta"), ""},
		{"gff without sequence, short reference", gffNoFasta, "gff", loadRef(t, "ref_short.fasta"), "(ref_short) is not the same length as the gff annotation"},
		{"gff without sequence or reference", gffNoFasta, "gff", fasta.EncodedRecord{}, "couldn't find a reference sequence in the gff, and none was given"},
		{"gff with two sequences, no reference", gffTwoFasta, "gff", fasta.EncodedRecord{}, "more than one sequence in gff ##FASTA section"},
		{"gff with two sequences, given reference", gffTwoFasta, "gff", loadRef(t, "ref.fasta"), ""},
	}
	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			ref, cds, _, err := LoadAnnotation(strings.NewReader(tc.anno), tc.format, tc.ref)
			if tc.wantErr != "" {
				if err == nil || !strings.Contains(err.Error(), tc.wantErr) {
					t.Fatalf("got error %v, want one containing %q", err, tc.wantErr)
				}
				return
			}
			if err != nil {
				t.Fatal(err)
			}
			if len(ref.Seq) != 120 || len(cds) == 0 {
				t.Errorf("got a reference of %d bases and %d CDS regions", len(ref.Seq), len(cds))
			}
			if len(tc.ref.Seq) == 0 && ref.ID != "annotation_fasta" {
				t.Errorf("reference ID %q, want annotation_fasta", ref.ID)
			}
		})
	}
}
