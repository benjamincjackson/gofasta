package variants

import (
	"bytes"
	"fmt"
	"strings"
	"testing"

	"github.com/virus-evolution/gofasta/pkg/encoding"
	"github.com/virus-evolution/gofasta/pkg/genbank"
	"github.com/virus-evolution/gofasta/pkg/gff"
)

// Spliced CDSs, joins that skip or reuse a base (ribosomal frameshifts, as in influenza PA-X, HIV-1 vpr
// and gag-pol) and a /translation that differs from the genome (HIV-1 HXB2 nef), read from small GenBank
// records. The expected calls are worked out by hand from the rules in the README.

// splicedRef has a spliced CDS g, join(1..8,20..29), coding ATG AAA GC|T TTT GGG TAA (M K A F G *): codon
// 3 is split across the junction (positions 7, 8 and 20). Positions 9-19 are the intron, 30-40 non-coding.
const splicedRef = "ATGAAAGC" + "ACACACACACA" + "TTTTGGGTAA" + "GGGGGGGGGGG"

// cdsSpec is one CDS of a test GenBank record
type cdsSpec struct {
	location, gene, translation string
}

// genbankText is a minimal GenBank record of seq with the CDSs given
func genbankText(seq string, cds ...cdsSpec) string {
	var b strings.Builder
	fmt.Fprintf(&b, "LOCUS       TEST %d bp    RNA     linear   VRL 01-JAN-2000\n", len(seq))
	b.WriteString("FEATURES             Location/Qualifiers\n")
	for _, c := range cds {
		fmt.Fprintf(&b, "     CDS             %s\n", c.location)
		fmt.Fprintf(&b, "                     /gene=\"%s\"\n", c.gene)
		b.WriteString("                     /codon_start=1\n")
		fmt.Fprintf(&b, "                     /translation=\"%s\"\n", c.translation)
	}
	b.WriteString("ORIGIN\n")
	for i := 0; i < len(seq); i += 60 {
		fmt.Fprintf(&b, "%9d %s\n", i+1, strings.ToLower(seq[i:min(i+60, len(seq))]))
	}
	b.WriteString("//\n")
	return b.String()
}

// callRegions runs GetVariantsPair on one pairwise alignment (rows with '-' for gaps) with the regions
// given, and returns the formatted variants
func callRegions(t *testing.T, cds []Region, intergenic []int, refRow, queRow string, appendSNP bool) string {
	t.Helper()
	EA := encoding.MakeEncodingArray()
	ref, que := []byte(refRow), []byte(queRow)
	for i := range ref {
		ref[i], que[i] = EA[ref[i]], EA[que[i]]
	}
	off, _ := GetMSAOffsets(ref)
	AS, err := GetVariantsPair(ref, que, "ref", "q", 0, cds, intergenic, off)
	if err != nil {
		t.Fatal(err)
	}
	out := make([]string, 0, len(AS.Vs))
	for _, v := range AS.Vs {
		s, err := FormatVariant(v, appendSNP)
		if err != nil {
			t.Fatal(err)
		}
		out = append(out, s)
	}
	return strings.Join(out, "|")
}

// genbankRegions reads a GenBank record and returns its CDS and intergenic regions
func genbankRegions(t *testing.T, text string, refLen int) ([]Region, []int) {
	t.Helper()
	gb, err := genbank.ReadGenBank(strings.NewReader(text))
	if err != nil {
		t.Fatal(err)
	}
	cds, intergenic, err := RegionsFromGenbank(gb, refLen)
	if err != nil {
		t.Fatal(err)
	}
	return cds, intergenic
}

// gffRegions reads a GFF3 annotation of seq and returns its CDS and intergenic regions
func gffRegions(t *testing.T, rows string, seq string) ([]Region, []int) {
	t.Helper()
	text := "##gff-version 3\n##sequence-region ref 1 " + fmt.Sprint(len(seq)) + "\n" + rows + "##FASTA\n>ref\n" + seq + "\n"
	g, err := gff.ReadGFF(strings.NewReader(text))
	if err != nil {
		t.Fatal(err)
	}
	cds, intergenic, err := RegionsFromGFF(g, seq)
	if err != nil {
		t.Fatal(err)
	}
	return cds, intergenic
}

type joinCase struct {
	name     string
	ref, que string // aligned rows ('-' for gaps)
	want     string // with --append-snps
}

func runJoinCases(t *testing.T, cds []Region, intergenic []int, cases []joinCase) {
	t.Helper()
	for _, c := range cases {
		if got := callRegions(t, cds, intergenic, c.ref, c.que, true); got != c.want {
			t.Errorf("%s:\n got  %s\n want %s", c.name, got, c.want)
		}
	}
}

// splicedCases are the calls for queries of splicedRef with the spliced CDS g (plus strand)
func splicedCases() []joinCase {
	r := splicedRef
	insRef, insQue := ins(r, 20, "AAA")            // between positions 20 and 21, both in g
	insRef2, insQue2 := ins(r, 23, "AAA")          // between codons 4 and 5
	insRefIntron, insQueIntron := ins(r, 8, "AAA") // after the last base of exon 1: in the intron
	return []joinCase{
		{"identical", r, r, ""},
		{"SNP in exon 1, codon split across the junction", r, snp(r, 8, "A"), "aa:g:A3D(nuc:C8A)"},
		{"SNP in exon 2, same codon, synonymous", r, snp(r, 20, "A"), "nuc:T20A"},
		{"SNPs either side of the junction in one codon", r, snp(snp(r, 8, "A"), 20, "G"), "aa:g:A3E(nuc:C8A;nuc:T20G)"},
		{"SNP in the intron", r, snp(r, 12, "G"), "nuc:C12G"},
		{"SNP after the CDS", r, snp(r, 35, "A"), "nuc:G35A"},
		{"in-frame deletion in exon 2", r, del(r, 21, 23), "aa:g:F4-(del:21:3)|del:21:3"},
		{"deletion of intron bases only", r, del(r, 12, 14), "del:12:3"},
		{"deletion from exon 1 into the intron, two coding bases (a frameshift)", r, del(r, 7, 12), "del:7:6"},
		// a deletion that removes a splice site gets no amino acid call, even if its coding bases are a
		// whole number of codons
		{"deletion from exon 1 into the intron, three coding bases", r, del(r, 6, 11), "del:6:6"},
		{"deletion from the intron into exon 2, three coding bases", r, del(r, 17, 22), "del:17:6"},
		// one that ends at the end of exon 1, leaving the splice site, does: ATG AA(T) TTT ... is M N F ...
		{"deletion ending at the end of exon 1", r, del(r, 6, 8), "aa:g:K2N(del:6:3)|del:6:3|aa:g:A3-(del:6:3)"},
		{"insertion inside the split codon's exon-2 side", insRef, insQue, "aa:g:3K(ins:20:3)|ins:20:3"},
		{"insertion between two codons of exon 2", insRef2, insQue2, "aa:g:4K(ins:23:3)|ins:23:3"},
		{"insertion at the junction (in the intron)", insRefIntron, insQueIntron, "ins:8:3"},
	}
}

func TestSplicedCDSPlusStrand(t *testing.T) {
	text := genbankText(splicedRef, cdsSpec{"join(1..8,20..29)", "g", "MKAFG"})
	cds, intergenic := genbankRegions(t, text, len(splicedRef))
	runJoinCases(t, cds, intergenic, splicedCases())
}

// The same CDS in a GFF3 annotation (two CDS rows with one ID) gives the same calls
func TestSplicedCDSPlusStrandGFF(t *testing.T) {
	rows := "ref\tt\tCDS\t1\t8\t.\t+\t0\tID=cds-g;Name=g\n" +
		"ref\tt\tCDS\t20\t29\t.\t+\t1\tID=cds-g;Name=g\n"
	cds, intergenic := gffRegions(t, rows, splicedRef)
	runJoinCases(t, cds, intergenic, splicedCases())
}

// mirror turns a plus-strand case on splicedRef into the same case on its reverse complement, where g
// is complement(join(12..21,33..40)): position p becomes 41-p, bases are complemented, and the residues
// are the same
func mirror(c joinCase) joinCase {
	const n = len(splicedRef) + 1
	comp := map[byte]byte{'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
	var out []string
	for _, e := range strings.Split(c.want, "|") {
		if e == "" {
			continue
		}
		out = append(out, mirrorEntry(e, n, comp))
	}
	return joinCase{c.name, revComp(c.ref), revComp(c.que), strings.Join(out, "|")}
}

// mirrorEntry rewrites one formatted variant (and its notes) for the reverse complement
func mirrorEntry(e string, n int, comp map[byte]byte) string {
	head, notes, hasNotes := strings.Cut(strings.TrimSuffix(e, ")"), "(")
	f := func(s string) string {
		switch {
		case strings.HasPrefix(s, "nuc:"):
			var a, b byte
			var p int
			fmt.Sscanf(s, "nuc:%c%d%c", &a, &p, &b)
			return fmt.Sprintf("nuc:%c%d%c", comp[a], n-p, comp[b])
		case strings.HasPrefix(s, "del:"):
			var p, l int
			fmt.Sscanf(s, "del:%d:%d", &p, &l)
			return fmt.Sprintf("del:%d:%d", n-(p+l-1), l)
		case strings.HasPrefix(s, "ins:"):
			var p, l int
			fmt.Sscanf(s, "ins:%d:%d", &p, &l)
			return fmt.Sprintf("ins:%d:%d", n-p-1, l)
		}
		return s
	}
	if !hasNotes {
		return f(head)
	}
	ns := strings.Split(notes, ";")
	for i := range ns {
		ns[i] = f(ns[i])
	}
	return f(head) + "(" + strings.Join(ns, ";") + ")"
}

func TestSplicedCDSMinusStrand(t *testing.T) {
	ref := revComp(splicedRef)
	text := genbankText(ref, cdsSpec{"complement(join(12..21,33..40))", "g", "MKAFG"})
	cds, intergenic := genbankRegions(t, text, len(ref))
	var cases []joinCase
	for _, c := range splicedCases() {
		cases = append(cases, mirror(c))
	}
	// the mirrored calls come out in the order of their new positions
	for i := range cases {
		got := callRegions(t, cds, intergenic, cases[i].ref, cases[i].que, true)
		if !sameEntries(got, cases[i].want) {
			t.Errorf("%s (minus strand):\n got  %s\n want %s", cases[i].name, got, cases[i].want)
		}
	}
}

// sameEntries reports whether two |-separated lists hold the same entries, in any order
func sameEntries(a, b string) bool {
	x, y := strings.Split(a, "|"), strings.Split(b, "|")
	if len(x) != len(y) {
		return false
	}
	seen := map[string]int{}
	for _, s := range x {
		seen[s]++
	}
	for _, s := range y {
		seen[s]--
	}
	for _, n := range seen {
		if n != 0 {
			return false
		}
	}
	return true
}

// A join that skips a base (a +1 ribosomal frameshift, as in influenza PA-X): base 10 is not in the CDS,
// so a change there is a nuc: entry, and the codons after it are read in the shifted frame.
func TestJoinSkippingABase(t *testing.T) {
	// join(1..9,11..19): ATG AAA CCC | (C) | GGG TTT TAA -> M K P G F *
	ref := "ATGAAACCC" + "C" + "GGGTTTTAA" + "AAAAA"
	text := genbankText(ref, cdsSpec{"join(1..9,11..19)", "g", "MKPGF"})
	cds, intergenic := genbankRegions(t, text, len(ref))
	runJoinCases(t, cds, intergenic, []joinCase{
		{"identical", ref, ref, ""},
		{"SNP at the skipped base", ref, snp(ref, 10, "A"), "nuc:C10A"},
		{"SNP before the skip", ref, snp(ref, 7, "G"), "aa:g:P3A(nuc:C7G)"},
		{"SNP after the skip", ref, snp(ref, 11, "C"), "aa:g:G4R(nuc:G11C)"},
	})
}

// A join that reuses a base (a -1 ribosomal frameshift, as in SARS-CoV-2 ORF1ab and HIV-1 gag-pol): base 9
// is the last base of codon 3 and the first of codon 4, so a change there can change both.
func TestJoinReusingABase(t *testing.T) {
	// join(1..9,9..14): ATG AAA CCG | GTT TAA -> M K P V *   (base 9 is the G of CCG and of GTT)
	ref := "ATGAAACCG" + "TTTAA" + "CCCCC"
	text := genbankText(ref, cdsSpec{"join(1..9,9..14)", "g", "MKPV"})
	cds, intergenic := genbankRegions(t, text, len(ref))
	runJoinCases(t, cds, intergenic, []joinCase{
		{"identical", ref, ref, ""},
		// CCG -> CCA is synonymous (P), GTT -> ATT is V4I
		{"SNP at the reused base", ref, snp(ref, 9, "A"), "aa:g:V4I(nuc:G9A)|nuc:G9A"},
		// CCG -> CCC is synonymous, GTT -> CTT is V4L
		{"another SNP at the reused base", ref, snp(ref, 9, "C"), "aa:g:V4L(nuc:G9C)|nuc:G9C"},
		{"SNP after the reused base", ref, snp(ref, 10, "C"), "aa:g:V4A(nuc:T10C)"},
	})
}

// Amino acids are compared with the annotation's /translation, also where it differs from the genome (as
// HIV-1 HXB2's nef, whose genome has a stop codon at residue 124 where /translation has W): a query with
// the reference's codon is reported as changed, and one with the codon /translation implies is not. The
// GFF path translates the genome itself, so it differs from GenBank at such a codon.
func TestTranslationDiffersFromGenome(t *testing.T) {
	// 1..12: ATG TAG AAA TAA, with /translation "MWK" (the genome's codon 2 is a stop)
	ref := "ATGTAGAAATAA" + "CCCCC"
	text := genbankText(ref, cdsSpec{"1..12", "g", "MWK"})
	cds, intergenic := genbankRegions(t, text, len(ref))
	runJoinCases(t, cds, intergenic, []joinCase{
		{"identical to the reference genome", ref, ref, "aa:g:W2*()"}, // (no SNPs to append)
		{"the codon /translation implies", ref, snp(ref, 5, "G"), "nuc:A5G"},
		{"another stop codon", ref, snp(ref, 6, "A"), "aa:g:W2*(nuc:G6A)"},
	})

	gcds, gintergenic := gffRegions(t, "ref\tt\tCDS\t1\t12\t.\t+\t0\tID=cds-g;Name=g\n", ref)
	runJoinCases(t, gcds, gintergenic, []joinCase{
		{"GFF: identical to the reference genome", ref, ref, ""},
		{"GFF: the codon /translation implies", ref, snp(ref, 5, "G"), "aa:g:*2W(nuc:A5G)"},
	})
}

// the rows used above must be the same width (a check on the helpers)
func TestJoinCaseRowsAreAligned(t *testing.T) {
	for _, c := range splicedCases() {
		if len(c.ref) != len(c.que) || !bytes.Equal([]byte(strings.ReplaceAll(c.ref, "-", "")), []byte(splicedRef)) {
			t.Errorf("%s: rows are not an alignment of splicedRef", c.name)
		}
	}
}

// A reference residue that is unknown (X in /translation, because the reference codon has an ambiguity code,
// as in some RefSeq influenza segments) gives no amino acid change, whatever the query has there; a change
// at a certainly different base is still a nuc: entry, and a deleted unknown residue is still reported. The
// GFF path translates such a codon as X too, so it gives the same calls.
func TestUnknownReferenceResidue(t *testing.T) {
	// 1..15: ATG GCT RAT AAA TAA (M A X K *): codon 3 is AAT (N) or GAT (D)
	ref := "ATGGCTRATAAATAA" + "CCCCC"
	insRef, insQue := ins(snp(ref, 7, "A"), 9, "GGG")
	cases := []joinCase{
		{"identical", ref, ref, ""},
		{"a residue the reference codon allows", ref, snp(ref, 7, "A"), ""},
		{"a residue it does not allow", ref, snp(ref, 7, "C"), "nuc:R7C"},
		{"an insertion next to the unknown residue", insRef, insQue, "aa:g:3G(ins:9:3)|ins:9:3"},
		{"the unknown residue deleted", ref, del(ref, 7, 9), "aa:g:X3-(del:7:3)|del:7:3"},
	}
	cds, intergenic := genbankRegions(t, genbankText(ref, cdsSpec{"1..15", "g", "MAXK"}), len(ref))
	runJoinCases(t, cds, intergenic, cases)

	gcds, gintergenic := gffRegions(t, "ref\tt\tCDS\t1\t15\t.\t+\t0\tID=cds-g;Name=g\n", ref)
	for i := range cases {
		cases[i].name = "GFF: " + cases[i].name
	}
	runJoinCases(t, gcds, gintergenic, cases)
}
