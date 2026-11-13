/*
Package variants implements functionality to annotate mutations relative
to a reference sequence for all records in a multiple sequence alignment
in fasta format.
*/
package variants

import (
	"bytes"
	"context"
	"errors"
	"io"
	"os"
	"slices"
	"sort"
	"strconv"
	"strings"
	"sync/atomic"

	"github.com/virus-evolution/gofasta/pkg/alphabet"
	"github.com/virus-evolution/gofasta/pkg/encoding"
	"github.com/virus-evolution/gofasta/pkg/fasta"
	"github.com/virus-evolution/gofasta/pkg/genbank"
	"github.com/virus-evolution/gofasta/pkg/gff"
	"github.com/virus-evolution/gofasta/pkg/pipeline"
)

type Region struct {
	Name        string // name of feature, if it has one
	Start       int    // 1-based 5'-most position of region on the forward strand, inclusive
	Stop        int    // 1-based 3'-most position of region on the forward strand, inclusive
	Translation string // amino acid sequence of this region if it is CDS
	Strand      int    // values in the set {-1, +1} only (and "0" for a mixture?!)
	Positions   []int  // all the (1-based, unadjusted) positions in order, on the reverse strand if needs be
}

// A Variant is a struct that contains information about one mutation (nuc, amino acid, indel) between
// reference and query
type Variant struct {
	Queryname      string
	RefAl          string
	QueAl          string
	Position       int    // (1-based) genomic location (for an amino acid change, the codon's first base in CDS order)
	Codon          [3]int // for an amino acid change, the reference positions of the codon, in CDS order
	Residue        int    // (1-based) amino acid location
	Changetype     string // one of {nuc,aa,ins,del}
	Feature        string // this should be, for example, the name of the CDS that the thing is in
	Length         int    // for indels
	SNPs           string // if this is an amino acid change, what are the snps
	Representation string
}

// AnnoStructs is for passing groups of Variant structs around with an index which is used to retain input
// order in the output
type AnnoStructs struct {
	Queryname string
	Vs        []Variant
	Idx       int
}

// Variants annotates the mutations in each record of a multiple sequence alignment in fasta format,
// relative to a reference sequence: the record refID in the alignment, else the one record of refIn (if
// not nil), else the sequence in the annotation. annoSuffix is the annotation's format, "gb" or "gff", or
// "" for no annotation, when only nucleotide changes, deletions and insertions are reported (and a
// reference must be given). A record that can't be processed (an invalid character, or a different
// width from the reference) is reported on stderr and skipped, unless strict, when it stops the run with
// an error.
// threads is the number of goroutines that process records (fewer than 1 means 1).
func Variants(msaIn io.Reader, stdin bool, refID string, refIn io.Reader, annoIn io.Reader, annoSuffix string, out io.Writer, start int, end int, aggregate bool, threshold float64, appendSNP bool, threads int, strict bool) error {

	var (
		ref fasta.EncodedRecord
		err error
	)

	// Find the reference
	// (Have to move the reader back to the beginning of the alignment, because we are scanning through it twice)
	if refID != "" {
		if x, ok := msaIn.(io.Seeker); ok && !stdin {
			ref, err = findReference(msaIn, refID)
			if err != nil {
				return err
			}
			_, err = x.Seek(0, io.SeekStart)
			if err != nil {
				return err
			}
		}
	}

	next := fasta.NewRecordSource(msaIn)

	// from stdin, the reference must be the first record; it is read here, and handed on again as the
	// first record of the alignment
	if stdin && refID != "" {
		first, err := next()
		if err != nil {
			return err
		}
		if first.ID != refID {
			return errors.New("--reference-id is not the first record in --query")
		}
		ref, err = first.Encode()
		if err != nil {
			return err
		}
		rest := next
		handedOn := false
		next = func() (fasta.Record, error) {
			if !handedOn {
				handedOn = true
				return first, nil
			}
			return rest()
		}
	}

	// or from its own file
	if refID == "" && refIn != nil {
		refs, err := fasta.LoadEncodeAlignment(refIn, false, false, false)
		if err != nil {
			return err
		}
		if len(refs) != 1 {
			return errors.New("--reference must have one record (it has " + strconv.Itoa(len(refs)) + ")")
		}
		ref = refs[0]
	}
	explicitRef := len(ref.Seq) > 0 // (otherwise it comes from the annotation)

	var (
		cdsregions []Region
		intregions []int
		refToMSA   []int
	)

	switch annoSuffix {
	case "":
		if !explicitRef {
			return errors.New("no reference sequence: give --reference, --reference-id or an annotation that contains the reference sequence")
		}

	case "gb":
		gb, err := genbank.ReadGenBank(annoIn)
		if err != nil {
			return err
		}

		// get the reference from the genbank source if required
		if len(ref.Seq) == 0 {
			ref = annotationReference(string(gb.ORIGIN))
		}

		refLenDegapped := len(ref.Decode().Degap().Seq)

		// check that the reference sequence is in the same coordinates as the annotation
		if refLenDegapped != len(gb.ORIGIN) {
			return errors.New("the degapped reference sequence (" + ref.ID + ", " + strconv.Itoa(refLenDegapped) + " bases) is not the same length as the genbank annotation's sequence (" + strconv.Itoa(len(gb.ORIGIN)) + " bases)")
		}

		// get a list of CDS + intergenic regions from the genbank file
		cdsregions, intregions, err = RegionsFromGenbank(gb, refLenDegapped)
		if err != nil {
			return err
		}

		// get the offsets accounting for insertions relative to the reference
		refToMSA, _ = GetMSAOffsets(ref.Seq)

	case "gff":
		gff, err := gff.ReadGFF(annoIn)
		if err != nil {
			return err
		}

		// get the reference from the gff FASTA if required
		if len(ref.Seq) == 0 {
			switch len(gff.FASTA) {
			case 0:
				return errors.New("couldn't find a reference sequence in the --query or the gff")
			case 1:
				for _, v := range gff.FASTA {
					ref = annotationReference(v.Seq)
				}
			default:
				return errors.New("more than one sequence in gff ##FASTA section")
			}
		}

		// get the offsets accounting for insertions relative to the reference
		refToMSA, _ = GetMSAOffsets(ref.Seq)

		refSeqDegapped := ref.Decode().Degap().Seq

		// a reference that was given must be in the same coordinates as the gff's own sequence, if it has one
		if explicitRef && len(gff.FASTA) == 1 {
			for _, v := range gff.FASTA {
				if len(refSeqDegapped) != len(v.Seq) {
					return errors.New("the degapped reference sequence (" + ref.ID + ", " + strconv.Itoa(len(refSeqDegapped)) + " bases) is not the same length as the gff annotation's sequence (" + strconv.Itoa(len(v.Seq)) + " bases)")
				}
			}
		}

		// check that the reference sequence is in the same coordinates as the annotation, if the gff
		// file has a ##sequence-region line
		if len(gff.SequenceRegions) > 1 {
			return errors.New("more than one sequence-region in gff header")
		}
		for _, region := range gff.SequenceRegions {
			if len(refSeqDegapped) != region.End {
				return errors.New("the degapped reference sequence (" + ref.ID + ") is not the same length as the gff annotation")
			}
		}

		// get a list of CDS + intergenic regions from the gff file
		cdsregions, intregions, err = RegionsFromGFF(gff, refSeqDegapped)
		if err != nil {
			return err
		}

	default:
		return errors.New("couldn't tell if --annotation was a .gb or a .gff file")
	}

	work := func(raw fasta.Record) (AnnoStructs, error) {
		record, err := raw.Encode()
		if err != nil {
			return AnnoStructs{}, pipeline.Named(raw.ID, err)
		}
		if len(record.Seq) != len(ref.Seq) {
			return AnnoStructs{}, pipeline.Named(raw.ID, errors.New("Gapped reference sequence and alignment are not the same width"))
		}
		if annoSuffix == "" {
			return AnnoStructs{Queryname: record.ID, Vs: getVariantsNoAnnotation(ref.Seq, record.Seq), Idx: record.Idx}, nil
		}
		AS, err := GetVariantsPair(ref.Seq, record.Seq, ref.ID, record.ID, record.Idx, cdsregions, intregions, refToMSA)
		return AS, pipeline.Named(raw.ID, err)
	}

	var write func(AnnoStructs) error
	var agg *Aggregator
	if aggregate {
		agg, err = NewAggregator(out, start, end, appendSNP, refID)
		if err != nil {
			return err
		}
		write = agg.Add
	} else {
		lw, err := NewLineWriter(out, start, end, appendSNP, refID)
		if err != nil {
			return err
		}
		write = lw.Write
	}

	h := pipeline.NewHandler(strict, os.Stderr, pipeline.Noun{One: "record", Many: "records"}, write)
	if err := pipeline.Run(context.Background(), threads, next, work, h.Sink); err != nil {
		return err
	}
	if err := h.Finish(); err != nil {
		return err
	}
	if aggregate {
		return agg.Write(threshold)
	}
	return nil
}

// annotationReference is the reference sequence taken from the annotation, encoded (a character that
// is not IUPAC is encoded as 0, as by the encoding array, not rejected), with a note on stderr
func annotationReference(seq string) fasta.EncodedRecord {
	EA := encoding.MakeEncodingArray()
	encoded := make([]byte, len(seq))
	for i := range seq {
		encoded[i] = EA[seq[i]]
	}
	os.Stderr.WriteString("using --annotation fasta as reference\n")
	return fasta.EncodedRecord{ID: "annotation_fasta", Seq: encoded}
}

// findReference gets the reference sequence from the msa if it is in there. Only the reference itself
// is checked: the other records are checked when they are processed.
func findReference(msaIn io.Reader, referenceID string) (fasta.EncodedRecord, error) {
	next := fasta.NewRecordSource(msaIn)
	for {
		record, err := next()
		if err == io.EOF || errors.Is(err, fasta.ErrEmptyFasta) {
			return fasta.EncodedRecord{}, errors.New("Couldn't find --reference-id (" + referenceID + ") in --query")
		}
		if err != nil {
			return fasta.EncodedRecord{}, err
		}
		if record.ID == referenceID {
			return record.Encode()
		}
	}
}

func RegionsFromGFF(anno gff.GFF, refSeqDegapped string) ([]Region, []int, error) {

	// rows sharing an ID are one CDS, a row without an ID is a CDS by itself; CDSs are kept in order of
	// their first row, so that CDSs with the same start are reported in file order
	groups := make([][]gff.Feature, 0)
	groupOfID := make(map[string]int)
	for _, f := range anno.Features {
		if !(f.Type == "CDS" || f.Type == "mature_protein_region_of_CDS") {
			continue
		}
		id := f.Attributes["ID"]
		if len(id) == 0 {
			groups = append(groups, []gff.Feature{f})
			continue
		}
		if g, ok := groupOfID[id[0]]; ok {
			groups[g] = append(groups[g], f)
			continue
		}
		groupOfID[id[0]] = len(groups)
		groups = append(groups, []gff.Feature{f})
	}

	tempcds := make([]Region, 0)
	for _, fs := range groups {
		r, err := CDSRegionfromGFF(fs, refSeqDegapped)
		if err != nil {
			return nil, nil, err
		}
		tempcds = append(tempcds, r)
	}

	// get a slide of positions that are not coding based on everything above
	inter := codes(tempcds, len(refSeqDegapped))

	// then make the final coding regions based on what has a name
	cds := make([]Region, 0)
	for _, r := range tempcds {
		if r.Name == "" {
			continue
		}
		cds = append(cds, r)
	}

	// sort by start position
	sort.SliceStable(cds, func(j, k int) bool {
		return cds[j].Start < cds[k].Start
	})

	return cds, inter, nil
}

func CDSRegionfromGFF(fs []gff.Feature, refSeqDegapped string) (Region, error) {
	r := Region{}
	r.Name = gffFeatureName(fs[0])
	// every feature must lie within the reference
	for _, f := range fs {
		if f.Start < 1 || f.End < f.Start || f.End > len(refSeqDegapped) {
			return r, errors.New("Error parsing gff: CDS feature " + strconv.Itoa(f.Start) + ".." + strconv.Itoa(f.End) + " is not within the reference sequence (length " + strconv.Itoa(len(refSeqDegapped)) + ")")
		}
	}
	// the positions in 5' to 3' order, and the phase of the 5' segment (the first on the plus strand, the
	// last on the minus strand); the phases of the other segments are not used
	pos := make([]int, 0)
	var phase int
	switch fs[0].Strand {
	case "+":
		for _, f := range fs {
			if f.Strand != "+" {
				return r, errors.New("Error parsing gff: mixed strands within a single ID")
			}
			for i := f.Start; i <= f.End; i++ {
				pos = append(pos, i)
			}
		}
		r.Strand = 1
		phase = fs[0].Phase
	case "-":
		for j := len(fs) - 1; j >= 0; j-- {
			f := fs[j]
			if f.Strand != "-" {
				return r, errors.New("Error parsing gff: mixed strands within a single ID")
			}
			for i := f.End; i >= f.Start; i-- {
				pos = append(pos, i)
			}
		}
		r.Strand = -1
		phase = fs[len(fs)-1].Phase
	default:
		return r, errors.New("Error parsing gff: protein coding feature needs a strand (+ or -), got \"" + fs[0].Strand + "\"")
	}

	var err error
	r.Positions, err = codingPositions(pos, phase)
	if err != nil {
		return r, err
	}
	r.Start = slices.Min(r.Positions)
	r.Stop = slices.Max(r.Positions)
	feat := make([]byte, len(r.Positions))
	for i, p := range r.Positions {
		feat[i] = refSeqDegapped[p-1]
	}
	refSeqFeat := strings.ToUpper(string(feat))
	if r.Strand == -1 {
		refSeqFeat = alphabet.Complement(refSeqFeat)
	}
	r.Translation, err = alphabet.Translate(refSeqFeat, true)
	if err != nil {
		return r, err
	}

	return r, nil
}

// codingPositions returns the positions of a CDS's complete codons, given all its positions in 5' to 3'
// order and the number of bases before its first complete codon (its GFF phase, or GenBank codon_start
// minus 1). With an offset, the CDS is partial at its 5' end, and a partial codon at its 3' end is left
// out too; without one, its length must be a whole number of codons.
func codingPositions(pos []int, offset int) ([]int, error) {
	if len(pos)-offset < 3 {
		return nil, errors.New("CDS of " + strconv.Itoa(len(pos)) + " bases has no complete codon after its first " + strconv.Itoa(offset))
	}
	coding := pos[offset:]
	if offset > 0 {
		coding = coding[:len(coding)-len(coding)%3]
	} else if len(coding)%3 != 0 {
		return nil, alphabet.ErrorCDSNotModThree
	}
	return coding, nil
}

// gffFeatureName is the name used to report a CDS: its Name attribute, or "" if it has none
// (a CDS without a Name is not reported). Only the first segment of a multi-segment CDS is consulted.
func gffFeatureName(f gff.Feature) string {
	if v := f.Attributes["Name"]; len(v) > 0 {
		return v[0]
	}
	return ""
}

// Parses a genbank flat format file of genome annotations to extract information about the
// the positions of CDS and intergenic regions, in order to annotate mutations within each
func RegionsFromGenbank(gb genbank.Genbank, refLength int) ([]Region, []int, error) {

	if refLength < 0 {
		return nil, nil, errors.New("negative reference length")
	}

	cds := make([]Region, 0)
	for _, f := range gb.FEATURES {
		if f.Feature == "CDS" {
			REGION, err := CDSRegionfromGenbank(f)
			if err != nil {
				return nil, nil, err
			}
			// every position must lie within the reference
			if len(REGION.Positions) > 0 && (slices.Min(REGION.Positions) < 1 || slices.Max(REGION.Positions) > refLength) {
				return nil, nil, errors.New("Genbank CDS feature \"" + REGION.Name + "\" is not within the reference sequence (length " + strconv.Itoa(refLength) + ")")
			}
			cds = append(cds, REGION)
		}
	}

	// Get a slice of the intergenic regions
	inter := codes(cds, refLength)

	return cds, inter, nil
}

func CDSRegionfromGenbank(f genbank.GenbankFeature) (Region, error) {

	if !f.HasAttribute("gene") {
		return Region{}, errors.New("No \"gene\" attibute in Genbank CDS feature")
	}
	if !f.HasAttribute("codon_start") {
		return Region{}, errors.New("No \"codon_start\" attibute in Genbank CDS feature")
	}

	if !f.HasAttribute("translation") {
		return Region{}, errors.New("No \"translation\" attibute in Genbank CDS feature")
	}

	r := Region{
		Name:        f.Value("gene"),
		Translation: f.Value("translation") + "*",
	}

	temp, err := f.Location.GetPositions()
	if err != nil {
		return Region{}, err
	}
	codon_start, err := strconv.Atoi(f.Value("codon_start"))
	if err != nil || codon_start < 1 || codon_start > 3 || codon_start > len(temp) {
		return Region{}, errors.New("Genbank CDS feature has an invalid \"codon_start\" attibute: " + f.Value("codon_start"))
	}
	r.Positions, err = codingPositions(temp, codon_start-1)
	if err != nil {
		return Region{}, err
	}
	if len(r.Translation) < len(r.Positions)/3 {
		return Region{}, errors.New("Genbank CDS feature \"" + r.Name + "\": the translation is shorter than the CDS")
	}

	r.Start = slices.Min(r.Positions)
	r.Stop = slices.Max(r.Positions)

	reverse, err := f.Location.IsReverse()
	if err != nil {
		return Region{}, err
	}
	r.Strand = 1
	if reverse {
		r.Strand = -1
	}

	return r, nil
}

// get a single slice of intergenic positions after parsing genbank or gff for protein-coding regions
// TO DO - return it in MSA coordinates?
// TO DO - just give it the length of the reference sequence?
func codes(proteincoding []Region, refLength int) []int {

	// true/false this site codes for a protein:
	codes := make([]bool, refLength) // initialises as falses
	for _, feature := range proteincoding {
		for _, pos := range feature.Positions {
			codes[pos-1] = true
		}
	}
	// then iterate over the genome and retain 1-based positions of sites that don't code for a protein
	intergenicregions := make([]int, 0)
	for i, b := range codes {
		if !b {
			intergenicregions = append(intergenicregions, i+1)
		}
	}

	return intergenicregions
}

// indelInSpan reports whether one of the indels (from getIndelsPair) touches start..stop: a deletion
// with a deleted base there, or an insertion after a base there. getIndelAAsPair can only find
// something in a CDS for which this is true, so it is skipped for the others.
func indelInSpan(indels []Variant, start, stop int) bool {
	for _, v := range indels {
		switch v.Changetype {
		case "del":
			if v.Position <= stop && v.Position+v.Length-1 >= start {
				return true
			}
		case "ins":
			if v.Position >= start && v.Position <= stop {
				return true
			}
		}
	}
	return false
}

// outsideRange reports whether a variant has no reference base in start..end: the bases of its codon for
// an amino acid change, the deleted bases for a deletion, and its position otherwise (for an insertion,
// the base it follows). A bound that is not > 0 is not applied, so either bound can be used on its own.
func outsideRange(v Variant, start, end int) bool {
	out := func(pos int) bool { return (start > 0 && pos < start) || (end > 0 && pos > end) }
	switch v.Changetype {
	case "aa":
		return out(v.Codon[0]) && out(v.Codon[1]) && out(v.Codon[2])
	case "del":
		return (start > 0 && v.Position+v.Length-1 < start) || (end > 0 && v.Position > end)
	}
	return out(v.Position)
}

// GetMSAOffsets returns two arrays which contain coordinate shifting information that can be used to
// convert reference to msa coordinates, and vice versa: refToMSA, the number of columns to add to each
// (0-based) reference position, and MSAToRef, the number to subtract from each column (0 for a column
// that is a gap in the reference)
func GetMSAOffsets(refseq []byte) ([]int, []int) {
	refToMSA := make([]int, 0, len(refseq))
	MSAToRef := make([]int, len(refseq))
	gaps := 0
	for i, nuc := range refseq {
		if nuc == 244 { // 244 is a gap ('-')
			gaps++
			continue
		}
		refToMSA = append(refToMSA, gaps)
		MSAToRef[i] = gaps
	}
	return refToMSA, MSAToRef
}

// zeroOffsets is a slice of zeros shared by RefToMSAOffsets
var zeroOffsets atomic.Pointer[[]int]

// RefToMSAOffsets returns GetMSAOffsets's first result, refToMSA, without the second. If refseq has no
// gaps every offset is 0, and the slice returned is shared: callers must not change it.
func RefToMSAOffsets(refseq []byte) []int {
	if bytes.IndexByte(refseq, 244) >= 0 {
		refToMSA, _ := GetMSAOffsets(refseq)
		return refToMSA
	}
	n := len(refseq)
	if z := zeroOffsets.Load(); z != nil && len(*z) >= n {
		return (*z)[:n:n]
	}
	z := make([]int, n)
	zeroOffsets.Store(&z)
	return z
}

func GetVariantsPair(ref, query []byte, refID, queryID string, idx int, cdsregions []Region, intregions []int, offsetRefCoord []int) (AnnoStructs, error) {

	indels := getIndelsPair(ref, query)
	nucs := getNucsPair(ref, query, intregions, offsetRefCoord)
	AAs := make([]Variant, 0)
	for _, r := range cdsregions {
		var handled map[int]bool
		if indelInSpan(indels, r.Start, r.Stop) {
			var indelAAs []Variant
			indelAAs, handled = getIndelAAsPair(ref, query, r, offsetRefCoord)
			AAs = append(AAs, indelAAs...)
		}
		AAs = append(AAs, getAAsPair(ref, query, r, offsetRefCoord, handled)...)
	}

	// (never nil, even with no variants)
	variants := append(append(append(make([]Variant, 0, len(indels)+len(nucs)+len(AAs)), indels...), nucs...), AAs...)
	sort.SliceStable(variants, func(i, j int) bool {
		return variants[i].Position < variants[j].Position || (variants[i].Position == variants[j].Position && variants[i].Changetype < variants[j].Changetype)
	})

	// there might be duplicates (adjacent, once sorted) if there was a snp in the region of a join()
	return AnnoStructs{Queryname: queryID, Vs: slices.Compact(variants), Idx: idx}, nil
}

// FormatVariant returns a string representation of a single mutation, the format
// of which varies given its type (aa/nuc/indel)
func FormatVariant(v Variant, appendSNP bool) (string, error) {
	var s string

	switch v.Changetype {
	case "del":
		s = "del:" + strconv.Itoa(v.Position) + ":" + strconv.Itoa(v.Length)
	case "ins":
		s = "ins:" + strconv.Itoa(v.Position) + ":" + strconv.Itoa(v.Length)
	case "nuc":
		s = "nuc:" + v.RefAl + strconv.Itoa(v.Position) + v.QueAl
	case "aa":
		s = "aa:" + v.Feature + ":" + v.RefAl + strconv.Itoa(v.Residue) + v.QueAl
		if appendSNP {
			s += "(" + v.SNPs + ")"
		}
	default:
		return "", errors.New("couldn't parse variant type")
	}

	return s, nil
}

// LineWriter writes each query's mutations as one line of CSV.
type LineWriter struct {
	w          io.Writer
	start, end int
	appendSNP  bool
	refID      string
}

// NewLineWriter writes the CSV header to w and returns a LineWriter that writes to w. Only mutations
// within start..end (either may be -1 for no limit) are written, and the reference itself (the
// record named refID) is not.
func NewLineWriter(w io.Writer, start, end int, appendSNP bool, refID string) (*LineWriter, error) {
	if _, err := w.Write([]byte("query,mutations\n")); err != nil {
		return nil, err
	}
	return &LineWriter{w: w, start: start, end: end, appendSNP: appendSNP, refID: refID}, nil
}

// Write writes one query's mutations.
func (lw *LineWriter) Write(AS AnnoStructs) error {
	if AS.Queryname == lw.refID {
		return nil
	}
	sa := make([]string, 0, len(AS.Vs))
	for _, v := range AS.Vs {
		if outsideRange(v, lw.start, lw.end) {
			continue
		}
		newVar, err := FormatVariant(v, lw.appendSNP)
		if err != nil {
			return err
		}
		sa = append(sa, newVar)
	}
	_, err := lw.w.Write([]byte(AS.Queryname + "," + strings.Join(sa, "|") + "\n"))
	return err
}

// Aggregator counts how many queries have each mutation, and writes their frequencies.
type Aggregator struct {
	w          io.Writer
	start, end int
	appendSNP  bool
	refID      string
	counts     map[Variant]float64
	n          float64
}

// NewAggregator writes the CSV header to w and returns an Aggregator that will write to w. Only
// mutations within start..end (either may be -1 for no limit) are counted, and the reference itself
// (the record named refID) is not counted.
func NewAggregator(w io.Writer, start, end int, appendSNP bool, refID string) (*Aggregator, error) {
	if _, err := w.Write([]byte("mutation,frequency\n")); err != nil {
		return nil, err
	}
	return &Aggregator{w: w, start: start, end: end, appendSNP: appendSNP, refID: refID, counts: make(map[Variant]float64)}, nil
}

// Add counts one query's mutations.
func (a *Aggregator) Add(AS AnnoStructs) error {
	if AS.Queryname == a.refID {
		return nil
	}
	a.n++
	for _, v := range AS.Vs {
		if outsideRange(v, a.start, a.end) {
			continue
		}
		rep, err := FormatVariant(v, a.appendSNP)
		if err != nil {
			return err
		}
		Vskinny := Variant{RefAl: v.RefAl, QueAl: v.QueAl, Position: v.Position, Residue: v.Residue, Changetype: v.Changetype, Feature: v.Feature, Length: v.Length, Representation: rep}
		a.counts[Vskinny]++
	}
	return nil
}

// Write writes the frequency of each mutation that is at least threshold, ordered by position and then
// by the other fields of the mutation.
func (a *Aggregator) Write(threshold float64) error {
	order := make([]Variant, 0, len(a.counts))
	for k := range a.counts {
		order = append(order, k)
	}

	// The entries came out of a map, so the order must be fully determined by the key: every field
	// that can differ between two entries has to be compared, otherwise ties come out in random order.
	sort.Slice(order, func(i, j int) bool {
		a, b := order[i], order[j]
		switch {
		case a.Position != b.Position:
			return a.Position < b.Position
		case a.Changetype != b.Changetype:
			return a.Changetype < b.Changetype
		case a.QueAl != b.QueAl:
			return a.QueAl < b.QueAl
		case a.Length != b.Length:
			return a.Length < b.Length
		case a.Feature != b.Feature:
			return a.Feature < b.Feature
		case a.Residue != b.Residue:
			return a.Residue < b.Residue
		case a.RefAl != b.RefAl:
			return a.RefAl < b.RefAl
		}
		return a.Representation < b.Representation
	})

	for _, V := range order {
		if a.counts[V]/a.n < threshold {
			continue
		}
		_, err := a.w.Write([]byte(V.Representation + "," + strconv.FormatFloat(a.counts[V]/a.n, 'f', 9, 64) + "\n"))
		if err != nil {
			return err
		}
	}
	return nil
}
