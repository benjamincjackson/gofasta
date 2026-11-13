package sam

import (
	"context"
	"errors"
	"fmt"
	"github.com/virus-evolution/gofasta/pkg/gfio"
	"io"
	"os"
	"path"
	"sort"
	"strconv"
	"strings"
	"unicode/utf8"

	"github.com/virus-evolution/gofasta/pkg/fasta"
	"github.com/virus-evolution/gofasta/pkg/pipeline"

	biogosam "github.com/biogo/hts/sam"
)

// alignPair is a struct for storing an aligned reference & query sequence pair
type alignPair struct {
	ref       []byte
	query     []byte
	refname   string
	queryname string
	idx       int // for retaining input order in the output
}

// insertionOccurence stores information about a single insertion from cigars in a block of SAM records,
// which is used to appropriately gap the other sam lines belonging to the same query sequence
type insertionOccurence struct {
	start     int
	length    int
	rowNumber int
}

// alignedBlockInfo stores arrays of sam fields for one queriy's sam lines -
// they should be equal length and the item at the same index in each array should
// originate from the same sam line
type alignedBlockInfo struct {
	seqpairArray []alignPair      // an array of ref + query seqs from primary + secondary mappings
	cigarArray   []biogosam.Cigar // the cigars for each mapping in seqpairArray
	posArray     []int            // the start POS from the sam file for each mapping in seqpairArray
}

// gapRowForInsertions returns a copy of one row (a reference/query pair from a single SAM record) that
// has been given the columns for every insertion in the block: the row's own insertions keep their
// query bases, and gaps are added for everyone else's. Rows are built from reference coordinates (the
// n-th reference base), so a row's own earlier insertions, and insertions added for other rows, cannot
// shift a later insertion to the wrong column.
//
// refRow and queRow are as returned by getOneLinePlusRef (the reference row has a '-' for each base the
// row has inserted), row is the index of this row in the block, insertions are all of the insertions in
// the block (sorted by start, with this row's own in cigar order) and insertionsAt groups them by start.
// Insertions that start after the end of the row are not added (the row is padded later).
func gapRowForInsertions(refRow, queRow []byte, row int, insertions []insertionOccurence, insertionsAt map[int][]insertionOccurence) ([]byte, []byte) {

	// split the row into columns that are reference bases, and chunks of its own inserted bases.
	// ownAt[c] holds the chunks inserted after the c-th reference base, in order
	baseRef := make([]byte, 0, len(refRow))
	baseQue := make([]byte, 0, len(refRow))
	ownAt := make(map[int][][]byte)

	col := 0
	shift := 0 // the number of this row's own inserted columns so far
	for _, ins := range insertions {
		if ins.rowNumber != row {
			continue
		}
		insCol := ins.start + shift // the column of the first inserted base in this row
		baseRef = append(baseRef, refRow[col:insCol]...)
		baseQue = append(baseQue, queRow[col:insCol]...)
		ownAt[ins.start] = append(ownAt[ins.start], queRow[insCol:insCol+ins.length])
		col = insCol + ins.length
		shift += ins.length
	}
	baseRef = append(baseRef, refRow[col:]...)
	baseQue = append(baseQue, queRow[col:]...)

	// put it back together, with a column for every insertion at each reference position
	outRef := make([]byte, 0, len(refRow))
	outQue := make([]byte, 0, len(queRow))
	used := make(map[int]int) // how many of this row's own chunks at each position have been placed

	for c := 0; c <= len(baseRef); c++ {
		for _, ins := range insertionsAt[c] {
			if ins.rowNumber == row {
				chunk := ownAt[c][used[c]]
				used[c]++
				for range chunk {
					outRef = append(outRef, '-')
				}
				outQue = append(outQue, chunk...)
			} else {
				for k := 0; k < ins.length; k++ {
					outRef = append(outRef, '-')
					outQue = append(outQue, '-')
				}
			}
		}
		if c < len(baseRef) {
			outRef = append(outRef, baseRef[c])
			outQue = append(outQue, baseQue[c])
		}
	}

	return outRef, outQue
}

// blockToSeqPair converts a block of sam records corresponding to the same query sequence to a pairwise
// alignment between reference and query, including insertions in the query relative to the reference
func blockToSeqPair(alignedBlock alignedBlockInfo, ref []byte) alignPair {
	// get all insertions in order of occurrence - earliest first
	//
	// and then collapse the sequences in the block down to one sequence
	//
	// The incoming sequences are all left-aligned. So for every insertion
	// in every cigar in cigarArray, need to gap all the OTHER imcoming sequences
	// except the one that it occurs in (because this has already been done)
	//
	// need to do this in order, left-most insertion first.

	insertions := make([]insertionOccurence, 0)

	for i, cigar := range alignedBlock.cigarArray {

		pos := alignedBlock.posArray[i]

		for _, op := range cigar {
			// populate orderedInsertions here...

			if op.Type().String() == "I" {

				I := insertionOccurence{start: pos, length: op.Len(), rowNumber: i}

				insertions = append(insertions, I)
			}

			if op.Type().Consumes().Reference == 1 {
				pos += op.Len()
			}
		}
	}

	// copy the sequences into new slices
	refSeqArray := make([][]byte, len(alignedBlock.seqpairArray))
	queSeqArray := make([][]byte, len(alignedBlock.seqpairArray))
	for i := range alignedBlock.seqpairArray {
		refSeqArray[i] = alignedBlock.seqpairArray[i].ref
		queSeqArray[i] = alignedBlock.seqpairArray[i].query
	}

	// if there are insertions, every row has to be rebuilt so that they all share the same columns
	if len(insertions) > 0 {
		// A stable sort keeps insertions that start at the same reference position in the order they
		// were found (row order, then cigar order), which is the order their columns are laid out in.
		sort.SliceStable(insertions, func(i, j int) bool { return insertions[i].start < insertions[j].start })

		// insertions grouped by the reference position they follow
		insertionsAt := make(map[int][]insertionOccurence)
		for _, insertion := range insertions {
			insertionsAt[insertion.start] = append(insertionsAt[insertion.start], insertion)
		}

		for j := range refSeqArray {
			refSeqArray[j], queSeqArray[j] = gapRowForInsertions(refSeqArray[j], queSeqArray[j], j, insertions, insertionsAt)
		}
	}

	max := 0

	for _, line := range refSeqArray {
		if len(line) > max {
			max = len(line)
		}
	}

	RBlock := make([][]byte, 0)
	for _, line := range refSeqArray {
		if len(line) < max {
			diff := max - len(line)
			stars := make([]byte, diff)
			for i, _ := range stars {
				stars[i] = '*'
			}
			line = append(line, stars...)
		}
		RBlock = append(RBlock, line)
	}

	R := checkAndGetFlattenedSeq(RBlock, "reference")

	QBlock := make([][]byte, 0)
	for _, line := range queSeqArray {
		if len(line) < max {
			diff := max - len(line)
			stars := make([]byte, diff)
			for i, _ := range stars {
				stars[i] = '*'
			}
			line = append(line, stars...)
		}
		QBlock = append(QBlock, line)
	}

	Q := checkAndGetFlattenedSeq(QBlock, alignedBlock.seqpairArray[0].queryname)

	// extend the alignment to the ref length + the total number of
	// insertions relative to the reference...
	totalInsertionLength := 0
	for _, I := range insertions {
		totalInsertionLength += I.length
	}

	if len(R) < totalInsertionLength+len(ref) {
		diff := (totalInsertionLength + len(ref)) - len(R)
		extendRight := make([]byte, diff)
		for i, _ := range extendRight {
			extendRight[i] = '*'
		}
		Q = append(Q, extendRight...)
		R = append(R, ref[len(ref)-diff:]...)
	}

	Q = swapInNs(Q)

	return alignPair{ref: R, query: Q}
}

// blockToPairwiseAlignment converts the SAM records for one query to a pairwise alignment between that
// query and the reference, with or without the query's insertions relative to the reference.
func blockToPairwiseAlignment(group samRecords, ref []byte, omitIns bool) (alignPair, error) {

	qname := group.records[0].Name
	refname := string(group.records[0].Ref.Name())

	if omitIns {
		Q := make([][]byte, 0, len(group.records))
		for _, line := range group.records {
			seq, _, err := alignRecord(line, ref, len(ref), false)
			if err != nil {
				return alignPair{}, err
			}
			Q = append(Q, seq)
		}
		Qflat := swapInNs(checkAndGetFlattenedSeq(Q, qname))
		return alignPair{ref: ref, query: Qflat, refname: refname, queryname: qname, idx: group.idx}, nil
	}

	infoBlock := alignedBlockInfo{}
	for _, line := range group.records {
		seq, refseq, err := alignRecord(line, ref, len(ref), true)
		if err != nil {
			return alignPair{}, err
		}
		infoBlock.seqpairArray = append(infoBlock.seqpairArray, alignPair{ref: refseq, query: seq, queryname: line.Name})
		infoBlock.cigarArray = append(infoBlock.cigarArray, line.Cigar)
		infoBlock.posArray = append(infoBlock.posArray, line.Pos)
	}

	pair := blockToSeqPair(infoBlock, ref)
	pair.refname = refname
	pair.queryname = qname
	pair.idx = group.idx
	return pair, nil
}

// get the number of bases to add to convert each reference position to MSA coordinates
func getRefOffset(refseq []byte) []int {

	// get the length of the reference sequence without any gaps
	degappedLen := 0
	for _, nuc := range refseq {
		// if there is no alignment gap at this site, ++
		if nuc != '-' {
			degappedLen++
		}
	}

	gapsum := 0
	refToMSA := make([]int, degappedLen)
	for i, nuc := range refseq {
		if nuc == '-' {
			gapsum++
			continue
		}
		refToMSA[i-gapsum] = gapsum
	}

	return refToMSA
}

// trimAlignment trims a pairwise alignment to coordinates in degapped reference space
func trimAlignment(pair alignPair, trimStart int, trimEnd int) alignPair {
	refToMSA := getRefOffset(pair.ref)

	adjTrimStart := trimStart + refToMSA[trimStart-1] - 1
	adjTrimEnd := trimEnd + refToMSA[trimEnd-1]

	pair.query = pair.query[adjTrimStart:adjTrimEnd]
	pair.ref = pair.ref[adjTrimStart:adjTrimEnd]
	return pair
}

func wrap(old string, wrap int) string {
	if wrap <= 0 {
		return old + "\n"
	}
	var (
		written int
		new     []byte
	)
	for {
		if written < len(old) {
			if written+wrap >= len(old) {
				new = append(new, []byte(old[written:]+"\n")...)
				written = written + wrap
			} else {
				new = append(new, []byte(old[written:written+wrap]+"\n")...)
				written = written + wrap
			}
		} else {
			break
		}
	}
	return string(new)
}

// writePair writes a pairwise alignment, with the reference first unless omitRef
func writePair(w io.Writer, AP alignPair, wrapWidth int, omitRef bool) error {
	if !omitRef {
		if _, err := io.WriteString(w, ">"+AP.refname+"\n"+wrap(string(AP.ref), wrapWidth)); err != nil {
			return err
		}
	}
	_, err := io.WriteString(w, ">"+AP.queryname+"\n"+wrap(string(AP.query), wrapWidth))
	return err
}

// pairFileName is the name of the file a query's alignment is written to: its name with '/' (not
// allowed in file names) replaced by '_', cut to 249 bytes (at a character boundary) so that, with
// ".fasta", it fits the usual 255-byte limit. It also reports whether the name was cut.
func pairFileName(queryname string) (string, bool) {
	des := strings.ReplaceAll(queryname, "/", "_")
	if len(des) <= 249 {
		return des + ".fasta", false
	}
	n := 249
	for n > 0 && !utf8.RuneStart(des[n]) {
		n--
	}
	return des[:n] + ".fasta", true
}

// writePairToFile writes a pairwise alignment to its own file in directory p, named after the query
func writePairToFile(p string, AP alignPair, wrapWidth int, omitRef bool) (err error) {
	name, cut := pairFileName(AP.queryname)
	if cut {
		fmt.Fprintf(os.Stderr, "Filename too long, truncating \"%s\" to: \"%s\"\n", strings.ReplaceAll(AP.queryname, "/", "_"), strings.TrimSuffix(name, ".fasta"))
	}
	f, err := os.Create(path.Join(p, name))
	if err != nil {
		return err
	}
	defer func() {
		if cerr := f.Close(); err == nil {
			err = cerr
		}
	}()
	return writePair(f, AP, wrapWidth, omitRef)
}

// oneFilePerName passes results to sink, but makes a query whose file name (see pairFileName) is that of
// an earlier query's a failed record, so it is skipped (or stops the run, if strict) rather than
// overwriting the earlier query's file. Results arrive in input order, so the first query keeps the name.
//
// Names are compared exactly, byte for byte. On a file system that ignores case (the default on macOS
// and Windows) or Unicode normalisation (macOS), two names that differ only in that way are different
// here but name the same file, so the later query's alignment overwrites the earlier one's. Detecting
// that would mean asking the file system (e.g. creating a probe file and looking for it under another
// case), which is not done.
func oneFilePerName(sink func(pipeline.Result[alignPair]) error) func(pipeline.Result[alignPair]) error {
	written := make(map[string]string) // file name -> the query written to it
	return func(r pipeline.Result[alignPair]) error {
		if r.Err == nil {
			name, _ := pairFileName(r.Val.queryname)
			if first, ok := written[name]; ok {
				r.Err = pipeline.Named(r.Val.queryname, fmt.Errorf("its output file %s is already used by %s", name, first))
			} else {
				written[name] = r.Val.queryname
			}
		}
		return sink(r)
	}
}

// ToPairAlign converts a SAM file containing pairwise alignments between assembled genomes into pairwise fasta-format alignments,
// optionally including the reference sequence and insertions relative to it, optionally trimmed to coordinates in (degapped-)reference space.
// The alignments are written to stdout if outpath is "stdout" or "-", otherwise to one file per query in the directory outpath.
// A query whose file name is the same as an earlier query's is skipped (see oneFilePerName); on a file system that ignores
// case, names that differ only in case are not caught, and the later query's file replaces the earlier one's.
// A query whose records can't be used is reported on stderr and skipped, unless strict, when it stops the run.
// threads is the number of goroutines that process records (fewer than 1 means 1).
func ToPairAlign(samIn, ref io.Reader, outpath string, wrap int, trimStart int, trimEnd int, omitRef bool, omitIns bool, threads int, strict bool) error {

	refs, err := fasta.LoadEncodeAlignment(ref, false, false, false)
	if err != nil {
		return err
	}
	if len(refs) != 1 {
		return errors.New("Need one record in --reference")
	}
	refSeq := refs[0].Decode().Seq

	trimStart, trimEnd, trim, err := checkArgs(len(refSeq), trimStart, trimEnd)
	if err != nil {
		return err
	}

	src, err := newSamSource(samIn)
	if err != nil {
		return err
	}
	header := src.header
	if header.Refs()[0].Len() != len(refSeq) {
		return errors.New("the length of the reference sequence (" + strconv.Itoa(len(refSeq)) + ") does not match the length in the sam header (" + strconv.Itoa(header.Refs()[0].Len()) + ")")
	}

	write := func(AP alignPair) error { return writePair(os.Stdout, AP, wrap, omitRef) }
	if !gfio.IsStdout(outpath) {
		if err := os.MkdirAll(outpath, 0755); err != nil {
			return err
		}
		write = func(AP alignPair) error { return writePairToFile(outpath, AP, wrap, omitRef) }
	}

	work := func(group samRecords) (alignPair, error) {
		group, err := prepareQuery(group, len(refSeq))
		if err != nil {
			return alignPair{}, err
		}
		pair, err := blockToPairwiseAlignment(group, []byte(refSeq), omitIns)
		if err != nil {
			return alignPair{}, pipeline.Named(group.records[0].Name, err)
		}
		if trim {
			pair = trimAlignment(pair, trimStart, trimEnd)
		}
		return pair, nil
	}

	h := pipeline.NewHandler(strict, os.Stderr, queryNoun, write)
	sink := h.Sink
	if !gfio.IsStdout(outpath) {
		sink = oneFilePerName(sink)
	}
	if err := pipeline.Run(context.Background(), threads, src.next, work, withNotes(src, sink)); err != nil {
		return err
	}
	src.writeTail()
	return h.Finish()
}
