package sam

import (
	"bufio"
	"bytes"
	"context"
	"errors"
	"fmt"
	"io"
	"os"
	"path"
	"sort"
	"strings"
	"unicode/utf8"

	"github.com/virus-evolution/gofasta/pkg/fasta"
	"github.com/virus-evolution/gofasta/pkg/gfio"
	"github.com/virus-evolution/gofasta/pkg/pipeline"
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

// gapRowForInsertions returns a copy of one row (a reference/query pair from a single SAM record) that
// has been given the columns for every insertion in the block: the row's own insertions keep their
// query bases, and gaps are added for everyone else's. Rows are built from reference coordinates (the
// n-th reference base), so a row's own earlier insertions, and insertions added for other rows, cannot
// shift a later insertion to the wrong column.
//
// refRow and queRow are as returned by alignRecord (the reference row has a '-' for each base the
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
				outRef = append(outRef, bytes.Repeat([]byte{'-'}, len(chunk))...)
				outQue = append(outQue, chunk...)
			} else {
				gaps := bytes.Repeat([]byte{'-'}, ins.length)
				outRef = append(outRef, gaps...)
				outQue = append(outQue, gaps...)
			}
		}
		if c < len(baseRef) {
			outRef = append(outRef, baseRef[c])
			outQue = append(outQue, baseQue[c])
		}
	}

	return outRef, outQue
}

// padTo pads each row to length n with '*' (positions it does not cover)
func padTo(rows [][]byte, n int) [][]byte {
	for i, row := range rows {
		if len(row) < n {
			rows[i] = append(row, bytes.Repeat([]byte{'*'}, n-len(row))...)
		}
	}
	return rows
}

// blockToSeqPair merges the aligned rows of one query's sam records (refRows[i] and queRows[i] from
// records[i], each with its own insertions) into one pairwise alignment between the reference and the
// query, including the query's insertions relative to the reference
func blockToSeqPair(refRows, queRows [][]byte, records []samRecord, ref []byte) alignPair {
	// every insertion in every record. Each row already has its own insertions; the other rows must be
	// gapped for them, in order, left-most first.
	insertions := make([]insertionOccurence, 0)
	for i, r := range records {
		pos := r.Pos
		for _, op := range r.Cigar {
			if op.Type() == cigarInsertion {
				insertions = append(insertions, insertionOccurence{start: pos, length: op.Len(), rowNumber: i})
			}
			if _, r := op.Type().consumes(); r == 1 {
				pos += op.Len()
			}
		}
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

		for j := range refRows {
			refRows[j], queRows[j] = gapRowForInsertions(refRows[j], queRows[j], j, insertions, insertionsAt)
		}
	}

	width := 0
	for _, row := range refRows {
		width = max(width, len(row))
	}
	R := checkAndGetFlattenedSeq(padTo(refRows, width), "reference")
	Q := checkAndGetFlattenedSeq(padTo(queRows, width), records[0].Name)

	// extend the alignment to the reference length plus the total length of the insertions
	totalInsertionLength := 0
	for _, I := range insertions {
		totalInsertionLength += I.length
	}
	if diff := totalInsertionLength + len(ref) - len(R); diff > 0 {
		Q = append(Q, bytes.Repeat([]byte{'*'}, diff)...)
		R = append(R, ref[len(ref)-diff:]...)
	}

	return alignPair{ref: R, query: swapInNs(Q)}
}

// blockToPairwiseAlignment converts the SAM records for one query to a pairwise alignment between that
// query and the reference, with or without the query's insertions relative to the reference.
func blockToPairwiseAlignment(group samRecords, ref []byte, omitIns bool) (alignPair, error) {

	qname := group.records[0].Name
	refname := group.refName

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

	refRows := make([][]byte, len(group.records))
	queRows := make([][]byte, len(group.records))
	for i, line := range group.records {
		var err error
		if queRows[i], refRows[i], err = alignRecord(line, ref, len(ref), true); err != nil {
			return alignPair{}, err
		}
	}

	pair := blockToSeqPair(refRows, queRows, group.records, ref)
	pair.refname = refname
	pair.queryname = qname
	pair.idx = group.idx
	return pair, nil
}

// getRefOffset returns, for each reference position (0-based), the number of columns to add to convert it
// to a column of the aligned reference row refseq (the gaps before it)
func getRefOffset(refseq []byte) []int {
	refToMSA := make([]int, 0, len(refseq))
	gaps := 0
	for _, nuc := range refseq {
		if nuc == '-' {
			gaps++
		} else {
			refToMSA = append(refToMSA, gaps)
		}
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

// wrap returns seq with a newline after every width characters and at the end (one line if width <= 0;
// nothing for an empty seq if width > 0)
func wrap(seq string, width int) string {
	if width <= 0 {
		return seq + "\n"
	}
	var b strings.Builder
	for i := 0; i < len(seq); i += width {
		b.WriteString(seq[i:min(i+width, len(seq))])
		b.WriteByte('\n')
	}
	return b.String()
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
	w := bufio.NewWriter(f)
	if err := writePair(w, AP, wrapWidth, omitRef); err != nil {
		return err
	}
	return w.Flush()
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
func ToPairAlign(samIn, ref io.Reader, outpath string, wrap int, trimStart int, trimEnd int, omitRef bool, omitIns bool, threads int, strict bool) (err error) {

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
	if err := src.checkRefLength(len(refSeq)); err != nil {
		return err
	}

	stdout := bufio.NewWriterSize(os.Stdout, 64<<10)
	defer func() {
		// written even after an error, so that what was written before it still reaches stdout
		if ferr := stdout.Flush(); err == nil {
			err = ferr
		}
	}()
	write := func(AP alignPair) error { return writePair(stdout, AP, wrap, omitRef) }
	if !gfio.IsStdout(outpath) {
		if err := os.MkdirAll(outpath, 0755); err != nil {
			return err
		}
		write = func(AP alignPair) error { return writePairToFile(outpath, AP, wrap, omitRef) }
	}

	work := func(group samRecords) (alignPair, error) {
		group, err := prepareQuery(group, []byte(refSeq), len(refSeq))
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
