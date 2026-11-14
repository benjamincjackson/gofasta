/*
Package sam provides functionality to handle sam format files.

It was written to process pairwise alignments between assembled SARS-CoV-2
genomes from minimap2.

Its main purpose is to generate fasta format alignments. It also has a routine
to list mutations between pairs of sequences.
*/
package sam

import (
	"bufio"
	"bytes"
	"errors"
	"fmt"
	"io"
	"os"
	"strconv"
	"sync"
	"sync/atomic"

	"github.com/virus-evolution/gofasta/pkg/pipeline"
)

// samRecords is a struct that carries a group of sam lines (belonging to the
// same sequence, probably) and a integer index which is used to keep track of
// the order of the input when we parallelise
type samRecords struct {
	name    string // the query's name (QNAME)
	records []samRecord
	refName string // the name of the reference (the first @SQ) they are mapped to
	// the query's primary lines (FLAG & 0x900 == 0, counting lines skipped as unmapped or mapped to another
	// reference): the SAM specification requires exactly one
	primaries int
	equals    bool // SEQ had "=" (the reference base), replaced by prepareQuery
	idx       int
	err       error // why one of the query's lines could not be parsed, if one could not
}

// RecordName returns the query's name (see pipeline.Namer)
func (g samRecords) RecordName() string { return g.name }

// isLetter reports whether b is an ASCII letter (a base or an ambiguity code, as opposed to '-' or '*')
func isLetter(b byte) bool { return 'A' <= b && b <= 'Z' || 'a' <= b && b <= 'z' }

// getNucFromSite flattens a site to a single nucleotide when a query sequence
// has secondary mappings (multiple records/lines) in the SAM file.
// * If there is more than one distinct alphabetic character at this site, an N is returned (with a
// message on stderr listing the distinct characters, in byte order).
// * Otherwise the largest byte at the site is returned, so alphabetic characters override '-'s and '*'s.
func getNucFromSite(s []byte, qname string, pos int) byte {
	var m byte
	var seen [256]bool
	letters := 0
	for _, e := range s {
		if e > m {
			m = e
		}
		if !seen[e] {
			seen[e] = true
			if isLetter(e) {
				letters++
			}
		}
	}

	if letters > 1 {
		distinct := make([]byte, 0, len(s))
		for b := range seen {
			if seen[b] {
				distinct = append(distinct, byte(b))
			}
		}
		os.Stderr.WriteString("ambiguous overlapping alignment: " + qname + ": " + strconv.Itoa(pos+1) + ": " + string(distinct) + "\n")
		return 'N'
	}

	return m
}

// checkAndGetFlattenedSeq applies getNucFromSite over all sites in a block
// of SAM records to get a single flattened sequence for one query
func checkAndGetFlattenedSeq(block [][]byte, qname string) []byte {

	seq := make([]byte, len(block[0]))
	if len(block) == 1 {
		copy(seq, block[0])
		return seq
	}
	site := make([]byte, len(block))

	for j := range block[0] {
		for i := range block {
			site[i] = block[i][j]
		}
		seq[j] = getNucFromSite(site, qname, j)
	}

	return seq
}

// getSeqFromBlock gets the aligned sequence of one query from its SAM records, flattening them into one
// sequence if there is more than one
func getSeqFromBlock(records []samRecord, refLen int, includeInsertions bool) ([]byte, error) {
	block := make([][]byte, len(records))
	for i, line := range records {
		row, _, err := alignRecord(line, nil, refLen, includeInsertions)
		if err != nil {
			return nil, err
		}
		block[i] = row
	}
	return checkAndGetFlattenedSeq(block, records[0].Name), nil
}

// swapInNs replaces unmapped positions with Ns
func swapInNs(seq []byte) []byte {
	for i, L := range seq {
		if L == '*' {
			seq[i] = 'N'
		}
	}
	return seq
}

// swapInGapsNs replaces internal unmapped positions with Ns and external
// unmapped positions with '-'s
func swapInGapsNs(seq []byte) []byte {
	// with no letters, every position is before the first letter and after the last one
	firstLetterIndx := len(seq)
	lastLetterIndx := -1

	for i, L := range seq {
		if isLetter(L) {
			if firstLetterIndx == len(seq) {
				firstLetterIndx = i
			}
			lastLetterIndx = i
		}
	}

	for i, L := range seq {
		if L != '*' {
			continue
		}
		if i > firstLetterIndx && i < lastLetterIndx {
			seq[i] = 'N'
		} else {
			seq[i] = '-'
		}
	}

	return seq
}

// validateRecord checks that a mapped sam record is internally consistent and fits the reference, so that
// the code that builds alignments from CIGAR strings can index its sequence without going out of range.
func validateRecord(rec *samRecord, refLen int) error {

	bad := func(msg string) error {
		return errors.New("bad sam record: " + msg)
	}

	if rec.Ref < 0 {
		return bad("no reference sequence (RNAME is *) but not flagged as unmapped")
	}
	if rec.Pos < 0 {
		return bad("negative position")
	}

	refSpan, queryLen := cigarLengths(rec.Cigar)
	refEnd := rec.Pos + refSpan

	if queryLen != len(rec.Seq) {
		return bad("the CIGAR string describes " + strconv.Itoa(queryLen) + " query bases but SEQ has " + strconv.Itoa(len(rec.Seq)))
	}
	if refEnd > refLen {
		return bad("the alignment extends beyond the end of the reference (" + strconv.Itoa(refEnd) + " > " + strconv.Itoa(refLen) + ")")
	}

	return nil
}

// samSource reads a sam file one query at a time (see newSamSource). The messages about skipped lines are
// kept with the query they were read with, so that they can be written in input order with the
// pipeline's own warnings (see withNotes).
type samSource struct {
	header samHeader
	next   func() (samRecords, error)

	mu    sync.Mutex
	notes map[int][]string // by the index of the query being read when the line was skipped
	tail  []string         // lines skipped after the last query

	otherSeq atomic.Int64 // queries with a record whose SEQ had characters the specification does not allow
}

func (s *samSource) note(idx int, msg string) {
	s.mu.Lock()
	defer s.mu.Unlock()
	s.notes[idx] = append(s.notes[idx], msg)
}

// withNotes wraps a sink so that the messages about the lines skipped while reading each query are written
// to stderr before the query's result is handled.
func withNotes[T any](s *samSource, sink func(pipeline.Result[T]) error) func(pipeline.Result[T]) error {
	return func(r pipeline.Result[T]) error {
		s.mu.Lock()
		msgs := s.notes[r.Idx]
		delete(s.notes, r.Idx)
		s.mu.Unlock()
		for _, m := range msgs {
			os.Stderr.WriteString(m)
		}
		return sink(r)
	}
}

// writeTail writes the messages about the lines skipped after the last query
func (s *samSource) writeTail() {
	s.mu.Lock()
	defer s.mu.Unlock()
	for _, m := range s.tail {
		os.Stderr.WriteString(m)
	}
	s.tail = nil
	if n := s.otherSeq.Load(); n > 0 {
		fmt.Fprintf(os.Stderr, "warning: %d %s had characters in SEQ that the SAM specification does not allow (it allows letters, '=' and '.'), read as N\n",
			n, map[bool]string{true: "query", false: "queries"}[n == 1])
	}
}

// errEmptySam is the error for a sam file with nothing in it
var errEmptySam = errors.New("the sam file is empty")

// newSamSource reads the header of a sam file and returns a samSource whose next function returns the
// records for the next query each time it is called, and io.EOF after the last query. The records for a query are
// consecutive lines with the same QNAME; unmapped and secondary records, and records mapped to any
// reference other than the first @SQ, are skipped (with a message on stderr), and blank lines are
// ignored. A line that can't be parsed is recorded as the error of its query (whose other lines are still
// read, so that the next query starts in the right place). The records are not checked here: see
// prepareQuery. An error from the returned function means the file could not be read.
func newSamSource(sam io.Reader) (*samSource, error) {

	br := bufio.NewReaderSize(sam, 1<<20)

	// readLine returns the next line without its line ending (a final line need not have one), or io.EOF
	readLine := func() ([]byte, error) {
		line, err := br.ReadBytes('\n')
		if err == io.EOF && len(line) > 0 {
			err = nil
		}
		if err != nil {
			return nil, err
		}
		line = bytes.TrimSuffix(line, []byte{'\n'})
		line = bytes.TrimSuffix(line, []byte{'\r'})
		return line, nil
	}

	if _, err := br.Peek(1); err == io.EOF {
		return nil, errEmptySam
	} else if err != nil {
		return nil, err
	}
	var headerText []byte
	for {
		p, err := br.Peek(1)
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, err
		}
		if p[0] != '@' {
			break
		}
		line, err := readLine()
		if err != nil {
			return nil, err
		}
		headerText = append(append(headerText, line...), '\n')
	}
	h, err := parseSamHeader(headerText)
	if err != nil {
		return nil, err
	}
	if len(h.refs) == 0 {
		return nil, errors.New("no reference sequence (@SQ line) in the sam header")
	}

	src := &samSource{header: *h, notes: make(map[int][]string)}
	counter := 0     // the index of the query being read
	started := false // whether a line of that query has been read
	// a skipped line is listed with the query that comes after it in the file: the one being read if
	// none of its lines has been read yet, otherwise (as far as we know yet) the next one
	noteIdx := func() int {
		if started {
			return counter + 1
		}
		return counter
	}

	// the primary lines read for each QNAME since the current query started (a query's lines are
	// consecutive, so only this query's and the next one's are needed)
	primaries := map[string]int{}

	// a parsed line: its QNAME, and the record or why it could not be parsed
	type samLine struct {
		name string
		rec  samRecord
		err  error
	}

	// readUsable returns the next line that is not skipped, or io.EOF
	readUsable := func() (samLine, error) {
		for {
			line, err := readLine()
			if err != nil {
				return samLine{}, err
			}
			if len(line) == 0 {
				continue
			}
			name, _, _ := bytes.Cut(line, []byte{'\t'})
			rec, err := h.parseRecord(line)
			if err != nil {
				return samLine{name: string(name), err: err}, nil
			}
			if rec.primary() {
				primaries[rec.Name]++
			}
			// the third bit (== 4) in the sam flag is set if the read is unmapped
			if rec.Flags&flagUnmapped != 0 {
				src.note(noteIdx(), "skipping unmapped read: "+rec.Name+"\n")
				continue
			}
			// the 9th bit (== 256) in the sam flag is set if the mapping is secondary
			if rec.Flags&flagSecondary != 0 {
				src.note(noteIdx(), "ignoring secondary mapping: "+rec.Name+"\n")
				continue
			}
			// only one reference sequence is supported: records mapped to any other @SQ are skipped
			if rec.Ref > 0 {
				src.note(noteIdx(), "skipping record mapped to a reference other than the first @SQ ("+h.refs[rec.Ref].Name+"): "+rec.Name+"\n")
				continue
			}
			return samLine{name: string(name), rec: rec}, nil
		}
	}

	var pending *samLine // the first line of the next query, once it has been read
	finished := false

	next := func() (samRecords, error) {
		group := samRecords{idx: counter, refName: h.refs[0].Name}
		started = false
		otherSeq := false
		for !finished {
			var l samLine
			if pending != nil {
				l = *pending
				pending = nil
			} else {
				l, err = readUsable()
				if err == io.EOF {
					finished = true
					break
				}
				if err != nil {
					return samRecords{}, err
				}
			}
			if started && l.name != group.name {
				pending = &l
				break
			}
			started = true
			group.name = l.name
			if l.err != nil {
				if group.err == nil {
					group.err = l.err
				}
				continue
			}
			group.records = append(group.records, l.rec)
			otherSeq = otherSeq || l.rec.otherSeq
		}
		if otherSeq {
			src.otherSeq.Add(1)
		}
		group.primaries = primaries[group.name]
		// (the lines of the next query read so far, its first usable line included, stay counted)
		next := 0
		if pending != nil {
			next = primaries[pending.name]
		}
		clear(primaries)
		if next > 0 {
			primaries[pending.name] = next
		}
		if !started {
			src.mu.Lock()
			src.tail = append(src.tail, src.notes[counter]...)
			delete(src.notes, counter)
			src.mu.Unlock()
			return samRecords{}, io.EOF
		}
		counter++
		return group, nil
	}

	src.next = next
	return src, nil
}

// checkRefLength checks that a reference sequence of length n is as long as the sam header's (first) @SQ
func (s *samSource) checkRefLength(n int) error {
	if l := s.header.refs[0].Length; l != n {
		return errors.New("the length of the reference sequence (" + strconv.Itoa(n) + ") does not match the length in the sam header (" + strconv.Itoa(l) + ")")
	}
	return nil
}

// prepareQuery checks the records for one query against a reference of length refLen (its bases, ref, may be
// nil), replaces each "=" in their SEQ by the reference base (or by N if ref is nil; group.equals says
// whether there were any), and resolves the overlaps between them. The error, if any, carries the query's
// name.
func prepareQuery(group samRecords, ref []byte, refLen int) (samRecords, error) {
	if group.err != nil {
		return samRecords{}, pipeline.Named(group.name, group.err)
	}
	if group.primaries != 1 {
		return samRecords{}, pipeline.Named(group.name, errors.New("bad sam records: the query has "+
			strconv.Itoa(group.primaries)+" primary lines (FLAG & 0x900 == 0), and the SAM specification requires exactly one"))
	}
	for i := range group.records {
		if err := validateRecord(&group.records[i], refLen); err != nil {
			return samRecords{}, pipeline.Named(group.name, err)
		}
		had, err := group.records[i].resolveEquals(ref)
		if err != nil {
			return samRecords{}, pipeline.Named(group.name, err)
		}
		group.equals = group.equals || had
	}
	group.records = resolveOverlaps(group.records)
	return group, nil
}

// queryNoun names the records of the sam pipelines (one query's records) in warnings
var queryNoun = pipeline.Noun{One: "query", Many: "queries"}
