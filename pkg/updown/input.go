package updown

import (
	"bufio"
	"bytes"
	"context"
	"encoding/binary"
	"encoding/csv"
	"errors"
	"io"
	"os"
	"slices"
	"strconv"
	"strings"

	"github.com/virus-evolution/gofasta/pkg/encoding"
	"github.com/virus-evolution/gofasta/pkg/fasta"
	"github.com/virus-evolution/gofasta/pkg/pipeline"
)

// getAmbArr parses the ambiguities field from one line of the output of gofasta
// updown list to an array of 1-based inclusive start-stop pairs of integers,
// which represent tracts of ambiguities. A tract is "start-stop", or "pos" for one site.
func getAmbArr(s string) ([]int, error) {
	A := make([]int, 0)
	if len(s) == 0 {
		return A, nil
	}
	for _, a := range strings.Split(s, "|") {
		bounds := strings.Split(a, "-")
		if len(bounds) > 2 {
			return nil, errors.New("error parsing ambiguity range from file")
		}
		start, err1 := strconv.Atoi(bounds[0])
		stop, err2 := strconv.Atoi(bounds[len(bounds)-1])
		if err1 != nil || err2 != nil {
			return nil, errors.New("error parsing ambiguity range from file")
		}
		A = append(A, start, stop)
	}
	return A, nil
}

// csvRow is one data line of a csv file in the format written by gofasta updown list, and its
// 0-based position among the data lines
type csvRow struct {
	fields []string
	idx    int
}

// RecordName returns the row's query name, its first field (see pipeline.Namer)
func (r csvRow) RecordName() string {
	if len(r.fields) == 0 {
		return ""
	}
	return r.fields[0]
}

// newCSVSource returns a function that reads one data line per call from a csv file in the format
// written by gofasta updown list, and returns io.EOF after the last one. It checks the header, returning
// badHeader if it is not updown list's, but not the data lines (see parseCSVRow).
func newCSVSource(in io.Reader, badHeader error) func() (csvRow, error) {
	r := csv.NewReader(in)
	header := true
	counter := 0
	return func() (csvRow, error) {
		record, err := r.Read()
		if err != nil {
			return csvRow{}, err
		}
		if header {
			if !slices.Equal(record, strings.Split(listHeader, ",")) {
				return csvRow{}, badHeader
			}
			header = false
			record, err = r.Read()
			if err != nil {
				return csvRow{}, err
			}
		}
		row := csvRow{fields: record, idx: counter}
		counter++
		return row, nil
	}
}

// parseCSVRow parses one data line of updown list's output
func parseCSVRow(row csvRow) (updownLine, error) {
	record := row.fields
	fail := func(err error) (updownLine, error) { return updownLine{}, pipeline.Named(record[0], err) }

	a, err := getAmbArr(record[2])
	if err != nil {
		return fail(err)
	}

	var snps []string
	var snpPos []int
	if len(record[1]) > 0 {
		snps = strings.Split(record[1], "|")
		snpPos = make([]int, len(snps))
		for i, snp := range snps {
			if len(snp) < 3 {
				return fail(errors.New("error parsing SNP \"" + snp + "\" from file"))
			}
			snpPos[i], err = strconv.Atoi(snp[1 : len(snp)-1])
			if err != nil {
				return fail(err)
			}
		}
	} else {
		snps = make([]string, 0)
		snpPos = make([]int, 0)
	}

	amb_count, err := strconv.Atoi(record[4])
	if err != nil {
		return fail(err)
	}

	snpsSorted := slices.Clone(snps)
	slices.Sort(snpsSorted)

	return updownLine{id: record[0], idx: row.idx, snps: snps, snpsSorted: snpsSorted, snpsPos: snpPos, ambs: a, ambCount: amb_count, ordered: inOrder(snpPos, a)}, nil
}

// collectQueries runs the query records through a pipeline, reporting the ones that fail on stderr and
// leaving them out (unless strict, when the first is an error), and returns the others in input order
func collectQueries[In any](next func() (In, error), work func(In) (updownLine, error), threads int, strict bool) ([]updownLine, error) {
	udla := make([]updownLine, 0)
	h := pipeline.NewHandler(strict, os.Stderr, pipeline.Noun{One: "query", Many: "queries"}, func(udL updownLine) error {
		udla = append(udla, udL)
		return nil
	})
	if err := pipeline.Run(context.Background(), threads, next, work, h.Sink); err != nil {
		return nil, err
	}
	if err := h.Finish(); err != nil {
		return nil, err
	}
	return udla, nil
}

// readCSVToUDLList reads a csv file in the format produced by gofasta updown
// list to an array of updownLine structs, each of which contains the snps and
// ambiguous positions for one query or target sequence
func readCSVToUDLList(in io.Reader, threads int, strict bool) ([]updownLine, error) {
	next := newCSVSource(in, errors.New("bad header when parsing --query csv: is this file the output of gofasta updown list?"))
	return collectQueries(next, parseCSVRow, threads, strict)
}

// fastaToUDLList converts a fasta format alignment to an array of updownLine
// structs, given a reference sequence, each of which contains the snps and
// ambiguous positions for one query or target sequence
func fastaToUDLList(in io.Reader, refSeq []byte, threads int, strict bool) ([]updownLine, error) {
	DA := encoding.MakeDecodingArray()
	work := func(record fasta.Record) (updownLine, error) { return lineFromRecord(refSeq, record, DA) }
	return collectQueries(fasta.NewRecordSource(in), work, threads, strict)
}

// lineFromRecord encodes a fasta record and gets its mutation and ambiguity lists
func lineFromRecord(refSeq []byte, record fasta.Record, DA [256]string) (updownLine, error) {
	FR, err := record.Encode()
	if err != nil {
		return updownLine{}, pipeline.Named(record.ID, err)
	}
	udL, err := getLine(refSeq, FR, DA)
	return udL, pipeline.Named(record.ID, err)
}

// getLine gets the mutation and ambiguity lists between the reference and one fasta record
func getLine(refSeq []byte, FR fasta.EncodedRecord, DA [256]string) (updownLine, error) {

	if len(FR.Seq) != len(refSeq) {
		return updownLine{}, errors.New("alignment and reference are not the same width")
	}

	var (
		snps      = make([]string, 0)
		snpPos    = make([]int, 0)
		ambs      = make([]int, 0) // [1,265,11083,11083,...]: each pair is the 1-based inclusive start/end of a tract of ambiguities
		ambCount  int
		amb_start int
		amb_stop  int
		cont      bool // in a tract of ambiguities
	)

	const known = 0x0808080808080808 // bit 3 (a known base) in every byte
	for i := 0; i < len(FR.Seq); i++ {
		// skip eight sites at a time where the query is the reference, every base known, outside a
		// tract of ambiguities
		for !cont && i+8 <= len(FR.Seq) {
			q := binary.LittleEndian.Uint64(FR.Seq[i:])
			if q != binary.LittleEndian.Uint64(refSeq[i:]) || q&known != known {
				break
			}
			i += 8
		}
		if i == len(FR.Seq) {
			break
		}
		que_nuc := FR.Seq[i]

		// if query nucleotide is a known base
		if que_nuc&8 == 8 {
			// if it is different from the reference:
			if (refSeq[i] & que_nuc) < 16 {
				snps = append(snps, DA[refSeq[i]]+strconv.Itoa(i+1)+DA[que_nuc])
				snpPos = append(snpPos, i+1)
			}
			// also need to finalise any previous ambiguity tract
			if cont {
				ambs = append(ambs, amb_start+1, amb_stop+1)
				cont = false
			}
			// otherwise update the ambiguities
		} else {
			ambCount++
			if cont {
				amb_stop = i
			} else {
				amb_start = i
				amb_stop = i
				cont = true
			}
		}
	}

	// need to finalise any ambiguity tract that reaches the edge of the sequence
	if cont {
		ambs = append(ambs, amb_start+1, amb_stop+1)
	}

	snpsSorted := slices.Clone(snps)
	slices.Sort(snpsSorted)

	return updownLine{
		id:         FR.ID,
		idx:        FR.Idx,
		snps:       snps,
		snpsPos:    snpPos,
		ambs:       ambs,
		ambCount:   ambCount,
		snpsSorted: snpsSorted,
		ordered:    inOrder(snpPos, ambs),
	}, nil
}

// listHeader is the header line of the CSV that gofasta updown list writes (see List)
const listHeader = "query,SNPs,ambiguities,SNPcount,ambcount"

// InputType tells whether r holds a fasta file or the CSV output of gofasta updown list, from its start: the
// first non-blank character of a fasta file is '>', and the CSV starts with its header. It returns "fasta"
// or "csv", and a reader that gives all of r (including what was looked at).
func InputType(r io.Reader) (io.Reader, string, error) {
	br := bufio.NewReaderSize(r, 64*1024)
	b, err := br.Peek(len(listHeader))
	if err != nil && err != io.EOF && err != bufio.ErrBufferFull {
		return br, "", err
	}
	if len(b) == len(listHeader) && string(b) == listHeader {
		return br, "csv", nil
	}
	// skip leading blank lines (up to the buffer's size) to find the first character of a fasta file
	b, _ = br.Peek(br.Size())
	switch text := bytes.TrimLeft(b, " \t\r\n"); {
	case len(text) > 0 && text[0] == '>':
		return br, "fasta", nil
	case len(text) == 0:
		return br, "", errors.New("it is empty")
	}
	return br, "", errors.New("it is neither a fasta file nor the CSV output of gofasta updown list (which starts with the header " + listHeader + ")")
}
