package fasta

import (
	"bufio"
	"bytes"
	"errors"
	"io"
)

// Errors from reading fasta files.
var (
	ErrBadlyFormedFasta = errors.New("badly formed fasta file")
	ErrEmptyFasta       = errors.New("empty fasta file")
	ErrDiffLenSeqs      = errors.New("different length sequences in input file: is this an alignment?")
)

// Reader reads fasta records one at a time
type Reader struct {
	br *bufio.Reader
}

// NewReader returns a Reader that reads from f
func NewReader(f io.Reader) *Reader {
	return &Reader{br: bufio.NewReader(f)}
}

// trimNewline removes a trailing "\n" or "\r\n" from a line
func trimNewline(line []byte) []byte {
	if l, ok := bytes.CutSuffix(line, []byte("\n")); ok {
		return bytes.TrimSuffix(l, []byte("\r"))
	}
	return line
}

// Read reads one fasta record from the underlying reader. The final record is
// returned with error = nil, and the next call to Read() returns an empty Record
// struct and error = io.EOF.
func (r *Reader) Read() (Record, error) {
	// the header line. One that ends the file without a newline is kept, as a record with an empty
	// sequence, the same as one that does end in a newline.
	line, err := r.br.ReadBytes('\n')
	if err != nil && !(err == io.EOF && len(line) > 0) {
		return Record{}, err
	}
	if line[0] != '>' {
		return Record{}, ErrBadlyFormedFasta
	}
	line = trimNewline(line)
	var FR Record
	if fields := bytes.Fields(line[1:]); len(fields) > 0 {
		FR.ID = string(fields[0]) // (empty if the header has nothing after the ">")
	}
	FR.Description = string(line[1:])

	// the sequence lines, up to the next header or the end of the file
	var buffer []byte
	for {
		peek, err := r.br.Peek(1)
		if err == io.EOF || (err == nil && peek[0] == '>') {
			break
		}
		if err != nil {
			return Record{}, err
		}
		// a final line without a newline gives io.EOF, which the next Peek sees
		line, err := r.br.ReadBytes('\n')
		if err != nil && err != io.EOF {
			return Record{}, err
		}
		buffer = append(buffer, trimNewline(line)...)
	}
	FR.Seq = string(buffer)
	return FR, nil
}

// NewRecordSource returns a function that reads the next record in fasta format each time it is called,
// setting its Idx, and returns io.EOF after the last one. It only splits the input into records (an error
// means the input could not be read as fasta); it does not check the records' sequences or lengths. It
// returns an error if the input has no records.
func NewRecordSource(f io.Reader) func() (Record, error) {
	r := NewReader(f)
	counter := 0
	return func() (Record, error) {
		record, err := r.Read()
		if err == io.EOF {
			if counter == 0 {
				return Record{}, ErrEmptyFasta
			}
			return Record{}, io.EOF
		}
		if err != nil {
			return Record{}, err
		}
		record.Idx = counter
		counter++
		return record, nil
	}
}

// EncodeRecord encodes a record (see Record.Encode; Record.EncodeHardGaps if hardGaps), optionally
// calculating its base content (atgc) and its completeness (score).
func EncodeRecord(record Record, hardGaps, atgc, score bool) (EncodedRecord, error) {
	encodedRecord, err := record.encode(hardGaps)
	if err != nil {
		return EncodedRecord{}, err
	}
	if score {
		encodedRecord.CalculateCompleteness()
	}
	if atgc {
		encodedRecord.CalculateBaseContent()
	}
	return encodedRecord, nil
}

// NewAlignmentSource returns a function that reads the next record of an alignment in fasta format each
// time it is called, setting its Idx. After the last record it returns io.EOF. It returns an error if the
// input has no records, or if a record's length differs from the first's.
func NewAlignmentSource(f io.Reader) func() (Record, error) {
	next := NewRecordSource(f)
	var width int
	return func() (Record, error) {
		record, err := next()
		if err != nil {
			return Record{}, err
		}
		if record.Idx == 0 {
			width = len(record.Seq)
		} else if width != len(record.Seq) {
			return Record{}, ErrDiffLenSeqs
		}
		return record, nil
	}
}

// NewEncodedAlignmentSource is as NewAlignmentSource, but encodes each record (see Record.Encode;
// Record.EncodeHardGaps if hardGaps), optionally calculating its base content (atgc) and its
// completeness (score).
func NewEncodedAlignmentSource(f io.Reader, hardGaps, atgc, score bool) func() (EncodedRecord, error) {
	next := NewAlignmentSource(f)
	return func() (EncodedRecord, error) {
		record, err := next()
		if err != nil {
			return EncodedRecord{}, err
		}
		return EncodeRecord(record, hardGaps, atgc, score)
	}
}

// LoadEncodeAlignment reads and encodes every record of an alignment in fasta format (see
// NewEncodedAlignmentSource) and returns them in input order
func LoadEncodeAlignment(f io.Reader, hardGaps, atgc, score bool) ([]EncodedRecord, error) {
	next := NewEncodedAlignmentSource(f, hardGaps, atgc, score)
	encodedRecords := make([]EncodedRecord, 0)
	for {
		encodedRecord, err := next()
		if err == io.EOF {
			return encodedRecords, nil
		}
		if err != nil {
			return make([]EncodedRecord, 0), err
		}
		encodedRecords = append(encodedRecords, encodedRecord)
	}
}

// WriteRecord writes one record in fasta format: its ID on the header line, then its sequence, in lines
// of wrap characters if wrap > 0, otherwise on one line. (With wrap > 0, an empty sequence gives no
// sequence lines.)
func WriteRecord(w io.Writer, record Record, wrap int) error {
	if _, err := io.WriteString(w, ">"+record.ID+"\n"); err != nil {
		return err
	}
	if wrap <= 0 {
		_, err := io.WriteString(w, record.Seq+"\n")
		return err
	}
	for written := 0; written < len(record.Seq); written += wrap {
		end := min(written+wrap, len(record.Seq))
		if _, err := io.WriteString(w, record.Seq[written:end]+"\n"); err != nil {
			return err
		}
	}
	return nil
}

// getAlignmentDims gets the dimensions (i.e. the number of sequences and the
// width of the alignment in number of nucleotides) of an alignment in fasta
// format
func getAlignmentDims(f io.Reader) (int, int, error) {
	n, l := 0, 0
	r := NewReader(f)
	for {
		record, err := r.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return 0, 0, err
		}
		if n == 0 {
			l = len(record.Seq)
		}
		n++
	}
	return n, l, nil
}
