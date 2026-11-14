// Package samfile reads files in the SAM text format, following the SAM specification (SAMv1, the SAM/BAM
// Format Specification Working Group): the header, and the mandatory fields of alignment records. Every
// mandatory field and every header line is checked against the specification; the optional fields of
// records, and the contents of @RG and @PG lines and of @SQ fields other than SN, LN and AN, are not. The
// one exception is SEQ, where a character the specification does not allow is read as N (see
// Record.OtherSeq).
package samfile

import (
	"bufio"
	"bytes"
	"errors"
	"fmt"
	"io"
	"strconv"
)

// ErrEmpty is the error NewReader returns for a file with nothing in it
var ErrEmpty = errors.New("the sam file is empty")

// Reader reads a SAM file: its header when it is made (see NewReader), then its alignment lines one at a
// time (see Next)
type Reader struct {
	Header *Header
	br     *bufio.Reader
}

// NewReader reads and checks the header of a SAM file (the lines at its start that begin with "@")
func NewReader(r io.Reader) (*Reader, error) {
	rd := &Reader{br: bufio.NewReaderSize(r, 1<<20)}
	if _, err := rd.br.Peek(1); err == io.EOF {
		return nil, ErrEmpty
	} else if err != nil {
		return nil, err
	}
	var text []byte
	for {
		p, err := rd.br.Peek(1)
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, err
		}
		if p[0] != '@' {
			break
		}
		line, err := rd.readLine()
		if err != nil {
			return nil, err
		}
		text = append(append(text, line...), '\n')
	}
	h, err := ParseHeader(text)
	if err != nil {
		return nil, err
	}
	rd.Header = h
	return rd, nil
}

// readLine returns the next line without its line ending (a final line need not have one), or io.EOF
func (r *Reader) readLine() ([]byte, error) {
	line, err := r.br.ReadBytes('\n')
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

// LineError is the error Next returns for an alignment line that does not follow the specification. Name
// is the line's first field (QNAME) as it is, so that the caller can tell which query the line belongs to.
type LineError struct {
	Name string
	Err  error
}

func (e *LineError) Error() string { return e.Err.Error() }
func (e *LineError) Unwrap() error { return e.Err }

// Next returns the next alignment line, skipping blank lines, or io.EOF after the last one. A line that
// does not follow the specification gives a *LineError, and the lines after it can still be read; any other
// error means the file could not be read.
func (r *Reader) Next() (Record, error) {
	for {
		line, err := r.readLine()
		if err != nil {
			return Record{}, err
		}
		if len(line) == 0 {
			continue
		}
		rec, err := r.Header.ParseRecord(line)
		if err != nil {
			name, _, _ := bytes.Cut(line, []byte{'\t'})
			return Record{}, &LineError{Name: string(name), Err: err}
		}
		return rec, nil
	}
}

// isLetter reports whether b is an ASCII letter
func isLetter(b byte) bool { return 'A' <= b && b <= 'Z' || 'a' <= b && b <= 'z' }

// CigarType is the type of a CIGAR operation
type CigarType uint8

const (
	CigarMatch CigarType = iota
	CigarInsertion
	CigarDeletion
	CigarSkipped
	CigarSoftClipped
	CigarHardClipped
	CigarPadded
	CigarEqual
	CigarMismatch
	cigarUnknown
)

var cigarLetters = "MIDNSHP=X?"

// String returns the CIGAR letter of t
func (t CigarType) String() string {
	if t > cigarUnknown {
		t = cigarUnknown
	}
	return cigarLetters[t : t+1]
}

// Consumes returns how many query and reference bases one base of an operation of type t uses
func (t CigarType) Consumes() (query, reference int) {
	switch t {
	case CigarMatch, CigarEqual, CigarMismatch:
		return 1, 1
	case CigarInsertion, CigarSoftClipped:
		return 1, 0
	case CigarDeletion, CigarSkipped:
		return 0, 1
	}
	return 0, 0
}

// CigarOp is one CIGAR operation: its length and type
type CigarOp uint64

func NewCigarOp(t CigarType, n int) CigarOp { return CigarOp(t) | CigarOp(n)<<4 }

func (o CigarOp) Type() CigarType { return CigarType(o & 0xf) }
func (o CigarOp) Len() int        { return int(o >> 4) }

// CigarLengths returns the number of reference bases a CIGAR spans and the number of query bases it uses
func CigarLengths(c []CigarOp) (ref, query int) {
	for _, o := range c {
		q, r := o.Type().Consumes()
		ref += o.Len() * r
		query += o.Len() * q
	}
	return ref, query
}

var cigarTypeOf = func() (t [256]CigarType) {
	for i := range t {
		t[i] = cigarUnknown
	}
	for i := range cigarLetters[:len(cigarLetters)-1] {
		t[cigarLetters[i]] = CigarType(i)
	}
	return t
}()

// maxInt31 is 2^31-1, the largest position, length or template length SAM allows
const maxInt31 = 1<<31 - 1

// parseCigar parses a CIGAR string: "*", or operations of a length (digits) and a type (one of MIDNSHP=X).
// Hard clips may only be first or last, and soft clips may only have hard clips between them and the ends.
func parseCigar(b []byte) ([]CigarOp, error) {
	if len(b) == 1 && b[0] == '*' {
		return nil, nil
	}
	bad := func(why string) ([]CigarOp, error) {
		return nil, fmt.Errorf("sam: failed to parse cigar string %q: %s", b, why)
	}
	if len(b) == 0 {
		return bad("empty")
	}
	nops := 0
	for _, c := range b {
		if c < '0' || c > '9' {
			nops++
		}
	}
	ops := make([]CigarOp, 0, nops)
	n, haveDigits := 0, false
	for _, c := range b {
		if '0' <= c && c <= '9' {
			if n = n*10 + int(c-'0'); n > maxInt31 {
				return bad("operation too long")
			}
			haveDigits = true
			continue
		}
		t := cigarTypeOf[c]
		switch {
		case t == cigarUnknown:
			return bad(fmt.Sprintf("unknown operation %q", c))
		case !haveDigits:
			return bad(fmt.Sprintf("operation %q has no length", c))
		}
		ops = append(ops, NewCigarOp(t, n))
		n, haveDigits = 0, false
	}
	if haveDigits {
		return bad("a length with no operation at the end")
	}
	for i, o := range ops {
		switch t := o.Type(); {
		case t == CigarHardClipped && i != 0 && i != len(ops)-1:
			return bad("hard clip (H) not at either end")
		case t == CigarSoftClipped && !clipToEnd(ops, i):
			return bad("soft clip (S) with something other than a hard clip between it and the end")
		}
	}
	return ops, nil
}

// clipToEnd reports whether operation i has only hard clips between it and one end of the CIGAR
func clipToEnd(ops []CigarOp, i int) bool {
	before, after := true, true
	for _, o := range ops[:i] {
		before = before && o.Type() == CigarHardClipped
	}
	for _, o := range ops[i+1:] {
		after = after && o.Type() == CigarHardClipped
	}
	return before || after
}

// Flags of a SAM record
const (
	FlagUnmapped      = 0x4
	FlagSecondary     = 0x100
	FlagSupplementary = 0x800
)

// Record holds the mandatory fields of a SAM alignment record that gofasta uses
type Record struct {
	Name  string
	Flags uint16
	Ref   int // the index of its @SQ line (RNAME), or -1 for "*"
	Pos   int // 0-based (POS - 1)
	Cigar []CigarOp
	Seq   []byte // upper case, "=" for "the reference base" (see ResolveEquals); nil if SEQ is "*"
	// SEQ had characters the SAM specification does not allow there (it allows letters, "=" and "."), such as
	// the "-" a FASTA alignment uses for a gap, which aligners copy from their input: read as N
	OtherSeq bool
}

// Primary reports whether r is its query's primary line (neither secondary nor supplementary)
func (r Record) Primary() bool { return r.Flags&(FlagSecondary|FlagSupplementary) == 0 }

// seqTable maps the characters SAM allows in SEQ (letters, "=" and ".") to what gofasta uses: IUPAC codes
// upper case, "=" as it is, and N for any other letter and for "."; 0 marks a character the specification
// does not allow in SEQ (ParseRecord reads it as N, and sets OtherSeq)
var seqTable = func() (t [256]byte) {
	for c := 'A'; c <= 'Z'; c++ {
		t[c], t[c+'a'-'A'] = 'N', 'N'
	}
	for _, c := range []byte("ACMGRSVTWYHKDBN") {
		t[c], t[c+'a'-'A'] = c, c
	}
	t['='], t['.'] = '=', 'N'
	return t
}()

// uintField parses a field that must be a decimal integer from 0 to max (leading zeros are allowed)
func uintField(b []byte, max int) (int, bool) {
	if len(b) == 0 {
		return 0, false
	}
	n := 0
	for _, c := range b {
		if c < '0' || c > '9' {
			return 0, false
		}
		if n = n*10 + int(c-'0'); n > max {
			return 0, false
		}
	}
	return n, true
}

// isRefName reports whether b is a valid reference name (SAM's [:rname:^*=][:rname:]*): printable ASCII
// other than \ , " ' ` ( ) [ ] { } < >, not starting with * or =
func isRefName(b []byte) bool {
	if len(b) == 0 || b[0] == '*' || b[0] == '=' {
		return false
	}
	for _, c := range b {
		if c < '!' || c > '~' || bytes.IndexByte([]byte("\\,\"'`()[]{}<>"), c) >= 0 {
			return false
		}
	}
	return true
}

// ParseRecord parses one alignment line of a SAM file, checking its mandatory fields against the SAM
// specification (the optional fields after them are not checked)
func (h *Header) ParseRecord(line []byte) (Record, error) {
	var f [11][]byte
	n := 0
	for ; n < 11; n++ {
		i := bytes.IndexByte(line, '\t')
		if i < 0 {
			f[n] = line
			n++
			break
		}
		f[n], line = line[:i], line[i+1:]
	}
	if n < 11 {
		return Record{}, errors.New("sam: missing SAM fields")
	}
	bad := func(field string, b []byte) (Record, error) {
		return Record{}, fmt.Errorf("sam: invalid %s: %q", field, b)
	}

	if len(f[0]) == 0 || len(f[0]) > 254 {
		return bad("QNAME (1 to 254 characters)", f[0])
	}
	for _, c := range f[0] {
		if c < '!' || c > '~' || c == '@' {
			return bad("QNAME (printable characters other than @)", f[0])
		}
	}
	r := Record{Name: string(f[0])}

	flags, ok := uintField(f[1], 0xffff)
	if !ok {
		return bad("FLAG (an integer from 0 to 65535)", f[1])
	}
	r.Flags = uint16(flags)

	var err error
	if r.Ref, err = h.refIndex(f[2]); err != nil {
		return Record{}, fmt.Errorf("sam: failed to assign reference: %w", err)
	}
	pos, ok := uintField(f[3], maxInt31)
	if !ok {
		return bad("POS (an integer from 0 to 2^31-1)", f[3])
	}
	r.Pos = pos - 1
	if _, ok := uintField(f[4], 255); !ok {
		return bad("MAPQ (an integer from 0 to 255)", f[4])
	}
	if r.Cigar, err = parseCigar(f[5]); err != nil {
		return Record{}, err
	}
	if !(len(f[6]) == 1 && (f[6][0] == '=' || f[6][0] == '*')) {
		if _, err := h.refIndex(f[6]); err != nil {
			return Record{}, fmt.Errorf("sam: failed to assign mate reference: %w", err)
		}
	}
	if _, ok := uintField(f[7], maxInt31); !ok {
		return bad("PNEXT (an integer from 0 to 2^31-1)", f[7])
	}
	if _, ok := uintField(bytes.TrimPrefix(f[8], []byte{'-'}), maxInt31); !ok {
		return bad("TLEN (an integer from -(2^31-1) to 2^31-1)", f[8])
	}
	if !(len(f[9]) == 1 && f[9][0] == '*') {
		if len(f[9]) == 0 {
			return bad("SEQ (empty)", f[9])
		}
		r.Seq = make([]byte, len(f[9]))
		for i, c := range f[9] {
			if r.Seq[i] = seqTable[c]; r.Seq[i] == 0 {
				r.Seq[i], r.OtherSeq = 'N', true
			}
		}
		if _, q := CigarLengths(r.Cigar); len(r.Cigar) != 0 && q != len(r.Seq) {
			return Record{}, errors.New("sam: sequence/CIGAR length mismatch")
		}
	}
	if !(len(f[10]) == 1 && f[10][0] == '*') {
		if len(f[10]) == 0 {
			return bad("QUAL (empty)", f[10])
		}
		for i, c := range f[10] {
			if c < '!' || c > '~' {
				return Record{}, fmt.Errorf("sam: invalid QUAL (printable characters): %q at query position %d", c, i+1)
			}
		}
		if len(f[10]) != len(r.Seq) {
			return Record{}, errors.New("sam: sequence/quality length mismatch")
		}
	}
	return r, nil
}

// ResolveEquals replaces each "=" in r's SEQ (which means "the reference base") by its reference base from
// ref, or by N if ref is nil (no reference is known), and reports whether there were any. A "=" in an
// insertion or a soft clip, which has no reference base, is an error, after which SEQ may have been
// partly rewritten. SEQ is rewritten in place: ParseRecord gives each record its own.
func (r *Record) ResolveEquals(ref []byte) (bool, error) {
	if bytes.IndexByte(r.Seq, '=') < 0 {
		return false, nil
	}
	seq := r.Seq
	q, p := 0, r.Pos
	for _, o := range r.Cigar {
		qn, rn := o.Type().Consumes()
		for k := 0; k < o.Len()*qn; k++ {
			if seq[q+k] != '=' {
				continue
			}
			switch {
			case rn == 0:
				return false, errors.New("bad sam record: '=' in SEQ in an insertion or soft clip, where there is no reference base")
			case ref == nil:
				seq[q+k] = 'N'
			default:
				seq[q+k] = ref[p+k]
			}
		}
		q += o.Len() * qn
		p += o.Len() * rn
	}
	return true, nil
}

// Reference is one @SQ line of a SAM header
type Reference struct {
	Name   string
	Length int
}

// Header holds what gofasta uses of a SAM header: its references, in order. Make one with ParseHeader.
type Header struct {
	Refs  []Reference
	index map[string]int // the index in Refs of each reference, by its name (SN)
}

// refIndex returns the index of the reference named name, or -1 for "*"
func (h *Header) refIndex(name []byte) (int, error) {
	if len(name) == 1 && name[0] == '*' {
		return -1, nil
	}
	if i, ok := h.index[string(name)]; ok {
		return i, nil
	}
	return -1, fmt.Errorf("no reference with name %q", name)
}

// Errors for a malformed SAM header (the line they are about is added)
var (
	errBadHeader    = errors.New("sam: malformed header line")
	errDupTag       = errors.New("sam: duplicate field")
	errDupReference = errors.New("sam: duplicate reference name")
	errDupReadGroup = errors.New("sam: duplicate read group name")
	errDupProgram   = errors.New("sam: duplicate program name")
	errBadLen       = errors.New("sam: reference length out of range")
	errHDNotFirst   = errors.New("sam: @HD must be the first line, and only one")
	errBadVersion   = errors.New("sam: @HD VN is not of the form 1.6")
	errBadRefName   = errors.New("sam: invalid reference name")
)

// ParseHeader parses the header lines of a SAM file (each starting with "@"), checking them against the
// SAM specification: each line is @HD, @SQ, @RG or @PG followed by TAG:VALUE fields (a letter, then a letter
// or digit, ':', and a value of printable characters), with no tag twice, or @CO and a tab; @HD, if there is
// one, is the first line, and has a VN of the form 1.6; @SQ has a valid SN and an LN from 1 to 2^31-1, and
// no reference name (SN, or one of the AN names) is used twice; @RG and @PG have an ID not used before.
func ParseHeader(text []byte) (*Header, error) {
	h := &Header{index: map[string]int{}}
	names := map[string]bool{} // the reference names (SN and AN) so far
	groups, progs := map[string]bool{}, map[string]bool{}
	first := true
	for i, l := range bytes.Split(text, []byte{'\n'}) {
		l = bytes.TrimSuffix(l, []byte{'\r'})
		if len(l) == 0 {
			continue
		}
		wrap := func(err error) error { return fmt.Errorf("%w: line %d: %q", err, i+1, l) }
		if len(l) < 4 || l[0] != '@' || l[3] != '\t' {
			return nil, wrap(errBadHeader)
		}
		kind := string(l[1:3])
		if kind == "CO" {
			first = false
			continue
		}
		tags, err := headerTags(l)
		if err != nil {
			return nil, wrap(err)
		}
		switch kind {
		case "HD":
			if !first {
				return nil, wrap(errHDNotFirst)
			}
			if !isVersion(tags["VN"]) {
				return nil, wrap(errBadVersion)
			}
		case "SQ":
			err = h.referenceLine(tags, names)
		case "RG":
			err = idLine(tags, groups, errDupReadGroup)
		case "PG":
			err = idLine(tags, progs, errDupProgram)
		default:
			err = errBadHeader
		}
		if err != nil {
			return nil, wrap(err)
		}
		first = false
	}
	return h, nil
}

// headerTags returns the TAG:VALUE fields of a header line (after the record type), checking their form:
// a tag of a letter then a letter or digit, ':', and a value of printable characters (spaces allowed), at
// least one field, and no tag twice
func headerTags(l []byte) (map[string][]byte, error) {
	fields := bytes.Split(l, []byte{'\t'})[1:]
	tags := make(map[string][]byte, len(fields))
	for _, f := range fields {
		if len(f) < 4 || f[2] != ':' || !isLetter(f[0]) || !(isLetter(f[1]) || '0' <= f[1] && f[1] <= '9') {
			return nil, errBadHeader
		}
		for _, c := range f[3:] {
			if c < ' ' || c > '~' {
				return nil, errBadHeader
			}
		}
		t := string(f[:2])
		if _, dup := tags[t]; dup {
			return nil, errDupTag
		}
		tags[t] = f[3:]
	}
	return tags, nil
}

// isVersion reports whether v is a SAM version number such as 1.6 (digits, '.', digits)
func isVersion(v []byte) bool {
	major, minor, ok := bytes.Cut(v, []byte{'.'})
	if !ok || len(major) == 0 || len(minor) == 0 {
		return false
	}
	for _, c := range append(append([]byte{}, major...), minor...) {
		if c < '0' || c > '9' {
			return false
		}
	}
	return true
}

// referenceLine reads an @SQ line: SN and LN are needed, the name and any alternative names (AN, comma
// separated) must be valid reference names not used before, and LN must be from 1 to 2^31-1
func (h *Header) referenceLine(tags map[string][]byte, names map[string]bool) error {
	name, hasName := tags["SN"]
	ln, hasLen := tags["LN"]
	if !hasName || !hasLen {
		return errBadHeader
	}
	length, ok := uintField(ln, maxInt31)
	if !ok {
		if _, err := strconv.Atoi(string(ln)); err != nil {
			return errBadHeader
		}
		return errBadLen
	}
	if length < 1 {
		return errBadLen
	}
	all := [][]byte{name}
	if an, ok := tags["AN"]; ok {
		all = append(all, bytes.Split(an, []byte{','})...)
	}
	for _, n := range all {
		if !isRefName(n) {
			return errBadRefName
		}
		if names[string(n)] {
			return errDupReference
		}
		names[string(n)] = true
	}
	h.index[string(name)] = len(h.Refs)
	h.Refs = append(h.Refs, Reference{Name: string(name), Length: length})
	return nil
}

// idLine checks an @RG or @PG line: it needs an ID not used before
func idLine(tags map[string][]byte, seen map[string]bool, errDup error) error {
	id, ok := tags["ID"]
	if !ok {
		return errBadHeader
	}
	if seen[string(id)] {
		return errDup
	}
	seen[string(id)] = true
	return nil
}
