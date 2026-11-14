package sam

import (
	"bytes"
	"errors"
	"fmt"
	"net/url"
	"strconv"
	"strings"
)

// This file parses the SAM text format: the header's @SQ lines and the fields of alignment records that
// gofasta uses. It checks the structure of what it reads (the eleven mandatory fields, numeric fields,
// CIGAR strings, sequence and quality lengths, references that exist, and well-formed header lines) but not
// the contents of fields gofasta does not use (optional tags, @RG and @PG details, @SQ M5 and UR values).

// cigarType is the type of a CIGAR operation
type cigarType uint8

const (
	cigarMatch cigarType = iota
	cigarInsertion
	cigarDeletion
	cigarSkipped
	cigarSoftClipped
	cigarHardClipped
	cigarPadded
	cigarEqual
	cigarMismatch
	cigarBack
	cigarUnknown
)

var cigarLetters = "MIDNSHP=XB?"

// String returns the CIGAR letter of t
func (t cigarType) String() string {
	if t > cigarUnknown {
		t = cigarUnknown
	}
	return cigarLetters[t : t+1]
}

// consumes returns how many query and reference bases one base of an operation of type t uses (a B
// operation steps back along the reference)
func (t cigarType) consumes() (query, reference int) {
	switch t {
	case cigarMatch, cigarEqual, cigarMismatch:
		return 1, 1
	case cigarInsertion, cigarSoftClipped:
		return 1, 0
	case cigarDeletion, cigarSkipped:
		return 0, 1
	case cigarBack:
		return 0, -1
	}
	return 0, 0
}

// cigarOp is one CIGAR operation: its length (at most 2^28-1) and type
type cigarOp uint32

func newCigarOp(t cigarType, n int) cigarOp { return cigarOp(t) | cigarOp(n)<<4 }

func (o cigarOp) Type() cigarType { return cigarType(o & 0xf) }
func (o cigarOp) Len() int        { return int(o >> 4) }

// maxCigarOpLen is the longest one operation can be; longer ones are split
const maxCigarOpLen = 1<<28 - 1

// cigarLengths returns the number of reference bases a CIGAR spans and the number of query bases it uses
func cigarLengths(c []cigarOp) (ref, query int) {
	for _, o := range c {
		q, r := o.Type().consumes()
		if o.Type() != cigarBack {
			ref += o.Len() * r
		}
		query += o.Len() * q
	}
	return ref, query
}

// cigarValid reports whether a CIGAR can describe a query sequence of the given length: hard clips only at
// the ends, soft clips only at the ends or next to a hard clip, and the query bases it uses add up
func cigarValid(c []cigarOp, length int) bool {
	pos := 0
	for i, o := range c {
		t := o.Type()
		if t == cigarHardClipped && i != 0 && i != len(c)-1 {
			return false
		}
		if t == cigarSoftClipped && i != 0 && i != len(c)-1 {
			if c[i-1].Type() != cigarHardClipped && c[i+1].Type() != cigarHardClipped {
				return false
			}
		}
		q, r := t.consumes()
		if pos < 0 && q != 0 {
			return false
		}
		length -= o.Len() * q
		pos += o.Len() * r
	}
	return length == 0
}

var cigarTypeOf = func() (t [256]cigarType) {
	for i := range t {
		t[i] = cigarUnknown
	}
	for i := range cigarLetters[:len(cigarLetters)-1] {
		t[cigarLetters[i]] = cigarType(i)
	}
	return t
}()

// parseCigar parses a CIGAR string ("*" for none). An operation with no number has length 0, a number at the
// end repeats the last operation (an error after one that was split, as every operation is), and an
// operation longer than maxCigarOpLen is split into several.
func parseCigar(b []byte) ([]cigarOp, error) {
	if len(b) == 1 && b[0] == '*' {
		return nil, nil
	}
	nops := 0
	for _, c := range b {
		if c < '0' || c > '9' {
			nops++
		}
	}
	ops := make([]cigarOp, 0, nops)
	t, n := cigarMatch, 0
	for i := 0; i < len(b); i++ {
		for j := i; j < len(b); j++ {
			if b[j] < '0' || '9' < b[j] {
				var err error
				if n, err = cigarNumber(b[i:j]); err != nil {
					return nil, err
				}
				t = cigarTypeOf[b[j]]
				i = j
				break
			}
		}
		if t == cigarUnknown {
			return nil, fmt.Errorf("sam: failed to parse cigar string %q: unknown operation %q", b, cigarUnknown.String())
		}
		for {
			if n < 0 { // (a number at the end, after a long operation was split)
				return nil, errors.New("error parsing a sam record: sam: illegal CIGAR op length")
			}
			ops = append(ops, newCigarOp(t, min(n, maxCigarOpLen)))
			if n -= maxCigarOpLen; n <= 0 {
				break
			}
		}
	}
	return ops, nil
}

// cigarNumber parses the length of a CIGAR operation (digits only; none is 0)
func cigarNumber(b []byte) (int, error) {
	if len(b) > 13 {
		return 0, fmt.Errorf("sam: integer overflow: %q", b)
	}
	n := 0
	for _, c := range b {
		n = n*10 + int(c-'0')
	}
	return n, nil
}

// Flags of a SAM record that gofasta uses
const (
	flagUnmapped      = 0x4
	flagSecondary     = 0x100
	flagSupplementary = 0x800
)

// samRecord holds the fields of a SAM alignment record that gofasta uses
type samRecord struct {
	Name  string
	Flags uint16
	Ref   int // the index of its @SQ line (RNAME), or -1 for "*"
	Pos   int // 0-based (POS - 1)
	Cigar []cigarOp
	Seq   []byte // upper case; nil if SEQ is "*"
}

// seqTable maps the characters of SEQ to what gofasta uses: IUPAC codes and "=" (upper case), the digits 0-3
// as A, C, G and T (as BAM conversion does), and N for anything else
var seqTable = func() (t [256]byte) {
	for i := range t {
		t[i] = 'N'
	}
	for _, c := range []byte("=ACMGRSVTWYHKDBN") {
		t[c] = c
		if 'A' <= c && c <= 'Z' {
			t[c+'a'-'A'] = c
		}
	}
	for i, c := range []byte("ACGT") {
		t['0'+i] = c
	}
	return t
}()

// digits parses a non-negative decimal number of at most 18 digits (the common case), reporting whether b was
// one; anything else is left to strconv, so that its rules and messages apply
func digits(b []byte) (int, bool) {
	if len(b) == 0 || len(b) > 18 {
		return 0, false
	}
	n := 0
	for _, c := range b {
		if c < '0' || c > '9' {
			return 0, false
		}
		n = n*10 + int(c-'0')
	}
	return n, true
}

// parseInt parses a decimal int as strconv.Atoi does
func parseInt(b []byte) (int, error) {
	if n, ok := digits(b); ok {
		return n, nil
	}
	return strconv.Atoi(string(b))
}

// parseRecord parses one alignment line of a SAM file
func (h *samHeader) parseRecord(line []byte) (samRecord, error) {
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
		return samRecord{}, errors.New("sam: missing SAM fields")
	}
	r := samRecord{Name: string(f[0])}

	// (FLAG is read as a Go integer literal, as biogo does: 0x10 is hex, and a leading 0 means octal)
	flags, ok := digits(f[1])
	if !ok || flags > 0xffff || (len(f[1]) > 1 && f[1][0] == '0') {
		v, err := strconv.ParseUint(string(f[1]), 0, 16)
		if err != nil {
			return samRecord{}, fmt.Errorf("sam: failed to parse flags: %w", err)
		}
		flags = int(v)
	}
	r.Flags = uint16(flags)

	var err error
	if r.Ref, err = h.refIndex(f[2]); err != nil {
		return samRecord{}, fmt.Errorf("sam: failed to assign reference: %w", err)
	}
	pos, err := parseInt(f[3])
	if err != nil {
		return samRecord{}, fmt.Errorf("sam: failed to parse position: %w", err)
	}
	r.Pos = pos - 1
	if q, ok := digits(f[4]); !ok || q > 255 {
		if _, err := strconv.ParseUint(string(f[4]), 10, 8); err != nil {
			return samRecord{}, fmt.Errorf("sam: failed to parse map quality: %w", err)
		}
	}
	if r.Cigar, err = parseCigar(f[5]); err != nil {
		if strings.HasPrefix(err.Error(), "error parsing") {
			return samRecord{}, err
		}
		return samRecord{}, fmt.Errorf("sam: failed to parse cigar string: %w", err)
	}
	if !bytes.Equal(f[2], f[6]) && !(len(f[6]) == 1 && f[6][0] == '=') {
		if _, err := h.refIndex(f[6]); err != nil {
			return samRecord{}, fmt.Errorf("sam: failed to assign mate reference: %w", err)
		}
	}
	if _, err := parseInt(f[7]); err != nil {
		return samRecord{}, fmt.Errorf("sam: failed to parse mate position: %w", err)
	}
	if _, err := parseInt(f[8]); err != nil {
		return samRecord{}, fmt.Errorf("sam: failed to parse template length: %w", err)
	}
	if !(len(f[9]) == 1 && f[9][0] == '*') {
		r.Seq = make([]byte, len(f[9]))
		for i, c := range f[9] {
			r.Seq[i] = seqTable[c]
		}
		if len(r.Cigar) != 0 && !cigarValid(r.Cigar, len(r.Seq)) {
			return samRecord{}, errors.New("sam: sequence/CIGAR length mismatch")
		}
	}
	if !(len(f[10]) == 1 && f[10][0] == '*') && len(f[10]) != 0 && len(f[10]) != len(r.Seq) {
		return samRecord{}, errors.New("sam: sequence/quality length mismatch")
	}
	return r, nil
}

// samRef is one @SQ line of a SAM header
type samRef struct {
	Name   string
	Length int
	tags   map[string]string // its other tags (not SN or LN), to compare repeated @SQ lines
}

// samHeader holds what gofasta uses of a SAM header: its references, in order
type samHeader struct {
	refs []samRef
}

// refIndex returns the index of the reference named name, or -1 for "*"
func (h *samHeader) refIndex(name []byte) (int, error) {
	if len(name) == 1 && name[0] == '*' {
		return -1, nil
	}
	for i := range h.refs {
		if h.refs[i].Name == string(name) {
			return i, nil
		}
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
)

// parseSamHeader parses the header lines of a SAM file (each starting with "@"). It checks the layout of
// every line (a known record type, then TAG:VALUE fields; a comment after a tab), that there is a VN by the
// end of each @HD line, that @SQ
// has SN and LN with LN in range, that @RG and @PG have an ID, and that no @SQ, @RG or @PG line repeats a tag
// and no names are repeated; the contents of other fields are not checked.
func parseSamHeader(text []byte) (*samHeader, error) {
	h := &samHeader{}
	seenRefs := map[string]int{}
	seenGroups := map[string]bool{}
	seenProgs := map[string]bool{}
	hasVersion := false // a VN has been seen in an @HD line
	for i, l := range bytes.Split(text, []byte{'\n'}) {
		l = bytes.TrimSuffix(l, []byte{'\r'})
		if len(l) == 0 {
			continue
		}
		if l[0] != '@' || len(l) < 3 {
			return nil, errBadHeader
		}
		var err error
		switch string(l[1:3]) {
		case "HD":
			err = headerLine(l, &hasVersion)
		case "SQ":
			err = h.referenceLine(l, seenRefs)
		case "RG":
			err = idLine(l, seenGroups, errDupReadGroup)
		case "PG":
			err = idLine(l, seenProgs, errDupProgram)
		case "CO":
			if bytes.IndexByte(l, '\t') < 0 {
				err = errBadHeader
			}
		default:
			return nil, errBadHeader
		}
		if err != nil {
			return nil, fmt.Errorf("%w: line %d: %q", err, i+1, l)
		}
	}
	return h, nil
}

// headerFields splits a header line into its TAG:VALUE fields (after the record type); the caller checks
// each with tagField, in order, so that the first problem on the line is the one reported
func headerFields(l []byte, atLeast int) ([][]byte, error) {
	fields := bytes.Split(l, []byte{'\t'})
	if len(fields) < atLeast {
		return nil, errBadHeader
	}
	return fields[1:], nil
}

// tagField checks that a header field is TAG:VALUE (at least three characters, the third a ':')
func tagField(f []byte) error {
	if len(f) < 3 || f[2] != ':' {
		return errBadHeader
	}
	return nil
}

// headerLine checks an @HD line: it, or an @HD line before it, must have a (non-empty) VN
func headerLine(l []byte, hasVersion *bool) error {
	fields, err := headerFields(l, 2)
	if err != nil {
		return err
	}
	for _, f := range fields {
		if err := tagField(f); err != nil {
			return err
		}
		*hasVersion = *hasVersion || (string(f[:2]) == "VN" && len(f) > 3)
	}
	if !*hasVersion {
		return errBadHeader
	}
	return nil
}

// referenceLine reads an @SQ line: SN and LN are needed, no tag may repeat, and a repeated name is
// accepted if it describes the same reference, or replaces the earlier line if that had no tags but SN, LN
// and knownRefTags
func (h *samHeader) referenceLine(l []byte, seen map[string]int) error {
	fields, err := headerFields(l, 3)
	if err != nil {
		return err
	}
	r := samRef{tags: map[string]string{}}
	var hasName, hasLen bool
	tags := map[string]bool{}
	for _, f := range fields {
		if err := tagField(f); err != nil {
			return err
		}
		t, v := string(f[:2]), string(f[3:])
		if tags[t] {
			return errDupTag
		}
		tags[t] = true
		switch t {
		case "SN":
			r.Name, hasName = v, true
		case "LN":
			n, err := strconv.Atoi(v)
			if err != nil {
				return errBadHeader
			}
			if n < 1 || n > 1<<31-1 {
				return errBadLen
			}
			r.Length, hasLen = n, true
		default:
			r.tags[t] = v
		}
	}
	if i, dup := seen[r.Name]; dup && hasName {
		old := h.refs[i]
		switch {
		case sameReference(old, r):
		case otherTags(old) == 0:
			h.refs[i] = r
		default:
			return errDupReference
		}
		return nil
	}
	if !hasName || !hasLen {
		return errBadHeader
	}
	seen[r.Name] = len(h.refs)
	h.refs = append(h.refs, r)
	return nil
}

// knownRefTags are the @SQ tags (besides SN and LN) that two lines for one reference may differ in by one
// giving it and the other not
var knownRefTags = map[string]func(string) string{
	"AS": func(s string) string { return s },
	"M5": strings.ToLower,
	"SP": func(s string) string { return s },
	"UR": func(s string) string {
		if u, err := url.Parse(s); err == nil {
			return u.String()
		}
		return s
	},
}

// otherTags returns how many tags of an @SQ line are neither SN, LN nor one of knownRefTags
func otherTags(r samRef) int {
	n := 0
	for t := range r.tags {
		if _, ok := knownRefTags[t]; !ok {
			n++
		}
	}
	return n
}

// sameReference reports whether two @SQ lines describe the same reference: the same name and length, the
// same assembly, MD5, species and URI where both give one, and the same other tags
func sameReference(a, b samRef) bool {
	if a.Name != b.Name || a.Length != b.Length {
		return false
	}
	known := knownRefTags
	others := 0
	for t, v := range a.tags {
		norm, ok := known[t]
		if !ok {
			others++
			if bv, ok := b.tags[t]; !ok || bv != v {
				return false
			}
			continue
		}
		if bv, ok := b.tags[t]; ok && v != "" && bv != "" && norm(v) != norm(bv) {
			return false
		}
	}
	for t := range b.tags {
		if _, ok := known[t]; !ok {
			others--
		}
	}
	return others == 0
}

// idLine checks an @RG or @PG line: it needs an ID not seen before, and no tag may repeat
func idLine(l []byte, seen map[string]bool, errDup error) error {
	fields, err := headerFields(l, 2)
	if err != nil {
		return err
	}
	tags := map[string]bool{}
	hasID := false
	for _, f := range fields {
		if err := tagField(f); err != nil {
			return err
		}
		t := string(f[:2])
		if tags[t] {
			return errDupTag
		}
		tags[t] = true
		if t == "ID" {
			if seen[string(f[3:])] {
				return errDup
			}
			seen[string(f[3:])] = true
			hasID = true
		}
	}
	if !hasID {
		return errBadHeader
	}
	return nil
}
