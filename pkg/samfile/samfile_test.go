package samfile

import (
	"bytes"
	"errors"
	"fmt"
	"io"
	"regexp"
	"strings"
	"testing"

	biogosam "github.com/biogo/hts/sam"
)

// The SAM parser follows the SAM specification (SAMv1). It is checked against the specification's own
// regular expressions, and against biogo/hts (a test-only dependency): whatever the parser accepts, biogo must
// accept too (the parser is stricter) and read the same way.

// The mandatory fields of an alignment line, as the specification gives them (section 1.4), with
// [:rname:^*=][:rname:]* written out (section 1.2.1)
var specFields = func() []*regexp.Regexp {
	rname := `[0-9A-Za-z!#$%&+./:;?@^_|~-][0-9A-Za-z!#$%&*+./:;=?@^_|~-]*`
	var out []*regexp.Regexp
	for _, e := range []string{
		`[!-?A-~]{1,254}`, `[0-9]+`, `\*|` + rname, `[0-9]+`, `[0-9]+`,
		`\*|([0-9]+[MIDNSHP=X])+`, `\*|=|` + rname, `[0-9]+`, `-?[0-9]+`, `\*|[A-Za-z=.]+`, `[!-~]+`,
	} {
		out = append(out, regexp.MustCompile(`^(`+e+`)$`))
	}
	return out
}()

// A header line, as the specification gives it (section 1.3)
var specHeaderLine = regexp.MustCompile(`^@(HD|SQ|RG|PG)(\t[A-Za-z][A-Za-z0-9]:[ -~]+)+$|^@CO\t.*$`)

func biogoHeader(text []byte) (h *biogosam.Header, err error) {
	defer func() {
		if p := recover(); p != nil {
			err = fmt.Errorf("panic: %v", p)
		}
	}()
	h, err = biogosam.NewHeader(nil, nil)
	if err != nil {
		return nil, err
	}
	return h, h.UnmarshalText(text)
}

func biogoRecord(h *biogosam.Header, line []byte) (r biogosam.Record, err error) {
	defer func() {
		if p := recover(); p != nil {
			err = fmt.Errorf("panic: %v", p)
		}
	}()
	err = r.UnmarshalSAM(h, line)
	return r, err
}

// sameRecord reports how our record differs from biogo's, if it does
func sameRecord(ours Record, theirs biogosam.Record) error {
	ref := -1
	if theirs.Ref != nil {
		ref = theirs.Ref.ID()
	}
	switch {
	case ours.Name != theirs.Name:
		return fmt.Errorf("name %q, biogo %q", ours.Name, theirs.Name)
	case ours.Flags != uint16(theirs.Flags):
		return fmt.Errorf("flags %d, biogo %d", ours.Flags, theirs.Flags)
	case ours.Ref != ref:
		return fmt.Errorf("reference %d, biogo %d", ours.Ref, ref)
	case ours.Pos != theirs.Pos:
		return fmt.Errorf("position %d, biogo %d", ours.Pos, theirs.Pos)
	case !bytes.Equal(ours.Seq, theirs.Seq.Expand()):
		return fmt.Errorf("SEQ %q, biogo %q", ours.Seq, theirs.Seq.Expand())
	}
	// biogo splits operations longer than 2^28-1 (a BAM limit) into several
	var theirOps []CigarOp
	for _, op := range theirs.Cigar {
		if n := len(theirOps); n > 0 && theirOps[n-1].Type() == CigarType(op.Type()) && theirOps[n-1].Len() > 0 && theirOps[n-1].Len()%(1<<28-1) == 0 {
			theirOps[n-1] = NewCigarOp(CigarType(op.Type()), theirOps[n-1].Len()+op.Len())
			continue
		}
		theirOps = append(theirOps, NewCigarOp(CigarType(op.Type()), op.Len()))
	}
	if fmt.Sprint(ours.Cigar) != fmt.Sprint(theirOps) {
		return fmt.Errorf("CIGAR %v, biogo %v", ours.Cigar, theirOps)
	}
	return nil
}

const testHeader = "@HD\tVN:1.6\n@SQ\tSN:ref\tLN:100\n@SQ\tSN:other\tLN:50\n"

func TestRecordsFollowTheSpecification(t *testing.T) {
	cases := []struct {
		line string
		ok   bool
	}{
		{"q\t0\tref\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\t*", true},
		{"q\t2048\tref\t5\t60\t3S10M2H\t*\t0\t0\tTTTACGTACGTAC\t*\tNM:i:0\tms:i:10", true},
		{"q\t0016\tother\t1\t60\t4M\t=\t1\t-4\tacgt\tIIII", true},  // leading zeros are decimal; lower case
		{"q\t0\tref\t1\t60\t4M\t*\t0\t0\tAC=.\t*", true},           // = and .
		{"q\t0\tref\t1\t60\t4M\t*\t0\t0\tACUX\t*", true},           // letters that are not IUPAC codes
		{"q\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*", true},                  // unmapped
		{"q\t0\tref\t1\t60\t4M\t*\t0\t0\tACGT\t*\tXX:i:abc", true}, // optional fields are not checked
		{"q\t0x10\tref\t1\t60\t4M\t*\t0\t0\tACGT\t*", false},       // FLAG is decimal
		{"q\t70000\tref\t1\t60\t4M\t*\t0\t0\tACGT\t*", false},
		{"q\t0\tref\t+1\t60\t4M\t*\t0\t0\tACGT\t*", false}, // no sign
		{"q\t0\tref\t2147483648\t60\t4M\t*\t0\t0\tACGT\t*", false},
		{"q\t0\tref\t1\t256\t4M\t*\t0\t0\tACGT\t*", false},
		{"q\t0\tref\t1\t60\t2M5\t*\t0\t0\tACGT\t*", false},    // a number after the last operation
		{"q\t0\tref\t1\t60\tM4M\t*\t0\t0\tACGT\t*", false},    // an operation with no number
		{"q\t0\tref\t1\t60\t4B\t*\t0\t0\tACGT\t*", false},     // B is not a SAM operation
		{"q\t0\tref\t1\t60\t2M2H2M\t*\t0\t0\tACGT\t*", false}, // a hard clip inside
		{"q\t0\tref\t1\t60\t1M2S1M\t*\t0\t0\tACGT\t*", false}, // a soft clip inside
		{"q\t0\tref\t1\t60\t2H2S2M\t*\t0\t0\tACGT\t*", true},
		{"q\t0\tref\t1\t60\t4M\t*\t0\t0\tACG\t*", false}, // CIGAR and SEQ disagree
		{"q\t0\tref\t1\t60\t4M\t*\t0\t0\tACGT\tIII", false},
		{"q\t0\tref\t1\t60\t4M\t*\t0\t0\tACGT\t", false},  // empty QUAL
		{"q\t0\tref\t1\t60\t4M\t*\t0\t0\t*\tIIII", false}, // QUAL without SEQ
		{"\t0\tref\t1\t60\t4M\t*\t0\t0\tACGT\t*", false},  // empty QNAME
		{"@q\t0\tref\t1\t60\t4M\t*\t0\t0\tACGT\t*", false},
		{"q r\t0\tref\t1\t60\t4M\t*\t0\t0\tACGT\t*", false},
		{"q\t0\tnowhere\t1\t60\t4M\t*\t0\t0\tACGT\t*", false},
		{"q\t0\tref\t1\t60\t4M\tnowhere\t0\t0\tACGT\t*", false},
		{"q\t0\tref\t1\t60\t4M\t*\t-1\t0\tACGT\t*", false},
		{"q\t0\tref\t1\t60\t4M\t*\t0\t+5\tACGT\t*", false},
		{"q\t0\tref\t1\t60\t4M\t*\t0\t0\tACGT", false}, // ten fields
	}
	h, err := ParseHeader([]byte(testHeader))
	if err != nil {
		t.Fatal(err)
	}
	for _, c := range cases {
		_, err := h.ParseRecord([]byte(c.line))
		if (err == nil) != c.ok {
			t.Errorf("%q: got error %v, want ok = %v", c.line, err, c.ok)
		}
	}
}

func TestHeadersFollowTheSpecification(t *testing.T) {
	cases := []struct {
		header string
		ok     bool
	}{
		{"@HD\tVN:1.6\tSO:unsorted\tGO:query\n@SQ\tSN:ref\tLN:29903\n@PG\tID:minimap2\tPN:minimap2\tVN:2.28\tCL:minimap2 -a\n", true},
		{"@SQ\tSN:ref\tLN:10\n@CO\tany text\n@RG\tID:a\tDT:yesterday\n", true},
		{"@SQ\tSN:ref\tLN:10\r\n", true},
		{"@SQ\tSN:chr1\tLN:10\tAN:1,chrA\n@SQ\tSN:chr2\tLN:5\n", true},
		{"@SQ\tSN:ref\tLN:10\n@SQ\tSN:ref\tLN:10\n", false},  // a name used twice
		{"@SQ\tSN:a\tLN:10\tAN:b\n@SQ\tSN:b\tLN:5\n", false}, // an alternative name used again
		{"@SQ\tSN:*ref\tLN:10\n", false},
		{"@SQ\tSN:re,f\tLN:10\n", false},
		{"@SQ\tSN:ref\n", false},
		{"@SQ\tLN:10\n", false},
		{"@SQ\tSN:ref\tLN:0\n", false},
		{"@SQ\tSN:ref\tLN:2147483648\n", false},
		{"@SQ\tSN:ref\tLN:10\tLN:10\n", false},
		{"@SQ\tSN:ref\tLN:10\txy:\n", false},  // an empty value
		{"@SQ\tSN:ref\tLN:10\t1x:a\n", false}, // a tag must start with a letter
		{"@HD\tSO:unsorted\n", false},         // no VN
		{"@HD\tVN:1.6\tVN:1.6\n", false},
		{"@HD\tVN:v1\n", false},
		{"@SQ\tSN:ref\tLN:10\n@HD\tVN:1.6\n", false}, // @HD not first
		{"@HD\tVN:1.6\n@HD\tVN:1.6\n", false},
		{"@CO\n", false},
		{"@XX\tAB:c\n", false},
		{"@RG\tID:a\n@RG\tID:a\n", false},
		{"@RG\tSM:x\n", false},
		{"@PG\tID:a\n@PG\tID:a\n", false},
	}
	for _, c := range cases {
		_, err := ParseHeader([]byte(c.header))
		if (err == nil) != c.ok {
			t.Errorf("%q: got error %v, want ok = %v", c.header, err, c.ok)
		}
	}
}

func TestSeqTableFollowsTheSpecification(t *testing.T) {
	for c := 0; c < 256; c++ {
		allowed := 'A' <= c && c <= 'Z' || 'a' <= c && c <= 'z' || c == '=' || c == '.'
		if (seqTable[c] != 0) != allowed {
			t.Errorf("SEQ character %q: allowed = %v, want %v", byte(c), seqTable[c] != 0, allowed)
		}
	}
}

// A SEQ character the specification does not allow (such as a FASTA alignment's "-") is read as N, and marked
func TestOtherSeqCharacters(t *testing.T) {
	h, err := ParseHeader([]byte(testHeader))
	if err != nil {
		t.Fatal(err)
	}
	for _, c := range []struct {
		seq, want string
		other     bool
	}{
		{"AC-T", "ACNT", true},
		{"AC0T", "ACNT", true},
		{"_c*\xe9", "NCNN", true},
		{"acg.", "ACGN", false},
	} {
		rec, err := h.ParseRecord([]byte("q\t0\tref\t1\t60\t4M\t*\t0\t0\t" + c.seq + "\t*"))
		if err != nil || string(rec.Seq) != c.want || rec.OtherSeq != c.other {
			t.Errorf("SEQ %q: got %q, OtherSeq %v, %v; want %q, %v", c.seq, rec.Seq, rec.OtherSeq, err, c.want, c.other)
		}
	}
}

func TestResolveEquals(t *testing.T) {
	h, err := ParseHeader([]byte(testHeader))
	if err != nil {
		t.Fatal(err)
	}
	ref := []byte(strings.Repeat("ACGT", 25))
	for _, c := range []struct {
		line, withRef, withoutRef string
		bad                       bool
	}{
		{line: "q\t0\tref\t3\t60\t6M\t*\t0\t0\tGT=CG=\t*", withRef: "GTACGT", withoutRef: "GTNCGN"},
		{line: "q\t0\tref\t1\t60\t2S2M2D2M\t*\t0\t0\tTT==C=\t*", withRef: "TTACCC", withoutRef: "TTNNCN"},
		{line: "q\t0\tref\t1\t60\t2M1I2M\t*\t0\t0\tAC=GT\t*", bad: true}, // in an insertion
		{line: "q\t0\tref\t1\t60\t1S3M\t*\t0\t0\t=CGT\t*", bad: true},    // in a soft clip
		{line: "q\t0\tref\t1\t60\t4M\t*\t0\t0\tACGT\t*", withRef: "ACGT", withoutRef: "ACGT"},
	} {
		for _, r := range [][]byte{ref, nil} {
			rec, err := h.ParseRecord([]byte(c.line))
			if err != nil {
				t.Fatal(err)
			}
			_, err = rec.ResolveEquals(r)
			want := c.withRef
			if r == nil {
				want = c.withoutRef
			}
			if c.bad {
				if err == nil {
					t.Errorf("%q: want an error", c.line)
				}
			} else if err != nil || string(rec.Seq) != want {
				t.Errorf("%q (reference given: %v): got %q, %v; want %q", c.line, r != nil, rec.Seq, err, want)
			}
		}
	}
}

// FuzzRecordsFollowTheSpecification: every line the parser accepts has mandatory fields that match the
// specification's expressions, and biogo (which is less strict) reads it the same way (its optional fields
// aside, which neither checks the same way)
func FuzzRecordsFollowTheSpecification(f *testing.F) {
	for _, s := range []string{
		"q\t0\tref\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\t*",
		"q\t2048\tref\t5\t60\t3S10M2H\t*\t0\t0\tTTTACGTACGTAC\t*\tNM:i:0",
		"q\t16\tother\t1\t60\t2M1I1D1M\t=\t1\t-5\tacgt\tIIII",
		"q\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*",
	} {
		f.Add(s)
	}
	bh, err := biogoHeader([]byte(testHeader))
	if err != nil {
		f.Fatal(err)
	}
	h, err := ParseHeader([]byte(testHeader))
	if err != nil {
		f.Fatal(err)
	}
	f.Fuzz(func(t *testing.T, line string) {
		if strings.Contains(line, "\n") {
			return
		}
		ours, err := h.ParseRecord([]byte(line))
		if err != nil {
			return
		}
		fields := strings.SplitN(line, "\t", 12)
		// SEQ may hold characters the specification does not allow (read as N, and marked): biogo reads some
		// of them as bases, so it is given N there
		if seqOK := specFields[9].MatchString(fields[9]); seqOK == ours.OtherSeq {
			t.Fatalf("%q: OtherSeq is %v, but SEQ %q matches the specification: %v", line, ours.OtherSeq, fields[9], seqOK)
		}
		if ours.OtherSeq {
			seq := []byte(fields[9])
			for i, c := range seq {
				if seqTable[c] == 0 {
					seq[i] = 'N'
				}
			}
			fields[9] = string(seq)
		}
		for i, re := range specFields {
			if !re.MatchString(fields[i]) {
				t.Fatalf("%q: field %d (%q) does not match the specification's %v", line, i+1, fields[i], re)
			}
		}
		if strings.HasPrefix(fields[1], "0") && fields[1] != "0" {
			fields[1] = strings.TrimLeft(fields[1], "0") // (biogo reads a leading 0 as octal)
			if fields[1] == "" {
				fields[1] = "0"
			}
		}
		theirs, err := biogoRecord(bh, []byte(strings.Join(fields[:11], "\t")))
		if err != nil {
			t.Fatalf("%q: accepted, but biogo rejects it: %v", line, err)
		}
		if err := sameRecord(ours, theirs); err != nil {
			t.Fatalf("%q: %v", line, err)
		}
	})
}

// FuzzHeadersFollowTheSpecification: every header the parser accepts has lines that match the
// specification's expression, and biogo reads the same references from it (where biogo accepts it at all:
// it also checks the contents of @RG, @PG and some @SQ fields)
func FuzzHeadersFollowTheSpecification(f *testing.F) {
	for _, s := range []string{testHeader, "@SQ\tSN:ref\tLN:10\tAN:r1,r2\n", "@CO\tx\n@RG\tID:a\n"} {
		f.Add(s)
	}
	f.Fuzz(func(t *testing.T, header string) {
		ours, err := ParseHeader([]byte(header))
		if err != nil {
			return
		}
		for _, l := range strings.Split(header, "\n") {
			if l = strings.TrimSuffix(l, "\r"); l != "" && !specHeaderLine.MatchString(l) {
				t.Fatalf("%q: line %q does not match the specification's expression", header, l)
			}
		}
		theirs, err := biogoHeader([]byte(header))
		if err != nil {
			return
		}
		refs := theirs.Refs()
		if len(refs) != len(ours.Refs) {
			t.Fatalf("%q: %d references, biogo %d", header, len(ours.Refs), len(refs))
		}
		for i, r := range refs {
			if ours.Refs[i].Name != r.Name() || ours.Refs[i].Length != r.Len() {
				t.Fatalf("%q: reference %d is %v, biogo %s:%d", header, i, ours.Refs[i], r.Name(), r.Len())
			}
		}
	})
}

// The Reader reads the header, then one alignment line at a time, skipping blank lines; a bad line gives a
// *LineError with its QNAME, and the lines after it can still be read
func TestReader(t *testing.T) {
	if _, err := NewReader(strings.NewReader("")); !errors.Is(err, ErrEmpty) {
		t.Errorf("empty file: got %v, want ErrEmpty", err)
	}
	if _, err := NewReader(strings.NewReader("@SQ\tSN:ref\n")); err == nil {
		t.Error("bad header: want an error")
	}
	text := testHeader +
		"a\t0\tref\t1\t60\t4M\t*\t0\t0\tACGT\t*\r\n" +
		"\n" +
		"b\t0\tref\t1\t60\t4X\t*\t0\t0\tACGT\t*\n" + // the next line is bad (RNAME is not an @SQ)
		"c\t0\tnowhere\t1\t60\t4M\t*\t0\t0\tACGT\t*\n" +
		"d\t16\tref\t5\t60\t2M\t*\t0\t0\tAC\t*" // no final newline
	r, err := NewReader(strings.NewReader(text))
	if err != nil {
		t.Fatal(err)
	}
	if len(r.Header.Refs) != 2 || r.Header.Refs[0] != (Reference{Name: "ref", Length: 100}) {
		t.Errorf("header references: %v", r.Header.Refs)
	}
	var got []string
	for {
		rec, err := r.Next()
		if err == io.EOF {
			break
		}
		var bad *LineError
		switch {
		case errors.As(err, &bad):
			got = append(got, "bad "+bad.Name)
		case err != nil:
			t.Fatal(err)
		default:
			got = append(got, fmt.Sprintf("%s %d %d %s", rec.Name, rec.Flags, rec.Pos, rec.Seq))
		}
	}
	want := []string{"a 0 0 ACGT", "b 0 0 ACGT", "bad c", "d 16 4 AC"}
	if fmt.Sprint(got) != fmt.Sprint(want) {
		t.Errorf("got %q, want %q", got, want)
	}
}
