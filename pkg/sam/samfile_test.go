package sam

import (
	"bytes"
	"fmt"
	"strings"
	"testing"

	biogosam "github.com/biogo/hts/sam"
)

// The SAM parser is checked against biogo/hts, which gofasta used before (a test-only dependency): whatever
// biogo accepts, the parser must accept and read the same way; whatever biogo rejects, it must reject too,
// except for what it deliberately does not check (optional tags, @RG and @PG contents, @SQ M5 and UR).

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

// biogoRecord parses a record with biogo, turning a panic into an error as gofasta used to
func biogoRecord(h *biogosam.Header, line []byte) (r biogosam.Record, err error) {
	defer func() {
		if p := recover(); p != nil {
			err = fmt.Errorf("error parsing a sam record: %v", p)
		}
	}()
	err = r.UnmarshalSAM(h, line)
	return r, err
}

// sameHeader reports how our header differs from biogo's, if it does
func sameHeader(ours *samHeader, theirs *biogosam.Header) error {
	refs := theirs.Refs()
	if len(ours.refs) != len(refs) {
		return fmt.Errorf("%d references, biogo %d", len(ours.refs), len(refs))
	}
	for i, r := range refs {
		if ours.refs[i].Name != r.Name() || ours.refs[i].Length != r.Len() {
			return fmt.Errorf("reference %d is %s:%d, biogo %s:%d", i, ours.refs[i].Name, ours.refs[i].Length, r.Name(), r.Len())
		}
	}
	return nil
}

// sameRecord reports how our record differs from biogo's, if it does
func sameRecord(ours samRecord, theirs biogosam.Record) error {
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
	case len(ours.Cigar) != len(theirs.Cigar):
		return fmt.Errorf("%d CIGAR operations, biogo %d", len(ours.Cigar), len(theirs.Cigar))
	}
	for i, op := range theirs.Cigar {
		if ours.Cigar[i].Len() != op.Len() || int(ours.Cigar[i].Type()) != int(op.Type()) {
			return fmt.Errorf("CIGAR operation %d is %d%v, biogo %v", i, ours.Cigar[i].Len(), ours.Cigar[i].Type(), op)
		}
	}
	return nil
}

const testHeader = "@HD\tVN:1.6\n@SQ\tSN:ref\tLN:100\n@SQ\tSN:other\tLN:50\n"

func TestRecordsMatchBiogo(t *testing.T) {
	cases := []struct {
		line             string
		onlyBiogoRejects bool // a check gofasta deliberately does not make
	}{
		{line: "q\t0\tref\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\t*"},
		{line: "q\t2048\tref\t5\t60\t3S10M2H\t*\t0\t0\tTTTACGTACGTAC\t*\tNM:i:0\tms:i:10"},
		{line: "q\t0x10\tref\t1\t60\t4M\t=\t1\t0\tacgt\tIIII"},                             // hex flags; lower case
		{line: "q\t16\tother\t1\t60\t4M\t*\t0\t0\tAC-X\t*"},                                // characters that are not IUPAC
		{line: "q\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*"},                                          // unmapped
		{line: "q\t0\tref\t+1\t60\t4M\t*\t0\t0\tACGT\t*"},                                  // a sign
		{line: "q\t0\tref\t1\t60\t2M5\t*\t0\t0\tACGT\t*"},                                  // a number after the last operation
		{line: "q\t0\tref\t1\t60\tM4M\t*\t0\t0\tACGT\t*"},                                  // an operation with no number
		{line: "q\t0\tref\t1\t60\t4M\t*\t0\t0\tACGT\t"},                                    // empty QUAL
		{line: "q\t0\tref\t1\t60\t4M\t*\t0\t0\tACGT\tIII"},                                 // QUAL of the wrong length
		{line: "q\t0\tref\t1\t60\t4M\t*\t0\t0\t*\tIIII"},                                   // QUAL without SEQ
		{line: "q\t0\tref\t1\t60\t4M\t*\t0\t0\tACG\t*"},                                    // CIGAR and SEQ disagree
		{line: "q\t0\tref\t1\t60\t2M2H2M\t*\t0\t0\tACGT\t*"},                               // a hard clip inside
		{line: "q\t0\tref\t1\t60\t4Q\t*\t0\t0\tACGT\t*"},                                   // an unknown operation
		{line: "q\t0\tref\t1\t60\t99999999999999M\t*\t0\t0\tACGT\t*"},                      // too long
		{line: "q\t0\tnowhere\t1\t60\t4M\t*\t0\t0\tACGT\t*"},                               // unknown RNAME
		{line: "q\t0\tref\t1\t60\t4M\tnowhere\t0\t0\tACGT\t*"},                             // unknown RNEXT
		{line: "q\t0\tref\t1\t60\t4M\tother\t0\t0\tACGT\t*"},                               // another RNEXT
		{line: "q\t70000\tref\t1\t60\t4M\t*\t0\t0\tACGT\t*"},                               // flags too large
		{line: "q\t0\tref\tx\t60\t4M\t*\t0\t0\tACGT\t*"},                                   // a bad POS
		{line: "q\t0\tref\t1\t256\t4M\t*\t0\t0\tACGT\t*"},                                  // MAPQ too large
		{line: "q\t0\tref\t1\t60\t4M\t*\tx\t0\tACGT\t*"},                                   // a bad PNEXT
		{line: "q\t0\tref\t1\t60\t4M\t*\t0\tx\tACGT\t*"},                                   // a bad TLEN
		{line: "q\t0\tref\t1\t60\t4M\t*\t0\t0\tACGT"},                                      // ten fields
		{line: "q\t0\tref\t1\t60\t4M\t*\t0\t0\tACGT\t*\tXX:i:abc", onlyBiogoRejects: true}, // a malformed optional tag
		{line: "q\t0\tref\t1\t60\t4M\t*\t0\t0\tACGT\t*\tbroken", onlyBiogoRejects: true},
	}
	bh, err := biogoHeader([]byte(testHeader))
	if err != nil {
		t.Fatal(err)
	}
	h, err := parseSamHeader([]byte(testHeader))
	if err != nil {
		t.Fatal(err)
	}
	for _, c := range cases {
		theirs, terr := biogoRecord(bh, []byte(c.line))
		ours, oerr := h.parseRecord([]byte(c.line))
		switch {
		case terr == nil && oerr != nil:
			t.Errorf("%q: biogo accepts it, gofasta does not: %v", c.line, oerr)
		case terr == nil:
			if err := sameRecord(ours, theirs); err != nil {
				t.Errorf("%q: %v", c.line, err)
			}
		case oerr == nil && !c.onlyBiogoRejects:
			t.Errorf("%q: biogo rejects it (%v), gofasta does not", c.line, terr)
		case oerr != nil && oerr.Error() != terr.Error():
			t.Errorf("%q: error %q, biogo %q", c.line, oerr, terr)
		}
	}
}

func TestHeadersMatchBiogo(t *testing.T) {
	cases := []struct {
		header           string
		onlyBiogoRejects bool
	}{
		{header: "@HD\tVN:1.6\tSO:unsorted\tGO:query\n@SQ\tSN:ref\tLN:29903\n@PG\tID:minimap2\tPN:minimap2\tVN:2.28\tCL:minimap2 -a\n"},
		{header: "@SQ\tSN:ref\tLN:10\n@SQ\tSN:ref\tLN:10\n"},       // the same @SQ twice
		{header: "@SQ\tSN:ref\tLN:10\n@SQ\tSN:ref\tLN:20\n"},       // a repeated name replaces a plain line
		{header: "@SQ\tSN:ref\tLN:10\tSP:x\n@SQ\tSN:ref\tLN:20\n"}, // but not one with more tags
		{header: "@SQ\tSN:ref\tLN:10\tSP:x\n@SQ\tSN:ref\tLN:10\n"}, // the same reference
		{header: "@SQ\tSN:ref\tLN:10\tM5:ABCDEF0123456789ABCDEF0123456789\n@SQ\tSN:ref\tLN:10\tM5:abcdef0123456789abcdef0123456789\n"},
		{header: "@SQ\tSN:ref\tLN:10\tXY:a\n@SQ\tSN:ref\tLN:10\tXY:b\n"}, // other tags differ
		{header: "@SQ\tSN:ref\n"},       // no LN
		{header: "@SQ\tLN:10\n"},        // no SN
		{header: "@SQ\tSN:ref\tLN:0\n"}, // LN out of range
		{header: "@SQ\tSN:ref\tLN:x\n"},
		{header: "@SQ\tSN:ref\tLN:10\tLN:10\n"}, // a repeated tag
		{header: "@SQ\tSN\tLN:10\n"},            // a field without a ':'
		{header: "@HD\tSO:unsorted\n"},          // no VN
		{header: "@HD\tVN:\tVN:1.6\n"},
		{header: "@HD\tVN:1.6\tVN:1.5\n", onlyBiogoRejects: true}, // @HD contents
		{header: "@HD\tVN:1.6\tSO:coordinate\tSO:queryname\n", onlyBiogoRejects: true},
		{header: "@CO\tany text\n"},
		{header: "@CO\n"},       // a comment needs a tab
		{header: "@XX\tAB:c\n"}, // an unknown record type
		{header: "SQ\tSN:ref\tLN:10\n"},
		{header: "@RG\tID:a\n@RG\tID:a\n"},                            // a repeated read group
		{header: "@RG\tSM:x\n"},                                       // no ID
		{header: "@RG\tID:a\tDT:yesterday\n", onlyBiogoRejects: true}, // @RG contents
		{header: "@PG\tID:a\n@PG\tID:a\n"},
		{header: "@SQ\tSN:ref\tLN:10\tM5:xyz\n", onlyBiogoRejects: true}, // @SQ M5
		{header: "@SQ\tSN:ref\tLN:10\r\n"},
		{header: "\n@SQ\tSN:ref\tLN:10\n\n"},
	}
	for _, c := range cases {
		theirs, terr := biogoHeader([]byte(c.header))
		ours, oerr := parseSamHeader([]byte(c.header))
		switch {
		case terr == nil && oerr != nil:
			t.Errorf("%q: biogo accepts it, gofasta does not: %v", c.header, oerr)
		case terr == nil:
			if err := sameHeader(ours, theirs); err != nil {
				t.Errorf("%q: %v", c.header, err)
			}
		case oerr == nil && !c.onlyBiogoRejects:
			t.Errorf("%q: biogo rejects it (%v), gofasta does not", c.header, terr)
		}
	}
}

// FuzzRecordsMatchBiogo checks the rule on random records: where biogo accepts a line, the parser reads it
// the same way; where only the parser accepts it, biogo objected only to its optional tags.
func FuzzRecordsMatchBiogo(f *testing.F) {
	for _, s := range []string{
		"q\t0\tref\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\t*",
		"q\t2048\tref\t5\t60\t3S10M2H\t*\t0\t0\tTTTACGTACGTAC\t*\tNM:i:0",
		"q\t16\tother\t1\t60\t2M1I1D1M\t=\t1\t0\tacgt\tIIII",
		"q\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*",
	} {
		f.Add(s)
	}
	bh, err := biogoHeader([]byte(testHeader))
	if err != nil {
		f.Fatal(err)
	}
	h, err := parseSamHeader([]byte(testHeader))
	if err != nil {
		f.Fatal(err)
	}
	f.Fuzz(func(t *testing.T, line string) {
		if strings.ContainsAny(line, "\n") {
			return
		}
		theirs, terr := biogoRecord(bh, []byte(line))
		ours, oerr := h.parseRecord([]byte(line))
		switch {
		case terr == nil && oerr != nil:
			t.Fatalf("%q: biogo accepts it, gofasta does not: %v", line, oerr)
		case terr == nil:
			if err := sameRecord(ours, theirs); err != nil {
				t.Fatalf("%q: %v", line, err)
			}
		case oerr != nil:
			// (the same message, as it reaches the user in gofasta's warnings)
			if oerr.Error() != terr.Error() {
				t.Fatalf("%q: error %q, biogo %q", line, oerr, terr)
			}
		case oerr == nil:
			// biogo must accept the line without its optional tags
			fields := strings.SplitN(line, "\t", 12)
			if len(fields) < 12 {
				t.Fatalf("%q: biogo rejects it (%v), gofasta does not", line, terr)
			}
			if _, err := biogoRecord(bh, []byte(strings.Join(fields[:11], "\t"))); err != nil {
				t.Fatalf("%q: biogo rejects it (%v), gofasta does not", line, terr)
			}
		}
	})
}

// FuzzHeadersMatchBiogo: where biogo accepts a header, the parser reads the same references
func FuzzHeadersMatchBiogo(f *testing.F) {
	for _, s := range []string{testHeader, "@SQ\tSN:ref\tLN:10\n@SQ\tSN:ref\tLN:20\n", "@CO\tx\n@RG\tID:a\n"} {
		f.Add(s)
	}
	f.Fuzz(func(t *testing.T, header string) {
		theirs, terr := biogoHeader([]byte(header))
		ours, oerr := parseSamHeader([]byte(header))
		if terr == nil {
			if oerr != nil {
				t.Fatalf("%q: biogo accepts it, gofasta does not: %v", header, oerr)
			}
			if err := sameHeader(ours, theirs); err != nil {
				t.Fatalf("%q: %v", header, err)
			}
			return
		}
		// the same message, where biogo's error is one of the checks the parser makes too (not one of the field
		// contents it leaves alone, or a panic on a field too short to hold a tag)
		shared := false
		for _, e := range []error{errBadHeader, errDupTag, errDupReference, errDupReadGroup, errDupProgram, errBadLen} {
			shared = shared || strings.HasPrefix(terr.Error(), e.Error()+":")
		}
		unchecked := false // a field whose contents the parser does not check, which biogo may object to first
		for _, x := range []string{"M5:", "UR:", "@HD", "@RG", "@PG"} {
			unchecked = unchecked || strings.Contains(header, x)
		}
		if oerr != nil && shared && !unchecked && oerr.Error() != terr.Error() {
			t.Fatalf("%q: error %q, biogo %q", header, oerr, terr)
		}
	})
}

func TestSeqTableMatchesBiogo(t *testing.T) {
	for c := 0; c < 256; c++ {
		if want := biogosam.NewSeq([]byte{byte(c)}).Expand()[0]; seqTable[c] != want {
			t.Errorf("SEQ character %q is %q, biogo %q", byte(c), seqTable[c], want)
		}
	}
}
