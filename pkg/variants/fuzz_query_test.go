package variants

import (
	"bytes"
	"os"
	"strconv"
	"strings"
	"testing"
	"time"
)

// fuzzRef is the reference of the variants fuzz targets, with fuzzGenbank and fuzzGFF its annotation: a
// forward CDS over the whole sequence and a reverse-strand one overlapping it in another frame
const fuzzRef = "ATGAAACCCGGGTTTAAACCCGGGTTTAAACCCGGGTTTAAATAA"

const fuzzGenbank = `LOCUS       X                 45 bp    RNA     linear   VRL 01-JAN-2000
FEATURES             Location/Qualifiers
     CDS             1..45
                     /gene="a"
                     /codon_start=1
                     /translation="MKPGFKPGFKPGFK"
     CDS             complement(3..32)
                     /gene="b"
                     /codon_start=1
                     /translation="GLNPGLNPGF"
ORIGIN
        1 atgaaacccg ggtttaaacc cgggtttaaa cccgggttta aataa
//
`

const fuzzGFF = "##gff-version 3\n##sequence-region r 1 45\n" +
	"r\t.\tCDS\t1\t45\t.\t+\t0\tID=a;Name=a\n" +
	"r\t.\tCDS\t3\t32\t.\t-\t0\tID=b;Name=b\n"

// fuzzOptions are the output flags of one Variants call, taken from the bits of mode
type fuzzOptions struct {
	annoSuffix, anno              string
	start, end                    int
	aggregate, appendSNP, groupAA bool
}

func optionsFromMode(mode uint8) fuzzOptions {
	o := fuzzOptions{start: -1, end: -1}
	switch mode % 3 {
	case 0:
		o.annoSuffix, o.anno = "gb", fuzzGenbank
	case 1:
		o.annoSuffix, o.anno = "gff", fuzzGFF
	}
	o.aggregate, o.appendSNP, o.groupAA = mode&4 != 0, mode&8 != 0, mode&16 != 0
	if mode&32 != 0 {
		o.start, o.end = 10, 40
	}
	return o
}

// runVariants runs Variants on msa with o, the reference refID in msa, or, if refID is "", fuzzRef as its own file
func runVariants(msa []byte, stdin bool, refID string, o fuzzOptions, threads int) (string, error) {
	var out bytes.Buffer
	var refIn *strings.Reader
	if refID == "" {
		refIn = strings.NewReader(">ref\n" + fuzzRef + "\n")
	}
	var annoIn *strings.Reader
	if o.anno != "" {
		annoIn = strings.NewReader(o.anno)
	}
	var err error
	switch {
	case refIn == nil && annoIn == nil:
		err = Variants(bytes.NewReader(msa), stdin, refID, nil, nil, "", &out, o.start, o.end, o.aggregate, 0, o.appendSNP, o.groupAA, threads, false)
	case refIn == nil:
		err = Variants(bytes.NewReader(msa), stdin, refID, nil, annoIn, o.annoSuffix, &out, o.start, o.end, o.aggregate, 0, o.appendSNP, o.groupAA, threads, false)
	case annoIn == nil:
		err = Variants(bytes.NewReader(msa), stdin, refID, refIn, nil, "", &out, o.start, o.end, o.aggregate, 0, o.appendSNP, o.groupAA, threads, false)
	default:
		err = Variants(bytes.NewReader(msa), stdin, refID, refIn, annoIn, o.annoSuffix, &out, o.start, o.end, o.aggregate, 0, o.appendSNP, o.groupAA, threads, false)
	}
	return out.String(), err
}

// FuzzVariantsQuery: arbitrary bytes as the --query alignment, with fuzzRef as the reference and any of the
// output flags, must give an error or output, never a panic or a hang.
func FuzzVariantsQuery(f *testing.F) {
	seeds := []string{
		">q1\n" + fuzzRef + "\n",
		">q1\nATGAAACCCGGGTTTAAACCCGGGTTTAAACCCGGGTTTAAATAA\n>q2\nATGAAAC--GGGTTTAAACCCGGGTTCAAACCCGGGTTTAAATAA\n",
		">q1\nATGNNNCCCGGGTTTAAACCCRGGTTTAAACCCGGGTTTAAATAG\n>q2 description\natgaaacccgggtttaaacc\ncgggtttaaacccgggtttaaataa\n",
		">q1\n---AAACCCGGGTTTAAACCCGGGTTTAAACCCGGGTTTAAA---\n>q2\nATGAAACCCGGG\n",
		">ref\n" + fuzzRef + "\n>q1\nATGAAACCCGGGTTTAAAGCCGGGTTTAAACCCGGGTTTAAATAA\n",
		"",
	}
	for i, s := range seeds {
		f.Add([]byte(s), uint8(i*13))
	}
	f.Fuzz(func(t *testing.T, data []byte, mode uint8) {
		quietStderr(t)
		refID := ""
		if mode&64 != 0 {
			refID = "ref"
		}
		collect(data, func() {
			runWithTimeout(t, func() { _, _ = runVariants(data, false, refID, optionsFromMode(mode), 2) })
		})
	})
}

// gappedAlignment builds an alignment from data: the reference (fuzzRef, "ref") with gap columns (insertions
// relative to it) after the reference positions in the first few bytes, and five queries, each the gapped
// reference changed at the (column, character) pairs of the rest of the bytes, taken in turn
func gappedAlignment(data []byte) []byte {
	const alphabet = "ACGTACGTACGTNRYK-?acgtn"
	n := 0
	if len(data) > 0 {
		n, data = int(data[0]%6), data[1:]
	}
	after := map[int]int{} // reference position (0 = before the first base) -> gap columns after it
	for i := 0; i < n && len(data) > 0; i++ {
		after[int(data[0])%(len(fuzzRef)+1)]++
		data = data[1:]
	}
	var ref strings.Builder
	ref.WriteString(strings.Repeat("-", after[0]))
	for i := 0; i < len(fuzzRef); i++ {
		ref.WriteByte(fuzzRef[i])
		ref.WriteString(strings.Repeat("-", after[i+1]))
	}
	rows := make([][]byte, 5)
	for i := range rows {
		rows[i] = []byte(ref.String())
	}
	for i := 0; i+1 < len(data); i += 2 {
		row := rows[(i/2)%len(rows)]
		row[int(data[i])%len(row)] = alphabet[int(data[i+1])%len(alphabet)]
	}
	var msa bytes.Buffer
	msa.WriteString(">ref\n" + ref.String() + "\n")
	for i, r := range rows {
		msa.WriteString(">q" + strconv.Itoa(i) + "\n" + string(r) + "\n")
	}
	return msa.Bytes()
}

// FuzzVariantsGappedReference: an alignment whose reference row (--reference-id) has gap columns, with queries
// changed anywhere, including in those columns (insertions). It must not panic or hang, and the output must
// not depend on the number of threads, nor on whether the alignment is read from a file or from stdin.
func FuzzVariantsGappedReference(f *testing.F) {
	f.Add([]byte{}, uint8(0))
	f.Add([]byte{2, 10, 30, 12, 5, 40, 7}, uint8(4))
	f.Add([]byte{3, 0, 45, 22, 23, 1, 24, 2, 25, 3, 0, 16}, uint8(25))
	f.Add([]byte{5, 3, 3, 3, 9, 9, 4, 1, 5, 2, 6, 3, 7, 16, 8, 16, 9, 16}, uint8(41))
	f.Fuzz(func(t *testing.T, data []byte, mode uint8) {
		quietStderr(t)
		msa := gappedAlignment(data)
		o := optionsFromMode(mode)
		collect(data, func() {
			runWithTimeout(t, func() {
				one, err1 := runVariants(msa, false, "ref", o, 1)
				four, err4 := runVariants(msa, false, "ref", o, 4)
				fromStdin, errS := runVariants(msa, true, "ref", o, 4)
				if (err1 == nil) != (err4 == nil) || (err1 == nil) != (errS == nil) {
					t.Fatalf("errors differ: 1 thread %v, 4 threads %v, stdin %v\n%s", err1, err4, errS, msa)
				}
				if one != four || one != fromStdin {
					t.Fatalf("output differs:\n1 thread:\n%s\n4 threads:\n%s\nstdin:\n%s\nalignment:\n%s", one, four, fromStdin, msa)
				}
			})
		})
	})
}

// runWithTimeout fails the test if fn does not return soon (a hang). The limit is 10 s, or
// GOFASTA_FUZZ_HANG_TIMEOUT seconds if that is set.
func runWithTimeout(t *testing.T, fn func()) {
	limit := 10 * time.Second
	if s, err := strconv.Atoi(os.Getenv("GOFASTA_FUZZ_HANG_TIMEOUT")); err == nil && s > 0 {
		limit = time.Duration(s) * time.Second
	}
	done := make(chan struct{})
	go func() { fn(); close(done) }()
	select {
	case <-done:
	case <-time.After(limit):
		t.Fatalf("hang: did not return within %v", limit)
	}
}

// quietStderr discards what is written to stderr (warnings about skipped records) until the test ends
func quietStderr(t *testing.T) {
	devnull, err := os.OpenFile(os.DevNull, os.O_WRONLY, 0)
	if err != nil {
		t.Fatal(err)
	}
	old := os.Stderr
	os.Stderr = devnull
	t.Cleanup(func() { os.Stderr = old; devnull.Close() })
}
