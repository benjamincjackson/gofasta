package sam

import (
	"bytes"
	"math/rand"
	"strconv"
	"strings"
	"testing"
)

// gapRowColumnByColumn is a direct statement of what gapRowForInsertions does: split the row into its
// reference bases and its own inserted chunks, then visit every reference position in turn, laying out
// the columns of every insertion there (the row's own chunk, or gaps for another row's) before the base.
func gapRowColumnByColumn(refRow, queRow []byte, row int, insertions []insertionOccurence) ([]byte, []byte) {
	var baseRef, baseQue []byte
	ownAt := map[int][][]byte{}
	col, shift := 0, 0
	for _, ins := range insertions {
		if ins.rowNumber != row {
			continue
		}
		insCol := ins.start + shift
		baseRef = append(baseRef, refRow[col:insCol]...)
		baseQue = append(baseQue, queRow[col:insCol]...)
		ownAt[ins.start] = append(ownAt[ins.start], queRow[insCol:insCol+ins.length])
		col = insCol + ins.length
		shift += ins.length
	}
	baseRef = append(baseRef, refRow[col:]...)
	baseQue = append(baseQue, queRow[col:]...)

	var outRef, outQue []byte
	used := map[int]int{}
	for c := 0; c <= len(baseRef); c++ {
		for _, ins := range insertions {
			if ins.start != c {
				continue
			}
			if ins.rowNumber == row {
				chunk := ownAt[c][used[c]]
				used[c]++
				outRef = append(outRef, bytes.Repeat([]byte{'-'}, len(chunk))...)
				outQue = append(outQue, chunk...)
			} else {
				outRef = append(outRef, bytes.Repeat([]byte{'-'}, ins.length)...)
				outQue = append(outQue, bytes.Repeat([]byte{'-'}, ins.length)...)
			}
		}
		if c < len(baseRef) {
			outRef = append(outRef, baseRef[c])
			outQue = append(outQue, baseQue[c])
		}
	}
	return outRef, outQue
}

// randomRecordLine returns a SAM line for query q mapped somewhere on ref, with a random CIGAR of
// matches, insertions, deletions and soft clips (insertions at the same position as other records', and
// two insertions in a row, included)
func randomRecordLine(rng *rand.Rand, ref string, flag int) string {
	pos := rng.Intn(len(ref) / 2)
	var cigar strings.Builder
	var seq strings.Builder
	addOp := func(n int, op byte) {
		cigar.WriteString(strconv.Itoa(n))
		cigar.WriteByte(op)
	}
	bases := func(n int) {
		for range n {
			seq.WriteByte("ACGT"[rng.Intn(4)])
		}
	}
	if rng.Intn(3) == 0 {
		n := 1 + rng.Intn(3)
		addOp(n, 'S')
		bases(n)
	}
	r := pos
	n := 1 + rng.Intn(4)
	addOp(n, 'M')
	bases(n)
	r += n
	for r < len(ref) && rng.Intn(6) != 0 {
		switch k := rng.Intn(5); {
		case k == 0:
			n = 1 + rng.Intn(3)
			addOp(n, 'I')
			bases(n)
			if rng.Intn(4) == 0 {
				m := 1 + rng.Intn(2)
				addOp(m, 'I')
				bases(m)
			}
		case k == 1 && r+3 < len(ref):
			n = 1 + rng.Intn(3)
			addOp(n, 'D')
			r += n
		}
		n = min(1+rng.Intn(6), len(ref)-r)
		if n == 0 {
			break
		}
		addOp(n, 'M')
		bases(n)
		r += n
	}
	if rng.Intn(3) == 0 {
		n := 1 + rng.Intn(3)
		addOp(n, 'S')
		bases(n)
	}
	return rec("q", flag, pos+1, cigar.String(), seq.String())
}

// gapRowForInsertions lays out every row of random blocks of one to four records exactly as the
// column-by-column statement of the same rule does, and a block of one record needs no rebuilding: its
// row is already in its final columns.
func TestGapRowForInsertionsMatchesColumnByColumn(t *testing.T) {
	rng := rand.New(rand.NewSource(1))
	ref := overlapTestRef()
	for trial := 0; trial < 3000; trial++ {
		lines := []string{randomRecordLine(rng, ref, 0)}
		for k := rng.Intn(4); k > 0; k-- {
			lines = append(lines, randomRecordLine(rng, ref, 2048))
		}
		_, recs := readTestSam(t, samLines(ref, lines...))
		insertions := blockInsertions(recs)
		for j, r := range recs {
			queRow, refRow, err := alignRecord(r, []byte(ref), len(ref), true)
			if err != nil {
				t.Fatal(err)
			}
			gotRef, gotQue := gapRowForInsertions(refRow, queRow, j, insertions)
			wantRef, wantQue := gapRowColumnByColumn(refRow, queRow, j, insertions)
			if string(gotRef) != string(wantRef) || string(gotQue) != string(wantQue) {
				t.Fatalf("block %q, row %d:\ngot  %s\n     %s\nwant %s\n     %s", lines, j, gotRef, gotQue, wantRef, wantQue)
			}
			if len(recs) == 1 && (string(gotRef) != string(refRow) || string(gotQue) != string(queRow)) {
				t.Fatalf("one record %q: rebuilding changed the row\nrow  %s\n     %s\ngot  %s\n     %s", lines, refRow, queRow, gotRef, gotQue)
			}
		}
	}
}
