package updown

import (
	"math/rand"
	"slices"
	"strconv"
	"testing"
)

// randomLine is an updownLine in updown list's order: SNPs (some sharing a position with another
// line's SNP but not its base) and ambiguity tracts, both ascending
func randomLine(rng *rand.Rand, length int) updownLine {
	var pos []int
	var snps []string
	for p := 1; p <= length; p++ {
		if rng.Intn(4) == 0 {
			pos = append(pos, p)
			snps = append(snps, "A"+strconv.Itoa(p)+string("CG"[rng.Intn(2)]))
		}
	}
	var ambs []int
	for p := 1; p <= length; {
		if rng.Intn(6) == 0 {
			end := p + rng.Intn(3)
			ambs = append(ambs, p, end)
			p = end + 2
		} else {
			p++
		}
	}
	sorted := slices.Clone(snps)
	slices.Sort(sorted)
	return updownLine{snps: snps, snpsPos: pos, snpsSorted: sorted, ambs: ambs, ordered: inOrder(pos, ambs)}
}

func TestCountMergeMatchesCountSearch(t *testing.T) {
	rng := rand.New(rand.NewSource(1))
	for i := 0; i < 100000; i++ {
		length := rng.Intn(30)
		q, tg := randomLine(rng, length), randomLine(rng, length)
		if !q.ordered || !tg.ordered {
			t.Fatal("randomLine gave a line that is not in order")
		}
		table1, d1 := countSearch(q, tg)
		table2, d2 := countMerge(q, tg)
		if table1 != table2 || d1 != d2 {
			t.Fatalf("query %+v\ntarget %+v\ncountSearch %v %d, countMerge %v %d", q, tg, table1, d1, table2, d2)
		}
	}
}

func TestInOrder(t *testing.T) {
	tests := []struct {
		pos, ambs []int
		want      bool
	}{
		{nil, nil, true},
		{[]int{1, 2, 9}, []int{3, 4, 6, 6}, true},
		{[]int{1, 1}, nil, false},       // a repeated position
		{[]int{2, 1}, nil, false},       // positions not ascending
		{nil, []int{5, 3}, false},       // a tract that ends before it starts
		{nil, []int{1, 4, 4, 6}, false}, // overlapping tracts
		{nil, []int{6, 7, 1, 2}, false}, // tracts not ascending
		{nil, []int{1, 2, 3}, false},    // an odd number of bounds
	}
	for _, tc := range tests {
		if got := inOrder(tc.pos, tc.ambs); got != tc.want {
			t.Errorf("inOrder(%v, %v) = %v, want %v", tc.pos, tc.ambs, got, tc.want)
		}
	}
}
