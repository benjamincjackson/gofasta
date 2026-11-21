package updown

import (
	"bytes"
	"strings"
	"testing"
)

// pushOptions are the --dist-* and --dist-push values of one TopRanking call
type pushOptions struct{ distAll, distUp, distDown, distSide, distPush int }

// topRankingCSV runs TopRanking with o on query and target lines in updown list's csv format (without
// the header), and returns its output
func topRankingCSV(t *testing.T, query, target string, o pushOptions) string {
	t.Helper()
	var out bytes.Buffer
	err := TopRanking(strings.NewReader(listHeader+"\n"+query), strings.NewReader(listHeader+"\n"+target), nil, &out,
		false, "csv", "csv", nil,
		0, 0, 0, 0, 0,
		o.distAll, o.distUp, o.distDown, o.distSide,
		0.1, 10000, false, o.distPush, 1, true)
	if err != nil {
		t.Fatal(err)
	}
	return out.String()
}

func TestPushSortsTheSameBinByAmbiguities(t *testing.T) {
	query := "Q,A1C,,1,0\n"
	// all on a polytomy with the query, with ambiguities away from its SNP
	target := "T5,A1C,10-14,1,5\nT1,A1C,20,1,1\nT3,A1C,30-32,1,3\n"
	got := topRankingCSV(t, query, target, pushOptions{distPush: 1})
	want := "query,closestsame,closestup,closestdown,closestside\nQ,T1;T3;T5,,,\n"
	if got != want {
		t.Errorf("got\n%swant\n%s", got, want)
	}
}

func TestPushFromTheDistBoundaries(t *testing.T) {
	query := "Q,A1C|A2C|A3C,,3,0\n"
	target := strings.Join([]string{
		// up: the query's SNPs are a superset, at distances 1, 2 and 3
		"U1,A1C|A2C,,2,0", "U2,A1C,,1,0", "U3,,,0,0",
		// down: the target's are a superset, at 1, 2 and 3
		"D1,A1C|A2C|A3C|A4C,,4,0", "D2,A1C|A2C|A3C|A4C|A5C,,5,0", "D3,A1C|A2C|A3C|A4C|A5C|A6C,,6,0",
		// side: both have SNPs of their own, at 2 and 3
		"S2,A1C|A2C|A7C,,3,0", "S3,A1C|A7C,,2,0",
	}, "\n") + "\n"
	header := "query,closestsame,closestup,closestdown,closestside\n"
	tests := []struct {
		name string
		o    pushOptions
		want string
	}{
		// no --dist-*: every bin from 0, so the closest distance(s)
		{"no dist", pushOptions{distPush: 1}, "Q,,U1,D1,S2\n"},
		{"no dist, push 2", pushOptions{distPush: 2}, "Q,,U1;U2,D1;D2,S2;S3\n"},
		// everything within the boundary, which already covers a distance
		{"dist-all 2", pushOptions{distAll: 2, distPush: 1}, "Q,,U1;U2,D1;D2,S2\n"},
		// pushed beyond the boundary to cover two distances
		{"dist-all 1, push 2", pushOptions{distAll: 1, distPush: 2}, "Q,,U1;U2,D1;D2,S2;S3\n"},
		// only the bins with a --dist-* value are pushed
		{"dist-up 3", pushOptions{distUp: 3, distPush: 1}, "Q,,U1;U2;U3,,\n"},
		{"dist-up 1, dist-side 2, push 3", pushOptions{distUp: 1, distSide: 2, distPush: 3}, "Q,,U1;U2;U3,,S2;S3\n"},
	}
	for _, tc := range tests {
		if got := topRankingCSV(t, query, target, tc.o); got != header+tc.want {
			t.Errorf("%s: got\n%swant\n%s%s", tc.name, got, header, tc.want)
		}
	}
}
