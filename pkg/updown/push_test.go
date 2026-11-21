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
