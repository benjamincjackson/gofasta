package variants

import (
	"reflect"
	"testing"
)

// CDSs with the same start are kept in file order, whether or not they have an ID
func TestRegionsFromGFFTiesKeepFileOrder(t *testing.T) {
	rows := []string{
		"ref1\ttest\tCDS\t1\t12\t.\t+\t0\tID=c1;Name=b",
		"ref1\ttest\tCDS\t1\t9\t.\t+\t0\tName=c",
		"ref1\ttest\tCDS\t1\t6\t.\t+\t0\tID=c2;Name=a",
		"ref1\ttest\tCDS\t1\t3\t.\t+\t0\tID=c3;Name=d",
	}
	want := []string{"b", "c", "a", "d"}
	// an order that depends on map iteration would only show up by chance, so check it many times
	for i := 0; i < 50; i++ {
		cds, _, err := regionsFromGFFRows(t, rows...)
		if err != nil {
			t.Fatal(err)
		}
		got := make([]string, len(cds))
		for j, r := range cds {
			got[j] = r.Name
		}
		if !reflect.DeepEqual(got, want) {
			t.Fatalf("got CDS order %v, want %v", got, want)
		}
	}
}
