package updown

import (
	"io"
	"testing"
)

func TestNegativeSizesAndDistancesAreErrors(t *testing.T) {
	if _, _, err := checkArgs(0, -2, 0, 0, 0, 0, 0, 0, 0, 0, io.Discard); err == nil {
		t.Errorf("negative --size-up: expected an error")
	}
	if _, _, err := checkArgs(0, 0, 0, 0, 0, 0, 0, -1, 0, 5, io.Discard); err == nil {
		t.Errorf("negative --dist-down: expected an error")
	}
}

func TestPushWithBothSizesAndDistancesIsAnError(t *testing.T) {
	// checkArgs(sizetotal, sizeup, sizedown, sizeside, sizesame, distall, distup, distdown, distside, distpush)
	tests := []struct {
		name    string
		args    [10]int
		wantErr bool
	}{
		{"--size-total, --dist-all, --dist-push", [10]int{10, 0, 0, 0, 0, 3, 0, 0, 0, 1}, true},
		{"--size-up, --dist-up, --dist-push", [10]int{0, 2, 0, 0, 0, 0, 3, 0, 0, 1}, true},
		{"--size-total, --size-same, --dist-side, --dist-push", [10]int{10, 0, 0, 0, 1, 0, 0, 0, 2, 2}, true},
		{"--dist-up, --dist-push", [10]int{0, 0, 0, 0, 0, 0, 3, 0, 0, 1}, false},
		{"--size-total, --dist-all", [10]int{10, 0, 0, 0, 0, 3, 0, 0, 0, 0}, false},
	}
	for _, tc := range tests {
		a := tc.args
		_, _, err := checkArgs(a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8], a[9], io.Discard)
		if (err != nil) != tc.wantErr {
			t.Errorf("%s: got error %v, want an error: %v", tc.name, err, tc.wantErr)
		}
	}
}
