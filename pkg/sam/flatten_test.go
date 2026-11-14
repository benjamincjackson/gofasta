package sam

import "testing"

// agreedNuc gives getNucFromSite's answer for every column of one to four rows over '*', '-', two bases,
// N and a byte below '*' whose rows agree (apart from '*'), and leaves every other column to it.
func TestAgreedNucMatchesGetNucFromSite(t *testing.T) {
	alphabet := []byte{'*', '-', 'A', 'C', 'N', '!'}
	for rows := 1; rows <= 4; rows++ {
		n := 1
		for range rows {
			n *= len(alphabet)
		}
		for k := 0; k < n; k++ {
			block := make([][]byte, rows)
			site := make([]byte, rows)
			x := k
			for i := range block {
				site[i] = alphabet[x%len(alphabet)]
				block[i] = []byte{site[i]}
				x /= len(alphabet)
			}
			// the rows agree if their bytes other than '*' are all the same
			agreed := true
			var first byte = '*'
			for _, b := range site {
				if b != '*' {
					if first == '*' {
						first = b
					} else if b != first {
						agreed = false
					}
				}
			}
			got, ok := agreedNuc(block, 0)
			if ok != agreed {
				t.Errorf("column %q: agreedNuc ok = %v, want %v", site, ok, agreed)
				continue
			}
			if ok {
				if want := getNucFromSite(site, "q", 0); got != want {
					t.Errorf("column %q: agreedNuc gives %q, getNucFromSite %q", site, got, want)
				}
			}
		}
	}
}
