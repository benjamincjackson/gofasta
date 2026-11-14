package variants

import (
	"encoding/binary"
	"slices"
	"strings"
)

// getVariantsNoAnnotation finds the changes in query relative to ref (an aligned pair, encoded) when
// there is no annotation: the same Variants as GetVariantsPair with every position intergenic (nuc:,
// del:, ins:), in one pass over the columns. Runs where the query equals the reference, with no gap in
// the reference, are skipped eight columns at a time, which makes this about as fast as comparing the
// two sequences.
func getVariantsNoAnnotation(ref, query []byte) []Variant {
	// counted first, so that the slice is made once at its final size (a query far from the reference
	// can have thousands of changes)
	n := 0
	scanNoAnnotation(ref, query, func(Variant) { n++ })
	variants := make([]Variant, 0, n)
	scanNoAnnotation(ref, query, func(v Variant) { variants = append(variants, v) })

	// the order GetVariantsPair gives: by position, then type. The variants are found in that order except
	// that an insertion after a SNP's position is found after the SNP, so they are rarely out of order.
	cmp := func(a, b Variant) int {
		if a.Position != b.Position {
			return a.Position - b.Position
		}
		return strings.Compare(a.Changetype, b.Changetype)
	}
	for i := 1; i < len(variants); i++ {
		if a, b := &variants[i-1], &variants[i]; b.Position < a.Position || (b.Position == a.Position && b.Changetype < a.Changetype) {
			slices.SortStableFunc(variants, cmp)
			break
		}
	}
	return variants
}

// scanNoAnnotation passes each change in query relative to ref to add, in the order found (see
// getVariantsNoAnnotation)
func scanNoAnnotation(ref, query []byte, add func(Variant)) {
	const gaps = 0xf4f4f4f4f4f4f4f4 // a gap (244) in every byte
	refPos := 0                     // the number of reference bases so far: the reference coordinate of the last one
	insLen, insAt := 0, 0
	delLen, delStart := 0, 0
	for i := 0; i < len(query); {
		if insLen == 0 && delLen == 0 && i+8 <= len(query) {
			r, q := binary.LittleEndian.Uint64(ref[i:]), binary.LittleEndian.Uint64(query[i:])
			if r == q && zeroBytes(r^gaps) == 0 {
				refPos += 8
				i += 8
				continue
			}
		}
		r, q := ref[i], query[i]
		i++
		if r == 244 { // a column with a gap in the reference: an insertion, if the query has a base
			if q != 244 {
				if insLen == 0 {
					insAt = refPos
				}
				insLen++
			}
			continue
		}
		if insLen > 0 {
			add(Variant{Changetype: "ins", Position: insAt, Length: insLen})
			insLen = 0
		}
		refPos++
		if q == 244 {
			if delLen == 0 {
				delStart = refPos
			}
			delLen++
			continue
		}
		if delLen > 0 {
			// a deletion at the start of the alignment is not reported (and one at the end is not closed)
			if delStart > 1 {
				add(Variant{Changetype: "del", Position: delStart, Length: delLen})
			}
			delLen = 0
		}
		if (q & r) < 16 {
			add(Variant{Changetype: "nuc", RefAl: decodingArray[r], QueAl: decodingArray[q], Position: refPos})
		}
	}
	if insLen > 0 {
		add(Variant{Changetype: "ins", Position: insAt, Length: insLen})
	}
}

// zeroBytes sets the top bit of each byte of x that is 0
func zeroBytes(x uint64) uint64 {
	const lo7, hi1 = 0x7f7f7f7f7f7f7f7f, 0x8080808080808080
	return ^(((x & lo7) + lo7) | x) & hi1
}
