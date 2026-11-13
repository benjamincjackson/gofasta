package variants

import (
	"encoding/binary"
	"sort"
)

// getVariantsNoAnnotation finds the changes in query relative to ref (an aligned pair, encoded) when
// there is no annotation: the same Variants as GetVariantsPair with every position intergenic (nuc:,
// del:, ins:), in one pass over the columns. Runs where the query equals the reference, with no gap in
// the reference, are skipped eight columns at a time, which makes this about as fast as comparing the
// two sequences.
func getVariantsNoAnnotation(ref, query []byte) []Variant {
	const gaps = 0xf4f4f4f4f4f4f4f4 // a gap (244) in every byte
	variants := make([]Variant, 0)
	refPos := 0 // the number of reference bases so far: the reference coordinate of the last one
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
			variants = append(variants, Variant{Changetype: "ins", Position: insAt, Length: insLen})
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
				variants = append(variants, Variant{Changetype: "del", Position: delStart, Length: delLen})
			}
			delLen = 0
		}
		if (q & r) < 16 {
			variants = append(variants, Variant{Changetype: "nuc", RefAl: decodingArray[r], QueAl: decodingArray[q], Position: refPos})
		}
	}
	if insLen > 0 {
		variants = append(variants, Variant{Changetype: "ins", Position: insAt, Length: insLen})
	}
	// the order GetVariantsPair gives: by position, then type
	sort.SliceStable(variants, func(i, j int) bool {
		return variants[i].Position < variants[j].Position || (variants[i].Position == variants[j].Position && variants[i].Changetype < variants[j].Changetype)
	})
	return variants
}

// zeroBytes sets the top bit of each byte of x that is 0
func zeroBytes(x uint64) uint64 {
	const lo7, hi1 = 0x7f7f7f7f7f7f7f7f, 0x8080808080808080
	return ^(((x & lo7) + lo7) | x) & hi1
}
