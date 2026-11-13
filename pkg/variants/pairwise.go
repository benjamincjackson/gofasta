package variants

import (
	"strconv"
	"strings"
)

func getIndelsPair(ref, query []byte) []Variant {

	var (
		insOpen   bool
		insLength int
		delOpen   bool
		delRefPos int
		delLength int

		// refPos is the number of (non-gap) reference bases seen so far, i.e. the
		// 1-based reference coordinate of the last reference base before this alignment column.
		// An insertion is reported at the reference base before it, a deletion at its first missing base.
		refPos    int
		insRefPos int
	)

	variants := make([]Variant, 0)

	for pos := range ref {
		if ref[pos] == 244 { // insertion relative to reference (somewhere in the alignment)
			if query[pos] == 244 { // insertion is not in this seq
				continue
			} else { // insertion is in this seq
				if insOpen { // not the first position of an insertion
					insLength++ // we increment the length counter
				} else { // the first position of an insertion
					insRefPos = refPos // and the number of reference bases that precede it (its reference coordinate)
					insLength = 1
					insOpen = true
				}
			}
		} else { // not an insertion relative to the reference at this position
			if insOpen { // first base after an insertion, so we need to log the insertion
				variants = append(variants, Variant{Changetype: "ins", Position: insRefPos, Length: insLength})
				insOpen = false
			}
			if query[pos] == 244 { // deletion in this seq
				if delOpen { // not the first position of a deletion
					delLength++ // we increment the length (there is not a deletion in the reference)
				} else { // the first position of a deletion
					delRefPos = refPos
					delLength = 1
					delOpen = true
				}
			} else { // no deletion in this seq
				if delOpen { // first base after a deletion, so we need to log the deletion
					if delRefPos != 0 { // (a deletion at the start of the alignment is not reported)
						variants = append(variants, Variant{Changetype: "del", Position: delRefPos + 1, Length: delLength})
					}
					delOpen = false // and reset things
				}
			}
			refPos++ // this alignment column is a reference base
		}
	}

	// (deletions that abut the end of the alignment are deliberately not reported)

	// catch insertions that abut the end of the alignment
	if insOpen {
		variants = append(variants, Variant{Changetype: "ins", Position: insRefPos, Length: insLength})
	}

	return variants
}

func getNucsPair(ref, query []byte, pos []int, offsetRefCoord []int) []Variant {
	variants := make([]Variant, 0)
	for _, p := range pos {
		alignPos := (p - 1) + offsetRefCoord[p-1]
		if (ref[alignPos] & query[alignPos]) < 16 { // check for SNPs
			variants = append(variants, Variant{Changetype: "nuc", RefAl: decodingArray[ref[alignPos]], QueAl: decodingArray[query[alignPos]], Position: p})
		}
	}
	return variants
}

// getAAsPair reports the amino acid changes in a CDS codon by codon, and the nucleotide changes in codons
// whose amino acid is unchanged, leaving out the codons in skip (0-based codon indices).
func getAAsPair(ref, query []byte, region Region, offsetRefCoord []int, skip map[int]bool) []Variant {

	variants := make([]Variant, 0)
	codonSNPs := make([]Variant, 0, 3) // reused for each codon
	var codon [3]byte                  // the query's codon, as IUPAC characters, in CDS order

	for k, refPos := range region.Positions {
		// here is the actual position in the msa:
		alignmentPos := (refPos - 1) + offsetRefCoord[refPos-1]

		// alignmentPos is the column of a reference base, so the reference is never a gap here: insertions
		// relative to the reference are not visited.
		if (query[alignmentPos] & ref[alignmentPos]) < 16 {
			codonSNPs = append(codonSNPs, Variant{Changetype: "nuc", RefAl: decodingArray[ref[alignmentPos]], QueAl: decodingArray[query[alignmentPos]], Position: refPos})
		}
		codon[k%3] = decodingByte[query[alignmentPos]]

		if k%3 != 2 {
			continue
		}
		aaCounter := k / 3

		if !skip[aaCounter] {
			// on the reverse strand the positions already run 3' to 5' on the forward strand, so the
			// complement (not the reverse complement) gives the codon
			if region.Strand == -1 {
				for i := range codon {
					codon[i] = compArray[codon[i]]
				}
			}
			aa, ok := codonDict[string(codon[:])]
			if !ok {
				aa = "X"
			}
			refaa := region.Translation[aaCounter : aaCounter+1]

			if aa != refaa && aa != "X" {
				temp := make([]string, 0, len(codonSNPs))
				for _, v := range codonSNPs {
					temp = append(temp, "nuc:"+v.RefAl+strconv.Itoa(v.Position)+v.QueAl)
				}
				variants = append(variants, Variant{Changetype: "aa", Feature: region.Name, RefAl: refaa, QueAl: aa, Position: region.Positions[3*aaCounter], Codon: codonPositions(region, aaCounter), Residue: aaCounter + 1, SNPs: strings.Join(temp, ";")})
			} else {
				variants = append(variants, codonSNPs...)
			}
		}
		codonSNPs = codonSNPs[:0]
	}

	return variants
}
