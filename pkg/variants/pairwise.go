package variants

import (
	"encoding/binary"
	"strconv"
	"strings"

	"github.com/virus-evolution/gofasta/pkg/alphabet"
)

// getIndelsPair finds the insertions and deletions in query relative to ref (an aligned pair, encoded). An
// insertion is reported at the reference base before it, a deletion at its first missing base; deletions
// at the start or end of the alignment are not reported.
func getIndelsPair(ref, query []byte) []Variant {
	var (
		// the number of reference bases seen so far: the 1-based reference coordinate of the last
		// reference base before the current alignment column
		refPos               int
		insOpen, delOpen     bool
		insRefPos, insLength int
		delRefPos, delLength int
	)
	variants := make([]Variant, 0)

	const gaps = 0xf4f4f4f4f4f4f4f4 // a gap (244) in every byte
	for pos := 0; pos < len(ref); pos++ {
		// skip eight columns at a time where neither sequence has a gap and no indel is open
		for !insOpen && !delOpen && pos+8 <= len(ref) &&
			zeroBytes(binary.LittleEndian.Uint64(ref[pos:])^gaps) == 0 && zeroBytes(binary.LittleEndian.Uint64(query[pos:])^gaps) == 0 {
			refPos += 8
			pos += 8
		}
		if pos == len(ref) {
			break
		}
		if ref[pos] == 244 { // a gap in the reference: an insertion, if the query has a base here
			if query[pos] != 244 {
				if !insOpen {
					insRefPos, insLength, insOpen = refPos, 0, true
				}
				insLength++
			}
			continue
		}
		if insOpen {
			variants = append(variants, Variant{Changetype: "ins", Position: insRefPos, Length: insLength})
			insOpen = false
		}
		if query[pos] == 244 { // a deletion
			if !delOpen {
				delRefPos, delLength, delOpen = refPos, 0, true
			}
			delLength++
		} else if delOpen {
			if delRefPos != 0 {
				variants = append(variants, Variant{Changetype: "del", Position: delRefPos + 1, Length: delLength})
			}
			delOpen = false
		}
		refPos++
	}

	// an insertion at the end of the alignment is reported (a deletion is not)
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
			aa := alphabet.TranslateCodon(codon[0], codon[1], codon[2]) // (0 if it has none)
			refaa := region.Translation[aaCounter : aaCounter+1]

			// (no change from an unknown reference residue: an X in /translation, or an ambiguous codon)
			if aa != 0 && aa != refaa[0] && refaa[0] != 'X' {
				temp := make([]string, 0, len(codonSNPs))
				for _, v := range codonSNPs {
					temp = append(temp, "nuc:"+v.RefAl+strconv.Itoa(v.Position)+v.QueAl)
				}
				variants = append(variants, Variant{Changetype: "aa", Feature: region.Name, RefAl: refaa, QueAl: string(rune(aa)), Position: region.Positions[3*aaCounter], Codon: codonPositions(region, aaCounter), Residue: aaCounter + 1, SNPs: strings.Join(temp, ";")})
			} else {
				variants = append(variants, codonSNPs...)
			}
		}
		codonSNPs = codonSNPs[:0]
	}

	return variants
}
