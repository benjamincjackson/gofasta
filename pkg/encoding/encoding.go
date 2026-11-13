/*
Package encoding provides mappings between character representations of nucleotides
and bitwise encodings, for fast sequence comparison.

The coding scheme is Emmanual Paradis' design, which is described here:
http://ape-package.ird.fr/misc/BitLevelCodingScheme.html
*/
package encoding

import (
	"math/bits"
	"strings"
	"unicode"
)

func DecodeToString(bs []byte) string {
	da := MakeDecodingArray()
	var s strings.Builder
	s.Grow(len(bs))
	for _, n := range bs {
		s.WriteString(da[n])
	}
	return s.String()
}

// codes gives each IUPAC nucleotide code (and '-', '?') with its Emmanuel Paradis encoding: bits 7-4
// are A, G, C and T, set for each base the code can stand for; bit 3 marks a code known to be one base;
// a gap ('-') is 244 (unknown, rather than EP's original 4, which buildEncodingArrayHardGaps uses)
var codes = []struct {
	char byte
	code byte
}{
	{'A', 136}, {'G', 72}, {'C', 40}, {'T', 24},
	{'R', 192}, {'M', 160}, {'W', 144}, {'S', 96}, {'K', 80}, {'Y', 48},
	{'V', 224}, {'H', 176}, {'D', 208}, {'B', 112}, {'N', 240},
	{'-', 244}, {'?', 242},
}

// hardGap is EP's original code for a gap: known, and none of A, C, G or T
const hardGap = 4

// buildEncodingArray returns an array whose indices are the byte representations
// of IUPAC codes and whose contents are Emmanual Paradis encodings.
// Lower case nucleotides are mapped to their upper case nucleotides's encoding
func buildEncodingArray() [256]byte {
	var byteArray [256]byte
	for _, c := range codes {
		byteArray[c.char] = c.code
		byteArray[unicode.ToLower(rune(c.char))] = c.code
	}
	return byteArray
}

var encodingArrayTable = buildEncodingArray()

// MakeEncodingArray returns a copy of the table buildEncodingArray builds, which is built only once
func MakeEncodingArray() [256]byte { return encodingArrayTable }

// buildEncodingArrayHardGaps is as MakeEncodingArray but with '-'
// set to EP's original code (4; character is not completely unknown,
// character represents NONE of ACTG)
func buildEncodingArrayHardGaps() [256]byte {
	byteArray := buildEncodingArray()
	byteArray['-'] = hardGap
	return byteArray
}

var encodingArrayHardGapsTable = buildEncodingArrayHardGaps()

// MakeEncodingArrayHardGaps returns a copy of the table buildEncodingArrayHardGaps builds, which is built only once
func MakeEncodingArrayHardGaps() [256]byte { return encodingArrayHardGapsTable }

// buildScoreArray makes an array that maps (the byte value of) IUPAC nucleotide
// characters (which are the arrays indices) to a score for how ambiguous that
// nucleotide is. The score is caculated as:
// 12 * 1/possible real nucleotides.
// E.g. an 'A' scores 12, but an 'N' (A, C, G or T) scores 3
func buildScoreArray() [256]int64 {
	var byteArray [256]int64
	for c, code := range encodingArrayTable {
		byteArray[c] = encodedScoreArrayTable[code]
	}
	return byteArray
}

var scoreArrayTable = buildScoreArray()

// MakeScoreArray returns a copy of the table buildScoreArray builds, which is built only once
func MakeScoreArray() [256]int64 { return scoreArrayTable }

// buildEncodedScoreArray is like MakeScoreArray() but the indices are EP bitwise representations of
// nucleotides (a gap and '?' count as N; a hard gap has no score)
func buildEncodedScoreArray() [256]int64 {
	var byteArray [256]int64
	for _, c := range codes {
		byteArray[c.code] = 12 / int64(bits.OnesCount8(c.code>>4))
	}
	return byteArray
}

var encodedScoreArrayTable = buildEncodedScoreArray()

// MakeEncodedScoreArray returns a copy of the table buildEncodedScoreArray builds, which is built only once
func MakeEncodedScoreArray() [256]int64 { return encodedScoreArrayTable }

// buildDecodingArray returns an array whose indices are Emmanual Paradis encodings
// of IUPAC codes and whose contents are IUPAC codes as strings (both gap codes decode to '-')
func buildDecodingArray() [256]string {
	var byteArray [256]string
	for _, c := range codes {
		byteArray[c.code] = string(c.char)
	}
	byteArray[hardGap] = "-"
	return byteArray
}

var decodingArrayTable = buildDecodingArray()

// MakeDecodingArray returns a copy of the table buildDecodingArray builds, which is built only once
func MakeDecodingArray() [256]string { return decodingArrayTable }
