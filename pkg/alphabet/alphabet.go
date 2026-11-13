// Package alphabet provides mappings between codons and
// amino acids
package alphabet

import (
	"errors"
	"slices"
	"strings"
	"unicode"

	"github.com/virus-evolution/gofasta/pkg/encoding"
)

var ErrorCDSNotModThree error = errors.New("Translation error: Nucleotide sequence not divisible by 3")

// codonDict is MakeCodonDict's table, built once and only read
var codonDict = MakeCodonDict()

// iupacIndex gives each upper-case IUPAC nucleotide code an index from 1 to 15; any other byte is 0
var iupacIndex = func() (idx [256]uint8) {
	for i, c := range []byte("ACGTRYSWKMBDHVN") {
		idx[c] = uint8(i + 1)
	}
	return idx
}()

// codonTable is codonDict indexed by codonKey (4 KiB): the amino acid of each codon, or 0 if it has none
var codonTable = func() (t [4096]byte) {
	for codon, aa := range codonDict {
		t[codonKey(codon[0], codon[1], codon[2])] = aa[0]
	}
	return t
}()

func codonKey(a, b, c byte) int {
	return int(iupacIndex[a])<<8 | int(iupacIndex[b])<<4 | int(iupacIndex[c])
}

// TranslateCodon returns the amino acid of a codon of upper-case IUPAC nucleotide codes, as in
// MakeCodonDict, or 0 if it has none (a codon with any other character, or with ambiguity codes that
// stand for more than one amino acid). It is a lookup in a small table, without the map's hashing.
func TranslateCodon(a, b, c byte) byte {
	return codonTable[codonKey(a, b, c)]
}

// Translate a nucleotide sequence (in upper or lower case) to a protein sequence
// Codons with ambiguous nucleotides are resolved if it can only possibly
// represent one amino acid
func Translate(nuc string, strict bool) (string, error) {
	if len(nuc)%3 != 0 {
		return "", ErrorCDSNotModThree
	}
	nuc = strings.ToUpper(nuc) // the codon table is upper case
	var translation strings.Builder
	translation.Grow(len(nuc) / 3)
	for i := 0; i+3 <= len(nuc); i += 3 {
		if aa := TranslateCodon(nuc[i], nuc[i+1], nuc[i+2]); aa != 0 {
			translation.WriteByte(aa)
		} else if strict {
			return "", errors.New("Translation error: Untranslatable codon: " + nuc[i:i+3])
		} else {
			translation.WriteByte('X')
		}
	}

	return translation.String(), nil
}

// Complement returns the complement of a nucleotide sequence (IUPAC codes, gaps and '?', in upper or
// lower case; any other character becomes 0)
func Complement(nuc string) string {
	ba := make([]byte, len(nuc))
	for i := 0; i < len(nuc); i++ {
		ba[i] = compArrayTable[nuc[i]]
	}
	return string(ba)
}

// ReverseComplement returns the reverse complement of a nucleotide sequence (see Complement)
func ReverseComplement(nuc string) string {
	ba := []byte(Complement(nuc))
	slices.Reverse(ba)
	return string(ba)
}

// complementPairs lists each IUPAC code (and '?' and '-') followed by its complement
const complementPairs = "ATTAGCCGRYYRSSWWKMMKBVVBDHHDNN??--"

func buildCompArray() [256]byte {
	var compArray [256]byte
	for i := 0; i < len(complementPairs); i += 2 {
		c, comp := complementPairs[i], complementPairs[i+1]
		compArray[c] = comp
		compArray[unicode.ToLower(rune(c))] = byte(unicode.ToLower(rune(comp)))
	}
	return compArray
}

var compArrayTable = buildCompArray()

// MakeCompArray returns a copy of the table buildCompArray builds, which is built only once
func MakeCompArray() [256]byte { return compArrayTable }

// buildEncodedCompArray is buildCompArray's table for nucleotides encoded by encoding.MakeEncodingArray
func buildEncodedCompArray() [256]byte {
	var compArray [256]byte
	EA := encoding.MakeEncodingArray()
	for c, comp := range compArrayTable {
		if comp != 0 {
			compArray[EA[c]] = EA[comp]
		}
	}
	return compArray
}

var encodedCompArrayTable = buildEncodedCompArray()

// MakeEncodedCompArray returns a copy of the table buildEncodedCompArray builds, which is built only once
func MakeEncodedCompArray() [256]byte { return encodedCompArrayTable }

// standardCode gives the amino acid of each codon of the standard genetic code: codon i has the bases
// codeBases[i/16], codeBases[i/4%4] and codeBases[i%4]
const (
	codeBases    = "TCAG"
	standardCode = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
)

// iupacBases gives the bases each IUPAC nucleotide code stands for
var iupacBases = map[byte]string{
	'A': "A", 'C': "C", 'G': "G", 'T': "T",
	'R': "AG", 'Y': "CT", 'S': "CG", 'W': "AT", 'K': "GT", 'M': "AC",
	'B': "CGT", 'D': "AGT", 'H': "ACT", 'V': "ACG", 'N': "ACGT",
}

// MakeCodonDict returns a map from codon (string, upper case) to amino acid code (string) for the
// standard genetic code. A codon with ambiguity codes is in it if every codon it can stand for gives
// the same amino acid.
func MakeCodonDict() map[string]string {
	plain := make(map[string]byte, 64)
	for i := range standardCode {
		plain[string([]byte{codeBases[i/16], codeBases[i/4%4], codeBases[i%4]})] = standardCode[i]
	}
	codonAA := make(map[string]string)
	for c1, b1 := range iupacBases {
		for c2, b2 := range iupacBases {
			for c3, b3 := range iupacBases {
				aa := byte(0)
				unique := true
				for _, x := range []byte(b1) {
					for _, y := range []byte(b2) {
						for _, z := range []byte(b3) {
							a := plain[string([]byte{x, y, z})]
							unique = unique && (aa == 0 || a == aa)
							aa = a
						}
					}
				}
				if unique {
					codonAA[string([]byte{c1, c2, c3})] = string(aa)
				}
			}
		}
	}
	return codonAA
}
