package variants

import (
	"github.com/virus-evolution/gofasta/pkg/alphabet"
	"github.com/virus-evolution/gofasta/pkg/encoding"
)

// Lookup tables, built once and only read (so they are safe to share between goroutines)
var (
	decodingArray = encoding.MakeDecodingArray() // encoded base -> IUPAC character, as a string
	decodingByte  = makeDecodingBytes()          // encoded base -> IUPAC character
	compArray     = alphabet.MakeCompArray()     // IUPAC character -> its complement
	codonDict     = alphabet.MakeCodonDict()     // codon (IUPAC characters) -> amino acid
)

func makeDecodingBytes() [256]byte {
	var b [256]byte
	for i, s := range decodingArray {
		if len(s) > 0 {
			b[i] = s[0]
		}
	}
	return b
}
