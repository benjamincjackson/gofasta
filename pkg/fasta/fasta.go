package fasta

import (
	"fmt"
	"slices"
	"strings"
	"unicode/utf8"

	"github.com/virus-evolution/gofasta/pkg/alphabet"
	"github.com/virus-evolution/gofasta/pkg/encoding"
)

// A struct for one Fasta record
type Record struct {
	ID          string
	Description string
	Seq         string
	Idx         int
}

// RecordName returns the record's ID (see pipeline.Namer)
func (r Record) RecordName() string { return r.ID }

// A struct for one Fasta record whose sequence is encoded using EP's scheme
type EncodedRecord struct {
	ID          string
	Description string
	Seq         []byte
	Idx         int
	Score       int64 // this is for genome completeness
	Count_A     int
	Count_T     int
	Count_G     int
	Count_C     int
}

func (FR Record) encode(hardGaps bool) (EncodedRecord, error) {
	EA := encoding.MakeEncodingArray()
	if hardGaps {
		EA = encoding.MakeEncodingArrayHardGaps()
	}
	seq := make([]byte, len(FR.Seq))
	for i := 0; i < len(FR.Seq); i++ {
		if seq[i] = EA[FR.Seq[i]]; seq[i] == 0 {
			// (every byte of a character that is not ASCII is 0 in the table)
			r, _ := utf8.DecodeRuneInString(FR.Seq[i:])
			return EncodedRecord{}, fmt.Errorf("invalid nucleotide in fasta record (\"%c\")", r)
		}
	}
	return EncodedRecord{ID: FR.ID, Description: FR.Description, Seq: seq, Idx: FR.Idx}, nil
}

// Convert a Record to an EncodedRecord
func (FR Record) Encode() (EncodedRecord, error) {
	return FR.encode(false)
}

// Convert a Record to an EncodedRecord with hard gaps
func (FR Record) EncodeHardGaps() (EncodedRecord, error) {
	return FR.encode(true)
}

// Strip the gaps from a Record's sequence, returning a new Record
func (FR Record) Degap() Record {
	FR.Seq = strings.ReplaceAll(FR.Seq, "-", "")
	return FR
}

// Complement a Record's sequence, returning a new Record
func (FR Record) Complement() Record {
	FR.Seq = alphabet.Complement(FR.Seq)
	return FR
}

// Reverse complement a Record's sequence, returning a new Record
func (FR Record) ReverseComplement() Record {
	FR.Seq = alphabet.ReverseComplement(FR.Seq)
	return FR
}

// Convert an EncodedRecord to a Record, returning a new Record
func (EFR EncodedRecord) Decode() Record {
	return Record{ID: EFR.ID, Description: EFR.Description, Seq: encoding.DecodeToString(EFR.Seq), Idx: EFR.Idx}
}

// Calculate the ATGC content of an EncodedRecord in place
func (EFR *EncodedRecord) CalculateBaseContent() {
	var counting [256]int
	for _, nuc := range EFR.Seq {
		counting[nuc]++
	}
	EFR.Count_A = counting[136]
	EFR.Count_T = counting[24]
	EFR.Count_G = counting[72]
	EFR.Count_C = counting[40]
}

// Score an EncodedRecord for completeness in place
func (EFR *EncodedRecord) CalculateCompleteness() {
	var score int64
	scoring := encoding.MakeEncodedScoreArray()
	for _, nuc := range EFR.Seq {
		score += scoring[nuc]
	}
	EFR.Score = score
}

// Complement an EncodedRecord's sequence, returning a new EncodedRecord (without base content or score)
func (EFR EncodedRecord) Complement() EncodedRecord {
	CA := alphabet.MakeEncodedCompArray()
	seq := make([]byte, len(EFR.Seq))
	for i, nuc := range EFR.Seq {
		seq[i] = CA[nuc]
	}
	return EncodedRecord{ID: EFR.ID, Description: EFR.Description, Seq: seq, Idx: EFR.Idx}
}

// Reverse complement an EncodedRecord's sequence, returning a new EncodedRecord (see Complement)
func (EFR EncodedRecord) ReverseComplement() EncodedRecord {
	NEFR := EFR.Complement()
	slices.Reverse(NEFR.Seq)
	return NEFR
}
