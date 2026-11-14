package sam

import (
	"bytes"
	"errors"
)

// alignRecord builds the aligned query row of a SAM record from its CIGAR, in reference coordinates:
// reference positions before the record and (if insertions are not kept) after it are '*', deletions
// are '-', and skipped regions (N) are '*'. If reference is not nil it also builds the aligned
// reference row, which starts with the reference bases before the record and has '-' opposite kept
// insertions; refLen is then len(reference). Insertions relative to the reference are kept if
// includeInsertions, and otherwise left out (so the query row has refLen columns).
func alignRecord(rec samRecord, reference []byte, refLen int, includeInsertions bool) ([]byte, []byte, error) {
	if rec.Pos < 0 {
		return nil, nil, errors.New("unmapped read")
	}
	withRef := reference != nil
	seq := rec.Seq

	query := bytes.Repeat([]byte{'*'}, rec.Pos)
	var ref []byte
	if withRef {
		ref = append(make([]byte, 0, len(reference)), reference[:rec.Pos]...)
	}

	q, r := 0, rec.Pos // the next query and reference positions
	for _, op := range rec.Cigar {
		n := op.Len()
		switch op.Type() {
		case cigarMatch, cigarEqual, cigarMismatch:
			query = append(query, seq[q:q+n]...)
			if withRef {
				ref = append(ref, reference[r:r+n]...)
			}
			q, r = q+n, r+n
		case cigarInsertion:
			if includeInsertions {
				query = append(query, seq[q:q+n]...)
				if withRef {
					ref = append(ref, bytes.Repeat([]byte{'-'}, n)...)
				}
			}
			q += n
		case cigarDeletion, cigarSkipped:
			fill := byte('-')
			if op.Type() == cigarSkipped {
				fill = '*'
			}
			query = append(query, bytes.Repeat([]byte{fill}, n)...)
			if withRef {
				ref = append(ref, reference[r:r+n]...)
			}
			r += n
		case cigarSoftClipped:
			q += n
		case cigarHardClipped, cigarPadded:
		default:
			// (validateRecord has already rejected these)
			return nil, nil, errors.New("unsupported CIGAR operation " + op.Type().String())
		}
	}

	if !includeInsertions {
		query = append(query, bytes.Repeat([]byte{'*'}, refLen-len(query))...)
	}
	return query, ref, nil
}
