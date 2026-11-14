package sam

import (
	"bytes"
	"os"
	"slices"
	"strconv"

	"github.com/virus-evolution/gofasta/pkg/encoding"
	"github.com/virus-evolution/gofasta/pkg/samfile"
)

// referenceView returns, for each reference position (0-based) from r.Pos onwards that r's alignment spans,
// the base that r aligns there, '-' for a deletion, and 0 for a skipped region.
func referenceView(r samfile.Record) []byte {
	seq := r.Seq
	view := make([]byte, 0)
	q := 0
	for _, op := range r.Cigar {
		n := op.Len()
		switch op.Type() {
		case samfile.CigarMatch, samfile.CigarEqual, samfile.CigarMismatch:
			view = append(view, seq[q:q+n]...)
			q += n
		case samfile.CigarInsertion, samfile.CigarSoftClipped:
			q += n
		case samfile.CigarDeletion:
			view = append(view, bytes.Repeat([]byte{'-'}, n)...)
		case samfile.CigarSkipped:
			view = append(view, make([]byte, n)...)
		}
	}
	return view
}

// interval is a half-open stretch [start, end) of reference positions (0-based)
type interval struct{ start, end int }

// recordSpan is the stretch of the reference that r's alignment spans
func recordSpan(r samfile.Record) interval {
	ref, _ := samfile.CigarLengths(r.Cigar)
	return interval{r.Pos, r.Pos + ref}
}

func anyContains(ivs []interval, c int) bool {
	for _, iv := range ivs {
		if iv.start <= c && c < iv.end {
			return true
		}
	}
	return false
}

// anyStrictlyInside reports whether the insertion that follows c reference bases (0-based boundary c, i.e. it
// sits between reference positions c-1 and c) lies in the interior of one of the intervals. An insertion at
// the very start or end of an interval is not inside it: the record that spans the interval cannot carry it.
func anyStrictlyInside(ivs []interval, c int) bool {
	for _, iv := range ivs {
		if iv.start < c && c < iv.end {
			return true
		}
	}
	return false
}

// resolveOverlaps makes the records that belong to one query (a primary alignment and its supplementary
// alignments) agree about who describes which stretch of the reference. The records are taken in priority
// order: the primary record first, then the supplementary records in the order they were read. Where
// records overlap on the reference, the alignment of the record that comes first in that order is used
// for the overlapping stretch and the later records only contribute outside it. (Overlaps come mostly from
// repeats at the genome's ends, such as poxvirus ITRs or retroviral LTRs: an aligner can place the
// query's copy at one end on the reference's copy at the other end too, and the primary alignment is
// the one in genome order.)
//
// Records are merged column by column along the reference. That works for substitutions and deletions, but an
// insertion has no reference column of its own, so nothing says that two records' insertions are the same
// bases of the query. So that an insertion carried by two overlapping records (or the two equally valid
// placements of an insertion in a repeat) is reported once, only the first record's alignment is used
// wherever records overlap.
//
// A later record is rewritten so that the stretches the earlier records span become skipped regions (and the
// query bases it aligned there, including any insertion inside them, are removed from its SEQ). A record
// that is left with nothing to contribute is dropped. Disagreements between the records in the overlap are
// reported on stderr.
func resolveOverlaps(records []samfile.Record) []samfile.Record {
	if len(records) < 2 {
		return records
	}
	if i := slices.IndexFunc(records, func(r samfile.Record) bool { return r.Flags&samfile.FlagSupplementary == 0 }); i > 0 {
		records = slices.Concat(records[i:i+1], records[:i], records[i+1:])
	}
	out := []samfile.Record{records[0]}
	covered := []interval{recordSpan(records[0])}

	for k := 1; k < len(records); k++ {
		r := records[k]
		span := recordSpan(r)
		overlaps := false
		for _, iv := range covered {
			if iv.start < span.end && span.start < iv.end {
				overlaps = true
				break
			}
		}
		if !overlaps {
			out = append(out, r)
			covered = append(covered, span)
			continue
		}
		if n := countDisagreements(records[:k], r); n > 0 {
			os.Stderr.WriteString("overlapping records disagree: " + r.Name + ": " + strconv.Itoa(n) +
				" reference positions, using the alignment of the first record at these\n")
		}
		if trimmed, ok := trimToUncovered(r, covered); ok {
			out = append(out, trimmed)
		}
		covered = append(covered, span)
	}
	return out
}

// countDisagreements counts the reference positions in the overlap of r with the earlier records where r and
// the first earlier record to align that position certainly differ (different bases, or a base against a deletion).
func countDisagreements(earlier []samfile.Record, r samfile.Record) int {
	enc := encoding.MakeEncodingArray()
	type view struct {
		start int
		bases []byte
	}
	views := make([]view, len(earlier))
	for i, e := range earlier {
		views[i] = view{e.Pos, referenceView(e)}
	}
	count := 0
	for i, b := range referenceView(r) {
		if b == 0 {
			continue
		}
		c := r.Pos + i
		for _, v := range views {
			j := c - v.start
			if j < 0 || j >= len(v.bases) || v.bases[j] == 0 {
				continue
			}
			a := v.bases[j]
			switch {
			case a == '-' && b == '-':
			case a == '-' || b == '-':
				count++
			case (enc[a] & enc[b]) < 16:
				count++
			}
			break
		}
	}
	return count
}

// trimToUncovered rewrites r so that it no longer aligns anything in the covered intervals. It returns false
// if nothing is left of it.
func trimToUncovered(r samfile.Record, covered []interval) (samfile.Record, bool) {
	type op struct {
		t samfile.CigarType
		n int
	}
	seq := r.Seq
	ops := make([]op, 0, len(r.Cigar))
	newSeq := make([]byte, 0, len(seq))
	push := func(t samfile.CigarType, n int) {
		if n == 0 {
			return
		}
		if len(ops) > 0 && ops[len(ops)-1].t == t {
			ops[len(ops)-1].n += n
			return
		}
		ops = append(ops, op{t, n})
	}

	ref, q := r.Pos, 0
	for _, o := range r.Cigar {
		n, t := o.Len(), o.Type()
		switch t {
		case samfile.CigarMatch, samfile.CigarEqual, samfile.CigarMismatch, samfile.CigarDeletion:
			hasQuery := t != samfile.CigarDeletion
			for i := 0; i < n; {
				cov := anyContains(covered, ref+i)
				j := i
				for j < n && anyContains(covered, ref+j) == cov {
					j++
				}
				if cov {
					push(samfile.CigarSkipped, j-i)
				} else {
					push(t, j-i)
					if hasQuery {
						newSeq = append(newSeq, seq[q+i:q+j]...)
					}
				}
				i = j
			}
			ref += n
			if hasQuery {
				q += n
			}
		case samfile.CigarSkipped:
			push(t, n)
			ref += n
		case samfile.CigarInsertion:
			if !anyStrictlyInside(covered, ref) {
				push(t, n)
				newSeq = append(newSeq, seq[q:q+n]...)
			}
			q += n
		case samfile.CigarSoftClipped:
			push(t, n)
			newSeq = append(newSeq, seq[q:q+n]...)
			q += n
		default: // hard clipping, padding
			push(t, n)
		}
	}

	clip := func(t samfile.CigarType) bool {
		return t == samfile.CigarHardClipped || t == samfile.CigarSoftClipped || t == samfile.CigarPadded
	}
	// a skipped region at either end is just a shorter alignment
	pos := r.Pos
	for i := 0; i < len(ops); i++ {
		if clip(ops[i].t) {
			continue
		}
		if ops[i].t == samfile.CigarSkipped {
			pos += ops[i].n
			ops = append(ops[:i], ops[i+1:]...)
		}
		break
	}
	for i := len(ops) - 1; i >= 0; i-- {
		if clip(ops[i].t) {
			continue
		}
		if ops[i].t == samfile.CigarSkipped {
			ops = append(ops[:i], ops[i+1:]...)
		}
		break
	}

	aligned := false
	for _, o := range ops {
		if o.t == samfile.CigarMatch || o.t == samfile.CigarEqual || o.t == samfile.CigarMismatch {
			aligned = true
		}
	}
	if !aligned {
		return r, false
	}

	cigar := make([]samfile.CigarOp, len(ops))
	for i, o := range ops {
		cigar[i] = samfile.NewCigarOp(o.t, o.n)
	}
	r.Pos = pos
	r.Cigar = cigar
	r.Seq = newSeq
	return r, true
}
