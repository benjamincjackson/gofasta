package variants

import (
	"sort"
	"strconv"
	"strings"
	"sync"

	"github.com/virus-evolution/gofasta/pkg/alphabet"
)

// cdsSlot is one reference base of a CDS, in CDS order, with what the query has there
type cdsSlot struct {
	refPos   int    // 1-based reference position
	refBase  string // the reference base
	queBase  string // the query base, "" if it is deleted or missing
	missing  bool   // the query has no data here: N, or a gap that runs to the start or end of the alignment
	snp      bool   // the query base certainly differs from the reference base
	inserted string // query bases inserted between this base and the next one of the CDS, in CDS order
}

// cdsSlots lines up the query with each base of a CDS. Insertions are recorded only between two bases
// of the CDS that are next to each other in the reference.
// Every slot is overwritten, so buf (whose length is ignored) may hold slots from an earlier call; it is
// grown if it is too small.
func cdsSlots(ref, query []byte, region Region, offsetRefCoord []int, buf []cdsSlot) []cdsSlot {
	// the query's first and last bases: gaps outside them are missing data, not deletions
	first, last := -1, -1
	for i, q := range query {
		if q != 244 {
			if first < 0 {
				first = i
			}
			last = i
		}
	}

	col := func(p int) int { return (p - 1) + offsetRefCoord[p-1] }
	slots := buf[:0]
	if cap(slots) < len(region.Positions) {
		slots = make([]cdsSlot, len(region.Positions))
	}
	slots = slots[:len(region.Positions)]
	for k, p := range region.Positions {
		c := col(p)
		s := cdsSlot{refPos: p, refBase: decodingArray[ref[c]]}
		switch {
		case decodingArray[query[c]] == "N":
			s.missing = true
		case query[c] != 244:
			s.queBase = decodingArray[query[c]]
			s.snp = (query[c] & ref[c]) < 16
		case c < first || c > last || first < 0:
			s.missing = true
		}
		if k+1 < len(region.Positions) && region.Positions[k+1] == p+region.Strand {
			next := col(region.Positions[k+1])
			var b strings.Builder
			for i := c + region.Strand; i != next; i += region.Strand {
				if ref[i] == 244 && query[i] != 244 {
					b.WriteString(decodingArray[query[i]])
				}
			}
			s.inserted = b.String()
		}
		slots[k] = s
	}
	return slots
}

// slotPool holds slot buffers for cdsSlots, so that a long CDS does not need a new one for each query
var slotPool = sync.Pool{New: func() any { return new([]cdsSlot) }}

// indelWindow is a group of indels in a CDS together with the codons they touch (0-based codon indices,
// inclusive)
type indelWindow struct {
	first, last int
}

// indelWindows groups the indels in a CDS by the codons they touch: a deletion touches the codons of its
// bases, an insertion the codons of the bases either side of it. Indels whose codons overlap are grouped.
func indelWindows(slots []cdsSlot) []indelWindow {
	windows := make([]indelWindow, 0)
	add := func(first, last int) {
		if n := len(windows); n > 0 && first <= windows[n-1].last {
			windows[n-1].last = max(windows[n-1].last, last)
			return
		}
		windows = append(windows, indelWindow{first, last})
	}
	for k := 0; k < len(slots); k++ {
		if slots[k].queBase == "" && !slots[k].missing {
			start := k
			for k+1 < len(slots) && slots[k+1].queBase == "" && !slots[k+1].missing {
				k++
			}
			add(start/3, k/3)
			if slots[k].inserted != "" {
				add(k/3, (k+1)/3)
			}
			continue
		}
		if slots[k].inserted != "" {
			add(k/3, (k+1)/3)
		}
	}
	return windows
}

// getIndelAAsPair finds the amino acid changes caused by in-frame indels in one CDS: deleted residues
// (aa:GENE:H69-), inserted residues (aa:GENE:214EPE, inserted after residue 214) and residues changed
// next to them (aa:GENE:E156G). It returns those changes, the nucleotide changes in the codons the indels
// touch that are not part of one of them, and the codons it has dealt with, which getAAsPair must skip.
// Residues after a stop codon in the query are not reported. Indels are left alone (no amino acid change is reported for them, and their codons are not dealt with
// here) if they shift the reading frame, come after one that does, run outside the CDS, touch a missing
// stretch of the query (N, or uncovered ends), or touch the start codon (a first residue M) or the stop codon (a last residue *).
func getIndelAAsPair(ref, query []byte, region Region, offsetRefCoord []int) ([]Variant, map[int]bool) {

	buf := slotPool.Get().(*[]cdsSlot)
	slots := cdsSlots(ref, query, region, offsetRefCoord, *buf)
	defer func() {
		*buf = slots
		slotPool.Put(buf)
	}()
	nCodons := len(slots) / 3
	handled := make(map[int]bool)
	variants := make([]Variant, 0)
	if nCodons == 0 || len(region.Translation) < nCodons {
		return variants, handled
	}
	refAA := region.Translation[:nCodons]

	// clean reports whether codon j is exactly as in the reference (so an indel may be shifted across it)
	clean := func(j int) bool {
		for k := 3 * j; k < 3*j+3; k++ {
			s := slots[k]
			if s.queBase == "" || s.missing || s.snp || (k < 3*j+2 && s.inserted != "") {
				return false
			}
		}
		return !handled[j]
	}

	// outside reports whether a deletion at the 5' or 3' end of the CDS carries on beyond it
	outside := func(w indelWindow) bool {
		ends := []int{}
		if w.first == 0 && slots[0].queBase == "" {
			ends = append(ends, region.Positions[0]-region.Strand)
		}
		if w.last == nCodons-1 && slots[len(slots)-1].queBase == "" {
			ends = append(ends, region.Positions[len(slots)-1]+region.Strand)
		}
		for _, p := range ends {
			if p < 1 || p > len(offsetRefCoord) {
				return true
			}
			if query[(p-1)+offsetRefCoord[p-1]] == 244 {
				return true
			}
		}
		return false
	}

	frameshifted := false
	for _, w := range indelWindows(slots) {
		if w.last >= nCodons {
			continue
		}
		deleted, inserted, missing := 0, 0, false
		for k := 3 * w.first; k < 3*w.last+3; k++ {
			if slots[k].missing {
				missing = true
			} else if slots[k].queBase == "" {
				deleted++
			}
			if k < 3*w.last+2 {
				inserted += len(slots[k].inserted)
			}
		}
		if (inserted-deleted)%3 != 0 {
			frameshifted = true
		}
		if frameshifted || missing || outside(w) ||
			(w.first == 0 && refAA[0] == 'M') || (w.last == nCodons-1 && refAA[nCodons-1] == '*') {
			continue
		}

		// the query's bases over the window's codons, translated in the CDS's frame
		var q strings.Builder
		for k := 3 * w.first; k < 3*w.last+3; k++ {
			q.WriteString(slots[k].queBase)
			if k < 3*w.last+2 {
				q.WriteString(slots[k].inserted)
			}
		}
		qNuc := q.String()
		if region.Strand == -1 {
			qNuc = alphabet.Complement(qNuc)
		}
		qAA := ""
		for i := 0; i+3 <= len(qNuc); i += 3 {
			if aa, ok := codonDict[qNuc[i:i+3]]; ok {
				qAA += aa
			} else {
				qAA += "X"
			}
		}
		// translation stops at a stop codon, so nothing after one is part of the change
		if i := strings.IndexByte(qAA, '*'); i >= 0 {
			qAA = qAA[:i+1]
		}
		rAA := refAA[w.first : w.last+1]

		// trim the residues they share, first at the N-terminal end, then at the C-terminal end
		a := 0
		for a < len(rAA) && a < len(qAA) && rAA[a] == qAA[a] {
			a++
		}
		b := 0
		for b < len(rAA)-a && b < len(qAA)-a && rAA[len(rAA)-1-b] == qAA[len(qAA)-1-b] {
			b++
		}
		rChanged := rAA[a : len(rAA)-b]
		qChanged := qAA[a : len(qAA)-b]
		firstChanged := w.first + a // codon index of the first reference residue in rChanged

		for j := w.first; j <= w.last; j++ {
			handled[j] = true
		}

		// the nucleotide indels in the window, for --append-snps
		indelNotes := windowIndels(slots, w, region.Strand)

		// residue firstChanged+i of the reference: changed to qChanged[i], or deleted if i >= len(qChanged)
		reported := make(map[int]bool) // codon indices with a reported change
		paired := min(len(rChanged), len(qChanged))
		newVariants := make([]Variant, 0)
		for i := 0; i < len(rChanged); i++ {
			j := firstChanged + i
			que := "-"
			if i < paired {
				que = string(qChanged[i])
				if que == "X" || que == string(rChanged[i]) {
					continue
				}
			}
			reported[j] = true
			newVariants = append(newVariants, Variant{Changetype: "aa", Feature: region.Name, RefAl: string(rChanged[i]), QueAl: que, Position: region.Positions[3*j], Codon: codonPositions(region, j), Residue: j + 1})
		}
		if len(qChanged) > paired {
			ins := qChanged[paired:]
			after := firstChanged + paired - 1 // codon index of the residue the insertion follows
			// a pure deletion or insertion in a run of residues is placed as far C-terminal as possible
			if len(rChanged) == 0 {
				for after+1 < nCodons && refAA[after+1] != '*' && ins[0] == refAA[after+1] && clean(after+1) {
					ins = ins[1:] + ins[:1]
					after++
				}
			}
			if after >= 0 { // (an insertion before the first residue can't be numbered)
				newVariants = append(newVariants, Variant{Changetype: "aa", Feature: region.Name, QueAl: ins, Position: region.Positions[3*after], Codon: codonPositions(region, after), Residue: after + 1})
			}
		} else if len(qChanged) == 0 && len(rChanged) > 0 {
			start, end := firstChanged, firstChanged+len(rChanged)-1
			for end+1 < nCodons && refAA[end+1] != '*' && refAA[start] == refAA[end+1] && clean(end+1) {
				start++
				end++
			}
			if start != firstChanged {
				newVariants = newVariants[:0]
				for j := start; j <= end; j++ {
					newVariants = append(newVariants, Variant{Changetype: "aa", Feature: region.Name, RefAl: string(refAA[j]), QueAl: "-", Position: region.Positions[3*j], Codon: codonPositions(region, j), Residue: j + 1})
				}
			}
		}

		// SNPs in codons with a reported change are part of it; the others are reported as nucleotide changes
		notes := append([]string{}, indelNotes...)
		for k := 3 * w.first; k < 3*w.last+3; k++ {
			s := slots[k]
			if !s.snp {
				continue
			}
			if reported[k/3] {
				notes = append(notes, "nuc:"+s.refBase+strconv.Itoa(s.refPos)+s.queBase)
			} else {
				variants = append(variants, Variant{Changetype: "nuc", RefAl: s.refBase, QueAl: s.queBase, Position: s.refPos})
			}
		}
		for i := range newVariants {
			newVariants[i].SNPs = strings.Join(notes, ";")
		}
		variants = append(variants, newVariants...)
	}
	return variants, handled
}

// codonPositions returns the reference positions of codon j (0-based) of a CDS, in CDS order
func codonPositions(region Region, j int) [3]int {
	return [3]int{region.Positions[3*j], region.Positions[3*j+1], region.Positions[3*j+2]}
}

// windowIndels describes the nucleotide deletions and insertions in a window, as they are reported
// (del:POS:LEN, ins:POS:LEN), in order of reference position
func windowIndels(slots []cdsSlot, w indelWindow, strand int) []string {
	type note struct {
		pos  int
		text string
	}
	notes := make([]note, 0)
	end := 3*w.last + 3
	for k := 3 * w.first; k < end; k++ {
		if slots[k].queBase == "" {
			start := k
			for k+1 < end && slots[k+1].queBase == "" {
				k++
			}
			lo := min(slots[start].refPos, slots[k].refPos)
			notes = append(notes, note{lo, "del:" + strconv.Itoa(lo) + ":" + strconv.Itoa(k-start+1)})
		}
		if k < end-1 && slots[k].inserted != "" {
			pos := slots[k].refPos
			if strand == -1 {
				pos = slots[k+1].refPos
			}
			notes = append(notes, note{pos, "ins:" + strconv.Itoa(pos) + ":" + strconv.Itoa(len(slots[k].inserted))})
		}
	}
	sort.SliceStable(notes, func(i, j int) bool { return notes[i].pos < notes[j].pos })
	texts := make([]string, len(notes))
	for i, n := range notes {
		texts[i] = n.text
	}
	return texts
}
