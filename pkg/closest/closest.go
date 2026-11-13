/*
Package closest provides routines to find the closest sequences
to a set of query sequences, by genetic distance
*/
package closest

import (
	"context"
	"encoding/binary"
	"errors"
	"fmt"
	"io"
	"math"
	"math/bits"
	"os"
	"runtime"
	"strconv"
	"strings"

	"github.com/virus-evolution/gofasta/pkg/encoding"
	"github.com/virus-evolution/gofasta/pkg/fasta"
	"github.com/virus-evolution/gofasta/pkg/pipeline"
)

// resultsStruct is a struct that contains information about a query sequence and its (current)
// single closest target by raw genetic distance, including the snps that distinguish them and
// the total snp-distance between them.
type resultsStruct struct {
	qname        string
	qidx         int
	tname        string
	completeness int64
	distance     float64
	snps         []string
}

// lessDistance reports whether distance a is strictly closer than distance b. An undefined
// distance (NaN, e.g. no comparable sites) is treated as further than every defined distance,
// and two undefined distances are never "less" than one another.
func lessDistance(a, b float64) bool {
	if math.IsNaN(a) {
		return false
	}
	if math.IsNaN(b) {
		return true
	}
	return a < b
}

// The distance measures classify sites without branching (which bases differ is unpredictable, so
// branches would often be mispredicted), eight sites at a time: the eight bytes are read as one uint64,
// each test sets the top bit of each byte (lane) where it holds, and the set bits are counted. The tests:
//   - zeroLanes(x): the lane is 0
//   - disjoint: the bases certainly differ (they share no base: the high nibble of a&b is 0)
//   - known: the base is A, C, G or T (bit 3, moved to bit 7)
// The last len%8 sites are done one at a time, with the same tests on a single byte.

const (
	lo7  = 0x7f7f7f7f7f7f7f7f
	hi1  = 0x8080808080808080
	hiN  = 0xf0f0f0f0f0f0f0f0
	bit3 = 0x0808080808080808
	allA = 0xc8c8c8c8c8c8c8c8 // A|G (200) in every lane
	allY = 0x3838383838383838 // C|T (56) in every lane
)

// zeroLanes sets the top bit of each byte of x that is 0
func zeroLanes(x uint64) uint64 { return ^(((x & lo7) + lo7) | x) & hi1 }

func disjointLanes(a, b uint64) uint64 { return zeroLanes(a & b & hiN) }

func knownLanes(a uint64) uint64 { return (a & bit3) << 4 }

// isZero, disjoint and known are the single-site tests, giving 0 or 1
func isZero(x uint32) uint32 { return (x - 1) >> 31 }

func disjoint(a, b byte) uint32 { return isZero(uint32(a&b) >> 4) }

func known(a byte) uint32 { return uint32(a>>3) & 1 }

// rawDistance is the number of sites where the bases certainly differ per site where they are certainly
// different or certainly the same
func rawDistance(query, target fasta.EncodedRecord) float64 {
	q, t := query.Seq, target.Seq[:len(query.Seq)]
	var n, same int
	i := 0
	for ; i+8 <= len(q); i += 8 {
		a, b := binary.LittleEndian.Uint64(q[i:]), binary.LittleEndian.Uint64(t[i:])
		n += bits.OnesCount64(disjointLanes(a, b))
		same += bits.OnesCount64(knownLanes(a) & zeroLanes(a^b))
	}
	for ; i < len(q); i++ {
		a, b := q[i], t[i]
		n += int(disjoint(a, b))
		same += int(known(a) & isZero(uint32(a^b)))
	}
	return float64(n) / float64(n+same)
}

// snpDistance is the number of sites where the bases certainly differ
func snpDistance(query, target fasta.EncodedRecord) float64 {
	q, t := query.Seq, target.Seq[:len(query.Seq)]
	n := 0
	i := 0
	for ; i+8 <= len(q); i += 8 {
		n += bits.OnesCount64(disjointLanes(binary.LittleEndian.Uint64(q[i:]), binary.LittleEndian.Uint64(t[i:])))
	}
	for ; i < len(q); i++ {
		n += int(disjoint(q[i], t[i]))
	}
	return float64(n)
}

// See equation (7) in Tamura K, Nei M. Estimation of the number of nucleotide substitutions in the control region of mitochondrial DNA in humans and chimpanzees.
// Mol Biol Evol. 1993 May;10(3):512-26. doi: 10.1093/oxfordjournals.molbev.a040023. PMID: 8336541.
// See also ape: https://github.com/cran/ape/blob/c2fd899f66d6493a80484033772a3418e5d706a4/src/dist_dna.c
func tn93Distance(query, target fasta.EncodedRecord) float64 {
	count_L, count_d, count_P1, count_P2 := tn93Counts(query.Seq, target.Seq)
	return tn93FromCounts(query, target, count_L, count_d, count_P1, count_P2)
}

// tn93Counts returns, over the sites where both bases are known for sure, the number of sites
// (count_L), of differences (count_d), of A <-> G differences (count_P1) and of C <-> T differences
// (count_P2)
func tn93Counts(query, target []byte) (count_L, count_d, count_P1, count_P2 int) {
	target = target[:len(query)]
	same := 0
	i := 0
	for ; i+8 <= len(query); i += 8 {
		a, b := binary.LittleEndian.Uint64(query[i:]), binary.LittleEndian.Uint64(target[i:])
		ka := knownLanes(a)
		diff := disjointLanes(a, b) & ka & knownLanes(b) // the bases are different (and known for sure)
		same += bits.OnesCount64(ka & zeroLanes(a^b))    // the bases are certainly the same
		count_d += bits.OnesCount64(diff)
		count_P1 += bits.OnesCount64(diff & zeroLanes((a|b)^allA)) // one is adenine and the other guanine
		count_P2 += bits.OnesCount64(diff & zeroLanes((a|b)^allY)) // one is cytosine and the other thymine
	}
	for ; i < len(query); i++ {
		a, b := query[i], target[i]
		ka := known(a)
		diff := disjoint(a, b) & ka & known(b)
		same += int(ka & isZero(uint32(a^b)))
		count_d += int(diff)
		count_P1 += int(diff & isZero(uint32(a|b)^200))
		count_P2 += int(diff & isZero(uint32(a|b)^56))
	}
	return same + count_d, count_d, count_P1, count_P2
}

// tn93FromCounts is the TN93 distance given the counts from tn93Counts. The base frequencies are
// those of the two whole sequences.
func tn93FromCounts(query, target fasta.EncodedRecord, count_L, count_d, count_P1, count_P2 int) float64 {

	// Total ATGC length of the two sequences
	L := float64(target.Count_A + target.Count_C + target.Count_G + target.Count_T + query.Count_A + query.Count_C + query.Count_G + query.Count_T)

	// estimates of the equilibrium base contents from the pair's sequence data
	g_A := float64(target.Count_A+query.Count_A) / L
	g_C := float64(target.Count_C+query.Count_C) / L
	g_G := float64(target.Count_G+query.Count_G) / L
	g_T := float64(target.Count_T+query.Count_T) / L

	g_R := float64(target.Count_A+query.Count_A+target.Count_G+query.Count_G) / L
	g_Y := float64(target.Count_C+query.Count_C+target.Count_T+query.Count_T) / L

	// tidies up the equations a bit, after ape
	k1 := 2.0 * g_A * g_G / g_R
	k2 := 2.0 * g_T * g_C / g_Y
	k3 := 2.0 * (g_R*g_Y - g_A*g_G*g_Y/g_R - g_T*g_C*g_R/g_Y)

	// estimated rates from this pairwise comparison
	P1 := float64(count_P1) / float64(count_L)                   // rate of changes which are transitional differences between purines (A ⇄ G)
	P2 := float64(count_P2) / float64(count_L)                   // rate of changes which are transitional differences between pyramidines (C ⇄ T)
	Q := float64(count_d-(count_P1+count_P2)) / float64(count_L) // rate of changes which are transversional differences  (A ⇄ C || A ⇄ T || G ⇄ T || C ⇄ G) (i.e. everything else)

	// tidies up the equations a bit, after ape
	w1 := 1.0 - P1/k1 - Q/(2*g_R)
	w2 := 1.0 - P2/k2 - Q/(2*g_Y)
	w3 := 1.0 - Q/(2*g_R*g_Y)

	// tn93 distance:
	d := -k1*math.Log(w1) - k2*math.Log(w2) - k3*math.Log(w3)

	if d == 0.0 {
		d = 0.0
	}

	return d
}

// scoredTarget is one target with its distance to each query
type scoredTarget struct {
	target    fasta.EncodedRecord
	distances []float64
}

// loadQueries reads and encodes the queries. A query that can't be used (an invalid character, or a
// different width from the first usable query) is reported on stderr and left out, unless strict, when
// it is an error. It writes the number of queries read to stderr.
func loadQueries(query io.Reader, strict bool) ([]fasta.EncodedRecord, error) {
	width := -1 // of the first usable query; only one worker, so the queries are seen in order
	work := func(record fasta.Record) (fasta.EncodedRecord, error) {
		q, err := fasta.EncodeRecord(record, false, true, false)
		if err != nil {
			return q, pipeline.Named(record.ID, err)
		}
		if width == -1 {
			width = len(q.Seq)
		} else if len(q.Seq) != width {
			return q, pipeline.Named(record.ID, errors.New("not the same width as the first query"))
		}
		return q, nil
	}
	queries := make([]fasta.EncodedRecord, 0)
	h := pipeline.NewHandler(strict, os.Stderr, pipeline.Noun{One: "query", Many: "queries"}, func(q fasta.EncodedRecord) error {
		queries = append(queries, q)
		return nil
	})
	if err := pipeline.Run(context.Background(), 1, fasta.NewRecordSource(query), work, h.Sink); err != nil {
		return nil, err
	}
	if err := h.Finish(); err != nil {
		return nil, err
	}
	fmt.Fprintf(os.Stderr, "number of sequences in query alignment: %d\n", h.Total())
	return queries, nil
}

// scanTargets reads the targets from an alignment, works out each one's distance to every query (on
// threads goroutines), and passes the targets with their distances to visit in input order. A target
// that can't be used (an invalid character, or a different width from the queries) is reported on
// stderr and skipped, unless strict, when it is an error. It writes the number of targets read to stderr.
func scanTargets(queries []fasta.EncodedRecord, target io.Reader, measure string, threads int, strict bool, visit func(scoredTarget)) error {
	var distance func(query, target fasta.EncodedRecord) float64
	switch measure {
	case "raw":
		distance = rawDistance
	case "snp":
		distance = snpDistance
	case "tn93":
		distance = tn93Distance
	default:
		return errors.New("unknown distance measure: " + measure)
	}

	work := func(record fasta.Record) (scoredTarget, error) {
		t, err := fasta.EncodeRecord(record, false, true, true)
		if err != nil {
			return scoredTarget{}, pipeline.Named(record.ID, err)
		}
		if len(queries) > 0 && len(t.Seq) != len(queries[0].Seq) {
			return scoredTarget{}, pipeline.Named(record.ID, errors.New("query and target alignments are not the same width"))
		}
		ds := make([]float64, len(queries))
		for i, q := range queries {
			ds[i] = distance(q, t)
		}
		return scoredTarget{target: t, distances: ds}, nil
	}
	h := pipeline.NewHandler(strict, os.Stderr, pipeline.Noun{One: "target", Many: "targets"}, func(st scoredTarget) error {
		visit(st)
		return nil
	})
	if err := pipeline.Run(context.Background(), threads, fasta.NewRecordSource(target), work, h.Sink); err != nil {
		return err
	}
	if err := h.Finish(); err != nil {
		return err
	}
	fmt.Fprintf(os.Stderr, "number of sequences in target alignment: %d\n", h.Total())
	return nil
}

// snpsBetween lists the sites at which query and target certainly differ, as position, query base, target base
func snpsBetween(query, target fasta.EncodedRecord, decoding [256]string) []string {
	snps := make([]string, 0)
	for i, tNuc := range target.Seq {
		if (query.Seq[i] & tNuc) < 16 {
			snps = append(snps, strconv.Itoa(i+1)+decoding[query.Seq[i]]+decoding[tNuc])
		}
	}
	return snps
}

// consider makes target the closest to query if it is closer than the current closest, or as close and
// more complete. The first target considered always becomes the closest.
func (closest *resultsStruct) consider(query, target fasta.EncodedRecord, distance float64, first bool, decoding [256]string) {
	if first || lessDistance(distance, closest.distance) || (distance == closest.distance && target.Score > closest.completeness) {
		*closest = resultsStruct{qname: query.ID, qidx: query.Idx, tname: target.ID, completeness: target.Score, distance: distance, snps: snpsBetween(query, target, decoding)}
	}
}

// writeClosest parses an array of resultsStructs in order to write them, usually to stdout or file
func writeClosest(results []resultsStruct, measure string, w io.Writer) error {

	var err error

	_, err = w.Write([]byte("query,closest,distance,SNPs\n"))
	if err != nil {
		return err
	}

	for _, result := range results {
		switch measure {
		case "raw":
			w.Write([]byte(result.qname + "," + result.tname + "," + strconv.FormatFloat(result.distance, 'f', 9, 64) + "," + strings.Join(result.snps, ";") + "\n"))
		case "snp":
			w.Write([]byte(result.qname + "," + result.tname + "," + strconv.Itoa(int(result.distance)) + "," + strings.Join(result.snps, ";") + "\n"))
		case "tn93":
			w.Write([]byte(result.qname + "," + result.tname + "," + strconv.FormatFloat(result.distance, 'f', 9, 64) + "," + strings.Join(result.snps, ";") + "\n"))
		}
	}

	return nil
}

// Closest finds the single closest sequence by genetic distance to a query/queries. It writes the results
// to stdout or to file. Ties for distance are broken by genome completeness
func Closest(query, target io.Reader, measure string, out io.Writer, threads int, strict bool) error {

	if threads == 0 {
		threads = runtime.NumCPU()
	} else if threads < runtime.NumCPU() {
		runtime.GOMAXPROCS(threads)
	}

	queries, err := loadQueries(query, strict)
	if err != nil {
		return err
	}

	nQ := len(queries)

	QResultsArray := make([]resultsStruct, nQ)
	for i, q := range queries {
		QResultsArray[i] = resultsStruct{qname: q.ID, qidx: q.Idx}
	}

	decoding := encoding.MakeDecodingArray()
	first := true
	err = scanTargets(queries, target, measure, threads, strict, func(st scoredTarget) {
		for i, q := range queries {
			QResultsArray[i].consider(q, st.target, st.distances[i], first, decoding)
		}
		first = false
	})
	if err != nil {
		return err
	}

	err = writeClosest(QResultsArray, measure, out)
	if err != nil {
		return err
	}

	return nil
}
