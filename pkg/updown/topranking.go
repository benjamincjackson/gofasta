/*
Package updown implements functions that calculate pseudo-tree aware
snp-distances between sequences. It relies on there being a suitable
sequence that one can use to represent the "root" (ancestor) of a
hypothetical phylogenetic tree, and for the distance between relevant
sequences to be small enough that multiple hits are unlikely
*/
package updown

import (
	"context"
	"errors"
	"io"
	"math"
	"os"
	"slices"
	"sort"
	"strconv"
	"strings"

	"github.com/virus-evolution/gofasta/pkg/encoding"
	"github.com/virus-evolution/gofasta/pkg/fasta"
	"github.com/virus-evolution/gofasta/pkg/pipeline"
)

/*
           Venn diagram of SNPs

Q: in query
T: in target
QT: in query and target

          _____             _____
   X‾‾‾‾‾‾     ‾‾‾‾xxxx‾‾‾‾‾     ‾‾‾‾X
  XX            XX     XXX           XXX
 XX           XXX        XXX            XX
X             X            X             |
|             |            |             |
|     Q       |     QT     |      T      |
|             |            |             |
|             X            X             |
XX            XXX         XX            XX
 XXX            XXX    XXX            XXX
   xx-----_____----XXXX----_____-----xx


Table of relationships

   Q   |  QT  |   T   |  interpretation           |   bucket for T
-------------------------------------------------------------------
   0       0      0      on a polytomy (with ref)     same
   0      >0      0      on a polytomy                same
  >0      >0      0      Q is child of T              up
  >0       0      0      Q is child of T              up
   0      >0     >0      T is child of Q              down
   0       0     >0      T is child of Q              down
  >0      >0     >0      Q and T are siblings         side
  >0       0     >0      Q and T are orphans          side?

*/

// isSiteAmb asks if a position is within an ambiguity tract, given an array of pairs of start-stop coordinates
// of such tracts. Co-ordinates are 1-based, inclusive
func isSiteAmb(pos int, a []int) bool {
	for i := 0; i < len(a); i += 2 {
		if pos >= a[i] && pos <= a[i+1] {
			return true
		}
	}
	return false
}

// inOrder reports whether SNP positions are strictly ascending and ambiguity tracts are ascending and
// do not overlap, as updown list writes them (so countMerge can be used)
func inOrder(pos, ambs []int) bool {
	for i := 1; i < len(pos); i++ {
		if pos[i] <= pos[i-1] {
			return false
		}
	}
	for i := 0; i+1 < len(ambs); i += 2 {
		if ambs[i] > ambs[i+1] || (i > 0 && ambs[i] <= ambs[i-1]) {
			return false
		}
	}
	return len(ambs)%2 == 0
}

// whichWay returns direction values of 0,1,2,3 = same,up,down,side respectively + SNP distance.
// distance is -1 if the pair fails the ambiguity threshold test.
func whichWay(q, t updownLine, thresh float32) (int, int) {
	// table has 4 items: number of Q SNPs; no. QT SNPs; no. T SNPs; no of consequential ambiguous sites for this pair
	var table [4]int
	var distance int
	if q.ordered && t.ordered {
		table, distance = countMerge(q, t)
	} else {
		table, distance = countSearch(q, t)
	}

	sum := 0
	for _, c := range table {
		sum += c
	}
	if (float32(table[3]) / float32(sum)) > thresh {
		return 0, -1
	}

	// Table of relationships
	//
	//  t[0]  | t[1] | t[2]  |  interpretation           |   bucket for T
	// -------------------------------------------------------------------
	//    0       0      0      on a polytomy (with ref)     same
	//    0      >0      0      on a polytomy                same
	//   >0      >0      0      Q is child of T              up
	//   >0       0      0      Q is child of T              up
	//    0      >0     >0      T is child of Q              down
	//    0       0     >0      T is child of Q              down
	//   >0      >0     >0      Q and T are siblings         side
	//   >0       0     >0      Q and T are orphans          side?
	var direction int
	switch {
	case table[0] == 0 && table[2] == 0:
		direction = 0
	case table[0] > 0 && table[2] == 0:
		direction = 1
	case table[0] == 0 && table[2] > 0:
		direction = 2
	case table[0] > 0 && table[2] > 0:
		direction = 3
	}

	// this is wrong (when there are multiple hits):
	// distance := table[0] + table[2]

	return direction, distance
}

// countMerge is countSearch for two lines that are inOrder: it walks the two SNP lists (and the
// ambiguity tracts) together once, by position. Each list has at most one SNP per position, and SNPs at
// the same position are the same SNP if their strings are equal.
func countMerge(q, t updownLine) (table [4]int, distance int) {
	// inAmb reports whether p is in one of the tracts a; k is where the last call left off (positions
	// are visited in ascending order)
	inAmb := func(a []int, k *int, p int) bool {
		for *k < len(a) && a[*k+1] < p {
			*k += 2
		}
		return *k < len(a) && a[*k] <= p
	}
	i, j, ki, kj := 0, 0, 0, 0
	for i < len(q.snps) || j < len(t.snps) {
		switch {
		case j == len(t.snps) || (i < len(q.snps) && q.snpsPos[i] < t.snpsPos[j]): // only in the query
			if inAmb(t.ambs, &ki, q.snpsPos[i]) {
				table[3]++
			} else {
				table[0]++
				distance++
			}
			i++
		case i == len(q.snps) || t.snpsPos[j] < q.snpsPos[i]: // only in the target
			if inAmb(q.ambs, &kj, t.snpsPos[j]) {
				table[3]++
			} else {
				table[2]++
				distance++
			}
			j++
		default: // both have a SNP at this position
			p, same := q.snpsPos[i], q.snps[i] == t.snps[j]
			counted := false // the position already counts towards the distance
			if inAmb(t.ambs, &ki, p) {
				table[3]++
			} else if same {
				table[1]++
			} else {
				table[0]++
				distance++
				counted = true
			}
			if inAmb(q.ambs, &kj, p) {
				table[3]++
			} else if !same {
				table[2]++
				if !counted {
					distance++
				}
			}
			i++
			j++
		}
	}
	return table, distance
}

// countSearch counts, for whichWay, the query's SNPs not in the target (table[0]), those in both
// (table[1]), the target's SNPs not in the query (table[2]) and SNPs at sites ambiguous in the other
// sequence (table[3]), and the number of positions at which they differ (distance)
func countSearch(q, t updownLine) (table [4]int, distance int) {

	// d is an array of unique SNP positions to get distance from
	d := make([]int, 0)

	for i, qsnp := range q.snps {
		// is the site ambiguous in the target
		if isSiteAmb(q.snpsPos[i], t.ambs) {
			table[3]++
			// if this exact snp (allele and position) is present in the target
		} else if _, found := slices.BinarySearch(t.snpsSorted, qsnp); found {
			table[1]++
		} else {
			// otherwise we +1 the Q snps bin and append the snp to d
			table[0]++
			d = append(d, q.snpsPos[i])
		}
	}
	d_plus := 0
	for i, tsnp := range t.snps {
		// is the site ambiguous in the query
		if isSiteAmb(t.snpsPos[i], q.ambs) {
			table[3]++
			// if there is no overlap between the exact snp (allele and position) and the query's snps
		} else if _, found := slices.BinarySearch(q.snpsSorted, tsnp); !found {
			// we +1 the T snps bin
			table[2]++
			// if this position has not already been added to the list of differences we use to calculate d, then add it here too (otherwise two mutations historically at the same position, but only counts 1 for distance)
			if _, found := slices.BinarySearch(d, t.snpsPos[i]); !found {
				d_plus++
			}
		}
	}
	return table, len(d) + d_plus
}

// a resultsStruct contains information about the relationship between one query and one target.
type resultsStruct struct {
	qname    string // name of the query
	qidx     int    // query's position in the input file
	tname    string // name of the target
	distance int    // snp distance between query and target
	ambCount int    // number of non-ATGC characters in the target
}

// an updownCatchmentSubStruct contains an array of resultsStructs which are the current closest neighbours of one query
// in one direction. The fields maxDist and minAmbig are used when deciding whether or not the next target should replace
// the last item in catchment
type updownCatchmentSubStruct struct {
	catchment []resultsStruct
	maxDist   int
	minAmbig  int
	// nDists    int
}

// an updownCatchmentStruct contains the current lists of targets in each direction for one query sequence
type updownCatchmentStruct struct {
	qname string
	qidx  int
	same  updownCatchmentSubStruct
	up    updownCatchmentSubStruct
	down  updownCatchmentSubStruct
	side  updownCatchmentSubStruct
	// minDistArray [4]int
}

// rearrangeCatchment sorts the catchment slice in an updownCatchmentSubStruct, pops the last item from it, and updates
// the distance and ambiguity fields in preparation for deciding what to do with the next target. It should only be
// called when the updownCatchmentSubStruct is over its capacity (catchmentSize)
func rearrangeCatchment(nS *updownCatchmentSubStruct, catchmentSize int) {
	byDistance(nS.catchment)
	nS.catchment = nS.catchment[0:catchmentSize]
	nS.maxDist = nS.catchment[catchmentSize-1].distance
	nS.minAmbig = nS.catchment[catchmentSize-1].ambCount
}

// a pushCatchmentSubStruct is used to keep a record of the targets at different snp distances away, in case there are no
// targets under/at the threshold distance defined on the command line and the user wants to expand the search to encompass
// the closest sequences in any particular direction. catchmentMap is a map with snp-distances as keys, and arrays of
// resultsStructs as values
type pushCatchmentSubStruct struct {
	catchmentMap map[int][]resultsStruct
	nDists       int // how many snp distances in total
	maxDist      int // the furthest snp distance
}

// getMaxKey gets the largest value key from a map whose keys are snp-distances
func getMaxKey(m map[int][]resultsStruct) int {
	max := 0
	for k := range m {
		if k > max {
			max = k
		}
	}
	return max
}

// refactorPushCatchment is like rearrangeCatchment but for pushCatchmentSubStructs. It either adds targets to snp distances that
// already exist, or adds new snp distances if they don't and the map isn't at capacity, or if it is but the new snp distance is
// lower than any that already exist. In the latter case, it deletes the key and value associated with the further snp distance
func refactorPushCatchment(pC *pushCatchmentSubStruct, rS resultsStruct, nodeDistance int) {
	if _, ok := pC.catchmentMap[rS.distance]; ok {
		// if this sequence's distance is already present, then add to it
		pC.catchmentMap[rS.distance] = append(pC.catchmentMap[rS.distance], rS)
		return
	}
	// otherwise add a new distance, first deleting the farthest one if the map is at capacity
	if len(pC.catchmentMap) == nodeDistance {
		delete(pC.catchmentMap, getMaxKey(pC.catchmentMap))
	}
	pC.catchmentMap[rS.distance] = []resultsStruct{rS}
	pC.maxDist = getMaxKey(pC.catchmentMap)
	pC.nDists = len(pC.catchmentMap)
}

// pushCatchment2Catchment converts a pushCatchmentSubStruct to an updownCatchmentSubStruct (in no particular
// order: see byDistance)
func pushCatchment2Catchment(pC pushCatchmentSubStruct) updownCatchmentSubStruct {
	var uC updownCatchmentSubStruct
	for _, v := range pC.catchmentMap {
		uC.catchment = append(uC.catchment, v...)
	}
	return uC
}

// byDistance sorts a catchment by snp distance, then by number of ambiguities, keeping ties in their order
func byDistance(c []resultsStruct) {
	sort.SliceStable(c, func(i, j int) bool {
		return c[i].distance < c[j].distance || (c[i].distance == c[j].distance && c[i].ambCount < c[j].ambCount)
	})
}

// consider adds rs to the catchment if it has fewer than size members, or if rs is closer than the
// furthest member (or as close, with fewer ambiguities)
func (nS *updownCatchmentSubStruct) consider(rs resultsStruct, size int) {
	switch {
	case len(nS.catchment) < size:
		nS.catchment = append(nS.catchment, rs)
		if len(nS.catchment) == size {
			rearrangeCatchment(nS, size)
		}
	case rs.distance < nS.maxDist,
		rs.distance == nS.maxDist && rs.ambCount < nS.minAmbig:
		nS.catchment = append(nS.catchment, rs)
		rearrangeCatchment(nS, size)
	}
}

// a search collects the neighbours of one query, one target at a time
type search interface {
	// add considers a target that is in direction (0,1,2,3 = same,up,down,side) at snp distance
	add(target updownLine, direction, distance int)
	// finish returns the neighbours
	finish() updownCatchmentStruct
}

// catchmentSearch assigns targets to the bins by direction, keeping the closest in each, within the
// distance limits for each bin. It can't do any pushing if any bins are empty after all targets have
// been processed.
type catchmentSearch struct {
	neighbours updownCatchmentStruct
	sizeArray  [4]int
	sizetotal  int // the maximum size of any bin (for filling in)
	distArray  [4]int
	nofill     bool
}

func newCatchmentSearch(q updownLine, sizeArray [4]int, nofill bool, distArray [4]int) *catchmentSearch {
	s := &catchmentSearch{sizeArray: sizeArray, distArray: distArray, nofill: nofill}
	s.neighbours = updownCatchmentStruct{qname: q.id, qidx: q.idx}
	s.sizetotal = sum4(sizeArray)
	for _, n := range sizeArray {
		if n == math.MaxInt32 {
			s.sizetotal = math.MaxInt32
		}
	}
	return s
}

func (s *catchmentSearch) bin(direction int) *updownCatchmentSubStruct {
	return [4]*updownCatchmentSubStruct{&s.neighbours.same, &s.neighbours.up, &s.neighbours.down, &s.neighbours.side}[direction]
}

func (s *catchmentSearch) add(target updownLine, direction, distance int) {
	// if the distance is too big, skip it
	if distance > s.distArray[direction] {
		return
	}
	s.bin(direction).consider(resultsStruct{tname: target.id, ambCount: target.ambCount, distance: distance}, s.sizetotal)
}

func (s *catchmentSearch) finish() updownCatchmentStruct {
	var sizeObserved [4]int
	for direction := range sizeObserved {
		b := s.bin(direction)
		// a bin that never filled up has not been sorted yet
		if n := len(b.catchment); n < s.sizetotal && n > 0 {
			rearrangeCatchment(b, n)
		}
		sizeObserved[direction] = len(b.catchment)
	}

	size := balance(s.sizetotal, s.sizeArray, sizeObserved, s.nofill)
	for direction := range size {
		b := s.bin(direction)
		b.catchment = b.catchment[0:size[direction]]
	}
	return s.neighbours
}

// pushSearch assigns targets to the bins by direction, keeping, for the up, down and side bins, all the
// targets at the pushDist smallest snp distances found in that direction (so that a bin is not empty
// just because nothing is within a fixed distance).
type pushSearch struct {
	q              updownLine
	pushDist       int
	same           updownCatchmentSubStruct // the sequences on the polytomy
	up, down, side pushCatchmentSubStruct
}

func newPushSearch(q updownLine, pushDist int) *pushSearch {
	return &pushSearch{
		q:        q,
		pushDist: pushDist,
		up:       pushCatchmentSubStruct{catchmentMap: make(map[int][]resultsStruct)},
		down:     pushCatchmentSubStruct{catchmentMap: make(map[int][]resultsStruct)},
		side:     pushCatchmentSubStruct{catchmentMap: make(map[int][]resultsStruct)},
	}
}

func (s *pushSearch) add(target updownLine, direction, distance int) {
	rs := resultsStruct{tname: target.id, ambCount: target.ambCount, distance: distance}
	if direction == 0 {
		s.same.catchment = append(s.same.catchment, rs)
		return
	}
	pC := [4]*pushCatchmentSubStruct{nil, &s.up, &s.down, &s.side}[direction]
	// is the distance lower, or are there not enough distances yet
	if distance <= pC.maxDist || pC.nDists < s.pushDist {
		refactorPushCatchment(pC, rs, s.pushDist)
	}
}

func (s *pushSearch) finish() updownCatchmentStruct {
	neighbours := updownCatchmentStruct{qname: s.q.id, qidx: s.q.idx, same: s.same}
	neighbours.up = pushCatchment2Catchment(s.up)
	byDistance(neighbours.up.catchment)
	neighbours.down = pushCatchment2Catchment(s.down)
	byDistance(neighbours.down.catchment)
	neighbours.side = pushCatchment2Catchment(s.side)
	byDistance(neighbours.side.catchment)
	// TO DO - balance the neighbours here if sizeArray is given? Would need a balance function specific to pushing
	return neighbours
}

// balance returns the size each bin's output should have, given the ideal sizes from the options and the
// observed sizes of the bins: a bin short of its ideal size keeps what it has, and (unless nofill) the
// others make up the shortfall from their surplus, one target at a time in bin order, up to sizetotal
func balance(sizetotal int, sizeIdeal [4]int, sizeObserved [4]int, nofill bool) [4]int {
	var size, surplus [4]int
	short := false
	for i := range sizeObserved {
		size[i] = min(sizeObserved[i], sizeIdeal[i])
		surplus[i] = max(sizeObserved[i]-sizeIdeal[i], 0)
		short = short || sizeObserved[i] < sizeIdeal[i]
	}
	if !short {
		return sizeIdeal
	}
	if nofill {
		return size
	}
	for {
		for i := range size {
			if sum4(surplus) == 0 {
				return size
			}
			if surplus[i] > 0 {
				size[i]++
				surplus[i]--
			}
			if sum4(size) == sizetotal {
				return size
			}
		}
	}
}

// binNames are the names of the bins, in direction order (0,1,2,3)
var binNames = [4]string{"same", "up", "down", "side"}

// bins returns a query's bins in direction order
func (r updownCatchmentStruct) bins() [4]updownCatchmentSubStruct {
	return [4]updownCatchmentSubStruct{r.same, r.up, r.down, r.side}
}

// writeUpDownCatchment writes the catchments for each query to a file/stdout
func writeUpDownCatchment(w io.Writer, results []updownCatchmentStruct) error {
	if _, err := w.Write([]byte("query,closestsame,closestup,closestdown,closestside\n")); err != nil {
		return err
	}
	for _, result := range results {
		fields := []string{result.qname}
		for _, bin := range result.bins() {
			names := make([]string, len(bin.catchment))
			for i, hit := range bin.catchment {
				names[i] = hit.tname
			}
			fields = append(fields, strings.Join(names, ";"))
		}
		if _, err := w.Write([]byte(strings.Join(fields, ",") + "\n")); err != nil {
			return err
		}
	}
	return nil
}

// writeUpdownTable writes the output in table format, including SNP-distances
func writeUpdownTable(w io.Writer, results []updownCatchmentStruct) error {
	if _, err := w.Write([]byte("query,direction,distance,target\n")); err != nil {
		return err
	}
	for _, result := range results {
		for direction, bin := range result.bins() {
			for _, neighbour := range bin.catchment {
				if _, err := w.Write([]byte(result.qname + "," + binNames[direction] + "," + strconv.Itoa(neighbour.distance) + "," + neighbour.tname + "\n")); err != nil {
					return err
				}
			}
		}
	}
	return nil
}

// true false every item in s == 0
func allZero(s []int) bool {
	for _, n := range s {
		if n != 0 {
			return false
		}
	}
	return true
}

// CheckTopRankingArgs checks the --size-* and --dist-* values given to TopRanking, which also checks
// them, so that a command can report a problem before it opens any input
func CheckTopRankingArgs(sizetotal int, sizeup int, sizedown int, sizeside int, sizesame int,
	distall int, distup int, distdown int, distside int, distpush int) error {
	_, _, err := checkArgs(sizetotal, sizeup, sizedown, sizeside, sizesame, distall, distup, distdown, distside, distpush, io.Discard)
	return err
}

// checkArgs sanity checks the command line options, writing warnings about overridden options to warn
func checkArgs(sizetotal int, sizeup int, sizedown int, sizeside int, sizesame int,
	distall int, distup int, distdown int, distside int, distpush int, warn io.Writer) ([4]int, [4]int, error) {

	for _, v := range []int{sizetotal, sizeup, sizedown, sizeside, sizesame, distall, distup, distdown, distside, distpush} {
		if v < 0 {
			return [4]int{}, [4]int{}, errors.New("--size-* and --dist-* values must not be negative")
		}
	}

	if sizetotal == 0 && allZero([]int{sizeup, sizedown, sizeside, sizesame}) && distpush == 0 && allZero([]int{distup, distdown, distside, distall}) {
		return [4]int{}, [4]int{}, errors.New(`
Please provide values to either --size-total,
or all or some of --size-up, --size-down, --size-side and --size-same,

                            * or *

to --dist-all,
or all or some of --dist-up, --dist-down, --dist-side
or to --dist-push`)
	}

	if (sizetotal != 0 && !allZero([]int{sizeup, sizedown, sizeside, sizesame})) && !allZero([]int{distup, distdown, distside, distall}) && (distpush > 0) {
		return [4]int{}, [4]int{}, errors.New("you can't combine --dist* and --size* flags if you also invoke --dist-push")
	}

	if sizetotal != 0 && !allZero([]int{sizeup, sizedown, sizeside, sizesame}) {
		io.WriteString(warn, "warning: setting --size-total overrides --size-up, --size-down, --size-side and --size-same\n")
	}

	if distall != 0 && !allZero([]int{distup, distdown, distside}) {
		io.WriteString(warn, "warning: setting --dist-all overrides --dist-up, --dist-down and --dist-side\n")
	}

	// bin sizes and maximum distances, for same,up,down,side respectively
	all := [4]int{math.MaxInt32, math.MaxInt32, math.MaxInt32, math.MaxInt32}

	sizeArray := all // get everything (sorted by distance, or filtered on distance, later)
	switch {
	case sizetotal > 0: // split evenly, any remainder going to same
		q := sizetotal / 4
		sizeArray = [4]int{sizetotal - 3*q, q, q, q}
	case !allZero([]int{sizeup, sizedown, sizeside, sizesame}):
		sizeArray = [4]int{sizesame, sizeup, sizedown, sizeside}
	}

	distArray := all
	switch {
	case distall > 0:
		distArray = [4]int{0, distall, distall, distall}
	case !allZero([]int{distup, distdown, distside}):
		distArray = [4]int{0, distup, distdown, distside}
	}

	return sizeArray, distArray, nil
}

// returns the sum of all the integers in an length-4 array of ints
func sum4(a [4]int) int {
	var t int
	for _, n := range a {
		t += n
	}
	return t
}

// targetIn is a target of topranking as read: a fasta record or a line of updown list's csv output
type targetIn struct {
	record *fasta.Record
	row    csvRow
}

// RecordName returns the target's name (see pipeline.Namer)
func (t targetIn) RecordName() string {
	if t.record != nil {
		return t.record.ID
	}
	return t.row.RecordName()
}

// TopRanking finds pseudo-tree-aware catchments for query sequences, given a large database of target sequences, the closest
// of which should be returned in the output. Targets are split into bins depending on whether they are likely direct ancestors of,
// direct descendants of, polyphyletic with, or exactly the same as, the query
// threads is the number of goroutines that process records (fewer than 1 means 1).
func TopRanking(query, target, reference io.Reader, out io.Writer, table bool,
	q_in_type, t_in_type string, ignoreArray []string,
	sizetotal int, sizeup int, sizedown int, sizeside int, sizesame int,
	distall int, distup int, distdown int, distside int,
	threshpair float32, threshtarg int, nofill bool, distpush int, threads int, strict bool) error {

	sizeArray, distArray, err := checkArgs(sizetotal, sizeup, sizedown, sizeside, sizesame, distall, distup, distdown, distside, distpush, os.Stderr)
	if err != nil {
		return err
	}

	var refSeq []byte
	if q_in_type == "fasta" || t_in_type == "fasta" {
		temp, err := fasta.LoadEncodeAlignment(reference, false, false, false)
		if err != nil {
			return err
		}
		if len(temp) > 1 {
			return errors.New("More than one record in --reference")
		}
		refSeq = temp[0].Seq
	}

	var queries []updownLine

	switch q_in_type {
	case "csv":
		queries, err = readCSVToUDLList(query, threads, strict)
		if err != nil {
			return err
		}
	case "fasta":
		queries, err = fastaToUDLList(query, refSeq, threads, strict)
		if err != nil {
			return err
		}
	}

	searches := make([]search, len(queries))
	for i, q := range queries {
		if distpush > 0 {
			searches[i] = newPushSearch(q, distpush)
		} else {
			searches[i] = newCatchmentSearch(q, sizeArray, nofill, distArray)
		}
	}

	// a target with its direction and snp distance from each query
	type rankedTarget struct {
		line                 updownLine
		skip                 bool
		directions, distance []int
	}

	var next func() (targetIn, error)
	switch t_in_type {
	case "csv":
		src := newCSVSource(target, errors.New("bad header when parsing --target csv: is this the output of gofasta updown list?"))
		next = func() (targetIn, error) {
			row, err := src()
			return targetIn{row: row}, err
		}
	case "fasta":
		src := fasta.NewRecordSource(target)
		next = func() (targetIn, error) {
			record, err := src()
			return targetIn{record: &record}, err
		}
	default:
		return errors.New("unknown --target file type: " + t_in_type)
	}

	DA := encoding.MakeDecodingArray()
	work := func(in targetIn) (rankedTarget, error) {
		var t updownLine
		var err error
		if in.record != nil {
			t, err = lineFromRecord(refSeq, *in.record, DA)
		} else {
			t, err = parseCSVRow(in.row)
		}
		if err != nil {
			return rankedTarget{}, err
		}
		if t.ambCount > threshtarg || slices.Contains(ignoreArray, t.id) {
			return rankedTarget{skip: true}, nil
		}
		rt := rankedTarget{line: t, directions: make([]int, len(queries)), distance: make([]int, len(queries))}
		for i, q := range queries {
			// direction values 0,1,2,3 = same,up,down,side; distance is -1 if the pair fails the ambiguity threshold
			rt.directions[i], rt.distance[i] = whichWay(q, t, threshpair)
		}
		return rt, nil
	}

	h := pipeline.NewHandler(strict, os.Stderr, pipeline.Noun{One: "target", Many: "targets"}, func(rt rankedTarget) error {
		if rt.skip {
			return nil
		}
		for i := range searches {
			if rt.distance[i] >= 0 {
				searches[i].add(rt.line, rt.directions[i], rt.distance[i])
			}
		}
		return nil
	})
	if err := pipeline.Run(context.Background(), threads, next, work, h.Sink); err != nil {
		return err
	}
	if err := h.Finish(); err != nil {
		return err
	}

	QResultsArray := make([]updownCatchmentStruct, len(queries))
	for i := range searches {
		QResultsArray[i] = searches[i].finish()
	}

	if table {
		return writeUpdownTable(out, QResultsArray)
	}
	return writeUpDownCatchment(out, QResultsArray)
}
