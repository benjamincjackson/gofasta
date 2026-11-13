package closest

import (
	"errors"
	"io"
	"math"
	"sort"
	"strings"

	"github.com/virus-evolution/gofasta/pkg/fasta"
)

// catchmentStruct contains information about the closest sequences to a particular query
type catchmentStruct struct {
	qname     string
	catchment []resultsStruct // once full, sorted (see rearrangeCatchment), so the last is the least close
}

// rearrangeCatchment sorts a catchmentStruct so that the sequences in its catchment field are in order of
// genetic distance, with ties broken by genome completeness. It is called before the catchment is written to
// output, or if a new sequence is added to a catchmentStruct which is already at capacity
func rearrangeCatchment(nS *catchmentStruct, catchmentSize int) {
	sort.SliceStable(nS.catchment, func(i, j int) bool {
		return lessDistance(nS.catchment[i].distance, nS.catchment[j].distance) || (nS.catchment[i].distance == nS.catchment[j].distance && nS.catchment[i].completeness > nS.catchment[j].completeness)
	})
	nS.catchment = nS.catchment[0:catchmentSize]
}

// consider adds target to the catchment if there is room, or if it is closer than the furthest member
// (or as close and more complete). A distance that is undefined (NaN) or greater than maxdist (unless
// maxdist is -1) never qualifies.
func (neighbours *catchmentStruct) consider(target fasta.EncodedRecord, distance float64, catchmentSize int, maxdist float64) {
	if maxdist != -1.0 {
		// an undefined distance (NaN) can't be said to be within maxdist
		if math.IsNaN(distance) || distance > maxdist {
			return
		}
	}
	rs := resultsStruct{tname: target.ID, completeness: target.Score, distance: distance}
	if len(neighbours.catchment) < catchmentSize {
		neighbours.catchment = append(neighbours.catchment, rs)
		if len(neighbours.catchment) == catchmentSize {
			rearrangeCatchment(neighbours, catchmentSize)
		}
		return
	}
	furthest := neighbours.catchment[catchmentSize-1]
	if lessDistance(distance, furthest.distance) || (distance == furthest.distance && target.Score > furthest.completeness) {
		neighbours.catchment = append(neighbours.catchment, rs)
		rearrangeCatchment(neighbours, catchmentSize)
	}
}

// writeClosestN writes each query's neighbours, closest first, on one line
func writeClosestN(results []catchmentStruct, w io.Writer) error {
	if _, err := io.WriteString(w, "query,closest\n"); err != nil {
		return err
	}
	for _, result := range results {
		names := make([]string, len(result.catchment))
		for i, hit := range result.catchment {
			names[i] = hit.tname
		}
		if _, err := io.WriteString(w, result.qname+","+strings.Join(names, ";")+"\n"); err != nil {
			return err
		}
	}
	return nil
}

// writeClosestNTable writes one line per query and neighbour, with the distance between them
func writeClosestNTable(results []catchmentStruct, w io.Writer, measure string) error {
	if _, err := io.WriteString(w, "query,target,distance\n"); err != nil {
		return err
	}
	for _, result := range results {
		for _, hit := range result.catchment {
			if _, err := io.WriteString(w, result.qname+","+hit.tname+","+formatDistance(hit.distance, measure)+"\n"); err != nil {
				return err
			}
		}
	}
	return nil
}

// ClosestN finds the closest sequence(s) by genetic distance to a query/queries. It writes the results
// to stdout or to file. Ties for distance are broken by genome completeness.
// threads is the number of goroutines that process records (fewer than 1 means 1).
func ClosestN(catchmentSize int, maxdist float64, query, target io.Reader, measure string, out io.Writer, table bool, threads int, strict bool) error {
	if catchmentSize < 0 {
		return errors.New("the number of neighbours to find must not be negative")
	}

	if maxdist != -1.0 && catchmentSize == 0 {
		catchmentSize = math.MaxInt
	}

	queries, err := loadQueries(query, strict)
	if err != nil {
		return err
	}

	QResultsArray := make([]catchmentStruct, len(queries))
	for i, q := range queries {
		QResultsArray[i] = catchmentStruct{qname: q.ID, catchment: make([]resultsStruct, 0)}
	}

	err = scanTargets(queries, target, measure, threads, strict, func(st scoredTarget) {
		for i := range QResultsArray {
			QResultsArray[i].consider(st.target, st.distances[i], catchmentSize, maxdist)
		}
	})
	if err != nil {
		return err
	}

	// a catchment that never filled up has not been sorted yet
	for i := range QResultsArray {
		if n := len(QResultsArray[i].catchment); n < catchmentSize && n > 0 {
			rearrangeCatchment(&QResultsArray[i], n)
		}
	}

	if table {
		return writeClosestNTable(QResultsArray, out, measure)
	}
	return writeClosestN(QResultsArray, out)
}
