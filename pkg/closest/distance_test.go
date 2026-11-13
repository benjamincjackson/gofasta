package closest

import (
	"math"
	"math/rand"
	"testing"

	"github.com/virus-evolution/gofasta/pkg/encoding"
	"github.com/virus-evolution/gofasta/pkg/fasta"
)

// Straightforward branching versions of the distance measures, to check the branch-free ones against

func refRawDistance(query, target fasta.EncodedRecord) float64 {
	n := 0
	d := 0
	for i, tNuc := range target.Seq {
		if (query.Seq[i] & tNuc) < 16 {
			n += 1
			d += 1
		}
		if (query.Seq[i]&8 == 8) && query.Seq[i] == tNuc {
			d += 1
		}
	}
	return float64(n) / float64(d)
}

func refSnpDistance(query, target fasta.EncodedRecord) float64 {
	n := 0
	for i, tNuc := range target.Seq {
		if (query.Seq[i] & tNuc) < 16 {
			n += 1
		}
	}
	return float64(n)
}

func refTn93Counts(query, target fasta.EncodedRecord) (count_L, count_d, count_P1, count_P2 int) {
	for i, tNuc := range target.Seq {
		if (query.Seq[i]&tNuc) < 16 && query.Seq[i]&8 == 8 && tNuc&8 == 8 {
			count_d++
			count_L++
			if (query.Seq[i] | tNuc) == 200 {
				count_P1++
			} else if (query.Seq[i] | tNuc) == 56 {
				count_P2++
			}
		} else if query.Seq[i]&8 == 8 && query.Seq[i] == tNuc {
			count_L++
		}
	}
	return
}

// allCodes returns every byte the encodings use, plus 0 (an unset entry)
func allCodes() []byte {
	seen := map[byte]bool{0: true}
	codes := []byte{0}
	for _, ea := range [][256]byte{encoding.MakeEncodingArray(), encoding.MakeEncodingArrayHardGaps()} {
		for _, c := range ea {
			if !seen[c] {
				seen[c] = true
				codes = append(codes, c)
			}
		}
	}
	return codes
}

func record(seq []byte) fasta.EncodedRecord {
	r := fasta.EncodedRecord{Seq: seq}
	for _, b := range seq {
		switch b {
		case 136:
			r.Count_A++
		case 40:
			r.Count_C++
		case 72:
			r.Count_G++
		case 24:
			r.Count_T++
		}
	}
	return r
}

// randomPair returns two related sequences of length n with Ns, ambiguity codes and gaps
func randomPair(rng *rand.Rand, n int) (fasta.EncodedRecord, fasta.EncodedRecord) {
	bases := []byte{136, 40, 72, 24}
	codes := allCodes()
	q, t := make([]byte, n), make([]byte, n)
	for i := range q {
		q[i] = bases[rng.Intn(4)]
		t[i] = q[i]
		switch r := rng.Float64(); {
		case r < 0.05:
			t[i] = bases[rng.Intn(4)]
		case r < 0.07:
			t[i] = codes[rng.Intn(len(codes))]
		case r < 0.09:
			q[i] = codes[rng.Intn(len(codes))]
		}
	}
	for _, s := range [][]byte{q, t} {
		start := rng.Intn(n / 2)
		for i := start; i < start+rng.Intn(n/10); i++ {
			s[i] = 240
		}
	}
	return record(q), record(t)
}

func sameFloat(a, b float64) bool {
	return a == b || (math.IsNaN(a) && math.IsNaN(b))
}

// Every pair of codes, one site at a time
func TestDistancesAllCodePairs(t *testing.T) {
	codes := allCodes()
	for _, a := range codes {
		for _, b := range codes {
			q, tg := record([]byte{a}), record([]byte{b})
			if got, want := rawDistance(q, tg), refRawDistance(q, tg); !sameFloat(got, want) {
				t.Errorf("raw %d %d: got %v, want %v", a, b, got, want)
			}
			if got, want := snpDistance(q, tg), refSnpDistance(q, tg); got != want {
				t.Errorf("snp %d %d: got %v, want %v", a, b, got, want)
			}
			L, d, p1, p2 := tn93Counts(q.Seq, tg.Seq)
			rL, rd, rp1, rp2 := refTn93Counts(q, tg)
			if L != rL || d != rd || p1 != rp1 || p2 != rp2 {
				t.Errorf("tn93 %d %d: got %d %d %d %d, want %d %d %d %d", a, b, L, d, p1, p2, rL, rd, rp1, rp2)
			}
		}
	}
}

// Every pair of codes in every lane of an eight-site block, with random codes in the other lanes
func TestDistancesAllCodePairsInBlocks(t *testing.T) {
	rng := rand.New(rand.NewSource(1))
	codes := allCodes()
	for _, a := range codes {
		for _, b := range codes {
			for lane := 0; lane < 8; lane++ {
				qs, ts := make([]byte, 16), make([]byte, 16)
				for i := range qs {
					qs[i], ts[i] = codes[rng.Intn(len(codes))], codes[rng.Intn(len(codes))]
				}
				qs[lane], ts[lane] = a, b
				q, tg := record(qs), record(ts)
				if got, want := rawDistance(q, tg), refRawDistance(q, tg); !sameFloat(got, want) {
					t.Errorf("raw %d %d lane %d: got %v, want %v", a, b, lane, got, want)
				}
				if got, want := snpDistance(q, tg), refSnpDistance(q, tg); got != want {
					t.Errorf("snp %d %d lane %d: got %v, want %v", a, b, lane, got, want)
				}
				L, d, p1, p2 := tn93Counts(q.Seq, tg.Seq)
				rL, rd, rp1, rp2 := refTn93Counts(q, tg)
				if L != rL || d != rd || p1 != rp1 || p2 != rp2 {
					t.Errorf("tn93 %d %d lane %d: got %d %d %d %d, want %d %d %d %d", a, b, lane, L, d, p1, p2, rL, rd, rp1, rp2)
				}
			}
		}
	}
}

// Long random sequences, including lengths that are not a multiple of any block size
func TestDistancesRandomSequences(t *testing.T) {
	rng := rand.New(rand.NewSource(1))
	for _, n := range []int{1, 7, 8, 9, 255, 256, 1000, 29903} {
		for rep := 0; rep < 5; rep++ {
			q, tg := randomPair(rng, max(n, 10))
			q.Seq, tg.Seq = q.Seq[:n], tg.Seq[:n]
			if got, want := rawDistance(q, tg), refRawDistance(q, tg); !sameFloat(got, want) {
				t.Errorf("raw, length %d: got %v, want %v", n, got, want)
			}
			if got, want := snpDistance(q, tg), refSnpDistance(q, tg); got != want {
				t.Errorf("snp, length %d: got %v, want %v", n, got, want)
			}
			L, d, p1, p2 := tn93Counts(q.Seq, tg.Seq)
			rL, rd, rp1, rp2 := refTn93Counts(q, tg)
			if L != rL || d != rd || p1 != rp1 || p2 != rp2 {
				t.Errorf("tn93, length %d: got %d %d %d %d, want %d %d %d %d", n, L, d, p1, p2, rL, rd, rp1, rp2)
			}
		}
	}
}

var benchSink float64

func benchmarkDistance(b *testing.B, f func(fasta.EncodedRecord, fasta.EncodedRecord) float64) {
	q, t := randomPair(rand.New(rand.NewSource(1)), 29903)
	b.SetBytes(29903)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		benchSink = f(q, t)
	}
}

func BenchmarkRawDistance(b *testing.B)     { benchmarkDistance(b, rawDistance) }
func BenchmarkSnpDistance(b *testing.B)     { benchmarkDistance(b, snpDistance) }
func BenchmarkTn93Distance(b *testing.B)    { benchmarkDistance(b, tn93Distance) }
func BenchmarkRefRawDistance(b *testing.B)  { benchmarkDistance(b, refRawDistance) }
func BenchmarkRefSnpDistance(b *testing.B)  { benchmarkDistance(b, refSnpDistance) }
func BenchmarkRefTn93Distance(b *testing.B) { benchmarkDistance(b, refTn93Distance) }

func refTn93Distance(query, target fasta.EncodedRecord) float64 {
	L, d, p1, p2 := refTn93Counts(query, target)
	return tn93FromCounts(query, target, L, d, p1, p2)
}
