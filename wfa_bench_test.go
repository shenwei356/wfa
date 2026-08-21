package wfa

import (
	"fmt"
	"math/rand"
	"testing"
)

func generateRandomSeq(length int, seed int64) []byte {
	bases := []byte{'A', 'C', 'G', 'T'}
	seq := make([]byte, length)
	r := rand.New(rand.NewSource(seed))
	for i := range seq {
		seq[i] = bases[r.Intn(len(bases))]
	}
	return seq
}

func mutateSubstitutions(src []byte, mismatchRate float64, seed int64) []byte {
	target := append([]byte(nil), src...)
	r := rand.New(rand.NewSource(seed))
	for i, base := range target {
		if r.Float64() >= mismatchRate {
			continue
		}
		for target[i] == base {
			target[i] = "ACGT"[r.Intn(4)]
		}
	}
	return target
}

func benchmarkAlign(b *testing.B, length int, mismatchRate float64, adaptive bool) {
	query := generateRandomSeq(length, 1)
	target := mutateSubstitutions(query, mismatchRate, 2)
	algn := New(DefaultPenalties, DefaultOptions)
	if adaptive {
		if err := algn.AdaptiveReduction(DefaultAdaptiveOption); err != nil {
			b.Fatal(err)
		}
	}
	defer RecycleAligner(algn)

	// Warm the aligner-owned arenas so allocations report steady-state use.
	result, err := algn.Align(query, target)
	if err != nil {
		b.Fatal(err)
	}
	RecycleAlignmentResult(result)

	b.ReportAllocs()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		result, err = algn.Align(query, target)
		if err != nil {
			b.Fatal(err)
		}
		RecycleAlignmentResult(result)
	}
}

func BenchmarkAlign(b *testing.B) {
	for _, length := range []int{100, 1000} {
		for _, mismatchRate := range []float64{0.05, 0.10, 0.20} {
			name := fmt.Sprintf("%dbp/mismatch-%02d", length, int(mismatchRate*100))
			b.Run(name, func(b *testing.B) {
				benchmarkAlign(b, length, mismatchRate, false)
			})
		}
	}
}

func BenchmarkAlignAdaptive(b *testing.B) {
	benchmarkAlign(b, 1000, 0.20, true)
}
