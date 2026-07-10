package wfa

import (
	"math/rand"
	"testing"
)

// generateRandomSeq generates a random DNA sequence
func generateRandomSeq(length int, seed int64) []byte {
	bases := []byte{'A', 'C', 'G', 'T'}
	seq := make([]byte, length)
	r := rand.New(rand.NewSource(seed))
	for i := 0; i < length; i++ {
		seq[i] = bases[r.Intn(4)]
	}
	return seq
}

func BenchmarkAlign_100bp(b *testing.B) {
	query := generateRandomSeq(100, 1)
	target := generateRandomSeq(100, 2)

	algn := New(DefaultPenalties, DefaultOptions)
	defer RecycleAligner(algn)

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_, err := algn.Align(query, target)
		if err != nil {
			b.Fatal(err)
		}
	}
}

func BenchmarkAlign_500bp(b *testing.B) {
	query := generateRandomSeq(500, 1)
	target := generateRandomSeq(500, 2)

	algn := New(DefaultPenalties, DefaultOptions)
	defer RecycleAligner(algn)

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_, err := algn.Align(query, target)
		if err != nil {
			b.Fatal(err)
		}
	}
}

func BenchmarkAlign_1000bp(b *testing.B) {
	query := generateRandomSeq(1000, 1)
	target := generateRandomSeq(1000, 2)

	algn := New(DefaultPenalties, DefaultOptions)
	defer RecycleAligner(algn)

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_, err := algn.Align(query, target)
		if err != nil {
			b.Fatal(err)
		}
	}
}

func BenchmarkAlign_WithAdaptive_1000bp(b *testing.B) {
	query := generateRandomSeq(1000, 1)
	target := generateRandomSeq(1000, 2)

	algn := New(DefaultPenalties, DefaultOptions)
	algn.AdaptiveReduction(DefaultAdaptiveOption)
	defer RecycleAligner(algn)

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_, err := algn.Align(query, target)
		if err != nil {
			b.Fatal(err)
		}
	}
}
