package wfa

import (
	"sync"
	"testing"
)

// TestConcurrentAlignments tests that multiple goroutines can use different Aligners safely
func TestConcurrentAlignments(t *testing.T) {
	numGoroutines := 100
	alignmentsPerGoroutine := 10

	var wg sync.WaitGroup
	errors := make(chan error, numGoroutines*alignmentsPerGoroutine)

	for i := 0; i < numGoroutines; i++ {
		wg.Add(1)
		go func(id int) {
			defer wg.Done()

			// Each goroutine gets its own Aligner
			algn := New(DefaultPenalties, DefaultOptions)
			defer RecycleAligner(algn)

			// Generate different sequences for each goroutine
			query := generateRandomSeq(100, int64(id*2))
			target := generateRandomSeq(100, int64(id*2+1))

			for j := 0; j < alignmentsPerGoroutine; j++ {
				result, err := algn.Align(query, target)
				if err != nil {
					errors <- err
					return
				}

				// Verify result is not nil
				if result == nil {
					errors <- err
					return
				}

				// Verify CIGAR can be generated
				cigar := result.CIGAR(false)
				if len(cigar) == 0 {
					errors <- err
					return
				}
			}
		}(i)
	}

	wg.Wait()
	close(errors)

	// Check for any errors
	for err := range errors {
		if err != nil {
			t.Fatalf("Concurrent alignment failed: %v", err)
		}
	}
}

// TestConcurrentAlignerPooling tests that Aligners from the pool work correctly under concurrency
func TestConcurrentAlignerPooling(t *testing.T) {
	numGoroutines := 50

	var wg sync.WaitGroup

	for i := 0; i < numGoroutines; i++ {
		wg.Add(1)
		go func(id int) {
			defer wg.Done()

			// Get Aligner from pool
			algn := New(DefaultPenalties, DefaultOptions)

			query := generateRandomSeq(50, int64(id))
			target := generateRandomSeq(50, int64(id+1000))

			result, err := algn.Align(query, target)
			if err != nil {
				t.Errorf("Goroutine %d: alignment failed: %v", id, err)
				return
			}

			if result == nil {
				t.Errorf("Goroutine %d: got nil result", id)
				return
			}

			// Recycle back to pool
			RecycleAligner(algn)
		}(i)
	}

	wg.Wait()
}

// TestConcurrentComponentOperations tests concurrent operations on different Components
func TestConcurrentComponentOperations(t *testing.T) {
	numGoroutines := 50

	var wg sync.WaitGroup

	for i := 0; i < numGoroutines; i++ {
		wg.Add(1)
		go func(id int) {
			defer wg.Done()

			cpt := NewComponent()
			defer RecycleComponent(cpt)

			// Each goroutine performs independent operations
			for s := uint32(0); s < 100; s++ {
				cpt.Set(s, id, uint32(id+int(s)), wfaMatch)
			}

			// Verify
			for s := uint32(0); s < 100; s++ {
				if !cpt.HasScore(s) {
					t.Errorf("Goroutine %d: score %d not found", id, s)
					return
				}
			}

			// Reset
			cpt.Reset()

			// Verify reset worked
			for s := uint32(0); s < 100; s++ {
				if cpt.HasScore(s) {
					t.Errorf("Goroutine %d: score %d still exists after reset", id, s)
					return
				}
			}
		}(i)
	}

	wg.Wait()
}

// TestConcurrentAlignmentWithAdaptive tests concurrent alignments with adaptive reduction
func TestConcurrentAlignmentWithAdaptive(t *testing.T) {
	numGoroutines := 20
	seqLength := 200

	var wg sync.WaitGroup
	errors := make(chan error, numGoroutines)

	for i := 0; i < numGoroutines; i++ {
		wg.Add(1)
		go func(id int) {
			defer wg.Done()

			algn := New(DefaultPenalties, DefaultOptions)
			err := algn.AdaptiveReduction(DefaultAdaptiveOption)
			if err != nil {
				errors <- err
				return
			}
			defer RecycleAligner(algn)

			query := generateRandomSeq(seqLength, int64(id*3))
			target := generateRandomSeq(seqLength, int64(id*3+1))

			result, err := algn.Align(query, target)
			if err != nil {
				errors <- err
				return
			}

			if result == nil {
				errors <- err
				return
			}

			// Verify alignment properties
			if result.AlignLen == 0 {
				errors <- err
				return
			}
		}(i)
	}

	wg.Wait()
	close(errors)

	for err := range errors {
		if err != nil {
			t.Fatalf("Concurrent alignment with adaptive reduction failed: %v", err)
		}
	}
}

// TestConcurrentAlignerReuse tests that Aligners can be safely reused in concurrent scenarios
func TestConcurrentAlignerReuse(t *testing.T) {
	numGoroutines := 30
	reusesPerGoroutine := 5

	var wg sync.WaitGroup

	for i := 0; i < numGoroutines; i++ {
		wg.Add(1)
		go func(id int) {
			defer wg.Done()

			algn := New(DefaultPenalties, DefaultOptions)
			defer RecycleAligner(algn)

			for j := 0; j < reusesPerGoroutine; j++ {
				query := generateRandomSeq(50, int64(id*100+j))
				target := generateRandomSeq(50, int64(id*100+j+50))

				result, err := algn.Align(query, target)
				if err != nil {
					t.Errorf("Goroutine %d, iteration %d: %v", id, j, err)
					return
				}

				if result == nil {
					t.Errorf("Goroutine %d, iteration %d: nil result", id, j)
					return
				}

				// Verify CIGAR
				cigar := result.CIGAR(false)
				if len(cigar) == 0 {
					t.Errorf("Goroutine %d, iteration %d: empty CIGAR", id, j)
					return
				}

				// Recycle the result
				RecycleAlignmentResult(result)
			}
		}(i)
	}

	wg.Wait()
}

// TestConcurrentAlignmentText tests concurrent generation of alignment text
func TestConcurrentAlignmentText(t *testing.T) {
	numGoroutines := 30

	var wg sync.WaitGroup

	for i := 0; i < numGoroutines; i++ {
		wg.Add(1)
		go func(id int) {
			defer wg.Done()

			algn := New(DefaultPenalties, DefaultOptions)
			defer RecycleAligner(algn)

			query := generateRandomSeq(50, int64(id))
			target := generateRandomSeq(50, int64(id+500))

			result, err := algn.Align(query, target)
			if err != nil {
				t.Errorf("Goroutine %d: alignment failed: %v", id, err)
				return
			}

			// Generate alignment text
			Q, A, T := result.AlignmentText(&query, &target, false)
			defer RecycleAlignmentText(Q, A, T)

			// Verify all three slices are not empty
			if len(*Q) == 0 || len(*A) == 0 || len(*T) == 0 {
				t.Errorf("Goroutine %d: empty alignment text", id)
				return
			}

			// Verify they have the same length
			if len(*Q) != len(*A) || len(*Q) != len(*T) {
				t.Errorf("Goroutine %d: alignment text length mismatch: Q=%d, A=%d, T=%d",
					id, len(*Q), len(*A), len(*T))
				return
			}
		}(i)
	}

	wg.Wait()
}

// TestConcurrentStressTest is a stress test with heavy concurrent load
func TestConcurrentStressTest(t *testing.T) {
	if testing.Short() {
		t.Skip("Skipping stress test in short mode")
	}

	numGoroutines := 100
	alignmentsPerGoroutine := 50

	var wg sync.WaitGroup
	successCount := make(chan int, numGoroutines)

	for i := 0; i < numGoroutines; i++ {
		wg.Add(1)
		go func(id int) {
			defer wg.Done()

			algn := New(DefaultPenalties, DefaultOptions)
			defer RecycleAligner(algn)

			successes := 0
			for j := 0; j < alignmentsPerGoroutine; j++ {
				// Vary sequence lengths
				length := 50 + (id+j)%150

				query := generateRandomSeq(length, int64(id*1000+j))
				target := generateRandomSeq(length, int64(id*1000+j+10000))

				result, err := algn.Align(query, target)
				if err != nil {
					t.Errorf("Goroutine %d, iteration %d: %v", id, j, err)
					continue
				}

				if result != nil {
					successes++
					RecycleAlignmentResult(result)
				}
			}

			successCount <- successes
		}(i)
	}

	wg.Wait()
	close(successCount)

	totalSuccesses := 0
	for count := range successCount {
		totalSuccesses += count
	}

	expectedSuccesses := numGoroutines * alignmentsPerGoroutine
	if totalSuccesses != expectedSuccesses {
		t.Errorf("Expected %d successful alignments, got %d", expectedSuccesses, totalSuccesses)
	}

	t.Logf("Stress test completed: %d alignments across %d goroutines", totalSuccesses, numGoroutines)
}
