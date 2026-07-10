package wfa

import (
	"testing"
)

// TestCIGARBuffer verifies the poolBytesBuffer fix
func TestCIGARBuffer(t *testing.T) {
	algn := New(DefaultPenalties, DefaultOptions)
	defer RecycleAligner(algn)

	query := []byte("ACCATACTCG")
	target := []byte("AGGATGCTCG")

	result, err := algn.Align(query, target)
	if err != nil {
		t.Fatal(err)
	}

	cigar1 := result.CIGAR(false)

	// Get CIGAR again to ensure buffer is properly reset
	cigar2 := result.CIGAR(false)

	if cigar1 != cigar2 {
		t.Errorf("CIGAR mismatch: first=%s, second=%s", cigar1, cigar2)
	}

	// Verify CIGAR format
	expected := "1M2X2M1X4M"
	if cigar1 != expected {
		t.Errorf("Expected CIGAR %s, got %s", expected, cigar1)
	}
}

// TestEmptySequences tests alignment with empty sequences
func TestEmptySequences(t *testing.T) {
	algn := New(DefaultPenalties, DefaultOptions)
	defer RecycleAligner(algn)

	tests := []struct {
		name   string
		query  []byte
		target []byte
		expect error
	}{
		{"empty query", []byte(""), []byte("ACGT"), ErrEmptySeq},
		{"empty target", []byte("ACGT"), []byte(""), ErrEmptySeq},
		{"both empty", []byte(""), []byte(""), ErrEmptySeq},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			_, err := algn.Align(tt.query, tt.target)
			if err != tt.expect {
				t.Errorf("Expected error %v, got %v", tt.expect, err)
			}
		})
	}
}

// TestSingleCharacter tests alignment with single character sequences
func TestSingleCharacter(t *testing.T) {
	algn := New(DefaultPenalties, DefaultOptions)
	defer RecycleAligner(algn)

	tests := []struct {
		name   string
		query  []byte
		target []byte
	}{
		{"same char", []byte("A"), []byte("A")},
		{"diff char", []byte("A"), []byte("C")},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			result, err := algn.Align(tt.query, tt.target)
			if err != nil {
				t.Fatal(err)
			}
			if result == nil {
				t.Error("Expected non-nil result")
			}
		})
	}
}

// TestIdenticalSequences tests alignment of identical sequences
func TestIdenticalSequences(t *testing.T) {
	algn := New(DefaultPenalties, DefaultOptions)
	defer RecycleAligner(algn)

	sequences := [][]byte{
		[]byte("A"),
		[]byte("ACGT"),
		[]byte("AAAAAAAA"),
		[]byte("ACGTACGTACGT"),
	}

	for _, seq := range sequences {
		t.Run(string(seq), func(t *testing.T) {
			result, err := algn.Align(seq, seq)
			if err != nil {
				t.Fatal(err)
			}
			// Identical sequences should have score 0 (no mismatches or gaps)
			if result.Score != 0 {
				t.Errorf("Expected score 0 for identical sequences, got %d", result.Score)
			}
			if result.Matches != uint32(len(seq)) {
				t.Errorf("Expected %d matches, got %d", len(seq), result.Matches)
			}
		})
	}
}

// TestCompletelyDifferent tests alignment of completely different sequences
func TestCompletelyDifferent(t *testing.T) {
	algn := New(DefaultPenalties, DefaultOptions)
	defer RecycleAligner(algn)

	query := []byte("AAAA")
	target := []byte("CCCC")

	result, err := algn.Align(query, target)
	if err != nil {
		t.Fatal(err)
	}

	// All positions should be mismatches
	if result.Matches != 0 {
		t.Errorf("Expected 0 matches, got %d", result.Matches)
	}
	expectedScore := DefaultPenalties.Mismatch * uint32(len(query))
	if result.Score != expectedScore {
		t.Errorf("Expected score %d, got %d", expectedScore, result.Score)
	}
}

// TestAlignerReuse tests that Aligner can be reused correctly
func TestAlignerReuse(t *testing.T) {
	algn := New(DefaultPenalties, DefaultOptions)
	defer RecycleAligner(algn)

	// First alignment
	result1, err := algn.Align([]byte("ACGT"), []byte("ACGT"))
	if err != nil {
		t.Fatal(err)
	}
	score1 := result1.Score

	// Second alignment with different sequences
	result2, err := algn.Align([]byte("AAAA"), []byte("CCCC"))
	if err != nil {
		t.Fatal(err)
	}
	score2 := result2.Score

	// Scores should be different
	if score1 == score2 {
		t.Error("Expected different scores for different alignments")
	}

	// Third alignment - same as first
	result3, err := algn.Align([]byte("ACGT"), []byte("ACGT"))
	if err != nil {
		t.Fatal(err)
	}
	score3 := result3.Score

	if score1 != score3 {
		t.Errorf("Expected same score for same sequences: %d vs %d", score1, score3)
	}
}

// TestLongSequences tests alignment with longer sequences
func TestLongSequences(t *testing.T) {
	algn := New(DefaultPenalties, DefaultOptions)
	defer RecycleAligner(algn)

	query := generateRandomSeq(200, 42)
	target := generateRandomSeq(200, 43)

	result, err := algn.Align(query, target)
	if err != nil {
		t.Fatal(err)
	}

	if result == nil {
		t.Error("Expected non-nil result")
	}

	// Verify CIGAR can be generated
	cigar := result.CIGAR(false)
	if len(cigar) == 0 {
		t.Error("Expected non-empty CIGAR string")
	}
}

// TestComponentReset verifies that Component.Reset() works correctly
func TestComponentReset(t *testing.T) {
	cpt := NewComponent()
	defer RecycleComponent(cpt)

	// Add some wavefronts
	cpt.Set(0, 0, 1, wfaMatch)
	cpt.Set(4, 2, 5, wfaMismatch)
	cpt.Set(10, -3, 10, wfaInsertOpen)

	// Verify they exist
	if !cpt.HasScore(0) || !cpt.HasScore(4) || !cpt.HasScore(10) {
		t.Error("Expected scores to exist before reset")
	}

	// Reset
	cpt.Reset()

	// Verify they're gone
	if cpt.HasScore(0) || cpt.HasScore(4) || cpt.HasScore(10) {
		t.Error("Expected scores to be cleared after reset")
	}

	// Verify we can reuse after reset
	cpt.Set(0, 0, 1, wfaMatch)
	if !cpt.HasScore(0) {
		t.Error("Expected to be able to add scores after reset")
	}
}

// TestWaveFrontArena verifies the arena allocation works correctly
func TestWaveFrontArena(t *testing.T) {
	cpt := NewComponent()
	defer RecycleComponent(cpt)

	// Allocate many wavefronts to trigger chunk growth
	for i := uint32(0); i < 1000; i++ {
		cpt.Set(i, 0, 1, wfaMatch)
	}

	// Verify all were allocated
	for i := uint32(0); i < 1000; i++ {
		if !cpt.HasScore(i) {
			t.Errorf("Score %d not found after allocation", i)
		}
	}

	// Reset and verify arena reuse
	cpt.Reset()

	// Allocate again - should reuse arena chunks
	for i := uint32(0); i < 500; i++ {
		cpt.Set(i, 0, 1, wfaMatch)
	}

	for i := uint32(0); i < 500; i++ {
		if !cpt.HasScore(i) {
			t.Errorf("Score %d not found after arena reuse", i)
		}
	}
}
