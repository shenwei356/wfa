package wfa

import "testing"

// TestTrimWithoutMatches keeps empty aligned regions safe for all renderers.
func TestTrimWithoutMatches(t *testing.T) {
	for _, global := range []bool{true, false} {
		a := New(DefaultPenalties, &Options{GlobalAlignment: global})
		q, s := []byte("A"), []byte("C")
		r, err := a.Align(q, s)
		if err != nil {
			t.Fatal(err)
		}
		if got := r.CIGAR(true); got != "" {
			t.Errorf("trimmed CIGAR=%q", got)
		}
		Q, A, T := r.AlignmentText(&q, &s, true)
		if len(*Q) != 0 || len(*A) != 0 || len(*T) != 0 {
			t.Error("nonempty trimmed alignment")
		}
		RecycleAlignmentText(Q, A, T)
		if r.CIGAR(false) != "1X" || r.QBegin != 0 || r.TBegin != 0 {
			t.Error("untrimmed result changed")
		}
		RecycleAlignmentResult(r)
		RecycleAligner(a)
	}
}

// TestCIGARLargeOperationCount must format uint32 counts without signed-int conversion.
func TestCIGARLargeOperationCount(t *testing.T) {
	r := NewAlignmentResult(true)
	defer RecycleAlignmentResult(r)
	r.AddN('M', ^uint32(0))
	if got := r.CIGAR(false); got != "4294967295M" {
		t.Fatalf("CIGAR = %s", got)
	}
}
