package wfa

import (
	"sync"
	"testing"
)

// TestComponentSparseScore checks setters beyond more than one growth block.
func TestComponentSparseScore(t *testing.T) {
	for _, raw := range []bool{false, true} {
		c := NewComponent()
		if raw {
			c.SetRaw(5000, -3, 7<<wfaTypeBits|wfaMatch)
		} else {
			c.Set(5000, -3, 7, wfaMatch)
		}
		offset, typ, ok := c.Get(5000, -3)
		if !ok || offset != 7 || typ != wfaMatch {
			t.Fatalf("got %d %d %v", offset, typ, ok)
		}
		RecycleComponent(c)
	}
}

// TestComponentDeleteMissingScore covers absent scores across allocation blocks.
func TestComponentDeleteMissingScore(t *testing.T) {
	c := NewComponent()
	defer RecycleComponent(c)
	c.Set(0, 0, 7, wfaMatch)
	for _, s := range []uint32{1, uint32(len(c.WaveFronts)), uint32(len(c.WaveFronts) * 2), ^uint32(0)} {
		c.Delete(s, 0)
	}
	if offset, typ, ok := c.Get(0, 0); !ok || offset != 7 || typ != wfaMatch {
		t.Fatalf("deleting absent scores changed score zero: %d %d %v", offset, typ, ok)
	}
	c.Delete(0, 0)
	if _, _, ok := c.Get(0, 0); ok {
		t.Fatal("existing offset was not deleted")
	}
}

// TestComponentReadsMissingScores checks both shifted readers at absent, sparse,
// and underflowing score positions without allocating additional wavefronts.
func TestComponentReadsMissingScores(t *testing.T) {
	c := NewComponent()
	defer RecycleComponent(c)
	c.Set(7, -2, 42, wfaMismatch)
	for _, tc := range []struct {
		s, diff uint32
		exists  bool
	}{
		{8, 1, true}, {7, 0, true}, {6, 0, false}, {6, 7, false},
		{uint32(len(c.WaveFronts)), 0, false}, {5000, 1, false}, {^uint32(0), 0, false},
	} {
		offset, typ, ok := c.GetAfterDiff(tc.s, tc.diff, -2)
		if ok != tc.exists || ok && (offset != 42 || typ != wfaMismatch) {
			t.Errorf("GetAfterDiff(%d,%d) = %d,%d,%v", tc.s, tc.diff, offset, typ, ok)
		}
		raw, ok := c.GetRawAfterDiff(tc.s, tc.diff, -2)
		if ok != tc.exists || ok && raw != 42<<wfaTypeBits|wfaMismatch {
			t.Errorf("GetRawAfterDiff(%d,%d) = %d,%v", tc.s, tc.diff, raw, ok)
		}
	}
}

// TestComponentIncreaseExistingOnly ensures absent scores and diagonal holes
// cannot panic or create offsets without a traceback type.
func TestComponentIncreaseExistingOnly(t *testing.T) {
	c := NewComponent()
	defer RecycleComponent(c)
	c.Set(7, -2, 42, wfaMismatch)
	c.Set(7, 2, 21, wfaMatch)
	length := len(c.WaveFronts)
	for _, s := range []uint32{0, 1, uint32(length), 5000, ^uint32(0)} {
		c.Increase(s, -2, 3)
		if c.HasScore(s) {
			t.Errorf("created missing score %d", s)
		}
	}
	for _, k := range []int{-1000, 0, 1000} {
		c.Increase(7, k, 3)
		if _, _, ok := c.Get(7, k); ok {
			t.Errorf("created missing diagonal %d", k)
		}
	}
	c.Increase(7, -2, 3)
	if offset, typ, ok := c.Get(7, -2); !ok || offset != 45 || typ != wfaMismatch {
		t.Fatalf("increased offset = %d,%d,%v", offset, typ, ok)
	}
	if len(c.WaveFronts) != length {
		t.Fatal("missing scores grew the component")
	}
}

// TestComponentPoolLargerBase simulates a component created before its initial
// allocation setting was increased. The injected pool also works across GC.
func TestComponentPoolLargerBase(t *testing.T) {
	oldBase, oldPool := WAVEFRONTS_BASE_SIZE, poolComponent
	defer func() { WAVEFRONTS_BASE_SIZE, poolComponent = oldBase, oldPool }()
	c := &Component{WaveFronts: make([]*WaveFront, 16)}
	c.Set(7, 0, 3, wfaMatch)
	poolComponent = &sync.Pool{New: func() interface{} { return c }}
	WAVEFRONTS_BASE_SIZE = 128
	got := NewComponent()
	if len(got.WaveFronts) != 128 || got.HasScore(7) {
		t.Fatal("pooled component was not resized and cleared")
	}
}

// TestComponentArenaLargerBase checks arena reuse after offset configuration
// changes and verifies that stale offsets stay absent when capacity is reused.
func TestComponentArenaLargerBase(t *testing.T) {
	oldBase := OFFSETS_BASE_SIZE
	defer func() { OFFSETS_BASE_SIZE = oldBase }()
	OFFSETS_BASE_SIZE = 16
	c := NewComponent()
	defer RecycleComponent(c)
	c.Set(0, 0, 1, wfaMatch)
	OFFSETS_BASE_SIZE = 512
	c.Reset()
	c.Set(0, 500, 8, wfaMatch)
	c.Reset()
	c.Set(0, 501, 9, wfaMismatch)
	if _, _, ok := c.Get(0, 500); ok {
		t.Fatal("stale offset survived arena reuse")
	}
	if offset, typ, ok := c.Get(0, 501); !ok || offset != 9 || typ != wfaMismatch {
		t.Fatalf("reused arena offset = %d,%d,%v", offset, typ, ok)
	}
}
