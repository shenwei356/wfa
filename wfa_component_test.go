package wfa

import "testing"

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
