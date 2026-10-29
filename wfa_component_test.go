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
