package wfa

import (
	"fmt"
	"sync"
	"testing"
)

// TestWaveFrontConfiguredGrowth exercises both setters with blocks larger and
// smaller than the default, including positive and negative diagonal indexes.
func TestWaveFrontConfiguredGrowth(t *testing.T) {
	oldBase := OFFSETS_BASE_SIZE
	defer func() { OFFSETS_BASE_SIZE = oldBase }()
	for _, base := range []int{16, 128, 512} {
		OFFSETS_BASE_SIZE = base
		for _, raw := range []bool{false, true} {
			t.Run(fmt.Sprintf("base=%d/raw=%t", base, raw), func(t *testing.T) {
				w := NewWaveFront()
				defer RecycleWaveFront(w)
				for _, k := range []int{-500, 0, 500} {
					if raw {
						w.SetRaw(k, 123<<wfaTypeBits|wfaMatch)
					} else {
						w.Set(k, 123, wfaMatch)
					}
					if offset, typ, ok := w.Get(k); !ok || offset != 123 || typ != wfaMatch {
						t.Fatalf("diagonal %d = %d,%d,%v", k, offset, typ, ok)
					}
				}
			})
		}
	}
}

// TestWaveFrontPoolLargerBase checks checkout of a smaller pooled wavefront.
func TestWaveFrontPoolLargerBase(t *testing.T) {
	oldBase, oldPool := OFFSETS_BASE_SIZE, poolWaveFront
	defer func() { OFFSETS_BASE_SIZE, poolWaveFront = oldBase, oldPool }()
	w := &WaveFront{Offsets: []uint32{42}}
	poolWaveFront = &sync.Pool{New: func() interface{} { return w }}
	OFFSETS_BASE_SIZE = 512
	got := NewWaveFront()
	if len(got.Offsets) != 512 || got.Offsets[0] != 0 {
		t.Fatal("pooled wavefront was not resized and cleared")
	}
}

// TestWaveFrontGrowthClearsReusedOffsets checks stale diagonal holes both when
// growth reuses capacity and when it allocates a larger backing array.
func TestWaveFrontGrowthClearsReusedOffsets(t *testing.T) {
	for _, k := range []int{200, 500} {
		for _, raw := range []bool{false, true} {
			t.Run(fmt.Sprintf("k=%d/raw=%t", k, raw), func(t *testing.T) {
				offsets := make([]uint32, 512)
				for i := 128; i < len(offsets); i++ {
					offsets[i] = 99<<wfaTypeBits | wfaMatch
				}
				w := &WaveFront{Offsets: offsets[:128]}
				w.Set(0, 7, wfaMismatch)
				if raw {
					w.SetRaw(k, 123<<wfaTypeBits|wfaMatch)
				} else {
					w.Set(k, 123, wfaMatch)
				}
				for diagonal := 1; diagonal < k; diagonal++ {
					if _, _, ok := w.Get(diagonal); ok {
						t.Fatalf("stale diagonal %d reappeared after growth", diagonal)
					}
				}
				if offset, typ, ok := w.Get(0); !ok || offset != 7 || typ != wfaMismatch {
					t.Fatal("growth changed the existing offset")
				}
				if offset, typ, ok := w.Get(k); !ok || offset != 123 || typ != wfaMatch {
					t.Fatal("growth lost the new offset")
				}
			})
		}
	}
}
