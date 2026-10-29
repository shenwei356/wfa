package wfa

import (
	"bytes"
	"io"
	"strings"
	"testing"
)

// TestPlotSparseComponents covers high M scores with shorter I/D storage.
func TestPlotSparseComponents(t *testing.T) {
	for _, global := range []bool{false, true} {
		p := &Penalties{Mismatch: 3000, GapOpen: 4000, GapExt: 1}
		a := New(p, &Options{GlobalAlignment: global})
		q, target := []byte("AC"), []byte("CC")
		r, err := a.Align(q, target)
		if err != nil {
			t.Fatal(err)
		}
		for _, c := range []*Component{a.M, a.I, a.D} {
			a.Plot(&q, &target, io.Discard, c, false, -1)
		}
		RecycleAlignmentResult(r)
		RecycleAligner(a)
	}
}

// TestPlotGapComponents keeps insertion/deletion tables free of match extension.
func TestPlotGapComponents(t *testing.T) {
	for _, swap := range []bool{false, true} {
		a := New(DefaultPenalties, DefaultOptions)
		q, target := []byte("AAAA"), []byte("AAACAAAA")
		if swap {
			q, target = target, q
		}
		r, err := a.Align(q, target)
		if err != nil {
			t.Fatal(err)
		}
		for _, c := range []*Component{a.I, a.D} {
			var out bytes.Buffer
			a.Plot(&q, &target, &out, c, false, -1)
			if strings.Contains(out.String(), "⬊") {
				t.Errorf("gap component rendered matches:\n%s", out.String())
			}
		}
		RecycleAlignmentResult(r)
		RecycleAligner(a)
	}
}
