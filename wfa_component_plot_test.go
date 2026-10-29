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

// TestPlotMatchExtensions distinguishes leading gaps outside the visible matrix
// from matches after those gaps, while retaining real terminal-gap symbols.
func TestPlotMatchExtensions(t *testing.T) {
	for _, tc := range []struct {
		q, target string
		row, col  int
		want      string
	}{
		{"A", "CA", 1, 2, "⬊ 8"},
		{"CA", "A", 2, 1, "⬊ 8"},
		{"T", "TG", 1, 1, "⬊ 0"}, // Reverse complements of A/CA.
		{"TG", "T", 1, 1, "⬊ 0"},
		{"AC", "CAC", 1, 2, "⬊ 8"},
		{"CAC", "AC", 2, 1, "⬊ 8"},
		{"GAAA", "CAAA", 1, 1, "⬂ 4"},
		{"GAAA", "CAAA", 2, 2, "⬊ 4"},
		{"AC", "ACCC", 2, 3, "⟼ 8"},
		{"AC", "ACCC", 2, 4, "🠦10"},
	} {
		t.Run(tc.q+"/"+tc.target+"/"+tc.want, func(t *testing.T) {
			a := New(DefaultPenalties, DefaultOptions)
			defer RecycleAligner(a)
			q, target := []byte(tc.q), []byte(tc.target)
			r, err := a.Align(q, target)
			if err != nil {
				t.Fatal(err)
			}
			defer RecycleAlignmentResult(r)
			var out bytes.Buffer
			a.Plot(&q, &target, &out, a.M, false, -1)
			lines := strings.Split(out.String(), "\n")
			if got := strings.Split(lines[tc.row+1], "\t")[tc.col+1]; got != tc.want {
				t.Fatalf("cell (%d,%d) = %q, want %q:\n%s", tc.row, tc.col, got, tc.want, out.String())
			}
		})
	}
}

// TestPlotTransitionSymbolsAndScoreLimit preserves the optional raw transition
// view and excludes wavefronts above the requested score limit.
func TestPlotTransitionSymbolsAndScoreLimit(t *testing.T) {
	a := New(DefaultPenalties, DefaultOptions)
	defer RecycleAligner(a)
	q, target := []byte("AC"), []byte("CAC")
	r, err := a.Align(q, target)
	if err != nil {
		t.Fatal(err)
	}
	defer RecycleAlignmentResult(r)
	for _, limit := range []int{-1, 4} {
		var out bytes.Buffer
		a.Plot(&q, &target, &out, a.M, true, limit)
		lines := strings.Split(out.String(), "\n")
		want := "⟼ 8"
		if limit == 4 {
			want = "  ."
		}
		for row, col := 1, 2; row <= 2; row, col = row+1, col+1 {
			if got := strings.Split(lines[row+1], "\t")[col+1]; got != want {
				t.Fatalf("limit=%d cell(%d,%d) = %q, want %q", limit, row, col, got, want)
			}
		}
	}
}
