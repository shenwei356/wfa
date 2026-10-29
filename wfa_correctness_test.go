package wfa

import (
	"context"
	"errors"
	"fmt"
	"math/rand"
	"os"
	"os/exec"
	"testing"
	"time"
)

const maxInt = int(^uint(0) >> 1)

// affineGlobalScore is a simple O(nm) Gotoh-style oracle for end-to-end global
// alignment. It intentionally shares no wavefront implementation code.
func affineGlobalScore(query, target []byte, p *Penalties) uint32 {
	cols := len(target) + 1
	m := make([]int, (len(query)+1)*cols)
	ins := make([]int, len(m))
	del := make([]int, len(m))
	for i := range m {
		m[i], ins[i], del[i] = maxInt/4, maxInt/4, maxInt/4
	}
	gapOpen := int(p.GapOpen + p.GapExt)
	gapExt := int(p.GapExt)
	m[0] = 0
	for j := 1; j <= len(target); j++ {
		idx := j
		ins[idx] = gapOpen + (j-1)*gapExt
		m[idx] = ins[idx]
	}
	for i := 1; i <= len(query); i++ {
		idx := i * cols
		del[idx] = gapOpen + (i-1)*gapExt
		m[idx] = del[idx]
	}

	for i := 1; i <= len(query); i++ {
		for j := 1; j <= len(target); j++ {
			idx := i*cols + j
			ins[idx] = min(m[idx-1]+gapOpen, ins[idx-1]+gapExt)
			del[idx] = min(m[idx-cols]+gapOpen, del[idx-cols]+gapExt)
			mismatch := 0
			if query[i-1] != target[j-1] {
				mismatch = int(p.Mismatch)
			}
			m[idx] = min(m[idx-cols-1]+mismatch, ins[idx], del[idx])
		}
	}
	return uint32(m[len(m)-1])
}

func validateGlobalResult(t *testing.T, query, target []byte, p *Penalties, result *AlignmentResult) {
	t.Helper()
	var qi, ti int
	var score uint32
	for _, packed := range result.Ops {
		op, n := Op(packed)
		switch op {
		case 'M':
			for j := uint32(0); j < n; j++ {
				if qi >= len(query) || ti >= len(target) || query[qi] != target[ti] {
					t.Fatalf("CIGAR %s: invalid M operation at query=%d target=%d", result.CIGAR(false), qi, ti)
				}
				qi++
				ti++
			}
		case 'X':
			for j := uint32(0); j < n; j++ {
				if qi >= len(query) || ti >= len(target) || query[qi] == target[ti] {
					t.Fatalf("invalid X operation at query=%d target=%d", qi, ti)
				}
				qi++
				ti++
				score += p.Mismatch
			}
		case 'I':
			ti += int(n)
			score += p.GapOpen + n*p.GapExt
		case 'D':
			qi += int(n)
			score += p.GapOpen + n*p.GapExt
		default:
			t.Fatalf("unexpected global CIGAR operation %q", op)
		}
	}
	if qi != len(query) || ti != len(target) {
		t.Fatalf("CIGAR consumed query=%d/%d target=%d/%d", qi, len(query), ti, len(target))
	}
	if score != result.Score {
		t.Fatalf("CIGAR score=%d, result score=%d", score, result.Score)
	}
}

func checkAgainstAffineOracle(t *testing.T, query, target []byte, p *Penalties) {
	t.Helper()
	algn := New(p, DefaultOptions)
	defer RecycleAligner(algn)
	result, err := algn.Align(query, target)
	if err != nil {
		t.Fatal(err)
	}
	defer RecycleAlignmentResult(result)

	want := affineGlobalScore(query, target, p)
	if result.Score != want {
		t.Fatalf("query=%q target=%q: score=%d, want %d", query, target, result.Score, want)
	}
	t.Run(fmt.Sprintf("q=%s/t=%s", query, target), func(t *testing.T) {
		validateGlobalResult(t, query, target, p, result)
	})
}

func allDNASequences(maxLen int) [][]byte {
	seqs := make([][]byte, 0)
	var add func([]byte, int)
	add = func(prefix []byte, remaining int) {
		if len(prefix) > 0 {
			seqs = append(seqs, append([]byte(nil), prefix...))
		}
		if remaining == 0 {
			return
		}
		for _, base := range []byte("AC") {
			add(append(prefix, base), remaining-1)
		}
	}
	add(nil, maxLen)
	return seqs
}

func TestGlobalAlignmentAgainstAffineOracle(t *testing.T) {
	penalties := []*Penalties{
		DefaultPenalties,
		{Mismatch: 2, GapOpen: 3, GapExt: 1},
		{Mismatch: 3, GapOpen: 7, GapExt: 1},
	}
	seqs := allDNASequences(5)
	for _, p := range penalties {
		t.Run(fmt.Sprintf("x%d-o%d-e%d", p.Mismatch, p.GapOpen, p.GapExt), func(t *testing.T) {
			for _, query := range seqs {
				for _, target := range seqs {
					checkAgainstAffineOracle(t, query, target, p)
				}
			}
		})
	}
}

func TestRandomGlobalAlignmentAgainstAffineOracle(t *testing.T) {
	r := rand.New(rand.NewSource(1))
	for i := 0; i < 1000; i++ {
		query := generateRandomSeq(1+r.Intn(20), r.Int63())
		target := generateRandomSeq(1+r.Intn(20), r.Int63())
		checkAgainstAffineOracle(t, query, target, DefaultPenalties)
	}
}

func TestGlobalAlignmentAllowsLeadingGaps(t *testing.T) {
	tests := []struct {
		query, target string
		cigar         string
	}{
		{query: "A", target: "CA", cigar: "1I1M"},
		{query: "CA", target: "A", cigar: "1D1M"},
	}
	for _, tt := range tests {
		t.Run(tt.query+"/"+tt.target, func(t *testing.T) {
			algn := New(DefaultPenalties, DefaultOptions)
			defer RecycleAligner(algn)
			result, err := algn.Align([]byte(tt.query), []byte(tt.target))
			if err != nil {
				t.Fatal(err)
			}
			defer RecycleAlignmentResult(result)
			if got := result.CIGAR(false); got != tt.cigar {
				t.Fatalf("CIGAR=%s, want %s", got, tt.cigar)
			}
		})
	}
}

func TestGlobalAlignmentWithTerminalGapTie(t *testing.T) {
	query := []byte("Bioinformatics helps Biology")
	target := []byte("We learn bioinformatics to help biologists")
	algn := New(DefaultPenalties, DefaultOptions)
	defer RecycleAligner(algn)
	result, err := algn.Align(query, target)
	if err != nil {
		t.Fatal(err)
	}
	defer RecycleAlignmentResult(result)

	const wantCIGAR = "9I1X14M3I4M1D1M1X5M3I1X"
	if got := result.CIGAR(false); got != wantCIGAR {
		t.Fatalf("CIGAR=%s, want %s", got, wantCIGAR)
	}
	if want := affineGlobalScore(query, target, DefaultPenalties); result.Score != want {
		t.Fatalf("score=%d, want %d", result.Score, want)
	}
	validateGlobalResult(t, query, target, DefaultPenalties, result)
}

func TestAlignerPoolResetsAdaptiveReduction(t *testing.T) {
	algn := New(DefaultPenalties, DefaultOptions)
	if err := algn.AdaptiveReduction(DefaultAdaptiveOption); err != nil {
		t.Fatal(err)
	}
	RecycleAligner(algn)

	algn = New(DefaultPenalties, DefaultOptions)
	defer RecycleAligner(algn)
	if algn.ad != nil {
		t.Fatal("pooled aligner retained adaptive-reduction settings")
	}
}

// TestAdaptiveGlobalTerminalGap checks the fragment that stalled genome compare.
// Isolate the alignment so a regression cannot leave an unbounded score loop.
func TestAdaptiveGlobalTerminalGap(t *testing.T) {
	if os.Getenv("WFA_TEST_TERMINAL_GAP") != "1" {
		exe, err := os.Executable()
		if err != nil {
			t.Fatal(err)
		}
		ctx, cancel := context.WithTimeout(context.Background(), 10*time.Second)
		defer cancel()
		cmd := exec.CommandContext(ctx, exe, "-test.run=^TestAdaptiveGlobalTerminalGap$", "-test.v")
		cmd.Env = append(os.Environ(), "WFA_TEST_TERMINAL_GAP=1")
		output, err := cmd.CombinedOutput()
		if ctx.Err() != nil {
			t.Fatalf("adaptive global alignment did not terminate: %v\n%s", ctx.Err(), output)
		}
		if err != nil {
			t.Fatalf("alignment failed: %v\n%s", err, output)
		}
		return
	}
	query := []byte("CCTGTTGCTCCTGTTGCTCCTGTTGGACCTGTTGCTCCATCTGCTCCTGTATTTCCTGTTGCTCCTGTTGGACCTGTTGCTCCATCTGCTCCTGTATTTCCTGTTGCTCCTGTTGCTCCTGTTGGACCTGCTACTCCATCTGCTCCTGTATTTCCTGTTGCTCCTGTTGGACCTGTTGCTCCATCTGCTCCTGTTGGACCTGTCGCTCCTGTTGCCCCTGTTGCTCCTGTTATTCCTACTCCTGTTGGGCCTGTTGCTCCTGTTGGACCTGTTGCTCCAGGTGTTCCTGCTGCTCCTGTTGCTCCTGTTGGGCCTACTAATCCATTTGCTCCTGTTGCCCCTGTTGCTCCTGTTGGACCTGTTACTCCTGTTGCTCCGTCTGCTCCTGTTGCCCCTGTTGCTCCTGTTGCTCCTGTT")
	target := []byte("CCTGTTGCTCCTGTTGCTCCTGTTGCTCCTGTTGGGCCTGTTGCTCCTGTATTTCCTGTTGCTCCTGTTGCTCCTGTTGGACCTGCTACTCCATCTGCTCCTGTATTTCCTGTTGCTCCTGTTGGACCTGTTGCTCCATCTGCTCCTGTTGGACCTGTCGCTCCTGTTGCCCCTGTTGCTCCTGTTATTCCTACTCCTGTTGGGCCTGTTGCTCCTGTTGGACCTGTTGCTCCAGGTGTTCCTGCTGCTCCTGTTGCTCCTGTTGGGCCTACTAATCCATTTGCTCCTGTTGCCCCTGTTGCTCCTGTTGGACCTGTTACTCCTGTTGCTCCGTCTGCTCCTGTTGCCCCTGTT")
	for _, pair := range [][2][]byte{{query, target}, {target, query}} {
		algn := New(DefaultPenalties, &Options{GlobalAlignment: true})
		if err := algn.AdaptiveReduction(DefaultAdaptiveOption); err != nil {
			t.Fatal(err)
		}
		result, err := algn.Align(pair[0], pair[1])
		if err != nil {
			t.Fatal(err)
		}
		validateGlobalResult(t, pair[0], pair[1], DefaultPenalties, result)
		RecycleAlignmentResult(result)
		RecycleAligner(algn)
	}
}

// TestAdaptiveReductionRetainsGlobalBoundary keeps terminal-gap paths eligible.
func TestAdaptiveReductionRetainsGlobalBoundary(t *testing.T) {
	algn := New(DefaultPenalties, &Options{GlobalAlignment: true})
	defer RecycleAligner(algn)
	if err := algn.AdaptiveReduction(&AdaptiveReductionOption{MinWFLen: 1, MaxDistDiff: 20}); err != nil {
		t.Fatal(err)
	}
	query, target := make([]byte, 100), make([]byte, 100)
	algn.initComponents(&query, &target)
	const score = 10
	algn.M.Set(score, -20, 80, wfaMatch) // Query end: 20 target bases remain.
	algn.M.Set(score, -10, 5, wfaMatch)  // A distant interior path triggers pruning.
	algn.M.Set(score, 0, 50, wfaMatch)
	algn.M.Set(score, 20, 100, wfaMatch) // Target end: 20 query bases remain.
	algn.reduce(&query, &target, score)
	for _, k := range []int{-20, 20} {
		if _, _, ok := algn.M.GetAfterDiff(score, 0, k); !ok {
			t.Errorf("pruned valid terminal-gap path at diagonal %d", k)
		}
	}
}

// TestAdaptiveScoreMatchesCIGAR covers adjacent gaps reopened after pruning.
func TestAdaptiveScoreMatchesCIGAR(t *testing.T) {
	query, target := []byte("TCAAGAGTCAAGGTCAAGGTCAAGGTCAAGGTCAAGGTCAAGGTCAAGGTCAAGGTCAAGGTCAAGGTCAAGGTCAAGGAG"), []byte("ATCAATCAAGGTCGGTCAAGGTCAACAAGGTCAAGGTCAAGGTCAAGGTCAAGGTCAAGGTCAAGGTCAAGGTCAAGATCAAGCAAGG")
	algn := New(DefaultPenalties, DefaultOptions)
	defer RecycleAligner(algn)
	if err := algn.AdaptiveReduction(DefaultAdaptiveOption); err != nil {
		t.Fatal(err)
	}
	result, err := algn.Align(query, target)
	if err != nil {
		t.Fatal(err)
	}
	defer RecycleAlignmentResult(result)
	validateGlobalResult(t, query, target, DefaultPenalties, result)
	if result.Score != 62 {
		t.Fatalf("score=%d, want 62", result.Score)
	}
}

// TestSemiGlobalLegacyExamples preserves the established endpoint constraints.
// A cheap corner overlap must not replace the README's longer alignment.
func TestSemiGlobalLegacyExamples(t *testing.T) {
	for _, adaptive := range []bool{false, true} {
		for _, tc := range []struct {
			q, target, cigar string
			score, matches   uint32
		}{
			{"Bioinformatics helps Biology", "We learn bioinformatics to help biologists", "9I1X14M3I4M1D1M1X5M1X3I", 32, 24},
			{"ACCATACTCG", "AGGATGCTCG", "1M2X2M1X4M", 12, 7},
			{"ACGATCTCG", "CAGGCTCCTCGG", "1I1M1X1M1X2M1I3M1I", 16, 7},
			{"AC", "CA", "2X", 8, 0},
			{"TTAC", "ACGG", "2H2M2I", 10, 2},
			{"AAAA", "AAAA", "4M", 0, 4},
			{"A", "CA", "1I1M", 0, 1},
			{"CA", "A", "1H1M", 0, 1},
		} {
			t.Run(fmt.Sprintf("adaptive=%v/%s/%s", adaptive, tc.q, tc.target), func(t *testing.T) {
				a := New(DefaultPenalties, &Options{GlobalAlignment: false})
				defer RecycleAligner(a)
				if adaptive {
					if err := a.AdaptiveReduction(DefaultAdaptiveOption); err != nil {
						t.Fatal(err)
					}
				}
				q, target := []byte(tc.q), []byte(tc.target)
				r, err := a.Align(q, target)
				if err != nil {
					t.Fatal(err)
				}
				defer RecycleAlignmentResult(r)
				if r.Score != tc.score || r.Matches != tc.matches || r.CIGAR(false) != tc.cigar {
					t.Fatalf("score=%d matches=%d cigar=%s; want %d, %d, %s", r.Score, r.Matches, r.CIGAR(false), tc.score, tc.matches, tc.cigar)
				}
				if tc.q == "Bioinformatics helps Biology" && (r.QBegin != 2 || r.QEnd != 27 || r.TBegin != 11 || r.TEnd != 38 || r.AlignLen != 29 || r.Gaps != 4) {
					t.Fatalf("README statistics changed: %+v", r)
				}
			})
		}
	}
}

// TestInvalidPenaltiesRejectBeforeAlignment prevents zero-cost recurrences.
func TestInvalidPenaltiesRejectBeforeAlignment(t *testing.T) {
	for _, p := range []*Penalties{{Mismatch: 0, GapOpen: 6, GapExt: 2}, {Mismatch: 4, GapOpen: 6, GapExt: 0}, {Mismatch: 4, GapOpen: ^uint32(0), GapExt: 1}} {
		algn := New(p, DefaultOptions)
		_, err := algn.Align([]byte("A"), []byte("C"))
		RecycleAligner(algn)
		if !errors.Is(err, ErrInvalidPenalties) {
			t.Errorf("penalties %+v: expected ErrInvalidPenalties, got %v", p, err)
		}
	}
}

// TestSemiGlobalLargeMismatch covers sparse initialization at a high score.
func TestSemiGlobalLargeMismatch(t *testing.T) {
	p := &Penalties{Mismatch: 5000, GapOpen: 6, GapExt: 2}
	algn := New(p, &Options{GlobalAlignment: false})
	defer RecycleAligner(algn)
	r, err := algn.Align([]byte("A"), []byte("C"))
	if err != nil {
		t.Fatal(err)
	}
	defer RecycleAlignmentResult(r)
	if r.Score != 5000 || r.CIGAR(false) != "1X" {
		t.Fatalf("score=%d cigar=%s", r.Score, r.CIGAR(false))
	}
}

// affineSemiGlobalScore is an independent DP oracle for the historical endpoint
// constraints: both consumed offsets must reach the shorter sequence length.
func affineSemiGlobalScore(q, target []byte, p *Penalties) uint32 {
	cols := len(target) + 1
	m, ins, del := make([]int, (len(q)+1)*cols), make([]int, (len(q)+1)*cols), make([]int, (len(q)+1)*cols)
	for i := range m {
		m[i], ins[i], del[i] = maxInt/4, maxInt/4, maxInt/4
	}
	best := maxInt / 4
	for i := 1; i <= len(q); i++ {
		for j := 1; j <= len(target); j++ {
			at := i*cols + j
			cost := 0
			if q[i-1] != target[j-1] {
				cost = int(p.Mismatch)
			}
			if i == 1 || j == 1 {
				m[at] = cost
			} else {
				m[at] = m[at-cols-1] + cost
			}
			if j > 1 {
				ins[at] = min(m[at-1]+int(p.GapOpen+p.GapExt), ins[at-1]+int(p.GapExt))
			}
			if i > 1 {
				del[at] = min(m[at-cols]+int(p.GapOpen+p.GapExt), del[at-cols]+int(p.GapExt))
			}
			m[at] = min(m[at], ins[at], del[at])
			if (i == len(q) && j >= len(q)) || (j == len(target) && i >= len(target)) {
				best = min(best, m[at])
			}
		}
	}
	return uint32(best)
}

// TestSemiGlobalAgainstAffineOracle covers edge selection, traceback, and reuse.
func TestSemiGlobalAgainstAffineOracle(t *testing.T) {
	for _, p := range []*Penalties{DefaultPenalties, {Mismatch: 3, GapOpen: 0, GapExt: 2}, {Mismatch: 9, GapOpen: 1, GapExt: 1}} {
		a := New(p, &Options{GlobalAlignment: false})
		sequences := allDNASequences(5)
		for _, q := range sequences {
			for _, s := range sequences {
				r, err := a.Align(q, s)
				if err != nil {
					t.Fatal(err)
				}
				if want := affineSemiGlobalScore(q, s, p); r.Score != want {
					t.Fatalf("q=%s t=%s p=%+v score=%d want=%d cigar=%s", q, s, p, r.Score, want, r.CIGAR(false))
				}
				qi, ti := 0, 0
				for _, packed := range r.Ops {
					op, n := Op(packed)
					for j := uint32(0); j < n; j++ {
						switch op {
						case 'M', 'X':
							if qi >= len(q) || ti >= len(s) || (op == 'M') != (q[qi] == s[ti]) {
								t.Fatalf("q=%s t=%s: invalid CIGAR %s", q, s, r.CIGAR(false))
							}
							qi++
							ti++
						case 'I':
							ti++
						case 'D', 'H':
							qi++
						default:
							t.Fatalf("unexpected CIGAR operation %c", op)
						}
					}
				}
				if qi != len(q) || ti != len(s) {
					t.Fatalf("CIGAR consumed %d/%d, want %d/%d", qi, ti, len(q), len(s))
				}
				RecycleAlignmentResult(r)
			}
		}
		RecycleAligner(a)
	}
}

// TestSemiGlobalPaidTerminalGaps distinguishes real free clipping from an end
// gap needed to satisfy the endpoint constraint, even after CIGAR runs merge.
func TestSemiGlobalPaidTerminalGaps(t *testing.T) {
	p := &Penalties{Mismatch: 9, GapOpen: 1, GapExt: 1}
	for _, tc := range []struct{ q, target, cigar string }{
		{"AAAC", "CCCCCCC", "3H1M6I"},
		{"CCCCCCC", "AAAC", "3I1M3D3H"},
	} {
		a := New(p, &Options{GlobalAlignment: false})
		r, err := a.Align([]byte(tc.q), []byte(tc.target))
		if err != nil {
			t.Fatal(err)
		}
		if r.Score != 4 || r.CIGAR(false) != tc.cigar {
			t.Fatalf("score=%d cigar=%s; want 4, %s", r.Score, r.CIGAR(false), tc.cigar)
		}
		RecycleAlignmentResult(r)
		RecycleAligner(a)
	}
}

// TestAdaptiveReductionUnequalComponentLengths checks pruning at score-block
// boundaries when insertion and deletion components have grown independently.
func TestAdaptiveReductionUnequalComponentLengths(t *testing.T) {
	for _, score := range []uint32{2048, 4096} {
		for _, grown := range []string{"neither", "insertion", "deletion"} {
			t.Run(fmt.Sprintf("score=%d/grown=%s", score, grown), func(t *testing.T) {
				algn := New(DefaultPenalties, &Options{GlobalAlignment: true})
				defer RecycleAligner(algn)
				if err := algn.AdaptiveReduction(DefaultAdaptiveOption); err != nil {
					t.Fatal(err)
				}
				query, target := make([]byte, 100), make([]byte, 100)
				algn.initComponents(&query, &target)
				algn.M.Set(score, -1, 1, wfaMatch)
				algn.M.Set(score, 0, 90, wfaMatch)
				algn.M.Set(score, 1, 1, wfaMatch)
				var gap *Component
				if grown == "insertion" {
					gap = algn.I
				} else if grown == "deletion" {
					gap = algn.D
				}
				if gap != nil {
					gap.Set(score, -1, 1, wfaMatch)
					gap.Set(score, 0, 90, wfaMatch)
					gap.Set(score, 1, 1, wfaMatch)
				}
				lengths := [3]int{len(algn.M.WaveFronts), len(algn.I.WaveFronts), len(algn.D.WaveFronts)}
				algn.reduce(&query, &target, score)
				if got := [3]int{len(algn.M.WaveFronts), len(algn.I.WaveFronts), len(algn.D.WaveFronts)}; got != lengths {
					t.Fatalf("pruning grew components: %v -> %v", lengths, got)
				}
				for _, c := range []*Component{algn.M, gap} {
					if c == nil {
						continue
					}
					for _, k := range []int{-1, 1} {
						if _, _, ok := c.Get(score, k); ok {
							t.Fatalf("distant diagonal %d survived pruning", k)
						}
					}
					if offset, _, ok := c.Get(score, 0); !ok || offset != 90 {
						t.Fatalf("best diagonal lost: offset=%d exists=%v", offset, ok)
					}
				}
			})
		}
	}
}

// validateAlignmentDetails independently checks CIGAR bases, consumption, and match-region statistics.
func validateAlignmentDetails(t *testing.T, q, target []byte, r *AlignmentResult, p *Penalties, global, exact bool) {
	t.Helper()
	qi, ti := 0, 0
	var score uint32
	first, last := -1, -1
	bounds := [4]int{}
	for i, packed := range r.Ops {
		op, n := Op(packed)
		if n == 0 {
			t.Fatalf("zero CIGAR op: %s", r.CIGAR(false))
		}
		if i > 0 && op == byte(r.Ops[i-1]>>32) {
			t.Fatalf("unmerged CIGAR: %s", r.CIGAR(false))
		}
		switch op {
		case 'M', 'X':
			for j := uint32(0); j < n; j++ {
				if qi >= len(q) || ti >= len(target) || (op == 'M') != (q[qi] == target[ti]) {
					t.Fatalf("invalid %c at %d,%d q=%q t=%q CIGAR=%s p=%+v global=%t exact=%t", op, qi, ti, q, target, r.CIGAR(false), p, global, exact)
				}
				qi++
				ti++
				if op == 'M' {
					if first < 0 {
						first = i
						bounds[0], bounds[2] = qi, ti
					}
					last = i
					bounds[1], bounds[3] = qi, ti
				}
			}
			if op == 'X' {
				score += n * p.Mismatch
			}
		case 'I':
			ti += int(n)
			score += p.GapOpen + n*p.GapExt
		case 'D':
			qi += int(n)
			score += p.GapOpen + n*p.GapExt
		case 'H':
			if global || i != 0 && i != len(r.Ops)-1 {
				t.Fatalf("invalid clipping: %s", r.CIGAR(false))
			}
			qi += int(n)
		default:
			t.Fatalf("invalid CIGAR op %c", op)
		}
	}
	if qi != len(q) || ti != len(target) {
		t.Fatalf("consumed %d/%d vs %d/%d: %s", qi, ti, len(q), len(target), r.CIGAR(false))
	}
	if got := [4]int{r.QBegin, r.QEnd, r.TBegin, r.TEnd}; got != bounds {
		t.Fatalf("bounds=%v want=%v q=%q t=%q CIGAR=%s global=%t exact=%t", got, bounds, q, target, r.CIGAR(false), global, exact)
	}
	var matches, gaps, gapRegions, alen uint32
	for i := first; i >= 0 && i <= last; i++ {
		op, n := Op(r.Ops[i])
		alen += n
		if op == 'M' {
			matches += n
		}
		if op == 'I' || op == 'D' {
			gaps += n
			gapRegions++
		}
	}
	if got, want := [4]uint32{r.Matches, r.Gaps, r.GapRegions, r.AlignLen}, [4]uint32{matches, gaps, gapRegions, alen}; got != want {
		t.Fatalf("stats=%v want=%v", got, want)
	}
	if global && score != r.Score {
		t.Fatalf("CIGAR score=%d reported=%d", score, r.Score)
	}
	want := affineGlobalScore(q, target, p)
	if !global {
		want = affineSemiGlobalScore(q, target, p)
	}
	if (exact && want != r.Score) || r.Score < want {
		t.Fatalf("q=%q t=%q p=%+v global=%t exact=%t score=%d optimal=%d cigar=%s", q, target, p, global, exact, r.Score, want, r.CIGAR(false))
	}
	for _, trim := range []bool{false, true} {
		Q, A, T := r.AlignmentText(&q, &target, trim)
		if len(*Q) != len(*T) || len(*Q) != len(*A) {
			t.Fatal("rendered lengths differ")
		}
		RecycleAlignmentText(Q, A, T)
	}
}

// TestAdaptiveAlignmentSequenceEnds prevents recurrence steps beyond either
// input, even when aggressive pruning makes a nonoptimal path win.
func TestAdaptiveAlignmentSequenceEnds(t *testing.T) {
	p := &Penalties{Mismatch: 3, GapOpen: 6, GapExt: 3}
	for _, reverse := range []bool{false, true} {
		for _, swap := range []bool{false, true} {
			t.Run(fmt.Sprintf("reverse=%t/swap=%t", reverse, swap), func(t *testing.T) {
				q, target := []byte("AC"), []byte("TTATA")
				if reverse {
					reverseComplement(q)
					reverseComplement(target)
				}
				if swap {
					q, target = target, q
				}
				a := New(p, DefaultOptions)
				defer RecycleAligner(a)
				if err := a.AdaptiveReduction(&AdaptiveReductionOption{MinWFLen: 1}); err != nil {
					t.Fatal(err)
				}
				for i := 0; i < 2; i++ {
					r, err := a.Align(q, target)
					if err != nil {
						t.Fatal(err)
					}
					validateAlignmentDetails(t, q, target, r, p, true, false)
					RecycleAlignmentResult(r)
				}
			})
		}
	}
}

// reverseComplement reverses DNA orientation for alignment regression cases.
func reverseComplement(s []byte) {
	for i, j := 0, len(s)-1; i <= j; i, j = i+1, j-1 {
		complement := func(b byte) byte {
			switch b {
			case 'A':
				return 'T'
			case 'T':
				return 'A'
			case 'C':
				return 'G'
			case 'G':
				return 'C'
			default:
				return b
			}
		}
		s[i], s[j] = complement(s[j]), complement(s[i])
	}
}

// FuzzAlignmentResult compares exact alignment with an independent affine-gap
// oracle. Heuristic alignment may cost more, but must return a valid full path,
// consistent match-region statistics, and renderable CIGAR coordinates.
func FuzzAlignmentResult(f *testing.F) {
	for _, pair := range [][2]string{{"AC", "TTATA"}, {"A", "C"}, {"ACGT", "ACGT"}, {"ACACACA", "CAC"}, {"AAAA", "AAACAAAA"}, {"aN-\x00", "N-a\x00"}} {
		for _, global := range []bool{true, false} {
			for _, adaptive := range []uint8{0, 1, 4} {
				f.Add([]byte(pair[0]), []byte(pair[1]), uint8(2), uint8(6), uint8(2), global, adaptive)
			}
		}
	}
	// Minimized fuzz regression: traceback must not select an invalid predecessor.
	f.Add([]byte("0CAC0C0"), []byte("CAC"), uint8(3), uint8(6), uint8(1), true, uint8(1))
	f.Fuzz(func(t *testing.T, q, target []byte, mismatch, open, extend uint8, global bool, adaptive uint8) {
		if len(q) == 0 || len(target) == 0 || len(q) > 96 || len(target) > 96 {
			return
		}
		p := &Penalties{Mismatch: 1 + uint32(mismatch)%15, GapOpen: uint32(open) % 10, GapExt: 1 + uint32(extend)%6}
		a := New(p, &Options{GlobalAlignment: global})
		defer RecycleAligner(a)
		mode := adaptive % 5
		if mode != 0 {
			ad := &AdaptiveReductionOption{MinWFLen: 1, MaxDistDiff: []uint32{0, 0, 1, 5, 50}[mode]}
			if err := a.AdaptiveReduction(ad); err != nil {
				t.Fatal(err)
			}
		}
		for reuse := 0; reuse < 2; reuse++ {
			r, err := a.Align(q, target)
			if err != nil {
				t.Fatal(err)
			}
			validateAlignmentDetails(t, q, target, r, p, global, mode == 0)
			RecycleAlignmentResult(r)
		}
	})
}

// TestAdaptiveReductionLargeThreshold checks unsigned options on 32-bit hosts.
func TestAdaptiveReductionLargeThreshold(t *testing.T) {
	a := New(DefaultPenalties, DefaultOptions)
	defer RecycleAligner(a)
	if err := a.AdaptiveReduction(&AdaptiveReductionOption{MinWFLen: 1, MaxDistDiff: ^uint32(0)}); err != nil {
		t.Fatal(err)
	}
	q, target := make([]byte, 100), make([]byte, 100)
	a.initComponents(&q, &target)
	for k := -1; k <= 1; k++ {
		a.M.Set(10, k, 90, wfaMatch)
	}
	a.reduce(&q, &target, 10)
	for k := -1; k <= 1; k++ {
		if _, _, ok := a.M.Get(10, k); !ok {
			t.Fatalf("large distance threshold pruned diagonal %d", k)
		}
	}
}

// TestAdaptiveReductionLargeMinimum keeps narrow wavefronts unpruned on 32-bit hosts.
func TestAdaptiveReductionLargeMinimum(t *testing.T) {
	p := &Penalties{Mismatch: 3, GapOpen: 6, GapExt: 3}
	a := New(p, DefaultOptions)
	defer RecycleAligner(a)
	if err := a.AdaptiveReduction(&AdaptiveReductionOption{MinWFLen: ^uint32(0)}); err != nil {
		t.Fatal(err)
	}
	q, target := []byte("AC"), []byte("TTATA")
	r, err := a.Align(q, target)
	if err != nil {
		t.Fatal(err)
	}
	defer RecycleAlignmentResult(r)
	validateAlignmentDetails(t, q, target, r, p, true, true)
}

// TestAdaptiveTracebackBoundary follows the stored predecessor instead of an
// inadmissible mismatch beyond a sequence end. Both DNA orientations are checked.
func TestAdaptiveTracebackBoundary(t *testing.T) {
	for _, reverse := range []bool{false, true} {
		for _, swap := range []bool{false, true} {
			q, target := []byte("GCACGCG"), []byte("CAC")
			if reverse {
				reverseComplement(q)
				reverseComplement(target)
			}
			if swap {
				q, target = target, q
			}
			a := New(DefaultPenalties, DefaultOptions)
			if err := a.AdaptiveReduction(&AdaptiveReductionOption{MinWFLen: 1}); err != nil {
				t.Fatal(err)
			}
			r, err := a.Align(q, target)
			if err != nil {
				t.Fatal(err)
			}
			validateAlignmentDetails(t, q, target, r, DefaultPenalties, true, false)
			RecycleAlignmentResult(r)
			RecycleAligner(a)
		}
	}
}

// TestInvalidAdaptiveOptions reports invalid settings without disturbing the
// previously configured heuristic or panicking on a nil options pointer.
func TestInvalidAdaptiveOptions(t *testing.T) {
	a := New(DefaultPenalties, DefaultOptions)
	defer RecycleAligner(a)
	if err := a.AdaptiveReduction(DefaultAdaptiveOption); err != nil {
		t.Fatal(err)
	}
	for _, ad := range []*AdaptiveReductionOption{nil, {MinWFLen: 0, CutoffStep: 1}} {
		if err := a.AdaptiveReduction(ad); err == nil {
			t.Fatal("invalid options were accepted")
		}
		if a.ad != DefaultAdaptiveOption {
			t.Fatal("invalid options changed the existing heuristic")
		}
	}
}

// TestRandomAlignmentModes keeps broad oracle and rendering checks in the
// ordinary suite, including default/aggressive heuristics and aligner reuse.
func TestRandomAlignmentModes(t *testing.T) {
	rng := rand.New(rand.NewSource(20261029))
	for i := 0; i < 1000; i++ {
		p := &Penalties{Mismatch: uint32(1 + rng.Intn(15)), GapOpen: uint32(rng.Intn(10)), GapExt: uint32(1 + rng.Intn(6))}
		global, mode := rng.Intn(2) == 0, rng.Intn(5)
		a := New(p, &Options{GlobalAlignment: global})
		if mode != 0 {
			if err := a.AdaptiveReduction(&AdaptiveReductionOption{MinWFLen: 1, MaxDistDiff: []uint32{0, 0, 1, 5, 50}[mode]}); err != nil {
				t.Fatal(err)
			}
		}
		q, target := generateRandomSeq(1+rng.Intn(150), rng.Int63()), generateRandomSeq(1+rng.Intn(150), rng.Int63())
		for orientation := 0; orientation < 2; orientation++ {
			r, err := a.Align(q, target)
			if err != nil {
				t.Fatal(err)
			}
			validateAlignmentDetails(t, q, target, r, p, global, mode == 0)
			RecycleAlignmentResult(r)
			reverseComplement(q)
			reverseComplement(target)
		}
		RecycleAligner(a)
	}
}
