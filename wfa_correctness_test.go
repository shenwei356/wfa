package wfa

import (
	"context"
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
