// Copyright © 2024 Wei Shen <shenwei356@gmail.com>
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

package wfa

import (
	"encoding/binary"
	"fmt"
	"math"
	"math/bits"
	"sync"
)

// Penalties contains the gap-affine penalties, Match is 0.
type Penalties struct {
	Mismatch uint32
	GapOpen  uint32
	GapExt   uint32
}

// DefaultPenalties is from the WFA paper.
var DefaultPenalties = &Penalties{
	Mismatch: 4,
	GapOpen:  6,
	GapExt:   2,
}

// AdaptiveReductionOption contains the parameters for adaptive reduction
type AdaptiveReductionOption struct {
	MinWFLen    uint32
	MaxDistDiff uint32
	CutoffStep  uint32 // not used yet.
}

// DefaultAdaptiveOption provides a default option with parameters from the official repo.
// attributes.heuristic.min_wavefront_length = 10;
// attributes.heuristic.max_distance_threshold = 50;
// attributes.heuristic.steps_between_cutoffs = 1;
var DefaultAdaptiveOption = &AdaptiveReductionOption{
	MinWFLen:    10,
	MaxDistDiff: 50,
	CutoffStep:  1,
}

// Options represents a list of options.
// Currently it only support global or semi-global alignment.
type Options struct {
	GlobalAlignment bool
}

// DefaultOptions is the default option
var DefaultOptions = &Options{
	GlobalAlignment: true,
}

// Aligner is the object for aligning,
// which can apply to multiple pairs of query and ref sequences.
// But it's not occurrence safe, which means you can't call Align() in multiple goroutines.
// Instead, you can create multiple aligners, one for each goroutine.
// Aligner objects are from a object pool, in case a large number of alignments are needed.
// Just remember to recyle it with RecycleAligner().
type Aligner struct {
	p *Penalties

	ad *AdaptiveReductionOption

	opt *Options

	M, I, D *Component
}

// object pool of aligners.
var poolAligner = &sync.Pool{New: func() interface{} {
	algn := Aligner{
		p: nil,
		M: NewComponent(),
		I: NewComponent(),
		D: NewComponent(),
	}
	algn.M.IsM = true
	return &algn
}}

// RecycleAligner recycles an Aligner object.
func RecycleAligner(algn *Aligner) {
	if algn != nil {
		// there's no need to recyle them, just leave them with the aligner.

		// RecycleComponent(algn.M)
		// RecycleComponent(algn.I)
		// RecycleComponent(algn.D)

		// algn.M = nil
		// algn.I = nil
		// algn.D = nil

		poolAligner.Put(algn)
	}
}

// New returns a new Aligner from the object pool.
// Do not forget to call RecycleAligner() after using it.
func New(p *Penalties, opt *Options) *Aligner {
	algn := poolAligner.Get().(*Aligner)
	algn.p = p
	algn.opt = opt

	// there's no need to recyle them, just leave them with the aligner.
	// algn.M = NewComponent()
	// algn.I = NewComponent()
	// algn.D = NewComponent()

	return algn
}

// AdaptiveReduction sets the adaptive reduction parameters
func (algn *Aligner) AdaptiveReduction(ad *AdaptiveReductionOption) error {
	if ad.MinWFLen == 0 {
		return fmt.Errorf("cutoff step should not be 0")
	}
	algn.ad = ad
	return nil
}

// initComponents resets the internal data before each alignment.
func (algn *Aligner) initComponents(q, t *[]byte) {
	// // clear all wavefronts
	// algn.M.Reset()
	// algn.I.Reset()
	// algn.D.Reset()

	m, n := len(*t), len(*q)
	M := algn.M

	var wfaType, score uint32

	// have to check the first bases
	if (*q)[0] == (*t)[0] { // M[0,0] = 0
		wfaType, score = wfaMatch, 0
	} else { // M[0,0] = 4
		wfaType, score = wfaMismatch, algn.p.Mismatch
	}
	M.Set(score, 0, 1, wfaType)

	// for semi-global alignment
	if !algn.opt.GlobalAlignment {
		for k := 1; k < m; k++ { // first row
			if (*q)[0] == (*t)[k] {
				wfaType, score = wfaMatch, 0
			} else {
				wfaType, score = wfaMismatch, algn.p.Mismatch
			}

			M.Set(score, k, uint32(k+1), wfaType)
		}

		for k := 1; k < n; k++ { // first column
			if (*q)[k] == (*t)[0] {
				wfaType, score = wfaMatch, 0
			} else {
				wfaType, score = wfaMismatch, algn.p.Mismatch
			}

			M.Set(score, -k, 1, wfaType)
		}
	}
}

// ErrEmptySeq means the query or target sequence is empty.
var ErrEmptySeq error = fmt.Errorf("wfa: invalid empty sequence")

// MaxSeqLen is the allowed longest sequence length.
const MaxSeqLen int = 1<<(32-wfaTypeBits) - 1

// ErrSeqTooLong means the sequence is too long.
var ErrSeqTooLong error = fmt.Errorf("wfa: sequences longer than %d are not supported", MaxSeqLen)

// Align performs alignment with two sequences.
func (algn *Aligner) Align(q, t []byte) (*AlignmentResult, error) {
	return algn.AlignPointers(&q, &t)
}

// AlignPointers performs alignment with two sequences. The arguments are pointers.
func (algn *Aligner) AlignPointers(q, t *[]byte) (*AlignmentResult, error) {
	m, n := len(*t), len(*q)

	if n == 0 || m == 0 {
		return nil, ErrEmptySeq
	}
	if n > MaxSeqLen || m > MaxSeqLen {
		return nil, ErrSeqTooLong
	}

	algn.initComponents(q, t)

	// -------------------------------------------------

	Ak := m - n
	Aoffset := uint32(m)
	var offset uint32

	M := algn.M

	var s uint32
	var lo, hi int
	reduce := algn.ad != nil
	var minWFLen int
	if reduce {
		minWFLen = int(algn.ad.MinWFLen)
	}
	for {
		// fmt.Printf("---------------------- s: %-3d ----------------------\n", s)
		if M.HasScore(s) {
			// fmt.Printf("extend:\n")
			lo, hi = algn.extend(q, t, s)
			// fmt.Printf("max offset: %d, Aoffset: %d\n", (*(*M)[s])[Ak], Aoffset)

			offset, _, _ = M.GetAfterDiff(s, 0, Ak)
			if offset >= Aoffset { // reached the end
				// fmt.Printf("reach end, s:%d, k:%d, offset:%d, Aoffset:%d\n", s, Ak, offset>>wfaTypeBits, Aoffset)
				break
			}

			// fmt.Printf("reduce:\n")
			if reduce && hi-lo+1 >= minWFLen {
				algn.reduce(q, t, s)
			}
		}

		s++

		// fmt.Printf("next:\n")
		algn.next(q, t, s)
	}

	// M.Print(os.Stdout, "M")

	// bottom right cell
	minS, lastK := s, Ak
	// fmt.Printf("min s:%d, k:%d\n", minS, lastK)
	if !algn.opt.GlobalAlignment { // find the minimum score on the last row/column
		minS, lastK = algn.backtraceStartPosistion(q, t, s)
		// fmt.Printf("new min s:%d, k:%d\n", minS, lastK)
	}
	// offset, _, _ = M.Get(minS, 0, lastK)
	// h := offset
	// v := h - uint32(lastK)
	// fmt.Printf("min s:%d, k:%d, h:%d, v:%d\n", minS, lastK, h, v)

	r := algn.backTrace(q, t, minS, lastK)

	// clear all wavefronts
	algn.M.Reset()
	algn.I.Reset()
	algn.D.Reset()

	return r, nil
}

func (algn *Aligner) backtraceStartPosistion(q, t *[]byte, s uint32) (uint32, int) {
	M := algn.M
	m, n := len(*t), len(*q)
	minS := s
	Ak := m - n
	lastK := Ak

	var offset uint32
	var ok bool
	var k int
	var lastRowOrCol bool
	var h, v int
	var lo, hi int

	// algn.Plot(q, t, os.Stdout, algn.M, true, -1)
	// fmt.Printf("m: %d, n: %d\n", m, n)

	for _s := s; _s >= 0; _s-- {
		if !M.HasScore(_s) {
			if _s == 0 {
				break
			}
			continue
		}

		lo, hi = M.KRange(_s, 0)
		// fmt.Printf("test s:%d, lo:%d, hi:%d\n", _s, lo, hi)

		lastRowOrCol = false
		k = Ak
		// fmt.Printf("a, s:%d, k:%d\n", _s, k)
		for {
			if k < lo {
				break
			}

			offset, _, ok = M.GetAfterDiff(_s, 0, k)
			if !ok {
				k--
				continue
			}
			h = int(offset)
			v = h - k

			if v <= 0 || v > n || h > m { // bound check
				break
			}

			// fmt.Printf("  s: %d, offset:%d, k:%d, h:%d, v:%d\n", _s, offset, k, h, v)
			if (v == n && h >= n) || (h == m && v >= m) {
				// fmt.Println("    ok", v == n && h >= n, h == m && v >= m)
				lastRowOrCol = true
				break
			}

			k--
		}

		if lastRowOrCol && _s <= minS {
			lastK = k
			minS = _s
		}

		lastRowOrCol = false
		k = Ak + 1
		// fmt.Printf("b, s:%d, k:%d\n", _s, k)
		for {
			if k > hi {
				break
			}

			offset, _, ok = M.GetAfterDiff(_s, 0, k)
			if !ok {
				k++
				continue
			}
			h = int(offset)
			v = h - k

			if v <= 0 || v > n || h > m { // bound check
				break
			}

			// fmt.Printf("  s: %d, offset:%d, k:%d, h:%d, v:%d\n", _s, offset, k, h, v)
			if (v == n && h >= n) || (h == m && v >= m) {
				// fmt.Println("    ok", v == n && h >= n, h == m && v >= m)
				lastRowOrCol = true
				break
			}

			k++
		}

		if lastRowOrCol && _s <= minS {
			lastK = k
			minS = _s
		}

		if _s == 0 {
			break
		}
	}

	// fmt.Printf("min s:%d, lastk:%d\n", minS, lastK)
	return minS, lastK
}

var be = binary.BigEndian

// extend refers to the WF_EXTEND method.
// The return bool value indicates whether the end of one sequence is reached.
func (algn *Aligner) extend(q, t *[]byte, s uint32) (int, int) {
	wf := algn.M.WaveFronts[s]
	lo, hi := wf.Lo, wf.Hi
	// fmt.Printf("  lo: %d, hi: %d, offsets: %d\n", lo, hi, *offsets)

	var offset uint32
	var v, h int
	lenQ := len(*q)
	lenT := len(*t)
	var q8, t8 uint64
	var n, N int

	var ok bool
	for k := hi; k >= lo; k-- {
		offset, _, ok = wf.Get(k)
		// fmt.Printf("    k:%d, ok:%v, offset:%d\n", k, ok, offset>>wfaTypeBits)

		if !ok {
			continue
		}

		h = int(offset)                       // x
		v = h - k                             // y
		if v <= 0 || v >= lenQ || h >= lenT { // bound check
			continue
		}

		// offset is 1-based, here it's checking the base in the next position.

		// compare every 8 bases, convert 8 bases to a uint64, xor them and count leading zeroes.
		if v+8 <= lenQ && h+8 <= lenT {
			N = 0
			// fmt.Printf("      block wise, start from: h: %d, v: %d\n", h, v)
			for {
				q8, t8 = be.Uint64((*q)[v:v+8]), be.Uint64((*t)[h:h+8])
				n = bits.LeadingZeros64(q8^t8) >> 3 // divide 8
				v += n
				h += n
				N += n
				if n < 8 || v+8 >= lenQ || h+8 >= lenT {
					break
				}
			}
			if N == 0 {
				continue
			}

			// fmt.Printf("        block wise, %d matches\n", N)
			wf.Increase(k, uint32(N))

			if !(n == 8 && v < lenQ && h < lenT) {
				continue
			}
			// fmt.Printf("         need to check left ------------------------------ \n")
		}

		// compare each base

		N = 0
		for (*q)[v] == (*t)[h] {
			v++
			h++
			N++

			if v == lenQ || h == lenT {
				break
			}
		}
		if N == 0 {
			continue
		}

		// fmt.Printf("      byte wise: k: %d, extend to h: %d, v: %d\n", k, h, v)
		wf.Increase(k, uint32(N))
	}

	return lo, hi
}

// adaptive reduction
func (algn *Aligner) reduce(q, t *[]byte, s uint32) {
	wf := algn.M.WaveFronts[s] // previously, we've checked. M.HasScore(s)
	lo, hi := wf.Lo, wf.Hi
	var offset uint32
	var v, h int
	lenQ := len(*q)
	lenT := len(*t)
	var ok bool

	var d, minDist int
	ds := poolDist.Get().(*[]int)
	*ds = (*ds)[:0]
	minDist = math.MaxInt
	for k := lo; k <= hi; k++ {
		offset, _, ok = wf.Get(k)
		if !ok {
			*ds = append(*ds, -1)
			continue
		}

		h = int(offset) // x
		v = h - k       // y
		if v < 0 || v >= lenQ || h >= lenT {
			*ds = append(*ds, -1)
			continue
		}

		d = max(lenT-h, lenQ-v)
		*ds = append(*ds, d)

		if d < minDist {
			minDist = d
		}
	}

	_lo := lo
	_hi := hi
	maxDistDiff := int(algn.ad.MaxDistDiff)
	updateLo := true
	found := false
	I := algn.I
	D := algn.D
	for i, d := range *ds {
		if d < 0 {
			continue
		}
		if d-minDist > maxDistDiff {
			found = true
			if updateLo {
				_lo = lo + i + 1
			}
			(*ds)[i] = -1 // mark it
		} else {
			updateLo = false
		}
	}
	if found { // found some distance where d-minDist > maxDistDiff
		for i := len(*ds) - 1; i >= 0; i-- {
			if (*ds)[i] >= 0 {
				_hi = lo + i
				break
			}
		}
	}

	for k := lo; k < _lo; k++ {
		wf.Delete(k)
		I.Delete(s, k)
		D.Delete(s, k)
	}
	for k := _hi + 1; k <= hi; k++ {
		wf.Delete(k)
		I.Delete(s, k)
		D.Delete(s, k)
	}

	wf.Lo, wf.Hi = _lo, _hi

	poolDist.Put(ds)
}

// poolDist is used in reduce()
var poolDist = &sync.Pool{New: func() interface{} {
	tmp := make([]int, 0, 128)
	return &tmp
}}

// next refers to the WF_NEXT method.
func (algn *Aligner) next(q, t *[]byte, s uint32) {
	M := algn.M
	I := algn.I
	D := algn.D
	p := algn.p
	lenQ := len(*q)
	lenT := len(*t)

	loMismatch, hiMismatch := M.KRange(s, p.Mismatch)       // M[s-x]
	loGapOpen, hiGapOpen := M.KRange(s, p.GapOpen+p.GapExt) // M[s-o-e]
	loInsert, hiInsert := I.KRange(s, p.GapExt)             // I[s-e]
	loDelete, hiDelete := D.KRange(s, p.GapExt)             // D[s-e]

	hi := min(int(lenT-1), max(hiMismatch, hiGapOpen, hiInsert, hiDelete)+1)
	lo := max(-int(lenQ-1), min(loMismatch, loGapOpen, loInsert, loDelete)-1)

	// fmt.Printf("s: %d, k: %d -> %d, lenQ: %d, lenT: %d\n", s, lo, hi, lenQ, lenT)

	var fromI, fromD, fromM bool
	var v1, v2 uint32
	var Isk, Dsk, Msk uint32
	var updatedI, updatedD bool
	var wfaTypeI, wfaTypeD, wfaTypeM uint32
	for k := lo; k <= hi; k++ {
		updatedI, updatedD = false, false
		wfaTypeI, wfaTypeD, wfaTypeM = 0, 0, 0
		// fmt.Printf(" k: %d\n", k)

		// --------------------------------------
		// insertion: 🠦
		v1, _, fromM = M.GetAfterDiff(s, p.GapOpen+p.GapExt, k-1)
		v2, _, fromI = I.GetAfterDiff(s, p.GapExt, k-1)
		if fromM && int(v1) > lenT {
			fromM = false
			v1 = 0
		}
		if fromI && int(v2) > lenT {
			fromI = false
			v2 = 0
		}
		Isk = max(v1, v2) + 1
		if fromM || fromI {
			if fromM && fromI {
				if v1 >= v2 {
					wfaTypeI = wfaInsertOpen
				} else {
					wfaTypeI = wfaInsertExt
				}
			} else if fromM {
				wfaTypeI = wfaInsertOpen
			} else {
				wfaTypeI = wfaInsertExt
			}

			updatedI = true
			I.Set(s, k, Isk, wfaTypeI)
			// fmt.Printf("  %d fromM:%v(%d), fromI:%v(%d), save I: s=%d, k=%d, offset:%d, type:%s\n",
			// 	Isk, fromM, v1, fromI, v2, s, k, Isk, wfaType2str(wfaTypeI))
		} else {
			Isk = 0
		}

		// --------------------------------------
		// deletion: 🠧

		v1, _, fromM = M.GetAfterDiff(s, p.GapOpen+p.GapExt, k+1)
		v2, _, fromD = D.GetAfterDiff(s, p.GapExt, k+1)
		if fromM && int(v1)-k > lenQ {
			fromM = false
			v1 = 0
		}
		if fromD && int(v2)-k > lenQ {
			fromD = false
			v2 = 0
		}

		Dsk = max(v1, v2)
		if fromM || fromD {
			if fromM && fromD {
				if v1 >= v2 {
					wfaTypeD = wfaDeleteOpen
				} else {
					wfaTypeD = wfaDeleteExt
				}
			} else if fromM {
				wfaTypeD = wfaDeleteOpen
			} else {
				wfaTypeD = wfaDeleteExt
			}

			updatedD = true
			D.Set(s, k, Dsk, wfaTypeD)
			// fmt.Printf("  %d fromM:%v(%d), fromD:%v(%d), save D: s=%d, k=%d, offset:%d, type:%s\n",
			// 	Dsk, fromM, v1, fromD, v2, s, k, Dsk, wfaType2str(wfaTypeD))
		} else {
			Dsk = 0
		}

		// --------------------------------------
		// mismatch: ⬂

		v1, _, fromM = M.GetAfterDiff(s, p.Mismatch, k)
		if fromM && (int(v1) > lenT || int(v1)-k > lenQ) { // it's the last column/row
			fromM = false
			v1 = 0
		}
		Msk = max(Isk, Dsk, v1+1)
		if updatedI || updatedD || fromM {
			if updatedI && updatedD && fromM {
				if Msk == v1+1 { // mismatch is prefered if it might come from 3 ways
					wfaTypeM = wfaMismatch
				} else if Msk == Isk {
					wfaTypeM = wfaTypeI
				} else {
					wfaTypeM = wfaTypeD
				}
			} else if updatedI {
				if updatedD { // updatedI && updatedD && !fromM
					if Msk == Isk {
						wfaTypeM = wfaTypeI
					} else {
						wfaTypeM = wfaTypeD
					}
				} else if fromM { // updatedI && !updatedD && fromM
					if Msk == v1+1 { // mismatch is prefered
						wfaTypeM = wfaMismatch
					} else {
						wfaTypeM = wfaTypeI
					}
				} else { // updatedI && !updatedD && !fromM
					wfaTypeM = wfaTypeI
				}
			} else if updatedD {
				if fromM { // !updatedI && updatedD && fromM
					if Msk == v1+1 { // mismatch is prefered
						wfaTypeM = wfaMismatch
					} else {
						wfaTypeM = wfaTypeD
					}
				} else { // !updatedI && updatedD && !fromM
					wfaTypeM = wfaTypeD
				}
			} else { // !updatedI && !updatedD && fromM
				wfaTypeM = wfaMismatch
			}

			M.Set(s, k, Msk, wfaTypeM)
			// fmt.Printf("  %d fromI:%v(%d), fromD:%v(%d), fromM:%v(%d), save M: s=%d, k=%d, offset:%d, type:%s\n",
			// 	Msk, updatedI, Isk, updatedD, Dsk, fromM, v1+1, s, k, Msk, wfaType2str(wfaTypeM))
		}
	}
}

// backTrace backtraces the alignment
func (algn *Aligner) backTrace(q, t *[]byte, s uint32, Ak int) *AlignmentResult {
	semiGlobal := !algn.opt.GlobalAlignment
	var M0 *Component
	M := algn.M
	I := algn.I
	D := algn.D
	p := algn.p
	lenQ := len(*q)
	lenT := len(*t)

	cigar := NewAlignmentResult()
	cigar.Score = s

	var ok bool
	var k, h, v int
	var offset, wfaType uint32
	var h0 int
	var op byte
	var qBegin, tBegin int

	var v1, v2, Isk, Dsk, offset0 uint32
	var fromMI, fromMD, fromItself bool
	var fromI, fromD, fromM bool
	var sMismatch, sGapOpen, sGapExt uint32
	var previousFromM bool
	var nMatches int

	k = Ak
	firstMatch := true

	// fmt.Printf("backtrace from M%d,%d, lenQ: %d, lenT: %d:\n", s, Ak, lenQ, lenT)

	// ------------------------------------------------
	// start point

	offset, _ = M.GetRaw(s, k)
	// fmt.Printf("------\nfirst s: %d, k: %d, offset: %d, existed: %v\n", s, k, offset>>wfaTypeBits, true)

	previousFromM = true
	wfaType = offset & wfaTypeMask
	h = int(offset >> wfaTypeBits)
	v = h - k

	if h < lenT {
		cigar.AddN(wfaOps[wfaInsertOpen], uint32(lenT)-uint32(h))
	} else if v < lenQ {
		cigar.AddN('H', uint32(lenQ)-uint32(v))
	}

LOOP:
	for v > 0 && h > 0 {
		// fmt.Printf("------\ncurrent s: %d, k: %d, h: %d, v: %d\n", s, k, h, v)

		// -----------------------------------------------------------------------------
		// compute the offset before extending

		// score of source
		sMismatch = s - p.Mismatch
		sGapOpen = s - p.GapOpen - p.GapExt
		sGapExt = s - p.GapExt

		// offset of the source
		fromMI, fromMD = false, false
		switch wfaType {
		case wfaInsertExt:
			v1, _, fromM = M.Get(sGapOpen, k-1)
			v2, _, fromI = I.Get(sGapExt, k-1)
			if fromM || fromI {
				fromMI = true
				offset0 = max(v1, v2) + 1
			} else {
				offset0 = 0
			}

			M0 = I // for get the wfaType of the next one
		case wfaDeleteExt:
			v1, _, fromM = M.Get(sGapOpen, k+1)
			v2, _, fromD = D.Get(sGapExt, k+1)
			if fromM || fromD {
				fromMD = true
				offset0 = max(v1, v2)
			} else {
				offset0 = 0
			}

			M0 = D
		default:
			v1, _, fromM = M.Get(sGapOpen, k-1)
			v2, _, fromI = I.Get(sGapExt, k-1)
			if fromM || fromI {
				fromMI = true
				Isk = max(v1, v2) + 1
			} else {
				Isk = 0
			}

			v1, _, fromM = M.Get(sGapOpen, k+1)
			v2, _, fromD = D.Get(sGapExt, k+1)
			if fromM || fromD {
				fromMD = true
				Dsk = max(v1, v2)
			} else {
				Dsk = 0
			}

			v1, _, fromM = M.Get(sMismatch, k)
			if fromMI || fromMD || fromM {
				offset0 = max(Isk, Dsk, v1+1)
				fromItself = false
			} else {
				fromItself = true
			}

			M0 = M
		}
		if fromItself {
			// fmt.Printf("  break as there's no valid source offset\n")
			break
		}
		if offset0 == 0 {
			// fmt.Printf("  break as there's no valid source offset\n")
			break
		}

		h0 = int(offset0)

		// fmt.Printf("  current type: %s, h0:%d, nMatches:%d\n",
		// 	wfaType2str(wfaType), h0, h-h0)

		// traceback matches
		if previousFromM {
			nMatches = h - h0
			// fmt.Printf("  fromM h0: %d, n:%d\n", h0, nMatches)

			// record matches
			if nMatches > 0 {
				if firstMatch { // record the end position of matched region
					firstMatch = false
					cigar.TEnd, cigar.QEnd = h, v
					// fmt.Printf("    == end position of matched region, t:%d, q:%d\n", h, v)
				}

				op = wfaOps[wfaMatch] // correct it as M
				cigar.AddN(op, uint32(nMatches))
				// fmt.Printf("    [ADD %s as match]: %d%c, h:%d, v:%d\n", wfaType2str(wfaType), nMatches, op, h, v)
			}

			// update coordinates with the offset before extention
			offset = offset0
			h = int(offset)
			v = h - k
			// fmt.Printf("  update h:%d, v:%d\n", h, v)

			// update the start position of matched region
			if wfaType == wfaMatch { // first line/row
				tBegin, qBegin = h, v
				// fmt.Printf("  -- update start position: h:%d, v:%d\n", h, v)
			} else if nMatches > 0 {
				tBegin, qBegin = h+1, v+1
				// fmt.Printf("  -- update start position: h:%d, v:%d\n", h+1, v+1)
			}

			if h <= 0 || v <= 0 {
				// fmt.Printf("  break as h<=0 || v <=0\n")
				break
			}
		}

		// record
		op = wfaOps[wfaType]
		cigar.AddN(op, 1)
		// fmt.Printf("  [ADD]: %s as %d%c, h:%d, v:%d\n", wfaType2str(wfaType), 1, op, h, v)

		if semiGlobal && (h == 1 || v == 1) {
			// fmt.Printf("  break as reached 1th row/col\n")
			break
		}

		// -----------------------------------------------------------------------------
		// for next one

		// update score, h, k according to wfaType of current one
		previousFromM = true
		switch wfaType {
		case wfaMismatch:
			s = sMismatch
			h--
		case wfaInsertOpen:
			s = sGapOpen
			k--
			h--
		case wfaInsertExt:
			s = sGapExt
			k--
			h--
			previousFromM = false
		case wfaDeleteOpen:
			s = sGapOpen
			k++
		case wfaDeleteExt:
			s = sGapExt
			k++
			previousFromM = false
		default:
			// fmt.Printf("  break as invalid wfa type\n")
			break LOOP
		}
		// update coordinates
		v = h - k
		// fmt.Printf("  %s from M: %v\n", wfaType2str(wfaType), previousFromM)

		// wfaType of the next one
		offset, ok = M0.GetRaw(s, k)
		if !ok {
			// fmt.Printf("  break as invalid wfa type. s: %d, k: %d\n", s, k)
			break
		}
		wfaType = offset & wfaTypeMask
		// fmt.Printf("\n  next type: %s, s:%d, k:%d\n", wfaType2str(wfaType), s, k)

		// fmt.Printf("  NEXT s: %d, k: %d, h:%d, v:%d\n", s, k, h, v)
	}

	// -----------------------------------------------------------------------------
	// the last one

	// fmt.Printf("------\nexit loop. h:%d, v:%d\n", h, v)
	if h > 0 && v > 0 {
		nMatches = min(h, v) - 1
		// fmt.Printf("nmatches: %d\n", nMatches)
		if nMatches > 0 {
			if firstMatch { // record the end position of matched region
				firstMatch = false
				cigar.TEnd, cigar.QEnd = h, v
				// fmt.Printf("    == end position of matched region, t:%d, q:%d\n", h, v)
			}

			op = wfaOps[wfaMatch] // correct it as M
			cigar.AddN(op, uint32(nMatches))
			// fmt.Printf("[ADD %s as match]: %d%c, h:%d, v:%d\n", wfaType2str(wfaType), nMatches, op, h, v)
			h -= nMatches
			v -= nMatches
			// fmt.Printf("h:%d, v:%d\n", h, v)

			// update the start position of matched region
			if wfaType == wfaMatch { // first line/row
				tBegin, qBegin = h, v
				// fmt.Printf("  -- update start position: h:%d, v:%d\n", h, v)
			} else if nMatches > 0 {
				tBegin, qBegin = h+1, v+1
				// fmt.Printf("  -- update start position: h:%d, v:%d\n", h+1, v+1)
			}
		} else if wfaType == wfaMatch { // first line/row
			tBegin, qBegin = h, v
			// fmt.Printf("  --b update start position: h:%d, v:%d\n", h, v)
			if firstMatch { // record the end position of matched region
				firstMatch = false
				cigar.TEnd, cigar.QEnd = h, v
				// fmt.Printf("    == end position of matched region, t:%d, q:%d\n", h, v)
			}
		}

		op = wfaOps[wfaType]
		cigar.AddN(op, 1)
		//	fmt.Printf("  final [ADD]: %s as %d%c, h:%d, v:%d\n", wfaType2str(wfaType), 1, op, h, v)
	}

	if v > 1 {
		cigar.AddN('H', uint32(v-1))
	}

	if h > 1 {
		cigar.AddN(wfaOps[wfaInsertOpen], uint32(h-1))
	}

	// fmt.Printf("  == start position of matched region, t:%d, q:%d\n", tBegin, qBegin)
	cigar.TBegin, cigar.QBegin = tBegin, qBegin

	cigar.process()
	return cigar
}
