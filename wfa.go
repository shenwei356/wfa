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

// DefaultPenalties is from the paper.
var DefaultPenalties = &Penalties{
	Mismatch: 4,
	GapOpen:  6,
	GapExt:   2,
}

// AdaptiveReductionOption contains the parameters for adaptive reduction
type AdaptiveReductionOption struct {
	MinWFLen    uint32
	MaxDistDiff uint32
	CutoffStep  uint32
}

// attributes.heuristic.min_wavefront_length = 10;
// attributes.heuristic.max_distance_threshold = 50;
// attributes.heuristic.steps_between_cutoffs = 1;
var DefaultAdaptiveOption = &AdaptiveReductionOption{
	MinWFLen:    10,
	MaxDistDiff: 50,
	CutoffStep:  1,
}

// Options represents a list of options
type Options struct {
	GlobalAlignment bool
}

// DefaultOptions is the default option
var DefaultOptions = &Options{
	GlobalAlignment: true,
}

// Aligner is the object for aligning,
// which can apply to multiple pairs of query and ref sequences.
// And it's from a object pool, in case a large number of alignment are needed.
type Aligner struct {
	p *Penalties

	ad *AdaptiveReductionOption

	opt *Options

	// The Wavefront component is a list of offsets for different scores.
	// So we have two indexes to save the offsets: s (score), k (k).
	//
	// To support fast access, we use a list to store offsets of different scores,
	// the nil data means there's no such a score.
	//
	// Since k might be negative and usually the values are symmetrical,
	// we store them like this:
	//	k value in []uint32: 0, -1, 1, -2, 2
	// if the value is 0, it means there's no records for that k.
	// The function kLowHigh is used to return the lowest and largest k values.
	M, I, D []*[]uint32

	SM, SI, SD []*[]uint32 // offsets before extending, will be used in backtrace
}

// object pool of aligners.
var poolAligner = &sync.Pool{New: func() interface{} {
	algn := Aligner{
		p:  nil,
		M:  make([]*[]uint32, 0, 1024),
		I:  make([]*[]uint32, 0, 1024),
		D:  make([]*[]uint32, 0, 1024),
		SM: make([]*[]uint32, 0, 1024),
		SI: make([]*[]uint32, 0, 1024),
		SD: make([]*[]uint32, 0, 1024),
	}
	return &algn
}}

// RecycleAligner recycles an Aligner object.
func RecycleAligner(algn *Aligner) {
	if algn != nil {
		poolAligner.Put(algn)
	}
}

// New returns a new Aligner from the object pool.
// Do not forget to call RecycleAligner after using it.
func New(p *Penalties, opt *Options) *Aligner {
	algn := poolAligner.Get().(*Aligner)
	algn.p = p
	algn.opt = opt
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
	m, n := len(*t), len(*q)

	// -------------------------------------------------
	// I, D, S

	algn.initComponent(&algn.I)
	algn.initComponent(&algn.D)
	algn.initComponent(&algn.SM)
	algn.initComponent(&algn.SI)
	algn.initComponent(&algn.SD)

	// -------------------------------------------------
	// M

	algn.initComponent(&algn.M) // M[0,0] = 0

	var wfaType, score uint32
	// have to check the first bases
	if (*q)[0] == (*t)[0] { // M[0,0] = 0
		wfaType, score = wfaMatch, 0
	} else { // M[0,0] = 4
		wfaType, score = wfaMismatch, algn.p.Mismatch
	}
	setOffset2(&algn.M, score, 0, (1<<wfaTypeBits)|wfaType)
	setOffset2(&algn.SM, score, 0, (1<<wfaTypeBits)|wfaType)

	if !algn.opt.GlobalAlignment { // for semi-global alignment
		for k := 1; k < m; k++ { // first row
			if (*q)[0] == (*t)[k] {
				wfaType, score = wfaMatch, 0
			} else {
				wfaType, score = wfaMismatch, algn.p.Mismatch
			}
			setOffset2(&algn.M, score, k, (uint32(k+1)<<wfaTypeBits)|wfaType)
			setOffset2(&algn.SM, score, k, (uint32(k+1)<<wfaTypeBits)|wfaType)
		}
		for k := 1; k < n; k++ { // first column
			if (*q)[k] == (*t)[0] {
				wfaType, score = wfaMatch, 0
			} else {
				wfaType, score = wfaMismatch, algn.p.Mismatch
			}
			setOffset2(&algn.M, score, -k, (1<<wfaTypeBits)|wfaType)
			setOffset2(&algn.SM, score, -k, (1<<wfaTypeBits)|wfaType)
		}
	}
}

// poolOffsets is a object pool for offsets.
var poolOffsets = &sync.Pool{New: func() interface{} {
	tmp := make([]uint32, 0, 128)
	return &tmp
}}

// initComponent inilializes a new WFA component.
func (algn *Aligner) initComponent(M *[]*[]uint32) {
	// resets a WFA component (a list of offsets).
	for _, v := range *M {
		if v != nil {
			*v = (*v)[:0]
			poolOffsets.Put(v)
		}
	}
	*M = (*M)[:0]

	// inilializes
	offsets := poolOffsets.Get().(*[]uint32)
	*M = append(*M, offsets)
}

// ErrEmptySeq means one of the query or target sequence is empty
var ErrEmptySeq error = fmt.Errorf("wfa: invalid empty sequence")

// MaxSeqLen is the allowed longest sequence length
const MaxSeqLen int = 1<<(32-wfaTypeBits) - 1

var ErrSeqTooLong error = fmt.Errorf("wfa: sequences longer than %d are not supported", MaxSeqLen)

// Align performs alignment with two sequences.
func (algn *Aligner) Align(q, t []byte) (*CIGAR, error) {
	return algn.AlignPointers(&q, &t)
}

// AlignPointers performs alignment with two sequences. The arguments are pointers.
func (algn *Aligner) AlignPointers(q, t *[]byte) (*CIGAR, error) {
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

	M := &algn.M

	var s uint32
	var lo, hi int
	reduce := algn.ad != nil
	var maxDistDiff int
	if reduce {
		maxDistDiff = int(algn.ad.MaxDistDiff)
	}
	for {
		// fmt.Printf("---------------------- s: %-3d ----------------------\n", s)
		if (*M)[s] != nil {
			// fmt.Printf("extend:\n")
			lo, hi = algn.extend((*M)[s], q, t, s)
			// fmt.Printf("max offset: %d, Aoffset: %d\n", (*(*M)[s])[Ak], Aoffset)

			offset, _, _ = getOffset((*M)[s], Ak)
			if offset>>wfaTypeBits >= Aoffset { // reached the end
				// fmt.Printf("reach end, s:%d, k:%d, offset:%d, Aoffset:%d\n", s, Ak, offset>>wfaTypeBits, Aoffset)
				break
			}

			// fmt.Printf("reduce:\n")
			if reduce && hi-lo+1 >= maxDistDiff {
				algn.reduce(q, t, s, lo, hi)
			}
		}

		s++
		// fmt.Printf("next:\n")
		algn.next(q, t, s)
	}

	minS, lastK := s, Ak
	// fmt.Printf("min s:%d, k:%d\n", minS, lastK)
	if !algn.opt.GlobalAlignment { // find the minimum score on the last line
		minS, lastK = algn.backtraceStartPosistion(q, t, s)
		// fmt.Printf("new min s:%d, k:%d\n", minS, lastK)
	}

	// offset, _, _ = getOffset((*M)[minS], lastK)
	// h := offset >> wfaTypeBits
	// v := h - uint32(lastK)
	// fmt.Printf("min s:%d, k:%d, h:%d, v:%d\n", minS, lastK, h, v)

	// algn.Plot(q, t, os.Stdout, algn.M, true, true, -1)
	// algn.Plot(q, t, os.Stdout, algn.I, false, true, -1)
	// algn.Plot(q, t, os.Stdout, algn.D, false, true, -1)

	return algn.backTrace(q, t, minS, lastK), nil
}

func (algn *Aligner) backtraceStartPosistion(q, t *[]byte, s uint32) (uint32, int) {
	M := &algn.M
	m, n := len(*t), len(*q)
	minS := s
	Ak := m - n
	lastK := Ak

	var offset uint32
	var ok, kInScope bool
	var k int
	var lastRowOrCol bool
	var h, v int
	for _s := s; _s >= 0; _s-- {
		lastRowOrCol = false

		k = Ak
		// fmt.Printf("a, s:%d, k:%d\n", _s, k)
		for {
			offset, ok, kInScope = getOffset2(M, _s, 0, k)
			// fmt.Printf("  s: %d, k: %d, kInScope: %v\n", _s, k, kInScope)
			if !kInScope {
				if k == 1 {
					k-- // try 0
					continue
				}
				break
			}
			if !ok {
				k--
				continue
			}
			h = int(offset >> wfaTypeBits)
			v = h - k

			if v <= 0 {
				break
			}

			// fmt.Printf("  s: %d, offset:%d, k:%d, h:%d, v:%d\n", _s, offset>>wfaTypeBits, k, h, v)
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
			offset, ok, kInScope = getOffset2(M, _s, 0, k)
			// fmt.Printf("  s: %d, k: %d, kInScope: %v\n", _s, k, kInScope)
			if !kInScope {
				if k == -1 {
					k++ // try 0
					continue
				}
				break
			}
			if !ok {
				k++
				continue
			}
			h = int(offset >> wfaTypeBits)
			v = h - k

			if v <= 0 {
				break
			}

			// fmt.Printf("  s: %d, offset:%d, k:%d, h:%d, v:%d\n", _s, offset>>wfaTypeBits, k, h, v)
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
func (algn *Aligner) extend(offsets *[]uint32, q, t *[]byte, s uint32) (int, int) {
	lo, hi := kLowHigh(offsets)
	// fmt.Printf("  lo: %d, hi: %d, offsets: %d\n", lo, hi, *offsets)

	var offset uint32
	var v, h int
	lenQ := len(*q)
	lenT := len(*t)
	var q8, t8 uint64
	var n, N int

	// processing in revsere order, just to reduce the append operations in setOffsetUpdate,
	// where offsets of k are saved like this:
	//   offset(k=0), offset(k=-1), offset(k=1), offset(k=-2), offset(k=2), ...
	var ok bool
	for k := hi; k >= lo; k-- {
		offset, ok, _ = getOffset(offsets, k)
		// fmt.Printf("    k:%d, ok:%v, offset:%d\n", k, ok, offset>>wfaTypeBits)

		if s > 0 && !ok {
			continue
		}

		h = int(offset >> wfaTypeBits) // x
		v = h - k                      // y
		if v < 0 || v >= lenQ || h >= lenT {
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
			// fmt.Printf("        block wise, %d matches\n", N)
			setOffsetUpdate(offsets, k, uint32(N)<<wfaTypeBits)

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
		// fmt.Printf("      byte wise: k: %d, extend to h: %d, v: %d\n", k, h, v)
		setOffsetUpdate(offsets, k, uint32(N)<<wfaTypeBits)
	}

	return lo, hi
}

// adaptive reduction
func (algn *Aligner) reduce(q, t *[]byte, s uint32, lo, hi int) {
	// fmt.Printf("  lo: %d, hi: %d, offsets: %d\n", lo, hi, *offsets
	offsets := algn.M[s]
	var offset uint32
	var v, h int
	lenQ := len(*q)
	lenT := len(*t)
	var ok bool

	// fmt.Printf("  lo: %d, hi: %d, offsets: %d\n", lo, hi, *offsets)
	var d, minDist int
	// var minDistK int
	ds := poolDist.Get().(*[]int)
	*ds = (*ds)[:0]
	minDist = math.MaxInt
	for k := lo; k < hi; k++ { // must in normal order, cause we apppend d later
		offset, ok, _ = getOffset(offsets, k)
		if s > 0 && !ok {
			*ds = append(*ds, -1)
			continue
		}

		h = int(offset >> wfaTypeBits) // x
		v = h - k                      // y
		if v < 0 || v >= lenQ || h >= lenT {
			continue
		}

		d = max(lenT-h, lenQ-v)
		*ds = append(*ds, d)

		if d < minDist {
			minDist = d
			// minDistK = k
		}
	}
	// offset, _ = getOffset(offsets, minDistK)
	// fmt.Printf("    minD: %d, from k: %d (offset: %d, %d), ds: %d\n", minDist, minDistK, offset, offset>>wfaTypeBits, *ds)

	_lo := lo
	_hi := hi
	maxDistDiff := int(algn.ad.MaxDistDiff)
	updateLo := true
	I := &algn.I
	D := &algn.D
	for i, d := range *ds {
		if d < 0 {
			continue
		}
		if d-minDist > maxDistDiff {
			if updateLo {
				_lo = lo + i + 1
				// fmt.Printf("    new lo: %d, i:%d/%d\n", _lo, i, len(*ds))
			}
			(*ds)[i] = -1 // mark it
			// fmt.Printf("    mark k: %d to delete\n", lo+i)
		} else {
			updateLo = false
		}
	}
	if updateLo { // found some distance where d-minDist > maxDistDiff
		for i := len(*ds) - 1; i >= 0; i-- {
			if (*ds)[i] >= 0 {
				_hi = lo + i
				// fmt.Printf("    new hi: %d, i:%d\n", _hi, i)
				break
			}
		}
	}

	for k := lo; k < _lo; k++ {
		// fmt.Printf("    remove s: %d, k: %d\n", s, k)
		removeOffset(offsets, k)
		removeOffset((*I)[s], k)
		removeOffset((*D)[s], k)
	}
	for k := _hi + 1; k <= hi; k++ {
		// fmt.Printf("    remove s: %d, k: %d\n", s, k)
		removeOffset(offsets, k)
		removeOffset((*I)[s], k)
		removeOffset((*D)[s], k)
	}
	// fmt.Printf("  new lo: %d, hi: %d, offsets: %d\n", _lo, _hi, *offsets)

	poolDist.Put(ds)
}

var poolDist = &sync.Pool{New: func() interface{} {
	tmp := make([]int, 0, 128)
	return &tmp
}}

const (
	// type of the 6 kinds of offsets, which will be saved as the lowest 3bits of the offset.
	wfaInsertOpen uint32 = iota + 1
	wfaInsertExt
	wfaDeleteOpen
	wfaDeleteExt
	wfaMismatch
	wfaMatch // only for backtrace, not saved in the component
)

var wfaOps []byte = []byte{'.', 'I', 'I', 'D', 'D', 'X', 'M', 'H'}
var wfaArrows []rune = []rune{'⊕', '⟼', '🠦', '↧', '🠧', '⬂', '⬊'} // ⬂

func wfaType2str(t uint32) string {
	switch t {
	case wfaInsertOpen:
		return "I.O"
	case wfaInsertExt:
		return "I.E"
	case wfaDeleteOpen:
		return "D.O"
	case wfaDeleteExt:
		return "D.E"
	case wfaMismatch:
		return "Mis"
	case wfaMatch:
		return "Mat"
	default:
		return "N/A"
	}
}

// number of bits to save the path.
const wfaTypeBits uint32 = 3
const wfaTypeMask uint32 = (1 << wfaTypeBits) - 1

// next refers to the WF_NEXT method.
func (algn *Aligner) next(q, t *[]byte, s uint32) {
	// PrintComponent(os.Stdout, algn.M, "M")
	M := &algn.M
	I := &algn.I
	D := &algn.D
	SM := &algn.SM
	SI := &algn.SI
	SD := &algn.SD
	p := algn.p
	lenQ := uint32(len(*q))
	lenT := uint32(len(*t))

	loMismatch, hiMismatch := kLowHigh2(M, s, p.Mismatch)       // M[s-x]
	loGapOpen, hiGapOpen := kLowHigh2(M, s, p.GapOpen+p.GapExt) // M[s-o-e]
	loInsert, hiInsert := kLowHigh2(I, s, p.GapExt)             // I[s-e]
	loDelete, hiDelete := kLowHigh2(D, s, p.GapExt)             // D[s-e]

	hi := max(hiMismatch, hiGapOpen, hiInsert, hiDelete) + 1
	lo := min(loMismatch, loGapOpen, loInsert, loDelete) - 1
	// fmt.Printf("s: %d, k: %d -> %d, lenQ: %d, lenT: %d\n", s, lo, hi, lenQ, lenT)

	var fromM bool
	var fromI, fromD bool
	var v1, v2 uint32
	var Isk, Dsk, Msk uint32
	var updatedI, updatedD, updatedM bool
	var wfaTypeI, wfaTypeD, wfaTypeM uint32
	for k := lo; k <= hi; k++ {
		updatedI, updatedD, updatedM = false, false, false
		wfaTypeI, wfaTypeD, wfaTypeM = 0, 0, 0
		// fmt.Printf(" k: %d\n", k)

		// --------------------------------------
		// insertion: 🠦
		v1, fromM, _ = getOffset2(M, s, p.GapOpen+p.GapExt, k-1)
		v2, fromI, _ = getOffset2(I, s, p.GapExt, k-1)
		v1 >>= wfaTypeBits
		v2 >>= wfaTypeBits
		if fromM && v1 > 1 && v1 == lenT { // it's the last column
			fromM = false
			v1 = 0
		}
		if fromI && v2 > 1 && v2 == lenT { // it's the last column
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
			setOffset2(I, s, k, Isk<<wfaTypeBits|wfaTypeI)
			setOffset2(SI, s, k, Isk<<wfaTypeBits|wfaTypeI) // another copy without extending
			// fmt.Printf("  I start, s:%d, k:%d, offset:%d\n", s, k, Isk)
			// fmt.Printf("  %d fromM:%v(%d), fromI:%v(%d), save I: s=%d, k=%d, offset:%d, type:%s\n",
			// 	Isk, fromM, v1, fromI, v2, s, k, Isk, wfaType2str(wfaTypeI))
		} else {
			Isk = 0
		}

		// --------------------------------------
		// deletion: 🠧

		v1, fromM, _ = getOffset2(M, s, p.GapOpen+p.GapExt, k+1)
		v2, fromD, _ = getOffset2(D, s, p.GapExt, k+1)
		v1 >>= wfaTypeBits
		v2 >>= wfaTypeBits
		if fromM && v1 > 1 && v1-uint32(k) == lenQ { // it's the last row
			fromM = false
			v1 = 0
		}
		if fromD && v2 > 1 && v2-uint32(k) == lenQ { // it's the last row
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
			setOffset2(D, s, k, Dsk<<wfaTypeBits|wfaTypeD)
			setOffset2(SD, s, k, Dsk<<wfaTypeBits|wfaTypeD) // another copy without extending
			// fmt.Printf("  D start, s:%d, k:%d, offset:%d\n", s, k, Dsk)
			// fmt.Printf("  %d fromM:%v(%d), fromD:%v(%d), save D: s=%d, k=%d, offset:%d, type:%s\n",
			// 	Dsk, fromM, v1, fromD, v2, s, k, Dsk, wfaType2str(wfaTypeD))
		} else {
			Dsk = 0
		}

		// --------------------------------------
		// mismatch: ⬂

		v1, fromM, _ = getOffset2(M, s, p.Mismatch, k)
		v1 >>= wfaTypeBits
		if fromM && v1 > 1 && (v1 == lenT || v1-uint32(k) == lenQ) { // it's the last column/row
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

			updatedM = true
			setOffset2(M, s, k, Msk<<wfaTypeBits|wfaTypeM)
			setOffset2(SM, s, k, Msk<<wfaTypeBits|wfaTypeM) // another copy without extending
			// fmt.Printf("  M start, s:%d, k:%d, offset:%d\n", s, k, Msk)
			// fmt.Printf("  %d fromI:%v(%d), fromD:%v(%d), fromM:%v(%d), save M: s=%d, k=%d, offset:%d, type:%s\n",
			// 	Msk, updatedI, Isk, updatedD, Dsk, fromM, v1+1, s, k, Msk, wfaType2str(wfaTypeM))
		}
	}

	// fill with nil, so the score would be the right index
	if !updatedM {
		*M = append(*M, nil)
	}
	if !updatedI {
		*I = append(*I, nil)
	}
	if !updatedD {
		*D = append(*D, nil)
	}
}

// backTrace backtraces the alignment
func (algn *Aligner) backTrace(q, t *[]byte, s uint32, Ak int) *CIGAR {
	// PrintComponent(os.Stdout, algn.M, "M")
	// PrintComponent(os.Stdout, algn.I, "I")
	// PrintComponent(os.Stdout, algn.D, "D")
	M := &algn.M
	I := &algn.I
	D := &algn.D
	SM := &algn.SM
	SI := &algn.SI
	SD := &algn.SD
	p := algn.p
	lenQ := len(*q)
	lenT := len(*t)
	cigar := NewCIGAR()
	cigar.Score = s
	var ok bool
	var k, h, v int
	var offset, wfaType uint32
	var offset0 uint32
	var h0 int
	var op byte
	var n uint32 // number of matched bases
	var qBegin, tBegin int

	// fmt.Printf("backtrace from M%d,%d, lenQ: %d, lenT: %d:\n", s, Ak, lenQ, lenT)

	k = Ak
	firstRecord := true
	firstMatch := true
	semiGlobal := !algn.opt.GlobalAlignment
LOOP:
	for {
		if firstRecord {
			firstRecord = false

			offset, _, _ = getOffset((*M)[s], k)
			// fmt.Printf("------\nfirst s: %d, k: %d, offset: %d, existed: %v\n", s, k, offset>>wfaTypeBits, ok)

			offset0, _, _ = getOffset((*SM)[s], k)
			// fmt.Printf("  get h/h0 from SM\n")

			h = int(offset >> wfaTypeBits)
			v = h - k

			if h < lenT {
				cigar.AddN(wfaOps[wfaInsertOpen], uint32(lenT)-uint32(h))
			} else if v < lenQ {
				cigar.AddN('H', uint32(lenQ)-uint32(v))
			}
		} else {
			// fmt.Printf("------\ncurret s: %d, k: %d\n", s, k)
			switch wfaType {
			case wfaInsertExt:
				offset, ok, _ = getOffset((*I)[s], k)
				offset0, _, _ = getOffset((*SI)[s], k)
				// fmt.Printf("  get h/h0 from SI\n")
			case wfaDeleteExt:
				offset, ok, _ = getOffset((*D)[s], k)
				offset0, _, _ = getOffset((*SD)[s], k)
				// fmt.Printf("  get h/h0 from SD\n")
			default:
				offset, ok, _ = getOffset((*M)[s], k)
				offset0, _, _ = getOffset((*SM)[s], k)
				// fmt.Printf("  get h/h0 from SM\n")
			}
			// fmt.Printf("  offset: %d, existed: %v\n", offset>>wfaTypeBits, ok)

			if !ok {
				// fmt.Printf("  break as there's no valid offset\n")
				break
			}
		}

		wfaType = offset & wfaTypeMask

		h0 = int(offset0 >> wfaTypeBits)
		n = uint32(h - h0)

		// fmt.Printf("  type: %s, h0:%d, h:%d, v:%d\n", wfaType2str(wfaType), h0, h, v)

		if v < 0 || h < 0 || v > lenQ {
			// fmt.Printf("  break as v(%d) < 0 || h(%d) < 0 || v(%d) > lenQ \n", h, v, v)
			break
		}

		if firstMatch { // record the end position of matched region
			if wfaType == wfaMatch || n > 0 {
				firstMatch = false
				cigar.TEnd, cigar.QEnd = h, v
				// fmt.Printf("  == end position of matched region, t:%d, q:%d\n", h, v)
			}
		}

		if n > 0 {
			tBegin, qBegin = h, v

			for v >= 1 && h >= 1 && (*q)[v-1] == (*t)[h-1] {
				// fmt.Printf("    back extend to h:%d, v:%d\n", h, v)
				if h == h0 {
					// fmt.Printf("    reached h0, h:%d\n", h)
					break
				}

				h--
				v--
			}

			op = wfaOps[wfaMatch] // correct it as M
			cigar.AddN(op, n)

			// fmt.Printf("  [ADD %s as match]: %d%c, h:%d, v:%d\n", wfaType2str(wfaType), n, op, h, v)
		} else if wfaType == wfaMatch {
			tBegin, qBegin = h, v
		}
		op = wfaOps[wfaType]
		cigar.AddN(op, 1)

		// fmt.Printf("  [ADD]: %s as %d%c, h:%d, v:%d\n", wfaType2str(wfaType), 1, op, h, v)

		if semiGlobal && (h == 1 || v == 1) {
			// fmt.Printf("  break as reached 1th row/col\n")
			break
		}

		switch wfaType {
		case wfaInsertOpen:
			s -= p.GapOpen + p.GapExt
			k--
			h--
		case wfaInsertExt:
			s -= p.GapExt
			k--
			h--
		case wfaDeleteOpen:
			s -= p.GapOpen + p.GapExt
			k++
			v--
		case wfaDeleteExt:
			s -= p.GapExt
			k++
			v--
		case wfaMismatch:
			s -= p.Mismatch
			h--
			v--
		default:
			// fmt.Printf("  break as invalid wfa type\n")
			break LOOP
		}

		// fmt.Printf("  NEXT s: %d, k: %d, h:%d, v:%d\n", s, k, h, v)
	}

	// fmt.Printf("------\nh:%d, v:%d\n", h, v)

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
