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
	"fmt"
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
}

// object pool of aligners.
var poolAligner = &sync.Pool{New: func() interface{} {
	algn := Aligner{
		p: nil,
		M: make([]*[]uint32, 0, 1024),
		I: make([]*[]uint32, 0, 1024),
		D: make([]*[]uint32, 0, 1024),
	}
	return &algn
}}

// RecycleAligner recycles an Aligner object.
func RecycleAligner(algn *Aligner) {
	if algn != nil {
		poolAligner.Put(algn)
	}
}

// New returns a new Aligner with default penalties from the object pool.
// Do not forget to call RecycleAligner after using it.
func New(p *Penalties, opt *Options) *Aligner {
	algn := poolAligner.Get().(*Aligner)
	algn.p = p
	algn.opt = opt
	return algn
}

// reset resets the internal data before each alignment.
func (algn *Aligner) reset() {
	algn.resetComponent(&algn.M)
	algn.resetComponent(&algn.I)
	algn.resetComponent(&algn.D)
}

// poolOffsets is a object pool for offsets.
var poolOffsets = &sync.Pool{New: func() interface{} {
	tmp := make([]uint32, 0, 128)
	return &tmp
}}

// initComponent inilializes a new WFA component.
func (algn *Aligner) initComponent(M *[]*[]uint32, lenQ, lenT int) {
	offsets := poolOffsets.Get().(*[]uint32)
	*offsets = append(*offsets, 0)
	*M = append(*M, offsets)
}

// resetComponent resets a WFA component (a list of offsets).
func (algn *Aligner) resetComponent(M *[]*[]uint32) {
	for _, v := range *M {
		if v != nil {
			*v = (*v)[:0]
			poolOffsets.Put(v)
		}
	}
	*M = (*M)[:0]
}

// ErrEmptySeq means one of the query or target sequence is empty
var ErrEmptySeq error = fmt.Errorf("wfa: invalid empty sequence")

// ErrLens means query seq is longer than target seq
var ErrLens error = fmt.Errorf("wfa: query seq longer than target seq")

// Align performs alignment for two sequence.
// The length of q should be <= that of t.
func (algn *Aligner) Align(q, t *[]byte) (*CIGAR, error) {
	m, n := len(*t), len(*q)

	if n == 0 || m == 0 {
		return nil, ErrEmptySeq
	}
	if n > m {
		return nil, ErrLens
	}

	// reset the stats
	algn.reset()

	algn.initComponent(&algn.M, n, m) // M[0,0] = 0
	if !algn.opt.GlobalAlignment {
		for k := 1; k < m; k++ {
			setOffset2(&algn.M, 0, k, uint32(k)<<wfaTypeBits)
		}
		for k := 1; k < n; k++ {
			setOffset2(&algn.M, 0, -k, 0)
		}
	}

	algn.initComponent(&algn.I, n, m) // I[0,0] = 0
	algn.initComponent(&algn.D, n, m) // D[0,0] = 0

	Ak := m - n
	Aoffset := uint32(m)
	var offset uint32

	M := &algn.M

	var s uint32
	for {
		// fmt.Printf("---------------------- s: %-3d ----------------------\n", s)
		if (*M)[s] != nil {
			// fmt.Printf("extend:\n")
			algn.extend((*M)[s], q, t, s)

			// fmt.Printf("max offset: %d, Aoffset: %d\n", (*(*M)[s])[Ak], Aoffset)

			offset, _ = getOffset((*M)[s], Ak)
			if offset>>wfaTypeBits >= Aoffset {
				break
			}
		}

		s++

		// fmt.Printf("next:\n")
		algn.next(q, t, s)
	}

	lastK := Ak
	minS := s
	if !algn.opt.GlobalAlignment { // find the minimum score on the last line
		var ok bool
		var k int
		var atTheLastLine bool
		for _s := s; _s > 0; _s-- {
			atTheLastLine = false
			for k = -m; k <= m; k++ {
				offset, ok = getOffset2(M, _s, 0, k)
				// fmt.Printf("  s: %d, test: offset:%d, k:%d\n", _s, offset>>wfaTypeBits, k)
				if ok && int(offset>>wfaTypeBits)-k == n {
					// fmt.Printf("  s: %d, ok: offset:%d, k:%d\n", _s, offset>>wfaTypeBits, k)
					atTheLastLine = true
					break
				}
			}

			if atTheLastLine && _s < minS {
				lastK = k
				minS = _s
			}
		}
	}
	// fmt.Printf("min s:%d, k:%d\n", minS, lastK)

	return algn.backTrace(q, t, minS, lastK), nil
}

// extend refers to the WF_EXTEND method.
// The return bool value indicates whether the end of one sequence is reached.
func (algn *Aligner) extend(offsets *[]uint32, q, t *[]byte, s uint32) {
	lo, hi := kLowHigh(offsets)
	// fmt.Printf("  lo: %d, hi: %d, offsets: %d\n", lo, hi, *offsets)

	var offset uint32
	var v, h int
	lenQ := len(*q)
	lenT := len(*t)

	// processing in revsere order, just to reduce the append operations in setOffsetUpdate,
	// where offsets of k are saved like this:
	//   offset(k=0), offset(k=-1), offset(k=1), offset(k=-2), offset(k=2), ...
	var ok bool
	for k := hi; k >= lo; k-- {
		offset, ok = getOffset(offsets, k)
		// fmt.Printf("    k:%d, ok:%v, offset:%d\n", k, ok, offset>>wfaTypeBits)
		if s > 0 && !ok {
			continue
		}

		h = int(offset >> wfaTypeBits) // x
		v = h - k                      // y
		if v < 0 {
			continue
		}
		if v >= lenQ || h >= lenT {
			continue
		}

		for (*q)[v] == (*t)[h] {
			setOffsetUpdate(offsets, k, 1<<wfaTypeBits)
			v++
			h++
			// fmt.Printf("      k: %d, extend to h: %d, v: %d\n", k, h, v)

			if v == lenQ || h == lenT {
				break
			}
		}
	}
}

const (
	// type of the 5 kinds of offsets, which will be saved as the lowest 3bits of the offset.
	wfaInsertOpen uint32 = iota + 1
	wfaInsertExt
	wfaDeleteOpen
	wfaDeleteExt
	wfaMismatch
	wfaMatch // only for backtrace
)

var wfaOps []byte = []byte{'.', 'I', 'I', 'D', 'D', 'X', 'M'}

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

const wfaTypeBits uint32 = 3
const wfaTypeMask uint32 = (1 << wfaTypeBits) - 1

// next refers to the WF_NEXT method.
func (algn *Aligner) next(q, t *[]byte, s uint32) {
	M := &algn.M
	I := &algn.I
	D := &algn.D
	p := algn.p

	loMismatch, hiMismatch := kLowHigh2(M, s, p.Mismatch)       // M[s-x]
	loGapOpen, hiGapOpen := kLowHigh2(M, s, p.GapOpen+p.GapExt) // M[s-o-e]
	loInsert, hiInsert := kLowHigh2(I, s, p.GapExt)             // I[s-e]
	loDelete, hiDelete := kLowHigh2(D, s, p.GapExt)             // D[s-e]

	hi := max(hiMismatch, hiGapOpen, hiInsert, hiDelete) + 1
	lo := min(loMismatch, loGapOpen, loInsert, loDelete) - 1
	// fmt.Printf("s: %d, k: %d -> %d\n", s, lo, hi)

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

		v1, fromM = getOffset2(M, s, p.GapOpen+p.GapExt, k-1)
		v2, fromI = getOffset2(I, s, p.GapExt, k-1)
		v1 >>= wfaTypeBits
		v2 >>= wfaTypeBits
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
			// fmt.Printf("  %d fromM:%v(%d), fromD:%v(%d), save I: s=%d, k=%d, offset:%d, type:%s\n",
			// 	Isk, fromM, v1, fromI, v2, s, k, Isk, wfaType2str(wfaTypeI))
		}

		v1, fromM = getOffset2(M, s, p.GapOpen+p.GapExt, k+1)
		v2, fromD = getOffset2(D, s, p.GapExt, k+1)
		v1 >>= wfaTypeBits
		v2 >>= wfaTypeBits
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
			// fmt.Printf("  %d fromM:%v(%d), fromD:%v(%d), save D: s=%d, k=%d, offset:%d, type:%s\n",
			// 	Dsk, fromM, v1, fromD, v2, s, k, Dsk, wfaType2str(wfaTypeD))
		}

		v1, fromM = getOffset2(M, s, p.Mismatch, k)
		v1 >>= wfaTypeBits
		Msk = max(Isk, Dsk, v1+1)
		if updatedI || updatedD || fromM {
			if updatedI && updatedD && fromM {
				if Msk == Isk {
					wfaTypeM = wfaTypeI
				} else if Msk == Dsk {
					wfaTypeM = wfaTypeD
				} else {
					wfaTypeM = wfaMismatch
				}
			} else if updatedI {
				if updatedD { // updatedI && updatedD && !fromM
					if Msk == Isk {
						wfaTypeM = wfaTypeI
					} else {
						wfaTypeM = wfaTypeD
					}
				} else if fromM { // updatedI && !updatedD && fromM
					if Msk == Isk {
						wfaTypeM = wfaTypeI
					} else {
						wfaTypeM = wfaMismatch
					}
				} else { // updatedI && !updatedD && !fromM
					wfaTypeM = wfaTypeI
				}
			} else if updatedD {
				if fromM { // !updatedI && updatedD && fromM
					if Msk == Dsk {
						wfaTypeM = wfaTypeD
					} else {
						wfaTypeM = wfaMismatch
					}
				} else { // !updatedI && updatedD && !fromM
					wfaTypeM = wfaTypeD
				}
			} else { // !updatedI && !updatedD && fromM
				wfaTypeM = wfaMismatch
			}

			updatedM = true
			setOffset2(M, s, k, Msk<<wfaTypeBits|wfaTypeM)
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
	M := &algn.M
	p := algn.p
	cigar := NewCIGAR()
	var ok bool
	var k, h, v int
	var offset, wfaType uint32
	// fmt.Printf("backtrace from M%d,%d:\n", s, Ak)

	k = Ak
	var op byte
	var n uint32
	var wfaTypePre uint32
LOOP:
	for {
		if wfaTypePre > 0 {
			op = wfaOps[wfaTypePre]
			cigar.AddN(op, 1)
			// fmt.Printf("  addPre. type: %s, op: %c, n=%d, h: %d, v: %d\n",
			// 	wfaType2str(wfaTypePre), op, 1, h+1, v+1)
		}

		offset, ok = getOffset((*M)[s], k)
		if !ok {
			break
		}
		h = int(offset>>wfaTypeBits) - 1
		wfaType = offset & wfaTypeMask
		op = wfaOps[wfaType] // set op temporally

		v = h - k

		n = 0
		for (*q)[v] == (*t)[h] {
			n++

			v--
			h--
			if v < 0 || h < 0 {
				break
			}
		}

		if n > 0 {
			op = wfaOps[wfaMatch] // correct it as M
			cigar.AddN(op, n)
			// fmt.Printf("  addMul. type: %s, op: %c, n=%d, h: %d, v: %d\n",
			// 	wfaType2str(wfaMatch), op, n, h+1, v+1)

			wfaTypePre = wfaType
		} else {
			cigar.AddN(op, 1)
			// fmt.Printf("  addSig. type: %s, op: %c, n=%d, h: %d, v: %d\n",
			// 	wfaType2str(wfaType), op, 1, h, v)

			wfaTypePre = 0
		}

		switch wfaType {
		case wfaInsertOpen:
			s -= p.GapOpen + p.GapExt
			k--
		case wfaInsertExt:
			s -= p.GapExt
			k--
		case wfaDeleteOpen:
			s -= p.GapOpen + p.GapExt
			k++
		case wfaDeleteExt:
			s -= p.GapExt
			k++
		case wfaMismatch:
			s -= p.Mismatch
		default:
			break LOOP
		}

		// fmt.Printf("  new s: %d, k: %d\n", s, k)
	}

	if wfaTypePre > 0 {
		op = wfaOps[wfaTypePre]
		cigar.AddN(op, 1)
		// fmt.Printf("  addPre. type: %s, op: %c, n=%d, h: %d, v: %d\n",
		// 	wfaType2str(wfaTypePre), op, 1, h+1, v+1)
	}

	return cigar
}
