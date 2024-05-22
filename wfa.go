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
	"sync"
)

// Penalties contains the penalties, Match is 0.
type Penalties struct {
	Mismatch uint32
	GapOpen  uint32
	GapExt   uint32
}

// DefaultPenalties is from the paper.
var DefaultPenalties = Penalties{
	Mismatch: 4,
	GapOpen:  6,
	GapExt:   2,
}

// Aligner is the object for aligning,
// which can apply to multiple pairs of query and ref sequences.
// And it's from a object pool, in case a large number of alignment are needed.
type Aligner struct {
	p *Penalties

	// The Wavefront component is a list of offsets for different scores.
	// So we have two indexes to save the offsets: s (score), k (k).
	// Since k might be negative and usually the values are symmetrical,
	// wo we store them like this:
	//	k value in []uint32: 0, -1, 1, -2, 2
	// if the value is 0, it means there's no records for that k.
	// The function kLowHigh is used to return the lowest and largest k values.
	M, I, D []*[]uint32
}

// object pool of aligners.
var poolAligner = &sync.Pool{New: func() interface{} {
	algn := Aligner{
		p: nil,
		M: make([]*[]uint32, 0, 128),
		I: make([]*[]uint32, 0, 128),
		D: make([]*[]uint32, 0, 128),
	}
	return &algn
}}

// RecycleAligner recycle a Aligner object.
func RecycleAligner(algn *Aligner) {
	poolAligner.Put(algn)
}

// New returns a new Aligner with default penalties from the object pool.
func New() *Aligner {
	algn := poolAligner.Get().(*Aligner)
	algn.p = &DefaultPenalties
	return algn
}

// NewWithPenality returns a new Aligner from the object pool.
func NewWithPenality(p *Penalties) *Aligner {
	algn := poolAligner.Get().(*Aligner)
	algn.p = p
	return algn
}

// reset resets the internal data before alignment.
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
func (algn *Aligner) initComponent(M *[]*[]uint32) {
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

// Align performs alignment for two sequence
func (algn *Aligner) Align(q, t *[]byte) error {
	// reset the stats
	algn.reset()

	algn.initComponent(&algn.M) // M[0,0] = 0
	algn.initComponent(&algn.I) // I[0,0] = 0
	algn.initComponent(&algn.D) // D[0,0] = 0

	m, n := len(*t), len(*q)
	Ak := m - n
	Aoffset := uint32(m)

	M := &algn.M

	var s uint32
	// var reachTheEnd bool
	for {
		// fmt.Printf("---------------------- s: %-3d ----------------------\n", s)
		if (*M)[s] != nil {
			// fmt.Printf("extend:\n")
			algn.extend((*M)[s], q, t, s)

			// fmt.Printf("max offset: %d, Aoffset: %d\n", (*(*M)[s])[Ak], Aoffset)

			// if reachTheEnd {
			// 	break
			// }
			if (*(*M)[s])[Ak]>>wfaTypeBits >= Aoffset {
				break
			}
		}

		s++

		// fmt.Printf("next:\n")
		algn.next(q, t, s)

	}

	return nil
}

// extend refers to the WF_EXTEND method.
func (algn *Aligner) extend(offsets *[]uint32, q, t *[]byte, s uint32) bool {
	lo, hi := kLowHigh(offsets)
	// fmt.Printf("  lo: %d, hi: %d, offsets: %d\n", lo, hi, *offsets)

	var offset uint32
	var v, h int
	lenQ := len(*q)
	lenT := len(*t)
	var reachTheEnd bool
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

		h = int(offset >> wfaTypeBits)
		v = h - k
		if v < 0 {
			continue
		}
		for (*q)[v] == (*t)[h] {
			setOffsetUpdate(offsets, k, 1<<wfaTypeBits)
			v++
			h++
			// fmt.Printf("      k: %d, extend to h: %d, v: %d\n", k, h, v)

			if v == lenQ || h == lenT {
				reachTheEnd = true
				break
			}
		}
	}
	return reachTheEnd
}

const (
	// type of the 5 kinds of offsets, which will be saved as the lowest 3bits of the offset.
	wfaInsertOpen uint32 = iota + 1
	wfaInsertExt
	wfaDeleteOpen
	wfaDeleteExt
	wfaMismatch
)

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
			updatedI = true
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
			setOffset2(I, s, k, Isk<<wfaTypeBits|wfaTypeI)
			// fmt.Printf("  fromM:%v, fromD:%v, save I: s=%d, k=%d, offset:%d, type:%s\n", fromM, fromI, s, k, Isk, wfaType2str(wfaTypeI))
		}

		v1, fromM = getOffset2(M, s, p.GapOpen+p.GapExt, k+1)
		v2, fromD = getOffset2(D, s, p.GapExt, k+1)
		v1 >>= wfaTypeBits
		v2 >>= wfaTypeBits
		Dsk = max(v1, v2)
		if fromM || fromD {
			updatedD = true
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
			setOffset2(D, s, k, Dsk<<wfaTypeBits|wfaTypeD)
			// fmt.Printf("  fromM:%v, fromD:%v, save D: s=%d, k=%d, offset:%d, type:%s\n", fromM, fromD, s, k, Dsk, wfaType2str(wfaTypeD))
		}

		v1, fromM = getOffset2(M, s, p.Mismatch, k)
		v1 >>= wfaTypeBits
		Msk = max(v1+1, Isk, Dsk)
		if updatedI || updatedD || fromM {
			updatedM = true
			if updatedI && updatedD && fromM {
				if Msk == v1+1 {
					wfaTypeM = wfaTypeI
				} else if Msk == Isk {
					wfaTypeM = wfaTypeD
				} else {
					wfaTypeM = wfaMismatch
				}
			} else if updatedI {
				if updatedD { // updatedI && updatedD && !fromM
					if Msk == v1+1 {
						wfaTypeM = wfaTypeI
					} else {
						wfaTypeM = wfaTypeD
					}
				} else if fromM { // updatedI && !updatedD && fromM
					if Msk == v1+1 {
						wfaTypeM = wfaTypeI
					} else {
						wfaTypeM = wfaMismatch
					}
				} else { // updatedI && !updatedD && !fromM
					wfaTypeM = wfaTypeI
				}
			} else if updatedD {
				if fromM { // !updatedI && updatedD && fromM
					if Msk == Isk {
						wfaTypeM = wfaTypeD
					} else {
						wfaTypeM = wfaMismatch
					}
				} else { // !updatedI && updatedD && !fromM
					wfaTypeM = wfaTypeD
				}
			} else { // !updatedI && !updatedD && !fromM
				wfaTypeM = wfaMismatch
			}

			setOffset2(M, s, k, Msk<<wfaTypeBits|wfaTypeM)
			// fmt.Printf("  fromI:%v, fromD:%v, fromM:%v, save M: s=%d, k=%d, offset:%d,type:%s\n", updatedI, updatedD, fromM, s, k, Msk, wfaType2str(wfaTypeM))
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
func (algn *Aligner) backTrace(q, t *[]byte, s uint32) {

}
