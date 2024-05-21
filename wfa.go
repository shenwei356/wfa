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

// New returns a new Aligner with default penalties from the object pool
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

// reset resets the internal data before alignment
func (algn *Aligner) reset() {
	algn.resetOffsets(&algn.M)
	algn.resetOffsets(&algn.I)
	algn.resetOffsets(&algn.D)
}

// poolOffsets is a object pool for offsets
var poolOffsets = &sync.Pool{New: func() interface{} {
	tmp := make([]uint32, 0, 128)
	return &tmp
}}

// resetOffsets resets offsets
func (algn *Aligner) resetOffsets(data *[]*[]uint32) {
	for _, v := range *data {
		if v != nil {
			*v = (*v)[:0]
			poolOffsets.Put(v)
		}
	}
	*data = (*data)[:0]
}

// Align performs alignment for two sequence
func (algn *Aligner) Align(q, t *[]byte) error {
	// reset the stats
	algn.reset()

	m, n := len(*t), len(*q)
	Ak := m - n
	Aoffset := uint32(m)

	// M[0,0] = 0
	M := &algn.M
	offsets := poolOffsets.Get().(*[]uint32)
	*offsets = append(*offsets, 0)
	*M = append(*M, offsets)

	// I[0,0] = 0
	I := &algn.I
	offsets = poolOffsets.Get().(*[]uint32)
	*offsets = append(*offsets, 0)
	*I = append(*I, offsets)

	// D[0,0] = 0
	D := &algn.D
	offsets = poolOffsets.Get().(*[]uint32)
	*offsets = append(*offsets, 0)
	*D = append(*D, offsets)

	var s uint32
	// var reachTheEnd bool
	for {
		// fmt.Printf("--------------------------------\n")
		// fmt.Printf("s: %d\n", s)
		// fmt.Printf("extend:\n")
		if (*M)[s] != nil {
			algn.extend((*M)[s], q, t)

			// fmt.Printf("\nM ----------------\n")
			// for _s, offsets := range *M {
			// 	if offsets != nil {
			// 		fmt.Printf("M%d: %d\n", _s, *offsets)
			// 	}
			// }
			// algn.Plot(q, t, os.Stdout, algn.M, true)
			// fmt.Printf("\nI ----------------\n")
			// for _s, offsets := range *I {
			// 	if offsets != nil {
			// 		fmt.Printf("I%d: %d\n", _s, *offsets)
			// 	}
			// }
			// algn.Plot(q, t, os.Stdout, algn.I, false)
			// fmt.Printf("\nD ----------------\n")
			// for _s, offsets := range *D {
			// 	if offsets != nil {
			// 		fmt.Printf("D%d: %d\n", _s, *offsets)
			// 	}
			// }
			// algn.Plot(q, t, os.Stdout, algn.D, false)

			// fmt.Printf("max offset: %d, Aoffset: %d\n", (*(*M)[s])[Ak], Aoffset)

			// if reachTheEnd {
			// 	break
			// }
			if (*(*M)[s])[Ak] >= Aoffset {
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
func (algn *Aligner) extend(offsets *[]uint32, q, t *[]byte) bool {
	lo, hi := kLowHigh(offsets)
	// fmt.Printf("  lo: %d, hi: %d, offsets: %d\n", lo, hi, *offsets)

	var offset uint32
	var v, h int
	lenQ := len(*q)
	lenT := len(*t)
	var reachTheEnd bool
	for k := lo; k <= hi; k++ {
		offset, _ = getOffset(offsets, k)
		h = int(offset)
		v = h - k
		// fmt.Printf("  for k: %d, v: %d, h: %d\n", k, v, h)

		for (*q)[v] == (*t)[h] {
			setOffsetUpdate(offsets, k, 1)
			v++
			h++

			if v == lenQ || h == lenT {
				reachTheEnd = true
				break
			}
		}
	}
	return reachTheEnd
}

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

	var ok1, ok2, ok3 bool
	var v1, v2 uint32
	var Isk, Dsk, Msk uint32
	var updatedI, updatedD, updatedM bool
	for k := lo; k <= hi; k++ {
		// fmt.Printf(" k: %d\n", k)

		v1, ok1 = getOffset2(M, s, p.GapOpen+p.GapExt, k-1)
		v2, _ = getOffset2(I, s, p.GapExt, k-1)
		Isk = max(v1, v2) + 1
		if ok1 {
			setOffset2(I, s, k, Isk)
			// fmt.Printf("  save I: s=%d, k=%d, offset:%d, %d\n", s, k, Isk, (*I)[s])
			updatedI = true
		}

		v1, ok2 = getOffset2(M, s, p.GapOpen+p.GapExt, k+1)
		v2, _ = getOffset2(D, s, p.GapExt, k+1)
		Dsk = max(v1, v2)
		if ok2 {
			setOffset2(D, s, k, Dsk)
			// fmt.Printf("  save D: s=%d, k=%d, offset:%d, %d\n", s, k, Dsk, (*D)[s])
			updatedD = true
		}

		v1, ok3 = getOffset2(M, s, p.Mismatch, k)
		Msk = max(v1+1, Isk, Dsk)
		if ok1 || ok2 || ok3 {
			setOffset2(M, s, k, Msk)
			// fmt.Printf("  save M: s=%d, k=%d, offset:%d, %d\n", s, k, Msk, (*M)[s])
			updatedM = true
		}

		// fmt.Printf("  Isk: %d, Dsk: %d, Msk: %d\n", Isk, Dsk, Msk)
	}
	if !updatedM {
		*M = append(*M, nil)
	}
	if !updatedI {
		*I = append(*I, nil)
	}
	if !updatedD {
		*D = append(*D, nil)
	}

	// fmt.Printf("  M: %v, Ms: %v\n", *M, (*M)[s])
	// fmt.Printf("  I: %v, Is: %v\n", *I, (*I)[s])
	// fmt.Printf("  D: %v, Ds: %v\n", *D, (*D)[s])
}

// backTrace backtraces the alignment
func (algn *Aligner) backTrace(q, t *[]byte, s uint32) {

}
