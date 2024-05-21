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
	"io"
	"sync"
)

// --------------------------------------------------------------

// kLowHigh is a function to return the lowest and largest k values.
func kLowHigh(offsets *[]uint32) (lo int, hi int) {
	if len(*offsets) == 1 {
		return 0, 0
	}

	l := len(*offsets)
	var v uint32
	// curently we simply read the list twice for simplicy.
	// k:  0, -1 , 1, -2, 2
	for i := 1; i < l; i += 2 {
		v = (*offsets)[i]
		if v != 0 {
			lo = -(i + 1) >> 1
		}
	}
	for i := 2; i < l; i += 2 {
		v = (*offsets)[i]
		if v != 0 {
			hi = i >> 1
		}
	}
	return
}

// kLowHigh2 check s first.
func kLowHigh2(M *[]*[]uint32, s uint32, d uint32) (lo int, hi int) {
	if d > s {
		return 0, 0
	}
	s = s - d
	if int(s) <= len(*M) || (*M)[s] == nil {
		return 0, 0
	}
	return kLowHigh((*M)[s])
}

func setOffset(offsets *[]uint32, k int, offset uint32) {
	// k:  0, -1 , 1, -2, 2
	if k == 0 {
		if len(*offsets) == 0 {
			*offsets = append(*offsets, 0)
		}
		(*offsets)[0] = offset
	} else if k < 0 {
		n := (-k)<<1 - len(*offsets)
		for i := 0; i < n; i++ {
			*offsets = append(*offsets, 0)
		}
		(*offsets)[((-k)<<1)-1] = offset
	} else {
		n := k<<1 - len(*offsets) + 1
		for i := 0; i < n; i++ {
			*offsets = append(*offsets, 0)
		}
		(*offsets)[k<<1] = offset
	}
	// fmt.Printf("s: %d, k: %d, offset: %d, len(offsets)=%d\n", s, k, offset, len(*offsets))
}

func setOffsetUpdate(offsets *[]uint32, k int, delta uint32) {
	// k:  0, -1 , 1, -2, 2
	if k == 0 {
		if len(*offsets) == 0 {
			*offsets = append(*offsets, 0)
		}
		(*offsets)[0] += delta
	} else if k < 0 {
		n := (-k)<<1 - len(*offsets)
		for i := 0; i < n; i++ {
			*offsets = append(*offsets, 0)
		}
		(*offsets)[((-k)<<1)-1] += delta
	} else {
		n := k<<1 - len(*offsets) + 1
		for i := 0; i < n; i++ {
			*offsets = append(*offsets, 0)
		}
		(*offsets)[k<<1] += delta
	}
}

func setOffset2(M *[]*[]uint32, s uint32, k int, offset uint32) {
	if int(s) >= len(*M) {
		n := int(s) + 1 - len(*M)
		for i := 0; i < n; i++ {
			*M = append(*M, nil)
		}
	}
	offsets := (*M)[s]
	if (*M)[s] == nil {
		offsets = poolOffsets.Get().(*[]uint32)
		(*M)[s] = offsets
	}

	setOffset(offsets, k, offset)
}

func getOffset(offsets *[]uint32, k int) (uint32, bool) {
	if offsets == nil {
		return 0, false
	}
	if k == 0 {
		if len(*offsets) < k+1 {
			return 0, false
		}
		return (*offsets)[0], true
	}
	var i int
	if k < 0 {
		i = ((-k) << 1) - 1
		if len(*offsets) < i+1 {
			return 0, false
		}
		return (*offsets)[i], true
	}
	i = k << 1
	if len(*offsets) < i+1 {
		return 0, false
	}
	return (*offsets)[i], true
}

func getOffset2(M *[]*[]uint32, s uint32, d uint32, k int) (uint32, bool) {
	if M == nil || d > s {
		return 0, false
	}
	s = s - d
	return getOffset((*M)[s], k)
}

// --------------------------------------------------------------

// Plot plots one WFA component of the alignment.
func (algn *Aligner) Plot(q, t *[]byte, wtr io.Writer, M []*[]uint32, isM bool) {
	// create the matrix
	m := poolMatrix.Get().(*[]*[]int32)
	for range *t {
		r := poolRow.Get().(*[]int32)
		for range *q {
			*r = append(*r, -1)
		}
		*m = append(*m, r)
	}
	// fill in scores
	var k, _k int
	var offset uint32
	var v, h int
	for s, offsets := range M {
		if offsets == nil {
			continue
		}
		// fmt.Println(*offsets)
		for _k, offset = range *offsets { // k:  0, -1 , 1, -2, 2
			if offset == 0 {
				continue
			}
			if _k == 0 { // 0
				k = 0
			} else if _k&1 == 1 { // 0, -1, 1, -2, 2
				k = -((_k + 1) >> 1)
			} else { // 1, 2
				k = _k >> 1
			}
			// fmt.Printf("  _k:%d, k:%d, offset:%d\n", _k, k, offset)

			if isM {
				for h = int(offset) - 1; h >= 0; h-- {
					v = h - k
					if v < 0 {
						break
					}
					if (*(*m)[v])[h] < 0 { // only set unsetted
						// fmt.Printf("     v:%d, h:%d\n", v, h)
						(*(*m)[v])[h] = int32(s)
					}
					if (*q)[v] != (*t)[h] {
						break
					}
				}
			} else {
				h = int(offset) - 1
				v = h - k
				(*(*m)[v])[h] = int32(s)
			}

		}
	}

	for _, b := range *q {
		fmt.Fprintf(wtr, "\t%c", b)
	}
	fmt.Fprintln(wtr)

	for v, b := range *t {
		fmt.Fprintf(wtr, "%c", b)
		for _, s := range *(*m)[v] {
			if s < 0 {
				fmt.Fprintf(wtr, "\t")
			} else {
				fmt.Fprintf(wtr, "\t%d", s)
			}
		}

		fmt.Fprintln(wtr)
	}

	recycleMatrix(m)
}

var poolMatrix = &sync.Pool{New: func() interface{} {
	tmp := make([]*[]int32, 0, 128)
	return &tmp
}}

var poolRow = &sync.Pool{New: func() interface{} {
	tmp := make([]int32, 0, 128)
	return &tmp
}}

func recycleMatrix(m *[]*[]int32) {
	for _, r := range *m {
		if r != nil {
			*r = (*r)[:0]
			poolRow.Put(r)
		}
	}
	*m = (*m)[:0]
}
