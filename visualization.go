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

// Plot plots one WFA component a text table.
//
// A table cell contains the alignment type symbol and the score.
// Symbols:
//
//	⊕    Unknown,
//	⟼    Gap open (Insertion)
//	🠦    Gap extension (Insertion)
//	↧    Gap open (Deletion)
//	🠧    Gap extension (Deletion)
//	⬂    Mismatch
//	⬊    Match
func (algn *Aligner) Plot(q, t *[]byte, wtr io.Writer, M []*[]uint32, isM bool, notChangeToMatch bool, maxScore int) {
	// create the matrix
	m := poolMatrix.Get().(*[]*[]int32)
	for range *q {
		r := poolRow.Get().(*[]int32)
		for range *t {
			*r = append(*r, -1)
		}
		*m = append(*m, r)
	}

	// fill in scores
	var k, _k int
	var offset uint32
	var wfaType uint32
	var v, h int
	var n, vp, hp, v0, h0 int

	S := &algn.SM
	var h00 int
	var offset00 uint32

	for s, offsets := range M {
		if offsets == nil {
			continue
		}

		if maxScore >= 0 && s > maxScore {
			break
		}

		// fmt.Printf("s: %d\n", s)
		for _k, offset = range *offsets { // k:  0, -1 , 1, -2, 2
			if offset == 0 {
				continue
			}
			if _k == 0 { // 0
				k = 0
			} else if _k&1 == 1 { // negative
				k = -((_k + 1) >> 1)
			} else { // 1, 2
				k = _k >> 1
			}
			// fmt.Printf("  _k:%d, k:%d, offset:%d\n", _k, k, offset>>wfaTypeBits)

			wfaType = offset & wfaTypeMask
			h = int(offset>>wfaTypeBits) - 1
			v = h - k

			// fmt.Printf("   v (0-based): %d, h (0-based): %d\n", v, h)
			if v < 0 || h < 0 || v >= len(*m) || h >= len(*(*m)[v]) {
				continue
			}

			if (*(*m)[v])[h] >= 0 { // recorded with a lower score.
				continue
			}

			// fmt.Printf("    fill h (1-based):%d, v (1-based):%d\n", h+1, v+1)
			(*(*m)[v])[h] = int32(s)<<wfaTypeBits | int32(wfaType)

			if !isM || (*q)[v] != (*t)[h] {
				continue
			}

			offset00, _, _ = getOffset2(S, uint32(s), 0, k)
			h00 = int(offset00>>wfaTypeBits) - 1
			// fmt.Printf("    start h: %d, current h: %d\n", h00+1, h+1)

			if h == h00 { // were not extended at all
				continue
			}

			// change it to match
			v0, h0 = v, h
			if !notChangeToMatch {
				// fmt.Printf("    change %s to match: h (1-based):%d, v (1-based):%d\n", wfaType2str(wfaType), h+1, v+1)
				(*(*m)[v0])[h0] = int32(s)<<wfaTypeBits | int32(wfaMatch)
			}
			n = 0
			for {
				h--
				v--
				if v < 0 || h < 0 {
					break
				}
				n++

				if (*(*m)[v])[h] >= 0 {
					continue
				}

				if !notChangeToMatch {
					(*(*m)[v])[h] = int32(s)<<wfaTypeBits | int32(wfaMatch) // mark as match
					// fmt.Printf("    change %s to match: h (1-based):%d, v (1-based):%d\n", wfaType2str(wfaType), h+1, v+1)
				} else {
					(*(*m)[v])[h] = int32(s)<<wfaTypeBits | int32(wfaType)
				}

				vp, hp = v, h // for the last one (or the original one in the normal order), we will restore it.

				if (*q)[v] != (*t)[h] || h == h00 {
					break
				}
			}
			if n == 0 { // just itself
				vp, hp = v0, h0
			}
			if !notChangeToMatch {
				(*(*m)[vp])[hp] = int32(s)<<wfaTypeBits | int32(wfaType) // set back to the original type
				// fmt.Printf("    change %s back: h (1-based):%d, v (1-based):%d\n", wfaType2str(wfaType), h+1, v+1)
			}
		}
	}

	// sequence q

	fmt.Fprintf(wtr, "   \t ")
	for h := range *t {
		fmt.Fprintf(wtr, "\t%3d", h+1)
	}
	fmt.Fprintln(wtr)
	fmt.Fprintf(wtr, "   \t ")
	for _, b := range *t {
		fmt.Fprintf(wtr, "\t%3c", b)
	}
	fmt.Fprintln(wtr)

	for v, b := range *q {
		fmt.Fprintf(wtr, "%3d\t%c", v+1, b) // a base in seq t
		for _, s := range *(*m)[v] {        // a row of the matrix
			if s < 0 {
				fmt.Fprintf(wtr, "\t  .")
			} else {
				fmt.Fprintf(wtr, "\t%c%2d", wfaArrows[s&int32(wfaTypeMask)], s>>int32(wfaTypeBits))
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
	poolMatrix.Put(m)
}
