// Copyright © 2024-2026 Wei Shen <shenwei356@gmail.com>
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

// Plot plots one WFA component as a tab-delimited text table.
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
func (algn *Aligner) Plot(q, t *[]byte, wtr io.Writer, _M *Component, notChangeToMatch bool, maxScore int) {
	M := algn.M
	I := algn.I
	D := algn.D
	p := algn.p

	lenQ := len(*q)
	lenT := len(*t)
	isM := _M.IsM

	// create the matrix
	m := poolMatrix.Get().(*[]*[]int32)
	for range *q {
		r := poolRow.Get().(*[]int32)
		for range *t {
			*r = append(*r, -1)
		}
		*m = append(*m, r)
	}

	// ----------------------------------------------------------------
	// fill in scores
	for s, wf := range _M.WaveFronts {
		if wf == nil {
			continue
		}
		if maxScore >= 0 && s > maxScore {
			break
		}
		for k := wf.Lo; k <= wf.Hi; k++ {
			offset, wfaType, ok := wf.Get(k)
			if !ok {
				continue
			}

			// Only M extends matches. Follow its recorded predecessor so rejected
			// transitions cannot change the first cell drawn for this wavefront.
			startOffset := offset
			if isM {
				var found bool
				switch wfaType {
				case wfaInsertOpen:
					startOffset, _, found = M.GetAfterDiff(uint32(s), p.GapOpen+p.GapExt, k-1)
					startOffset++
				case wfaInsertExt:
					startOffset, _, found = I.GetAfterDiff(uint32(s), p.GapExt, k-1)
					startOffset++
				case wfaDeleteOpen:
					startOffset, _, found = M.GetAfterDiff(uint32(s), p.GapOpen+p.GapExt, k+1)
				case wfaDeleteExt:
					startOffset, _, found = D.GetAfterDiff(uint32(s), p.GapExt, k+1)
				case wfaMismatch:
					startOffset, _, found = M.GetAfterDiff(uint32(s), p.Mismatch, k)
					startOffset++
				case wfaMatch:
					startOffset, found = 0, true // The global score-zero origin.
					if !algn.opt.GlobalAlignment {
						startOffset = uint32(max(k+1, 1)) // First-row/column seed.
					}
				}
				if !found {
					startOffset = offset
					if !algn.opt.GlobalAlignment && wfaType == wfaMismatch {
						startOffset = uint32(max(k+1, 1)) // A mismatch seed has no predecessor.
					}
				}
			}

			// Draw the initial operation and its extension separately. A leading
			// gap may start outside the matrix; do not move its symbol onto a match.
			startH := int(startOffset) - 1
			for h := int(offset) - 1; h >= startH; h-- {
				v := h - k
				if v < 0 || h < 0 {
					break
				}
				if v >= lenQ || h >= lenT || (*(*m)[v])[h] >= 0 {
					continue // Keep cells already reached at a lower score.
				}
				typ := wfaType
				if isM && !notChangeToMatch && h > startH {
					typ = wfaMatch
				}
				(*(*m)[v])[h] = int32(s)<<wfaTypeBits | int32(typ)
			}
		}
	}

	// ----------------------------------------------------------------
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
