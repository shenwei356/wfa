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
)

// --------------------------------------------------------------

// kLowHigh is a function to return the lowest and highest k values from a list of offsets,
// where the offsets are saved like this:
//
//	offset(k=0), offset(k=-1), offset(k=1), offset(k=-2), offset(k=2), ...
//
// If the offset is 0, that means there's no record.
func kLowHigh(offsets *[]uint32) (lo int, hi int) {
	if len(*offsets) == 1 {
		return 0, 0
	}

	l := len(*offsets) - 1
	var gotLo, gotHi bool
	// k:  0, -1 , 1, -2, 2
	for i := l; i >= 1; i-- { // start from the end
		if (*offsets)[i] == 0 {
			continue
		}
		if i&1 == 1 {
			if !gotLo {
				lo = -(i + 1) >> 1
				gotLo = true

				if gotHi {
					break
				}
			}
		} else if !gotHi {
			hi = i >> 1
			gotHi = true

			if gotLo {
				break
			}
		}
	}

	return
}

// kLowHigh2 is similar with kLowHigh, but the input is offsets of all scores.
// it returns result for the score 's-d'. So it checks the scores first.
// d
func kLowHigh2(M *[]*[]uint32, s uint32, d uint32) (lo int, hi int) {
	if d > s {
		return 0, 0
	}

	s = s - d
	if int(s) >= len(*M) || (*M)[s] == nil {
		return 0, 0
	}

	return kLowHigh((*M)[s])
}

// setOffset adds an offset of a k.
func setOffset(offsets *[]uint32, k int, offset uint32) {
	// k:  0, -1 , 1, -2, 2
	if k == 0 {
		if len(*offsets) == 0 {
			*offsets = append(*offsets, offset)
		} else {
			(*offsets)[0] = offset
		}
	} else if k > 0 {
		n := k<<1 - len(*offsets) + 1
		for i := 0; i < n; i++ {
			*offsets = append(*offsets, 0)
		}
		(*offsets)[k<<1] = offset
	} else {
		n := (-k)<<1 - len(*offsets)
		for i := 0; i < n; i++ {
			*offsets = append(*offsets, 0)
		}
		(*offsets)[((-k)<<1)-1] = offset
	}
}

// setOffset updates the offset of a k.
func setOffsetUpdate(offsets *[]uint32, k int, delta uint32) {
	// k:  0, -1 , 1, -2, 2
	if k == 0 {
		if len(*offsets) == 0 {
			*offsets = append(*offsets, delta)
		} else {
			(*offsets)[0] += delta
		}
	} else if k > 0 {
		n := k<<1 - len(*offsets) + 1
		for i := 0; i < n; i++ {
			*offsets = append(*offsets, 0)
		}
		(*offsets)[k<<1] += delta
	} else {
		n := (-k)<<1 - len(*offsets)
		for i := 0; i < n; i++ {
			*offsets = append(*offsets, 0)
		}
		(*offsets)[((-k)<<1)-1] += delta
	}
}

// setOffset2 is similar with setOffset, but the input is offsets of all scores.
func setOffset2(M *[]*[]uint32, s uint32, k int, offset uint32) {
	if int(s) >= len(*M) { // fill the list of all offsets
		n := int(s) + 1 - len(*M)
		for i := 0; i < n; i++ {
			*M = append(*M, nil)
		}
	}

	offsets := (*M)[s]
	if (*M)[s] == nil { // creates a list of offsets
		offsets = poolOffsets.Get().(*[]uint32)
		(*M)[s] = offsets
	}

	setOffset(offsets, k, offset)
}

// getOffset returns the offset of a k and if it exists.
func getOffset(offsets *[]uint32, k int) (uint32, bool, bool) {
	if offsets == nil {
		return 0, false, false
	}

	if k == 0 {
		if len(*offsets) < k+1 {
			return 0, false, false
		}
		return (*offsets)[0], (*offsets)[0] > 0, true
	}

	var i int
	if k < 0 {
		i = ((-k) << 1) - 1
		if len(*offsets) < i+1 {
			return 0, false, false
		}
		return (*offsets)[i], (*offsets)[i] > 0, true
	}

	// k > 0
	i = k << 1
	if len(*offsets) < i+1 {
		return 0, false, false
	}
	return (*offsets)[i], (*offsets)[i] > 0, true
}

// removeOffset remove the offset of a k.
func removeOffset(offsets *[]uint32, k int) (uint32, bool) {
	if offsets == nil {
		return 0, false
	}

	var old uint32
	var ok bool
	if k == 0 {
		old, ok = (*offsets)[0], (*offsets)[0] > 0
		(*offsets)[0] = 0
		return old, ok
	}

	var i int
	if k < 0 {
		i = ((-k) << 1) - 1
		if len(*offsets) < i+1 {
			return 0, false
		}

		old, ok = (*offsets)[i], (*offsets)[i] > 0
		(*offsets)[i] = 0
		return old, ok
	}

	// k > 0
	i = k << 1
	if len(*offsets) < i+1 {
		return 0, false
	}
	old, ok = (*offsets)[i], (*offsets)[i] > 0
	(*offsets)[i] = 0
	return old, ok
}

// getOffset2 is similar with getOffset.
func getOffset2(M *[]*[]uint32, s uint32, d uint32, k int) (uint32, bool, bool) {
	if M == nil || d > s {
		return 0, false, false
	}
	s = s - d
	return getOffset((*M)[s], k)
}

// --------------------------------------------------------------

// PrintComponent lists the component details.
func PrintComponent(wtr io.Writer, M []*[]uint32, name string) {
	var _k, k int
	var offset uint32
	for _s, offsets := range M {
		if offsets != nil {
			fmt.Fprintf(wtr, "%s%d:", name, _s)
			for _k, offset = range *offsets {
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

				fmt.Fprintf(wtr, " k(%d):%d(%s)", k, offset>>wfaTypeBits, wfaType2str(offset&wfaTypeMask))
			}
			fmt.Fprintln(wtr)
		}
	}
}
