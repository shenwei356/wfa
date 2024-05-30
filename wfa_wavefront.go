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
	"bytes"
	"fmt"
	"math"
	"sync"
)

const OFFSET_BASE_SIZE = 2048

var OFFSET_GROW_SLICE = make([]uint32, OFFSET_BASE_SIZE)

// The wavefront is a list of offsets for different k values.
// We also use the low 3bits to store the backtrace type.
//
// Since k might be negative and usually the values are symmetrical,
// we store them like this:
//
//	index: 0,  1,  2,  3,  4,  5,  6
//	k:     0, -1,  1, -2,  2, -3,  3
//
// if the value is 0, it means there's no records for that k.
type WaveFront struct {
	Lo, Hi  int
	Offsets []uint32
}

func NewWaveFront() *WaveFront {
	wf := poolWaveFront.Get().(*WaveFront)
	wf.Lo = math.MaxInt
	wf.Hi = math.MinInt
	wf.Offsets = wf.Offsets[:OFFSET_BASE_SIZE]
	clear(wf.Offsets) // reset all values as 0's.

	return wf
}

var poolWaveFront = &sync.Pool{New: func() interface{} {
	wf := WaveFront{
		Offsets: make([]uint32, OFFSET_BASE_SIZE),
	}
	return &wf
}}

func RecycleWaveFront(wf *WaveFront) {
	if wf != nil {
		poolWaveFront.Put(wf)
	}
}

func k2i(k int) int {
	if k >= 0 {
		return k << 1
	}
	return ((-k) << 1) - 1
}

func (wf *WaveFront) Set(k int, offset uint32, _type uint32) {
	i := k2i(k)
	if i >= len(wf.Offsets) { // grow the slice
		wf.Offsets = append(wf.Offsets, OFFSET_GROW_SLICE...)
	}
	wf.Offsets[i] = offset<<wfaTypeBits | _type

	// update k range
	wf.Lo = min(wf.Lo, k)
	wf.Hi = max(wf.Hi, k)
}

func (wf *WaveFront) SetRaw(k int, offset uint32) {
	i := k2i(k)
	if i >= len(wf.Offsets) { // grow the slice
		wf.Offsets = append(wf.Offsets, OFFSET_GROW_SLICE...)
	}
	wf.Offsets[i] = offset

	// update k range
	wf.Lo = min(wf.Lo, k)
	wf.Hi = max(wf.Hi, k)
}

func (wf *WaveFront) Add(k int, delta uint32) {
	i := k2i(k)
	if i >= len(wf.Offsets) { // grow the slice
		wf.Offsets = append(wf.Offsets, OFFSET_GROW_SLICE...)
	}
	wf.Offsets[i] += delta << wfaTypeBits

	// update k range
	wf.Lo = min(wf.Lo, k)
	wf.Hi = max(wf.Hi, k)
}

func (wf *WaveFront) Get(k int) (uint32, uint32, bool) {
	if !(k >= wf.Lo && k <= wf.Hi) { // check k range
		return 0, 0, false
	}
	offset := wf.Offsets[k2i(k)]
	return offset >> wfaTypeBits, offset & wfaTypeMask, offset > 0
}

func (wf *WaveFront) GetRaw(k int) (uint32, bool) {
	if !(k >= wf.Lo && k <= wf.Hi) { // check k range
		return 0, false
	}
	offset := wf.Offsets[k2i(k)]
	return offset, offset > 0
}

func (wf *WaveFront) Delete(k int) {
	if !(k >= wf.Lo && k <= wf.Hi) { // check k range
		return
	}
	wf.Offsets[k2i(k)] = 0

	// update k range
	if k == wf.Hi {
		wf.Hi--
	} else if k == wf.Lo {
		wf.Lo++
	}
}

func (wf *WaveFront) String() string {
	var buf bytes.Buffer
	buf.WriteString(fmt.Sprintf("k range: [%d, %d].", wf.Lo, wf.Hi))
	var ok bool
	var offset, _type uint32
	for k := wf.Lo; k <= wf.Hi; k++ {
		offset, _type, ok = wf.Get(k)
		if ok {
			buf.WriteString(fmt.Sprintf(" k(%d):%d(%s)", k, offset, wfaType2str(_type)))
		}
	}
	return buf.String()
}
