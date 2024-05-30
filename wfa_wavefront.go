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

// OFFSETS_BASE_SIZE is the base size of the offset slice.
const OFFSETS_BASE_SIZE = 2048

var _OFFSETS_GROW_SLICE = make([]uint32, OFFSETS_BASE_SIZE)

// WaveFront is a list of offsets for different k values.
// We also use the low 3 bits to store the backtrace type.
//
// Since k might be negative and usually the values are symmetrical,
// we store them like this:
//
//	index: 0,  1,  2,  3,  4,  5,  6
//	k:     0, -1,  1, -2,  2, -3,  3
//
// if the value is 0, it means there's no records for that k.
type WaveFront struct {
	Lo, Hi  int      // Lowest and Highest k.
	Offsets []uint32 // offset data. We preset 2048 values to avoid frequent append operation.
}

// NewWaveFront creates a new WaveFront object.
// If you do not need it, do not remember to use RecycleWaveFront() to recycle it.
func NewWaveFront() *WaveFront {
	wf := poolWaveFront.Get().(*WaveFront)
	wf.Lo = math.MaxInt
	wf.Hi = math.MinInt
	wf.Offsets = wf.Offsets[:OFFSETS_BASE_SIZE] // reset the length to the base size
	clear(wf.Offsets)                           // reset all values as 0's.

	return wf
}

var poolWaveFront = &sync.Pool{New: func() interface{} {
	wf := WaveFront{
		Offsets: make([]uint32, OFFSETS_BASE_SIZE), // preset 2048 values.
	}
	return &wf
}}

// RecycleWaveFront recycles a WaveFront.
func RecycleWaveFront(wf *WaveFront) {
	if wf != nil {
		poolWaveFront.Put(wf)
	}
}

// convert k to slice index.
func k2i(k int) int {
	if k >= 0 {
		return k << 1
	}
	return ((-k) << 1) - 1
}

// Set sets an offset with a given backtrace type.
func (wf *WaveFront) Set(k int, offset uint32, wfaType uint32) {
	i := k2i(k)
	if i >= len(wf.Offsets) { // grow the slice
		n := (i - len(wf.Offsets) + OFFSETS_BASE_SIZE) / OFFSETS_BASE_SIZE
		for j := 0; j < n; j++ {
			wf.Offsets = append(wf.Offsets, _OFFSETS_GROW_SLICE...)
		}
	}
	wf.Offsets[i] = offset<<wfaTypeBits | wfaType

	// update k range
	wf.Lo = min(wf.Lo, k)
	wf.Hi = max(wf.Hi, k)
}

// Set sets an offset which has already contain a backtrace type.
// Here, offsetWithType = offset<<wfaTypeBits | wfaType.
func (wf *WaveFront) SetRaw(k int, offsetWithType uint32) {
	i := k2i(k)
	if i >= len(wf.Offsets) { // grow the slice
		n := (i - len(wf.Offsets) + OFFSETS_BASE_SIZE) / OFFSETS_BASE_SIZE
		for j := 0; j < n; j++ {
			wf.Offsets = append(wf.Offsets, _OFFSETS_GROW_SLICE...)
		}
	}
	wf.Offsets[i] = offsetWithType

	// update k range
	wf.Lo = min(wf.Lo, k)
	wf.Hi = max(wf.Hi, k)
}

// Increase increases the offset by delta.
// Here delta does not contain the backtrace type.
func (wf *WaveFront) Increase(k int, delta uint32) {
	i := k2i(k)
	if i >= len(wf.Offsets) { // grow the slice
		n := (i - len(wf.Offsets) + OFFSETS_BASE_SIZE) / OFFSETS_BASE_SIZE
		for j := 0; j < n; j++ {
			wf.Offsets = append(wf.Offsets, _OFFSETS_GROW_SLICE...)
		}
	}
	wf.Offsets[i] += delta << wfaTypeBits

	// update k range
	wf.Lo = min(wf.Lo, k)
	wf.Hi = max(wf.Hi, k)
}

// Get returns offset, wfaType, existed.
func (wf *WaveFront) Get(k int) (uint32, uint32, bool) {
	if !(k >= wf.Lo && k <= wf.Hi) { // check k range
		return 0, 0, false
	}
	offset := wf.Offsets[k2i(k)]
	return offset >> wfaTypeBits, offset & wfaTypeMask, offset > 0
}

// GetRaw returns "offset<<wfaTypeBits |  wfaType", existed.
func (wf *WaveFront) GetRaw(k int) (uint32, bool) {
	if !(k >= wf.Lo && k <= wf.Hi) { // check k range
		return 0, false
	}
	offset := wf.Offsets[k2i(k)]
	return offset, offset > 0
}

// Delete delete an offset of a k.
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

// String lists all the offsets.
func (wf *WaveFront) String() string {
	var buf bytes.Buffer
	buf.WriteString(fmt.Sprintf("k range: [%d, %d].", wf.Lo, wf.Hi))
	var ok bool
	var offset, wfaType uint32
	for k := wf.Lo; k <= wf.Hi; k++ {
		offset, wfaType, ok = wf.Get(k)
		if ok {
			buf.WriteString(fmt.Sprintf(" k(%d):%d(%s)", k, offset, wfaType2str(wfaType)))
		}
	}
	return buf.String()
}
