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
	"bytes"
	"fmt"
	"math"
	"slices"
	"sync"
)

// OFFSETS_BASE_SIZE is the base size of the offset slice.
var OFFSETS_BASE_SIZE = 128

// WaveFront is a list of offsets for different k values.
// We also use the low 3 bits to store the backtrace type.
//
// Since k might be negative and usually the values are symmetrical,
// we store them like this (ZigZag encoding):
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
	// Honor a larger base size even when the pool returns an older wavefront.
	if cap(wf.Offsets) < OFFSETS_BASE_SIZE {
		wf.Offsets = make([]uint32, OFFSETS_BASE_SIZE)
	} else {
		wf.Offsets = wf.Offsets[:OFFSETS_BASE_SIZE]
		clear(wf.Offsets)
	}

	return wf
}

var poolWaveFront = &sync.Pool{New: func() interface{} {
	wf := WaveFront{
		Offsets: make([]uint32, OFFSETS_BASE_SIZE),
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
// ZigZag encoding: https://gist.github.com/mfuerstenau/ba870a29e16536fdbaba
func k2i(k int64) int64 {
	// if k >= 0 {
	// 	return k << 1
	// }
	// return ((-k) << 1) - 1
	return (k << 1) ^ (k >> 63)
}

// k2i_inline is an inline version for performance-critical paths
func k2i_inline(k int) int {
	k64 := int64(k)
	return int((k64 << 1) ^ (k64 >> 63))
}

// Set sets an offset with a given backtrace type.
func (wf *WaveFront) Set(k int, offset uint32, wfaType uint32) {
	i := k2i_inline(k)
	wf.growOffsets(i)
	wf.Offsets[i] = offset<<wfaTypeBits | wfaType

	// update k range
	// wf.Lo = min(wf.Lo, k)
	// wf.Hi = max(wf.Hi, k)
	if k < wf.Lo {
		wf.Lo = k
	}
	if k > wf.Hi {
		wf.Hi = k
	}
}

// SetRaw sets an offset which has already contain a backtrace type.
// Here, offsetWithType = offset<<wfaTypeBits | wfaType.
func (wf *WaveFront) SetRaw(k int, offsetWithType uint32) {
	i := k2i_inline(k)
	wf.growOffsets(i)
	wf.Offsets[i] = offsetWithType

	// update k range
	// wf.Lo = min(wf.Lo, k)
	// wf.Hi = max(wf.Hi, k)
	if k < wf.Lo {
		wf.Lo = k
	}
	if k > wf.Hi {
		wf.Hi = k
	}
}

// growOffsets ensures index is writable, using the current configured block size.
// Clear reused capacity so old alignment data cannot reappear after a reset.
func (wf *WaveFront) growOffsets(index int) {
	if index < len(wf.Offsets) {
		return
	}
	oldLen := len(wf.Offsets)
	newLen := (index/OFFSETS_BASE_SIZE + 1) * OFFSETS_BASE_SIZE
	wf.Offsets = slices.Grow(wf.Offsets, newLen-oldLen)[:newLen]
	// Grow preserves old capacity contents, including offsets from prior alignments.
	clear(wf.Offsets[oldLen:])
}

// Increase increases an existing offset by delta, preserving its backtrace type.
// Missing diagonals are left unchanged.
func (wf *WaveFront) Increase(k int, delta uint32) {
	if k < wf.Lo || k > wf.Hi {
		return
	}
	i := k2i_inline(k)
	if wf.Offsets[i] != 0 {
		wf.Offsets[i] += delta << wfaTypeBits
	}
}

// Get returns offset, wfaType, existed.
func (wf *WaveFront) Get(k int) (uint32, uint32, bool) {
	if k < wf.Lo || k > wf.Hi { // check k range
		return 0, 0, false
	}
	offset := wf.Offsets[k2i_inline(k)]
	return offset >> wfaTypeBits, offset & wfaTypeMask, offset > 0
}

// GetRaw returns "offset<<wfaTypeBits |  wfaType", existed.
func (wf *WaveFront) GetRaw(k int) (uint32, bool) {
	if k < wf.Lo || k > wf.Hi { // check k range
		return 0, false
	}
	offset := wf.Offsets[k2i_inline(k)]
	return offset, offset > 0
}

// Delete deletes an offset of a k.
func (wf *WaveFront) Delete(k int) {
	if k < wf.Lo || k > wf.Hi { // check k range
		return
	}
	wf.Offsets[k2i_inline(k)] = 0

	// update k range
	switch k {
	case wf.Hi:
		wf.Hi--
	case wf.Lo:
		wf.Lo++
	}
	// same speed
	//	if k == wf.Hi {
	//		wf.Hi--
	//	} else if k == wf.Lo {
	//		wf.Lo++
	//	}
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
