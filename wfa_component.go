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

// WAVEFRONTS_BASE_SIZE is the base size of the wavefront slice.
var WAVEFRONTS_BASE_SIZE = 2048

var _WAVEFRONTS_GROW_SLICE = make([]*WaveFront, WAVEFRONTS_BASE_SIZE)

// Component is the wavefront component, it's a list of wavefronts for different scores.
// To support fast access, we use a list to them,
// the nil data means there's no such wavefront for a given score.
type Component struct {
	IsM bool // if it is true, the visualization is slightly different

	WaveFronts []*WaveFront // WaveFronts
}

// NewComponent returns a new Component object.
// If you do not need it, do not remember to use RecycleComponent() to recycle it.
func NewComponent() *Component {
	cpt := poolComponent.Get().(*Component)
	cpt.IsM = false

	cpt.Reset()

	cpt.WaveFronts = cpt.WaveFronts[:WAVEFRONTS_BASE_SIZE]

	return cpt
}

// Reset clears all existing wavefronts for new using.
func (cpt *Component) Reset() {
	for i, wf := range cpt.WaveFronts {
		if wf != nil {
			poolWaveFront.Put(wf)
			cpt.WaveFronts[i] = nil // reset it to nil
		}
	}
}

var poolComponent = &sync.Pool{New: func() interface{} {
	cpt := Component{
		WaveFronts: make([]*WaveFront, WAVEFRONTS_BASE_SIZE), // preset 2048 values.
	}
	return &cpt
}}

// RecycleComponent recycles a Component.
func RecycleComponent(cpt *Component) {
	if cpt != nil {
		poolComponent.Put(cpt)
	}
}

// HasScore tells if a score exists.
func (cpt *Component) HasScore(s uint32) bool {
	if s >= uint32(len(cpt.WaveFronts)) {
		return false
	}
	return cpt.WaveFronts[s] != nil
}

// KRange returns the lowest and highest values of k for score s-diff.
// Since scores are saved in uint32, if diff > s, s-diff would be a large value.
// So we have to check the two values first.
func (cpt *Component) KRange(s, diff uint32) (int, int) {
	if diff > s {
		return 0, 0
	}
	s -= diff
	if s >= uint32(len(cpt.WaveFronts)) || cpt.WaveFronts[s] == nil {
		return 0, 0
	}
	wf := cpt.WaveFronts[s]
	return wf.Lo, wf.Hi
}

// Set sets an offset with a given backtrace type for a score.
func (cpt *Component) Set(s uint32, k int, offset uint32, wfaType uint32) {
	if s >= uint32(len(cpt.WaveFronts)) {
		cpt.WaveFronts = append(cpt.WaveFronts, _WAVEFRONTS_GROW_SLICE...)
	}
	wf := cpt.WaveFronts[s]
	if wf == nil {
		wf = NewWaveFront()
		cpt.WaveFronts[s] = wf
	}

	wf.Set(k, offset, wfaType)
}

// Set sets an offset which has already contain a backtrace type for a score.
// Here, offsetWithType = offset<<wfaTypeBits | wfaType.
func (cpt *Component) SetRaw(s uint32, k int, offset uint32) {
	if s >= uint32(len(cpt.WaveFronts)) {
		cpt.WaveFronts = append(cpt.WaveFronts, _WAVEFRONTS_GROW_SLICE...)
	}
	wf := cpt.WaveFronts[s]
	if wf == nil {
		wf = NewWaveFront()
		cpt.WaveFronts[s] = wf
	}

	wf.SetRaw(k, offset)
}

// Increase increases the offset by delta.
// Here delta does not contain the backtrace type
func (cpt *Component) Increase(s uint32, k int, delta uint32) {
	if s >= uint32(len(cpt.WaveFronts)) {
		cpt.WaveFronts = append(cpt.WaveFronts, _WAVEFRONTS_GROW_SLICE...)
	}
	cpt.WaveFronts[s].Increase(k, delta)
}

// Get returns offset, wfaType, existed.
func (cpt *Component) Get(s uint32, k int) (uint32, uint32, bool) {
	if s >= uint32(len(cpt.WaveFronts)) || cpt.WaveFronts[s] == nil {
		return 0, 0, false
	}
	return cpt.WaveFronts[s].Get(k)
}

// GetRaw returns "offset<<wfaTypeBits |  wfaType", existed.
func (cpt *Component) GetRaw(s uint32, k int) (uint32, bool) {
	if s >= uint32(len(cpt.WaveFronts)) || cpt.WaveFronts[s] == nil {
		return 0, false
	}
	return cpt.WaveFronts[s].GetRaw(k)
}

// GetAfterDiff returns offset, wfaType, existed for s-diff and k.
func (cpt *Component) GetAfterDiff(s uint32, diff uint32, k int) (uint32, uint32, bool) {
	if diff > s {
		return 0, 0, false
	}
	s -= diff
	if s >= uint32(len(cpt.WaveFronts)) || cpt.WaveFronts[s] == nil {
		return 0, 0, false
	}
	return cpt.WaveFronts[s].Get(k)
}

// GetRaw returns "offset<<wfaTypeBits |  wfaType", existed for s-diff and k.
func (cpt *Component) GetRawAfterDiff(s uint32, diff uint32, k int) (uint32, bool) {
	if diff > s {
		return 0, false
	}
	s -= diff
	if s >= uint32(len(cpt.WaveFronts)) || cpt.WaveFronts[s] == nil {
		return 0, false
	}
	return cpt.WaveFronts[s].GetRaw(k)
}

// Delete delete an offset of a s and k.
func (cpt *Component) Delete(s uint32, k int) {
	if s >= uint32(len(cpt.WaveFronts)) || cpt.WaveFronts[s] == nil {
		return
	}
	cpt.WaveFronts[s].Delete(k)
}

// Print lists all offsets for all scores and k values.
func (cpt *Component) Print(wtr io.Writer, name string) {
	var ok bool
	var offset, wfaType uint32

	for _s, wf := range cpt.WaveFronts {
		if wf == nil {
			continue
		}

		fmt.Fprintf(wtr, "%s%d: k[%d, %d]: ", name, _s, wf.Lo, wf.Hi)
		for k := wf.Lo; k <= wf.Hi; k++ {
			offset, wfaType, ok = wf.Get(k)
			if ok {
				fmt.Fprintf(wtr, " k(%d):%d(%s)", k, offset, wfaType2str(wfaType))
			}
		}
		fmt.Fprintln(wtr)
	}
}
