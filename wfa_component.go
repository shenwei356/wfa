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

const COMPONENT_BASE_SIZE = 1024

var COMPONENT_SLICE = make([]*WaveFront, COMPONENT_BASE_SIZE)

// The wavefront component is a list of wavefronts for different scores.
// To support fast access, we use a list to them,
// the nil data means there's no such wavefront for a given score.
type Component struct {
	IsM bool // if it is, the visualization is slightly different

	WaveFronts []*WaveFront
}

func NewComponent() *Component {
	cpt := poolComponent.Get().(*Component)
	cpt.IsM = false
	for i, wf := range cpt.WaveFronts {
		if wf != nil {
			RecycleWaveFront(wf)
			cpt.WaveFronts[i] = nil // reset it to nil
		}
	}

	cpt.WaveFronts = cpt.WaveFronts[:COMPONENT_BASE_SIZE]

	return cpt
}

func ClearWaveFronts(cpt *Component) {
	cpt.IsM = false
	for i, wf := range cpt.WaveFronts {
		if wf != nil {
			RecycleWaveFront(wf)
			cpt.WaveFronts[i] = nil // reset it to nil
		}
	}
}

var poolComponent = &sync.Pool{New: func() interface{} {
	cpt := Component{
		WaveFronts: make([]*WaveFront, COMPONENT_BASE_SIZE),
	}
	return &cpt
}}

func RecycleComponent(cpt *Component) {
	if cpt != nil {
		poolComponent.Put(cpt)
	}
}

func (cpt *Component) HasScore(s uint32) bool {
	if s >= uint32(len(cpt.WaveFronts)) {
		return false
	}
	return cpt.WaveFronts[s] != nil
}

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

func (cpt *Component) Set(s uint32, k int, offset uint32, _type uint32) {
	if s >= uint32(len(cpt.WaveFronts)) {
		cpt.WaveFronts = append(cpt.WaveFronts, COMPONENT_SLICE...)
	}
	wf := cpt.WaveFronts[s]
	if wf == nil {
		wf = NewWaveFront()
		cpt.WaveFronts[s] = wf
	}

	wf.Set(k, offset, _type)
}

func (cpt *Component) SetRaw(s uint32, k int, offset uint32) {
	if s >= uint32(len(cpt.WaveFronts)) {
		cpt.WaveFronts = append(cpt.WaveFronts, COMPONENT_SLICE...)
	}
	wf := cpt.WaveFronts[s]
	if wf == nil {
		wf = NewWaveFront()
		cpt.WaveFronts[s] = wf
	}

	wf.SetRaw(k, offset)
}

func (cpt *Component) Add(s uint32, k int, delta uint32) {
	if s >= uint32(len(cpt.WaveFronts)) {
		cpt.WaveFronts = append(cpt.WaveFronts, COMPONENT_SLICE...)
	}
	cpt.WaveFronts[s].Add(k, delta)
}

func (cpt *Component) Get(s uint32, k int) (uint32, uint32, bool) {
	if s >= uint32(len(cpt.WaveFronts)) || cpt.WaveFronts[s] == nil {
		return 0, 0, false
	}
	return cpt.WaveFronts[s].Get(k)
}

func (cpt *Component) GetRaw(s uint32, k int) (uint32, bool) {
	if s >= uint32(len(cpt.WaveFronts)) || cpt.WaveFronts[s] == nil {
		return 0, false
	}
	return cpt.WaveFronts[s].GetRaw(k)
}

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

func (cpt *Component) Delete(s uint32, k int) {
	if s >= uint32(len(cpt.WaveFronts)) {
		return
	}
	cpt.WaveFronts[s].Delete(k)
}

// PrintComponent lists the component details.
func (cpt *Component) Print(wtr io.Writer, name string) {
	var ok bool
	var offset, _type uint32

	for _s, wf := range cpt.WaveFronts {
		if wf == nil {
			continue
		}

		fmt.Fprintf(wtr, "%s%d: k[%d, %d]: ", name, _s, wf.Lo, wf.Hi)
		for k := wf.Lo; k <= wf.Hi; k++ {
			offset, _type, ok = wf.Get(k)
			if ok {
				fmt.Fprintf(wtr, " k(%d):%d(%s)", k, offset, wfaType2str(_type))
			}
		}
		fmt.Fprintln(wtr)
	}
}
