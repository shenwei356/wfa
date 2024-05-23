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
	"strconv"
	"sync"
)

// CIGAR represent a CIGAR structure.
type CIGAR struct {
	Ops   []*CIGARRecord
	Score uint32

	reversed bool
}

// CIGARRecord records the operation and the number.
type CIGARRecord struct {
	N  uint32
	Op byte
}

// NewCIGAR returns a new CIGAR from the object pool.
func NewCIGAR() *CIGAR {
	cigar := poolCIGAR.Get().(*CIGAR)
	cigar.reset()
	return cigar
}

// reset resets a CIGAR.
func (cigar *CIGAR) reset() {
	for _, r := range cigar.Ops {
		poolCIGARRecord.Put(r)
	}
	cigar.Ops = cigar.Ops[:0]
	cigar.Score = 0
	cigar.reversed = false
}

// RecycleCIGAR recycle a CIGAR object.
func RecycleCIGAR(cigar *CIGAR) {
	if cigar != nil {
		poolCIGAR.Put(cigar)
	}
}

// object pool of a CIGAR.
var poolCIGAR = &sync.Pool{New: func() interface{} {
	cigar := CIGAR{
		Ops: make([]*CIGARRecord, 0, 128),
	}
	return &cigar
}}

// object pool of CIGARRecord.
var poolCIGARRecord = &sync.Pool{New: func() interface{} {
	return &CIGARRecord{}
}}

// Add adds a new record in backtrace.
func (cigar *CIGAR) Add(op byte) {
	cigar.AddN(op, 1)
}

// Add adds a new record in backtrace and set its number as n.
func (cigar *CIGAR) AddN(op byte, n uint32) {
	r := poolCIGARRecord.Get().(*CIGARRecord)
	r.Op = op
	r.N = n
	cigar.Ops = append(cigar.Ops, r)
}

// Update updates the last record.
func (cigar *CIGAR) Update(n uint32) {
	l := len(cigar.Ops)
	if l > 0 {
		cigar.Ops[l-1].N += n
	}
}

// reverse just reverses the order of all operations.
func (cigar *CIGAR) reverse() {
	if cigar.reversed {
		return
	}
	s := cigar.Ops
	for i, j := 0, len(s)-1; i < j; i, j = i+1, j-1 {
		s[i], s[j] = s[j], s[i]
	}
	cigar.reversed = true
}

// String returns the CIGAR
func (cigar *CIGAR) String() string {
	cigar.reverse()
	var buf bytes.Buffer
	for _, op := range cigar.Ops {
		buf.WriteString(strconv.Itoa(int(op.N)))
		buf.WriteByte(op.Op)
	}
	return buf.String()
}
