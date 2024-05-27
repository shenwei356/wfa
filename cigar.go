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
	Score uint32 // Alignment score

	TBegin, TEnd int // 0-based location of the alignment in target seq, no including flanking clipping/insertion sequences
	QBegin, QEnd int // 0-based location of the alignment in query seq, no including flanking clipping/insertion sequences

	// Stats of the aligned region, no including flanking clipping/insertion sequences
	AlignLen   uint32
	Matches    uint32
	Gaps       uint32
	GapRegions uint32

	proccessed bool
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
	cigar.proccessed = false

	cigar.AlignLen = 0
	cigar.Matches = 0
	cigar.Gaps = 0
	cigar.GapRegions = 0
}

// RecycleCIGAR recycles a CIGAR object.
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

// process processes the data
func (cigar *CIGAR) process() {
	if cigar.proccessed {
		return
	}
	s := &cigar.Ops

	// reverse the order of all operations.
	var i, j int
	for i, j = 0, len(*s)-1; i < j; i, j = i+1, j-1 {
		(*s)[i], (*s)[j] = (*s)[j], (*s)[i]
	}

	// merge operations of the same type.
	var opPre, op *CIGARRecord
	var newOp bool
	i, j = 0, 0
	opPre = (*s)[0]
	for i = 1; i < len(*s); i++ {
		op = (*s)[i]
		if op.Op == opPre.Op {
			opPre.N += op.N // update count

			if !newOp {
				j = i // mark insert position
				newOp = true
			}
			continue
		}

		if newOp {
			(*s)[j] = op
			j++
		}

		opPre = op
	}
	if j > 0 {
		*s = (*s)[:j]
	}

	// count matches, gaps
	var begin, end int
	for i, op = range *s {
		if op.Op == 'M' {
			begin = i
			break
		}
	}
	for i = len(*s) - 1; i >= 0; i-- {
		op = (*s)[i]
		if op.Op == 'M' {
			end = i
			break
		}
	}
	var alen uint32
	var matches uint32
	var gaps uint32
	var gapRegions uint32

	for i = begin; i <= end; i++ {
		op = (*s)[i]
		alen += op.N
		switch op.Op {
		case 'M':
			matches += op.N
		case 'I', 'D':
			gaps += op.N
			gapRegions++
		}
	}
	cigar.AlignLen = alen
	cigar.Matches = matches
	cigar.Gaps = gaps
	cigar.GapRegions = gapRegions

	cigar.proccessed = true
}

// CIGAR returns the CIGAR string
func (cigar *CIGAR) CIGAR() string {
	cigar.process()
	buf := poolBytesBuffer.Get().(*bytes.Buffer)
	buf.Reset()

	for _, op := range cigar.Ops {
		buf.WriteString(strconv.Itoa(int(op.N)))
		buf.WriteByte(op.Op)
	}

	text := buf.String()
	poolBytesBuffer.Put(buf)
	return text
}

// CIGAR returns the formated alignment strings for Query, Alignment, and Target.
// Do not forget to recycle them with RecycleAlignment().
func (cigar *CIGAR) Alignment(q, t *[]byte) (*[]byte, *[]byte, *[]byte) {
	cigar.process()

	Q := poolBytes.Get().(*[]byte)
	A := poolBytes.Get().(*[]byte)
	T := poolBytes.Get().(*[]byte)

	// var n int
	var h, v int

	v, h = 0, 0
	var i uint32
	for _, op := range cigar.Ops {
		switch op.Op {
		case 'M':
			for i = 0; i < op.N; i++ {
				*Q = append(*Q, (*q)[v])
				*A = append(*A, '|')
				*T = append(*T, (*t)[h])
				v++
				h++
			}
		case 'X':
			for i = 0; i < op.N; i++ {
				*Q = append(*Q, (*q)[v])
				*A = append(*A, ' ')
				*T = append(*T, (*t)[h])
				v++
				h++
			}
		case 'I':
			for i = 0; i < op.N; i++ {
				*Q = append(*Q, '-')
				*A = append(*A, ' ')
				*T = append(*T, (*t)[h])
				h++
			}
		case 'D', 'H':
			for i = 0; i < op.N; i++ {
				*Q = append(*Q, (*q)[v])
				*A = append(*A, ' ')
				*T = append(*T, '-')
				v++
			}
		}
	}

	return Q, A, T
}

// object pool of aligners.
var poolBytesBuffer = &sync.Pool{New: func() interface{} {
	buf := make([]byte, 1024)
	return bytes.NewBuffer(buf)
}}

var poolBytes = &sync.Pool{New: func() interface{} {
	buf := make([]byte, 0, 1024)
	return &buf
}}

// RecycleAlignment recycle alignment strings
func RecycleAlignment(Q, A, T *[]byte) {
	*Q = (*Q)[:0]
	*A = (*A)[:0]
	*T = (*T)[:0]
	poolBytes.Put(Q)
	poolBytes.Put(A)
	poolBytes.Put(T)
}
