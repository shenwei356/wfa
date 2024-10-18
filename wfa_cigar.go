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

// AlignmentResult represent a AlignmentResult structure.
type AlignmentResult struct {
	// Ops   []*CIGARRecord
	Ops []uint64 // 24-bit null + 8 bit Op + 32-bit N

	Score uint32 // Alignment score

	TBegin, TEnd int // 1-based location of the alignment in target seq, no including flanking clipping/insertion sequences
	QBegin, QEnd int // 1-based location of the alignment in query seq, no including flanking clipping/insertion sequences

	// Stats of the aligned region, no including flanking clipping/insertion sequences
	AlignLen   uint32
	Matches    uint32
	Gaps       uint32
	GapRegions uint32

	proccessed bool
}

// // CIGARRecord records the operation and the number.
// type CIGARRecord struct {
// 	N  uint32
// 	Op byte
// }

// Op extracts operation type and count.
func Op(op uint64) (byte, uint32) {
	return byte(op >> 32), uint32(op & MaskLower32)
}

const OpM = uint64('M')
const OpD = uint64('D')
const OpI = uint64('I')
const OpX = uint64('X')
const OpH = uint64('H')
const MaskLower32 = 4294967295

// NewAlignmentResult returns a new CIGAR from the object pool.
func NewAlignmentResult() *AlignmentResult {
	cigar := poolCIGAR.Get().(*AlignmentResult)
	cigar.reset()
	return cigar
}

// reset resets a CIGAR.
func (cigar *AlignmentResult) reset() {
	// for _, r := range cigar.Ops {
	// 	poolCIGARRecord.Put(r)
	// }
	cigar.Ops = cigar.Ops[:0]
	cigar.Score = 0
	cigar.proccessed = false

	cigar.AlignLen = 0
	cigar.Matches = 0
	cigar.Gaps = 0
	cigar.GapRegions = 0
}

// RecycleAlignmentResult recycles a CIGAR object.
func RecycleAlignmentResult(cigar *AlignmentResult) {
	if cigar != nil {
		poolCIGAR.Put(cigar)
	}
}

// object pool of a CIGAR.
var poolCIGAR = &sync.Pool{New: func() interface{} {
	cigar := AlignmentResult{
		// Ops: make([]*CIGARRecord, 0, 128),
		Ops: make([]uint64, 0, 1024),
	}
	return &cigar
}}

// // object pool of CIGARRecord.
// var poolCIGARRecord = &sync.Pool{New: func() interface{} {
// 	return &CIGARRecord{}
// }}

// Add adds a new record in backtrace.
func (cigar *AlignmentResult) Add(op byte) {
	cigar.AddN(op, 1)
}

// Add adds a new record in backtrace and set its number as n.
func (cigar *AlignmentResult) AddN(op byte, n uint32) {
	// r := poolCIGARRecord.Get().(*CIGARRecord)
	// r.Op = op
	// r.N = n
	// cigar.Ops = append(cigar.Ops, r)
	cigar.Ops = append(cigar.Ops, uint64(op)<<32|uint64(n))
}

// Update updates the last record.
func (cigar *AlignmentResult) Update(n uint32) {
	l := len(cigar.Ops)
	if l > 0 {
		// cigar.Ops[l-1].N += n
		cigar.Ops[l-1] += uint64(n)
	}
}

// process processes the data.
func (cigar *AlignmentResult) process() {
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
	// var opPre, op *CIGARRecord
	var opPre, op uint64
	var iPre int
	var newOp bool
	i, j = 0, 0
	opPre = (*s)[0]
	iPre = 0
	for i = 1; i < len(*s); i++ {
		op = (*s)[i]
		// if op.Op == opPre.Op {
		if op>>32 == opPre>>32 {
			// opPre.N += op.N // update count
			(*s)[iPre] = opPre + op&MaskLower32 // update count

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
		iPre = i
	}
	if j > 0 {
		*s = (*s)[:j]
	}

	// count matches, gaps
	var begin, end int
	for i, op = range *s {
		// if op.Op == 'M' {
		if op>>32 == OpM {
			begin = i
			break
		}
	}
	for i = len(*s) - 1; i >= 0; i-- {
		op = (*s)[i]
		// if op.Op == 'M' {
		if op>>32 == OpM {
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
		// alen += op.N
		alen += uint32(op & MaskLower32)
		// switch op.Op {
		switch op >> 32 {
		// case 'M':
		case OpM:
			matches += uint32(op & MaskLower32)
		// case 'I', 'D':
		case OpI, OpD:
			gaps += uint32(op & MaskLower32)
			gapRegions++
		}
	}
	cigar.AlignLen = alen
	cigar.Matches = matches
	cigar.Gaps = gaps
	cigar.GapRegions = gapRegions

	cigar.proccessed = true
}

// CIGAR returns the CIGAR string.
func (cigar *AlignmentResult) CIGAR() string {
	cigar.process()
	buf := poolBytesBuffer.Get().(*bytes.Buffer)
	buf.Reset()

	for _, op := range cigar.Ops {
		// buf.WriteString(strconv.Itoa(int(op.N)))
		buf.WriteString(strconv.Itoa(int(op & MaskLower32)))
		// buf.WriteByte(op.Op)
		buf.WriteByte(byte(op >> 32))
	}

	text := buf.String()
	poolBytesBuffer.Put(buf)
	return text
}

// AlignmentText returns the formated alignment text for Query, Alignment, and Target.
// Do not forget to recycle them with RecycleAlignmentText().
func (cigar *AlignmentResult) AlignmentText(q, t *[]byte) (*[]byte, *[]byte, *[]byte) {
	cigar.process()

	Q := poolBytes.Get().(*[]byte)
	A := poolBytes.Get().(*[]byte)
	T := poolBytes.Get().(*[]byte)

	// var n int
	var h, v int
	// lenQ := len(*q)
	// lenT := len(*t)

	v, h = 0, 0
	// var i uint32
	var i, n uint64
	for _, op := range cigar.Ops {
		n = op & MaskLower32

		// switch op.Op {
		switch op >> 32 {
		// case 'M':
		case OpM:
			// for i = 0; i < op.N; i++ {
			for i = 0; i < n; i++ {
				*Q = append(*Q, (*q)[v])
				*A = append(*A, '|')
				*T = append(*T, (*t)[h])
				v++
				h++
			}
		// case 'X':
		case OpX:
			// for i = 0; i < op.N; i++ {
			for i = 0; i < n; i++ {
				*Q = append(*Q, (*q)[v])
				*A = append(*A, ' ')
				*T = append(*T, (*t)[h])
				v++
				h++
			}
		// case 'I':
		case OpI:
			// for i = 0; i < op.N; i++ {
			for i = 0; i < n; i++ {
				*Q = append(*Q, '-')
				*A = append(*A, ' ')
				*T = append(*T, (*t)[h])
				h++
			}
		// case 'D', 'H':
		case OpD, OpH:
			// for i = 0; i < op.N; i++ {
			for i = 0; i < n; i++ {
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

// RecycleAlignmentText recycle alignment text.
func RecycleAlignmentText(Q, A, T *[]byte) {
	if Q != nil {
		*Q = (*Q)[:0]
		poolBytes.Put(Q)
	}
	if A != nil {
		*A = (*A)[:0]
		poolBytes.Put(A)
	}
	if T != nil {
		*T = (*T)[:0]
		poolBytes.Put(T)
	}
}
