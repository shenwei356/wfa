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

// the number of bits to save the path.
const wfaTypeBits uint32 = 3
const wfaTypeMask uint32 = (1 << wfaTypeBits) - 1

const (
	// type of the 6 kinds of offsets, which will be saved as the lowest 3bits of the offset.
	wfaInsertOpen uint32 = iota + 1
	wfaInsertExt
	wfaDeleteOpen
	wfaDeleteExt
	wfaMismatch
	wfaMatch // only for first row/column
)

var wfaOps []byte = []byte{'.', 'I', 'I', 'D', 'D', 'X', 'M', 'H'} // for backtrace

var wfaArrows []rune = []rune{'⊕', '⟼', '🠦', '↧', '🠧', '⬂', '⬊'} // for visualization

// for showing offsets.
func wfaType2str(t uint32) string {
	switch t {
	case wfaInsertOpen:
		return "I.O"
	case wfaInsertExt:
		return "I.E"
	case wfaDeleteOpen:
		return "D.O"
	case wfaDeleteExt:
		return "D.E"
	case wfaMismatch:
		return "Mis"
	case wfaMatch:
		return "Mat"
	default:
		return "N/A"
	}
}
