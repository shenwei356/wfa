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
	"os"
	"testing"
)

func TestWFA(_t *testing.T) {
	algn := New()

	// from the paper Bioinformatics, 37(4), 2021, 456–463
	// q := []byte("GAGATA")
	// t := []byte("GATACA")

	// from https://aacbb-workshop.github.io/slides/2022/WFA.ISCA.v6.pdf page15.
	q := []byte("AGGATGCTCG")
	t := []byte("ACCATACTCG")

	algn.Align(&q, &t)

	// wtr := os.Stdout

	fmt.Printf("\n---------------- M ----------------\n")
	for _s, offsets := range algn.M {
		if offsets != nil {
			fmt.Printf("M%d: %d\n", _s, *offsets)
		}
	}
	algn.Plot(&q, &t, os.Stdout, algn.M, true)

	fmt.Printf("\n---------------- I ----------------\n")
	for _s, offsets := range algn.I {
		if offsets != nil {
			fmt.Printf("I%d: %d\n", _s, *offsets)
		}
	}
	algn.Plot(&q, &t, os.Stdout, algn.I, false)

	fmt.Printf("\n---------------- D ----------------\n")
	for _s, offsets := range algn.D {
		if offsets != nil {
			fmt.Printf("D%d: %d\n", _s, *offsets)
		}
	}
	algn.Plot(&q, &t, os.Stdout, algn.D, false)

	RecycleAligner(algn)
}
