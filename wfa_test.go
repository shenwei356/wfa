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
	// q := []byte("GATACA")
	// t := []byte("GAGATA")

	// from https://aacbb-workshop.github.io/slides/2022/WFA.ISCA.v6.pdf page15.
	q := []byte("ACCATACTCG")
	t := []byte("AGGATGCTCG")

	// from https://github.com/smarco/WFA2-lib
	//    PATTERN    A-GCTA-GTGTC--AATGGCTACT-T-T-TCAGGTCCT
	//                |  ||| |||||    |||||||| | | |||||||||
	//    TEXT       AA-CTAAGTGTCGG--TGGCTACTATATATCAGGTCCT
	//    ALIGNMENT  1M1I1D3M1I5M2I2D8M1I1M1I1M1I9M
	// q := []byte("AGCTAGTGTCAATGGCTACTTTTCAGGTCCT")
	// t := []byte("AACTAAGTGTCGGTGGCTACTATATATCAGGTCCT")

	cigar, err := algn.Align(&q, &t)
	if err != nil {
		_t.Error(err)
		return
	}

	// wtr := os.Stdout

	fmt.Printf("\n---------------- M ----------------\n")
	PrintComponent(os.Stdout, algn.M, "M")
	algn.Plot(&q, &t, os.Stdout, algn.M, true)

	fmt.Printf("\n---------------- I ----------------\n")
	PrintComponent(os.Stdout, algn.I, "I")
	algn.Plot(&q, &t, os.Stdout, algn.I, false)

	fmt.Printf("\n---------------- D ----------------\n")
	PrintComponent(os.Stdout, algn.D, "D")
	algn.Plot(&q, &t, os.Stdout, algn.D, false)

	if cigar != nil {
		fmt.Printf("\nAlignment: %s\n", cigar.String())
		RecycleCIGAR(cigar)
	}
	RecycleAligner(algn)
}
