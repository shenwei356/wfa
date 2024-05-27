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
	"os"
	"testing"
)

func TestWFA(_t *testing.T) {
	algn := New(
		&Penalties{
			Mismatch: 4,
			GapOpen:  6,
			GapExt:   2,
		},
		&Options{
			GlobalAlignment: true,
		})
	algn.AdaptiveReduction(&AdaptiveReductionOption{
		MinWFLen:    10,
		MaxDistDiff: 50,
		CutoffStep:  1,
	})

	// from the paper Bioinformatics, 37(4), 2021, 456–463
	// q := []byte("GATACA")
	// t := []byte("GAGATA")

	// from https://aacbb-workshop.github.io/slides/2022/WFA.ISCA.v6.pdf page15.
	q := []byte("ACCATACTCG")
	t := []byte("AGGATGCTCG")

	// t := []byte("ccccccccccccccccccccccaACCATACTCGaaaaaaaaaaaaaaaaaaaaa")
	// q := []byte("AGGATGCTCG")

	// demo in readme
	// q := []byte("ACGATCTCG")
	// t := []byte("CAGGCTCCTCGG")

	// from https://aacbb-workshop.github.io/slides/2022/WFA.ISCA.v6.pdf page15.
	// q := []byte("TCTTTACTCGCGCGTTGGAGAAATACAATAGT")
	// t := []byte("TCTATACTGCGCGTTTGGAGAAATAAAATAGT")

	// q := []byte("TGGTGAATACCACAATTAGGACAAGTCCAGTTTCTATCATCTAATGTTAGTTTGTCCGTCCCATTAGTGCCCATTACAAAGCCACAATCGCGACATGTTTCGGTAGTATTTCTTGGGTTGATTGTGATAAATTGACGCCCATAAAGCTTTGCTTTGTAAGCCAACATGCCAAGAAATGTTCGCCAGCCAACGTCAGAAATACTAAGTGCCAAAGCATGATTTTTAAGCAAATTCTTGCTACGCAACTCCTCGGCTACTACTAAATCGTGGTTCTTGATTAATGCGGTAGAAATTTGTTGGAGAAAATTATGCCTTTGGTTCATTCCTTGGCATGAAGTTAGCGACTAACAAGCGCTGTTTTTGATAATTTTTACTATCTCGTAAAGAACGATGTTCTTTTTTGCACGTTGTTGCCGTCTAGATAAAATGCGCTGTTCTTTGGCTAATTTGCCTTTAATAGTGCGGTAATATCGTGGATTAGGAACTATGTTGCCTTCACTGTCGGTTAAGAAGTTATCAGTATTAAGATCAATTCCAACATGTCCATGAGTAGCTTTGGACACTTTAACAAAAGATTCATCTGAAGCTAACTGCATTGATAAAAAGAAGCGATCCGCTGAATCTTTAGTCAACGTCACGGTACCAATTCTAGTCTCGCATATTCTTTTCAAAAGCCGTGCTTG")
	// t := []byte("ATGGTGAATACCACAGTTAGGACAAGTCCATTTTCGGTCATCTAACATTAGTTTGTCCGTCCCATTAGTGCCCATCACAAAGCCACAATCGCAACATGTTTGGGTAGTGTTTCTTGGATTGATTGTGATAAATTGACGCCCATAAAGCTTTGCTTTATAAGCCAACATGCCAAGAAATGTTCGCCAACCAACGTCAGAAATACTAAGTGCCAAAGCATGATTTTTAAGCATATTCTTGCTACGCAACTCTTCGGCTACTACTAAATCGTGGTTCTTGATTAATGCAGTAGAGATTTGTTGGAGAAAATTATGTCTTTGGTTCATTACTTTGGCATGAAGTTTGGCAACTAGCAAGCGTTGCTTTTGGTAATTTTTACTATCGCGTAAAGAACGATGTTCTTTTTTTGCACGTTGTTGCCGTCTAGATAAAATGCGCTGTTCTTTGGCTAATTTACCTTTAATAGTGCGATAATATCGTGGATTAGGAACTATGTTACCTCCACTGTCGGTTAAGAAGTTATCAGTATTAAGATCAATTCCAACATGCCCATGAGTAGCTTTGGCCACTTTAACAAAAGATTCATCTGAAGCTAACTGCATTGATAAGAAGAAGCGATCCGCTGAATCTTTAGTCAACGTCACAGTACCAATTCTAGTCTCGCATATTCTTTTCAAAAGCCGTGCTTG")

	// q := []byte("GGATTAGCTAGTAGGTGGGGTAACGGCTCACCTAGGCGACGATCCCTAGCTGGTCTGAGAGGATGACCAGCCACACTGGAACTGAGACACGGTCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGCACAATGGGCGCAAGCCTGATGCAGCCATGCCGCGTGTATGAAGAAGGCCTTCGGGTTGTAAAGTACTTTCAGCGGGGAGGAAGGGAGTAAAGTTAATACCTTTGCTCATTGACGTTACCCGCAGAAGAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTTTGTTAAGTCAGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATCTGATACTGGCAAGCTTGAGTCTCGTAGAGGGGGGTAGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAATACCGGTGGCGAAGGCGGCCCCCTGGACGAAGACTGACGCTCAGGTGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGTCGACTTGGAGGTTGTGCCCTTGAGGCGTGGCTTCCGGAGCTAACGCGTTAAGTCGACCGCCTGGGGAGTACGGCCGCAAGGTTAAAACTCAAATGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGATGCAACGCGAAGAACCTTACCTGGTCTTGACATCCACGGAAGTTTTCAGAGATGAGAATGTGCCTTCGGGAACCGTGAGACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTTGTGAAATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTATCCTTTGTTGCCAGCGGTCCGGCCGGGAACTCAAAGGAGACTGCCAGTGATAAACTGGAGGAAGGTGGGGATGACGTCAAGTCATCATGGCCCTTACGACCAGGGCTACACACGTGCTACAATGGCGCATACAAAGAGAAGCGACCTCGCGAGAGCAAGCGGACCTCATAAAGTGCGTCGTAGTCCGGATTGGAGTCTGCAACTCGACTCCATGAAGTCGGAATCGCTAGTAATCGTGGATCAGAATGCCACGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCATGGGAGT")
	// t := []byte("GGATTAGCTAGTTGGTGAGGTAACGGCTCACCAAGGCGACGATCCATAGCTGGTCTGAGAGGATGATCAGCCACACTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGGACAATGGGGGCAACCCTGATCCAGCCATGCCGCGTGAGTGATGAAGGCCCTAGGGTTGTAAAGCTCTTTCACCGGTGAAGATAATGACGGTAACCGGAGAAGAAGCCCCGGCTAACTTCGTGCCAGCAGCCGCGGTAATACGAAGGGGGCTAGCGTTGTTCGGATTTACTGGGCGTAAAGCGCACGTAGGCGGATATTTAAGTCAGAGGTGAAATCCCAGGGCTCAACCCTGGAACTGCCTTTGATACTGGATATCTTGAGTGTGGAAGAGGTGAGTGGAATTCCGAGTGTAGAGGTAAAATTCGTAGATATTCGGAGGAACACCAGTGGCGAAGGCGGCTCACTGGTCCATTACTGACGCTGAGGTGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGAATGTTAGCCGTCGGGCGGTTTACTGCTCGGTGGCGCAGCTAACGCATTAAACATTCCGCCTGGGGAGTACGGTCGCAAGATTAAAACTCAAAGGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGCAGAACCTTACCAGCCCTTGACATCCCGATCGCGGAGAGTGGAGACACCCTCCTTCAGTTAGGCTGGATCGGAGACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTCGCCCTTAGTTGCCAGCATTAAGTTGGGCACTCTAAGGGGACTGCCGGTGATAAGCCGAGAGGAAGGTGGGGATGACGTCAAGTCCTCATGGCCCTTACGGGCTGGGCTACACACGTGCTACAATGGTGGTGACAGTGGGCAGCGAGACCGCGAGGTCGAGCTAATCTCCAAAAGCCATCTCAGTTCGGATTGCACTCTGCAACTCGAGTGCATGAAGTTGGAATCGCTAGTAATCGTGGATCAGCATGCCACGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCATGGGAGT")

	// q := []byte("ATGTCAAAAATTATTGAACTGCCAGAAGTTTTAGCCAATCAAATTGCTGCTGGAGAGGTTATTGAAAGACCAGCTAGTGTCGTAAAAGAACTAGTTGAGAATGCTATTGACGCTGATAGCAGCCAAATTACAATAGAAATTGAAGAATCAGGACTGAAAAAAATTCAAGTGACTGATAATGGTCAGGGAATAGAACAAGCAGATGTAATTATGAGCTTACGTCGTCATGCAACTAGTAAGATTAAAAAACAATCCGATCTTTTTAGAATTAGGACCTTAGGTTTTAGAGGTGAGGCCCTTCCATCAATTGCTTCTATTAGCAGACTAACATTAAAAACAGCAACAAAAGGGGAAATACATGGCACTTTAGTAATCGCAAATGGAGGTAAAATTGAAAAAGAAGAAGCCATTAGTACCCCTATTGGAACCAAAGTGAGTGTTGAAAATCTTTTCTTCAATACCCCTGCTCGCCTTAAATATATGAAAAGTTTGCAGGCAGAGCTTGCTCATATTATTGATGTCATTAACCGTCTTAGCTTGGCGCATCCCGAAATTGCTTTCACGCTTATTAATGATGGCCGTGAATTAACAAAGACGGCAGGTAAAGGCGATCTGCGTCAAGCTTTAGCGGGTATTTATGGGGTAACCACTGCTAAAAAAATGGTCGAAATTTCAAATGCAGATCTTGATTTTGAGGTATCTGGTTATATTAGTTTGCCTGAATTAACACGAGCTAATCGCAATTATATTACCATTCTTATTAATGGGCGTTATATCAAAAATTTTTTGCTCAATCGTGCTATTCTTGATGGTTATGGTTCAAAATTAATGGTGGGACGTTTTCCAATTGCTGTTATTGATATTCAGATTGACCCTTATTTAGCAGATGTTAATGTTCATCCAACTAAGCAGGAAGTTCGCATTTCTAAAGAAAAAGAGTTAATGAATCTCATTAGTTCGGCAATTGCAGATAGTTTAAGAGAACAAGAGTTAATTCCTGATGCTTTGGAAAATCTTGCTAAAACTAGCACAAACAAGAAACAAAAGTTTGAACAAACACAGCTTCCTCTTAAACAATCAAATCTTTACTATGATAAGACTCAAAATGACTTCTTTCTCAAAGAAACGACAGTTTCAGAAGCATCAAGCGAGTTTACAGTAGTTGATAAAGCTGTAAATATTACTGATAATGTATCAGGACATAGTAGTGTAAAATATGCTCAACGTCAAATGAGAACTTGTGAAAACAAAGAGCATCCAAGTCTTACTATGACAGATAAGCAACAAAAGCGACAGCTTAGTAGAATTGTTGAAAGTCTAGAAAATGAAGAAAAGTCAACTTTTCCTGAGTTAGAATATTTTGGTCAAATGCATGGAACTTACCTCTTTGCTCAAGGTGAAGGGGGGCTTTATATTATAGATCAACATGCTGCTCAAGAACGTGTTAAATATGAATATTATCGTGATAAAATTGGAGAAGTAGACAATAGTTTACAGCAGCTTTTAGTACCTTATCTTTTTGAATTTCCTGCTAATGATTTTATGACCTTACAAGAAAAAATGAGCATTTTAAATGAGGTGGGAATTCACTTAGAAAACTATGGAGAGAATACTTTTATTTTGAGGGAACATCCTATTTGGTTACAGGAAAAAGAAATTGAATCAGCCGTTTATGAGATGTGTGACATGTTGCTGTTGACAAATGAAGTCTCCATTAAGAAATATCGAGCGGAATTAGCTATTATGATGTCTTGCAAGCGTTCTATCAAGGCTAATCATACTTTGGATGATTATTCTGCCAGAAATCTCTTAATTCAATTATCACAGTGCAAAAATCCTTACAATTGTCCTCATGGACGTCCTGTGCTTGTAAATTTCACAAAAGCTGATATGGAAAAGATGTTTCGTCGCATTCAGGAAAACCATACTAGTTTGAGAGAATTAGG")
	// t := []byte("ATGTCAAAAATTATTGAACTGCCAGAAGTTTTAGCCAATCAAATTGCTGCTGGAGAGGTTATTGAAAGACCAGCTAGTGTCGTAAAAGAACTAGTTGAGAATGCTATTGACGCTGATAGCAGCCAAATTACAATAGAAATTGAAGAATCAGGACTGAAAAAAATTCAAGTGACTGATAATGGTCAAGGAATAGAACAAGCGGATGTAATTATGAGCTTACGTCGTCATGCAACTAGTAAGATTAAAAAACAATCCGATCTTTTTAGAATTAGGACCTTAGGTTTTAGAGGTGAGGCCCTTCCATCAATTGCTTCTATTAGCAGACTAACATTAAAAACAGCAACAAAAGGGGAAATATATGGCACTTTATTAATCGCAAATGGAGGTAAAATTGAAAAAGAAGAAGCCATTAGTACCCCTATTGGAACCAAAGTGAGTGTTGAAAATCTTTTCTTCAATACCCCTGCTCGCCTTAAATATATGAAAAGTTTGCAGGCAGAGCTTGCTCATATTATTGATGTCATTAACCGTCTTAGCTTGGCGCATCCCGAAATTGCTTTCACGCTTATTAATGATGGCCGTGAATTAACAAAGACGTCAGGTAAAGGCGATCTGCGTCAAGCTTTAGCGGGTATTTATGGGGTAACCACTGCTAAAAAAATGGTCGAAATTTCAAATGCAGATCTTGATTTTGAGGTATCTGGTTATATTAGTTTGCCTGAATTAACACGAGCTAATCGCAATTATATTACCATTCTTATTAATGGGCGTTATATCAAAAATTTTTTGCTCAATCGCGCTATTCTTGATGGTTATGGTTCAAAATTAATGGTGGGACGTTTTCCAATTGCTGTTATTGATATTCAGATTGACCCTTATTTAGCAGATGTTAATGTTCATCCAACTAAGCAGGAAGTTCGCATTTCTAAAGAAAAAGAGTTAATGAATCTCATTAGTTCGGCAATTGCAGATAGTTTAAGAGAACAAGAGTTAATTCCTGATGCTTTGGAAAATCTTGCTAAAACTAGCACAAACAAGAAACAAAAGTTTGAACAAACACAACTTCCTCTTAAACAATCAAATCTTTACTATGATAAGACTCAAAATGACTTCTTTCTCAAAGAAACGACAGTTTCAGAAGCATCAAGCGAGTTTACAGTAGTTGATAAAGCTGTAAATATTACTGATAATGTATCAGGACATAGCAGTGTAAAATATGCTCAACGTCAAATGAGAACTTGTGAAAACAAAGAGCATCCAAGTCTTACTATGACAAATAAGCAACAAAAGCGACAGCTTAGTAGAATTGTTGAAAGTCTAGAAAATGAAGAAAAGTCAACTTTTCCTGAGTTAGAATATTTTGGTCAAATGCATGGAACTTACCTCTTTGCTCAAGGTGAAGGGGGGCTTTATATCATAGATCAACATGCTGCTCAAGAACGTGTTAAATATGAATATTATCGTGATAAAATTGGAGAAGTAGACAATAGTTTACAGCAGCTTTTAGTACCTTATCTTTTTGAATTTCCTGCTAATGATTTTATGACCTTACAAGAAAAAATGAGCATTTTAAATGAGGTTGGAGTTCACTTAGAAAACTATGGAGAGAATACTTTCATTTTGAGAGAACATCCTATTTGGTTACAGGAAAAAGAAATTGAATCAGCCGTTTATGAGATGTGTGACATGTTGCTGTTGACAAATGAAGTCTCCATTAAGAAATATCGAGCGGAATTAGCTATTATGATGTCTTGCAAGCGTTCTATCAAGGCTAATCATACTTTGGATGATTATTCTGCCAGAAATCTCTTAATTCAATTATCACAGTGCAAAAATCCTTACAACTGTCCTCATGGACGTCCTGTGCTTGTGAATTTCACAAAAGCTGATATGGAAAAGATGTTTCGTCGCATTCAGGAAAACCATACTAGTTTGAGAGAATTAGG")

	// from https://github.com/smarco/WFA2-lib
	//    PATTERN    AGCTA-GTGTCAATGGCTACT---TTTCAGGTCCT
	//               | ||| |||||  ||||||||   | |||||||||
	//    TEXT       AACTAAGTGTCGGTGGCTACTATATATCAGGTCCT
	//    ALIGNMENT  1M1X3M1I5M2X8M3I1M1X9M
	// q := []byte("AGCTAGTGTCAATGGCTACTTTTCAGGTCCT")
	// t := []byte("AACTAAGTGTCGGTGGCTACTATATATCAGGTCCT")

	// q := []byte("ATTAAAATTATCATAACCAGTTGTCTGACATTTATAAAAATAAAATACCGCTGACCGCAGACTACGAGTCGGGGTTGTATGCTAGAGAGTGGACTGAAC")
	// t := []byte("GAATTAAAATTATCATAACCAGTTGTCTGATCATATAAAAATAAAATACCGCTGACCGCAGACTACGAGTCGGGGTTGTATGCTAGAGAGTGGACTGAAC")

	// q := []byte("ACTATAAGCGTCCTCTGCGAGACCGGATGCGTTGATGACAGCGAATTGAGTTGAACTCCCTAAGGACACTCAATAATATTGGTCTATGCAAAAAGTCATT")
	// t := []byte("CCTATAAGCGTCCTCTGCGAGACCGGATGCGTTGATGACAGCGAATTGAGTTGAACTCCCTAAGGACACTCAATAATATTGGTCTATGCAAAAAGTCATT")

	// q := []byte("GC")
	// t := []byte("AC")

	// q := []byte("C")
	// t := []byte("C")

	q = bytes.ToUpper(q)
	t = bytes.ToUpper(t)
	cigar, err := algn.Align(q, t)
	if err != nil {
		_t.Error(err)
		return
	}

	// fmt.Printf("\n---------------- M ----------------\n")
	// PrintComponent(os.Stdout, algn.M, "M")
	algn.Plot(&q, &t, os.Stdout, algn.M, true, false)

	// fmt.Printf("\n---------------- I ----------------\n")
	// PrintComponent(os.Stdout, algn.I, "I")
	// algn.Plot(&q, &t, os.Stdout, algn.I, false, false)

	// fmt.Printf("\n---------------- D ----------------\n")
	// PrintComponent(os.Stdout, algn.D, "D")
	// algn.Plot(&q, &t, os.Stdout, algn.D, false, false)

	// fmt.Printf("\n---------------- S ----------------\n")
	// PrintComponent(os.Stdout, algn.S, "S")
	// algn.Plot(&q, &t, os.Stdout, algn.S, false, false)

	if cigar != nil {
		fmt.Println()
		fmt.Printf("CIGAR:  %s\n", cigar.CIGAR())
		Q, A, T := cigar.Alignment(&q, &t)
		fmt.Printf("query   %s\n", *Q)
		fmt.Printf("        %s\n", *A)
		fmt.Printf("target  %s\n", *T)

		fmt.Println()
		fmt.Printf("align-score : %d\n", cigar.Score)
		fmt.Printf("align-region: q[%d, %d] vs t[%d, %d]\n",
			cigar.QBegin+1, cigar.QEnd+1, cigar.TBegin+1, cigar.TEnd+1)
		fmt.Printf("align-length: %d, matches: %d (%.2f%%), gaps: %d, gapRegions: %d\n",
			cigar.AlignLen, cigar.Matches, float64(cigar.Matches)/float64(cigar.AlignLen)*100,
			cigar.Gaps, cigar.GapRegions)
		fmt.Println()

		RecycleAlignment(Q, A, T)
		RecycleCIGAR(cigar)
	}
	RecycleAligner(algn)
}
