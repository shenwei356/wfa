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

package main

import (
	"bufio"
	"flag"
	"fmt"
	"os"
	"path/filepath"

	"github.com/pkg/profile"
	"github.com/shenwei356/wfa"
)

var version = "0.1.0"

func main() {
	app := filepath.Base(os.Args[0])
	usage := fmt.Sprintf(`
WFA alignment in Golang

 Author: Wei Shen <shenwei356@gmail.com>
   Code: https://github.com/shenwei356/wfa
Version: v%s

Input file format:
  see https://github.com/smarco/WFA-paper?tab=readme-ov-file#41-introduction-to-benchmarking-wfa-simple-tests
  Example:
  >ATTGGAAAATAGGATTGGGGTTTGTTTATATTTGGGTTGAGGGATGTCCCACCTTCGTCGTCCTTACGTTTCCGGAAGGGAGTGGTTAGCTCGAAGCCCA
  <GATTGGAAAATAGGATGGGGTTTGTTTATATTTGGGTTGAGGGATGTCCCACCTTGTCGTCCTTACGTTTCCGGAAGGGAGTGGTTGCTCGAAGCCCA
  >CCGTAGAGTTAGACACTCGACCGTGGTGAATCCGCGACCACCGCTTTGACGGGCGCTCTACGGTATCCCGCGATTTGTGTACGTGAAGCAGTGATTAAAC
  <CCTAGAGTTAGACACTCGACCGTGGTGAATCCGCGATCTACCGCTTTGACGGGCGCTCTACGGTATCCCGCGATTTGTGTACGTGAAGCGAGTGATTAAAC

Usage: 
  1. Align two sequences from the positional arguments.

        %s [options] <query seq> <target seq>

  2. Align sequence pairs from the input file (described above).

        %s [options] -i input.txt

Options/Flags:
`, version, app, app)

	flag.Usage = func() {
		fmt.Fprint(os.Stderr, usage)
		flag.PrintDefaults()
	}

	help := flag.Bool("h", false, "print help message")
	infile := flag.String("i", "", "input file. ")
	noGlobal := flag.Bool("g", false, "do not use global alignment")
	noAdaptive := flag.Bool("a", false, "do not use adaptive reduction")
	noOutput := flag.Bool("N", false, "do not output alignment (for benchmark)")

	pprofCPU := flag.Bool("p", false, "cpu pprof. go tool pprof -http=:8080 cpu.pprof")
	pprofMem := flag.Bool("m", false, "mem pprof. go tool pprof -http=:8080 mem.pprof")

	flag.Parse()

	if *help {
		flag.Usage()
		return
	}

	// go tool pprof -http=:8080 cpu.pprof
	if *pprofCPU {
		defer profile.Start(profile.CPUProfile, profile.ProfilePath(".")).Stop()
	} else if *pprofMem {
		defer profile.Start(profile.MemProfile, profile.ProfilePath(".")).Stop()
	}

	outfh := bufio.NewWriter(os.Stdout)

	algn := wfa.New(wfa.DefaultPenalties, &wfa.Options{
		GlobalAlignment: !*noGlobal,
	})

	if !*noAdaptive {
		algn.AdaptiveReduction(wfa.DefaultAdaptiveOption)
	}

	defer func() {
		wfa.RecycleAligner(algn)
		outfh.Flush()
	}()

	falign2Seq := func(q, t string) {

		_q, _t := []byte(q), []byte(t)
		cigar, err := algn.Align(_q, _t)
		if err != nil {
			checkError(err)
		}

		if !*noOutput {
			Q, A, T := cigar.Alignment(&_q, &_t)

			// fmt.Fprintln(outfh, q, t)
			fmt.Fprintf(outfh, "query   %s\n", *Q)
			fmt.Fprintf(outfh, "        %s\n", *A)
			fmt.Fprintf(outfh, "target  %s\n", *T)
			fmt.Fprintf(outfh, "cigar   %s\n", cigar.CIGAR())
			fmt.Fprintf(outfh, "length: %d, matches: %d (%.2f%%), gaps: %d, gap regions: %d\n",
				cigar.AlignLen, cigar.Matches, float64(cigar.Matches)/float64(cigar.AlignLen)*100,
				cigar.Gaps, cigar.GapRegions)
			fmt.Fprintln(outfh)

			wfa.RecycleAlignment(Q, A, T)
		}
		wfa.RecycleCIGAR(cigar)
	}

	var q, t string

	// two sequences from positional arguments

	if *infile == "" {
		if flag.NArg() != 2 {
			checkError(fmt.Errorf("if flag -i not given, please give me two sequences"))
		}
		q = flag.Arg(0)
		t = flag.Arg(1)

		falign2Seq(q, t)

		return
	}

	// sequence pairs from a file

	fh, err := os.Open(*infile)
	if err != nil {
		checkError(fmt.Errorf("failed to read file: %s", *infile))
	}

	scanner := bufio.NewScanner(fh)
	var flag bool
	for scanner.Scan() {
		q = scanner.Text()
		flag = scanner.Scan()
		if !flag {
			break
		}

		t = scanner.Text()

		falign2Seq(q[1:], t[1:])
	}
	if err = scanner.Err(); err != nil {
		checkError(fmt.Errorf("something wrong in reading file: %s", *infile))
	}

}

func checkError(err error) {
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
}
