// Copyright © 2024-2026 Wei Shen <shenwei356@gmail.com>
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
	"errors"
	"flag"
	"fmt"
	"os"
	"path/filepath"

	"github.com/pkg/profile"
	"github.com/shenwei356/wfa"
)

var version = "0.5.3"

func main() {
	checkError(run())
}

// run returns errors after flushing output and recycling alignment resources.
func run() (err error) {
	app := filepath.Base(os.Args[0])
	usage := fmt.Sprintf(`
WFA alignment in Golang

 Author: Wei Shen <shenwei356@gmail.com>
   Code: https://github.com/shenwei356/wfa
Version: v%s

Input file format:
  Each pair must contain a nonempty >query line followed by a nonempty <target line.
  Blank lines, invalid prefixes, and incomplete pairs are errors.
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
	noGlobal := flag.Bool("g", false, "use semi-global alignment with free ends")
	noAdaptive := flag.Bool("a", false, "do not use adaptive reduction")
	noOutput := flag.Bool("N", false, "do not output alignment (for benchmark)")
	trim := flag.Bool("t", false, "only show the region from the first to last exact match (empty if none)")

	pprofCPU := flag.Bool("p", false, "cpu pprof. go tool pprof -http=:8080 cpu.pprof")
	pprofMem := flag.Bool("m", false, "mem pprof. go tool pprof -http=:8080 mem.pprof")

	flag.Parse()

	if *help {
		flag.Usage()
		return nil
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
		algn.AdaptiveReduction(&wfa.AdaptiveReductionOption{
			MinWFLen:    10,
			MaxDistDiff: 50,
			CutoffStep:  1,
		})
	}

	defer func() {
		wfa.RecycleAligner(algn)
		err = errors.Join(err, outfh.Flush())
	}()

	falign2Seq := func(q, t string) error {

		_q, _t := []byte(q), []byte(t)
		result, err := algn.Align(_q, _t)
		if err != nil {
			return err
		}
		defer wfa.RecycleAlignmentResult(result)

		if !*noOutput {
			Q, A, T := result.AlignmentText(&_q, &_t, *trim)
			defer wfa.RecycleAlignmentText(Q, A, T)
			identity := 0.0
			if result.AlignLen > 0 {
				identity = float64(result.Matches) / float64(result.AlignLen) * 100
			}

			// fmt.Fprintln(outfh, q, t)
			fmt.Fprintf(outfh, "query   %s\n", *Q)
			fmt.Fprintf(outfh, "        %s\n", *A)
			fmt.Fprintf(outfh, "target  %s\n", *T)
			fmt.Fprintf(outfh, "cigar   %s\n", result.CIGAR(*trim))
			fmt.Fprintln(outfh)
			fmt.Fprintf(outfh, "align-score : %d\n", result.Score)
			fmt.Fprintf(outfh, "match-region: q[%d, %d]/%d vs t[%d, %d]/%d\n",
				result.QBegin, result.QEnd, len(_q), result.TBegin, result.TEnd, len(_t))
			fmt.Fprintf(outfh, "align-length: %d, matches: %d (%.2f%%), gaps: %d, gap regions: %d\n",
				result.AlignLen, result.Matches, identity,
				result.Gaps, result.GapRegions)
			// bufio.Writer retains earlier write errors; check the final write too.
			_, err := fmt.Fprintln(outfh)
			return err
		}
		return nil
	}

	var q, t string

	// two sequences from positional arguments

	if *infile == "" {
		if flag.NArg() != 2 {
			return fmt.Errorf("if flag -i not given, please give me two sequences. type \"%s -h\" for help.", app)
		}
		q = flag.Arg(0)
		t = flag.Arg(1)

		return falign2Seq(q, t)
	}

	// sequence pairs from a file

	fh, err := os.Open(*infile)
	if err != nil {
		return fmt.Errorf("failed to read file %s: %w", *infile, err)
	}
	defer func() { err = errors.Join(err, fh.Close()) }()

	scanner := bufio.NewScanner(fh)
	// Permit the library's supported sequence length, plus prefix and CRLF.
	scanner.Buffer(make([]byte, 4096), wfa.MaxSeqLen+3)
	line := 0
	for scanner.Scan() {
		line++
		q = scanner.Text()
		if len(q) < 2 || q[0] != '>' {
			return fmt.Errorf("%s:%d: expected '>' followed by a nonempty query sequence", *infile, line)
		}
		if !scanner.Scan() {
			if err := scanner.Err(); err != nil {
				return fmt.Errorf("read %s:%d: %w", *infile, line+1, err)
			}
			return fmt.Errorf("%s:%d: missing target sequence for query on line %d", *infile, line+1, line)
		}
		line++
		t = scanner.Text()
		if len(t) < 2 || t[0] != '<' {
			return fmt.Errorf("%s:%d: expected '<' followed by a nonempty target sequence", *infile, line)
		}
		if err := falign2Seq(q[1:], t[1:]); err != nil {
			return fmt.Errorf("%s:%d: %w", *infile, line, err)
		}
	}
	if err := scanner.Err(); err != nil {
		return fmt.Errorf("read %s:%d: %w", *infile, line+1, err)
	}
	return nil
}

func checkError(err error) {
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
}
