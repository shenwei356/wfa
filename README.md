# Wavefront alignment algorithm (WFA) in Golang (WIP)

[![Go Reference](https://pkg.go.dev/badge/github.com/shenwei356/wfa.svg)](https://pkg.go.dev/github.com/shenwei356/wfa)

This golang packages implements Wavefront alignment algorithm (WFA), not BiWFA (maybe in the future).

## Table of Contents

+ [Details](#details)
+ [Examples](#examples)
+ [Usages](#usages)
+ [CLI](#cli)
+ [Benchmark](#benchmark)
+ [Reference](#reference)

## Details

- A WFA component is saved with a 2-dimension (`s`, `k`) slice `[]*[]uint32`.
    - Score `s` is the index of the outer slice `[]*[]uint32`. 
      For non-existent scores, the pointer `*[]uint32` is `nil`.
    - Diagonal `k` is the index of the inner slice `[]uint32`.
      To support negative `k` values, we use this layout:

            index: 0,  1,  2,  3,  4,  5,  6
            k:     0, -1,  1, -2,  2, -3,  3

      Value `0` means the `k` does not exist.

    - Offsets are saved with `uint32` integers, with the lower 3 bits for
      saving 5 possible paths which are used for backtrace.

            wfaInsertOpen uint32 = iota + 1
            wfaInsertExt
            wfaDeleteOpen
            wfaDeleteExt
            wfaMismatch
            wfaMatch // only for backtrace, not saved in the component

- Maximum sequence length: `536,870,911` (1<<(32-3) - 1)
- All objects are saved in object pool for computation efficiency.
  Just don't forget to recycle them.

## Examples

Each WFA component can be visualized as a table.
A table cell contains the alignment type symbol and the score.

    ⊕    Unknown,
    ⟼    Gap open (Insertion)
    🠦    Gap extension (Insertion)
    ↧    Gap open (Deletion)
    🠧    Gap extension (Deletion)
    ⬂    Mismatch
    ⬊    Match

Global alignment

|   |   |  1   |  2  |  3   |  4   |  5  |  6  |  7  |  8  |  9  | 10  |
|:-:|:-:|:----:|:---:|:----:|:----:|:---:|:---:|:---:|:---:|:---:|:---:|
|   |   |  A   |  G  |  G   |  A   |  T  |  G  |  C  |  T  |  C  |  G  |
|1  |A  |⬊ 0   |⟼ 8  |🠦10   |🠦12   |  .  |  .  |  .  |  .  |  .  |  .  |
|2  |C  |↧ 8   |⬂ 4  |⬂12   |  .   |  .  |  .  |  .  |  .  |  .  |  .  |
|3  |C  |🠧10   |⬂12  |⬂ 8   |  .   |  .  |  .  |  .  |  .  |  .  |  .  |
|4  |A  |🠧12   |  .  |  .   |⬊ 8   |  .  |  .  |  .  |  .  |  .  |  .  |
|5  |T  |  .   |  .  |  .   |  .   |⬊ 8  |  .  |  .  |  .  |  .  |  .  |
|6  |A  |  .   |  .  |  .   |  .   |  .  |⬂12  |  .  |  .  |  .  |  .  |
|7  |C  |  .   |  .  |  .   |  .   |  .  |  .  |⬊12  |  .  |  .  |  .  |
|8  |T  |  .   |  .  |  .   |  .   |  .  |  .  |  .  |⬊12  |  .  |  .  |
|9  |C  |  .   |  .  |  .   |  .   |  .  |  .  |  .  |  .  |⬊12  |  .  |
|10 |G  |  .   |  .  |  .   |  .   |  .  |  .  |  .  |  .  |  .  |⬊12  |

```
CIGAR:  1M2X2M1X4M
query   ACCATACTCG
        |  || ||||
target  AGGATGCTCG

align-score : 12
align-region: q[1, 10] vs t[1, 10]
align-length: 10, matches: 7 (70.00%), gaps: 0, gapRegions: 0
```

Semi-global alignment

|   |   |  1  |  2  |  3   |  4   |  5   |  6   |  7  |  8   |  9   | 10  | 11   | 12  |
|:-:|:-:|:---:|:---:|:----:|:----:|:----:|:----:|:---:|:----:|:----:|:---:|:----:|:---:|
|   |   |  C  |  A  |  G   |  G   |  C   |  T   |  C  |  C   |  T   |  C  |  G   |  G  |
|1  |A  |⬂ 4  |⬊ 0  |⬂ 4   |⬂ 4   |⬂ 4   |⬂ 4   |⬂ 4  |⬂ 4   |⬂ 4   |⬂ 4  |⬂ 4   |⬂ 4  |
|2  |C  |⬊ 0  |⬂ 8  |⬂ 4   |⬂ 8   |⬊ 4   |⬂ 8   |⬊ 4  |⬊ 4   |⬂ 8   |⬊ 4  |⬂ 8   |⬂ 8  |
|3  |G  |⬂ 4  |⬂ 4  |⬊ 8   |⬊ 4   |⬂12   |⬂ 8   |⬂12  |⬂ 8   |⬂ 8   |  .  |⬊ 4   |⬊ 8  |
|4  |A  |⬂ 4  |⬊ 4  |⬂ 8   |⬂12   |⬂ 8   |🠧14   |⬂12  |🠧14   |⬂12   |⬂12  |↧12   |⬂ 8  |
|5  |T  |⬂ 4  |⬂ 8  |⬂ 8   |⬂12   |  .   |⬊ 8   |  .  |⬂16   |⬊14   |  .  |🠧14   |⬂16  |
|6  |C  |⬊ 0  |⬂ 8  |🠦10   |⬂12   |⬊12   |🠦18   |⬊ 8  |⟼16   |🠦18   |⬊14  |🠧16   |⬂18  |
|7  |T  |⬂ 4  |⬂ 4  |  .   |  .   |  .   |⬊12   |↧16  |⬂12   |⬊16   |  .  |⬂18   |⬂20  |
|8  |C  |⬊ 0  |⬂ 8  |⬂ 8   |🠦12   |🠦14   |🠦16   |⬊12  |⬊16   |⬂16   |⬊16  |🠧20   |  .  |
|9  |G  |⬂ 4  |⬂ 4  |⬊ 8   |⬊ 8   |⬂16   |⬂18   |⬂20  |⬂16   |⬂20   |⬂20  |⬊16   |⬊20  |
```
CIGAR:  1I1M1X1M1X1M1I4M1I
query   -ACGAT-CTCG-
         | | | ||||
target  CAGGCTCCTCGG

align-score : 16
align-region: q[1, 9] vs t[2, 11]
align-length: 10, matches: 7 (70.00%), gaps: 1, gapRegions: 1
```

## Usages

```
import "github.com/shenwei356/wfa"

// aligner
algn := wfa.New(
    &wfa.Penalties{
        Mismatch: 4,
        GapOpen:  6,
        GapExt:   2,
    },
    &wfa.Options{
        GlobalAlignment: false,
    })

// set adaptive reduction parameters
algn.AdaptiveReduction(&wfa.AdaptiveReductionOption{
    MinWFLen:    10,
    MaxDistDiff: 50,
    CutoffStep:  1,
})

q := []byte("ACCATACTCG")
t := []byte("AGGATGCTCG")

// align
cigar, err := algn.Align(&q, &t)
checkErr(err)

// score table of M
algn.Plot(&q, &t, os.Stdout, algn.M, true)

if outputAlignment {
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

    // !! important, recycle objects
    wfa.RecycleAlignment(Q, A, T)
}

// !! important, recycle objects
wfa.RecycleCIGAR(cigar)
wfa.RecycleAligner(algn)
```

## CLI

A [CLI](https://github.com/shenwei356/wfa/blob/main/benchmark/wfa-go.go) is available to
align two sequences from either positional arguments or an input file
([format](https://github.com/smarco/WFA-paper?tab=readme-ov-file#41-introduction-to-benchmarking-wfa-simple-tests)).

<details>
<summary>Usage</summary>

```
WFA alignment in Golang

 Author: Wei Shen <shenwei356@gmail.com>
   Code: https://github.com/shenwei356/wfa
Version: v0.1.0

Input file format:
  see https://github.com/smarco/WFA-paper?tab=readme-ov-file#41-introduction-to-benchmarking-wfa-simple-tests
  Example:
  >ATTGGAAAATAGGATTGGGGTTTGTTTATATTTGGGTTGAGGGATGTCCCACCTTCGTCGTCCTTACGTTTCCGGAAGGGAGTGGTTAGCTCGAAGCCCA
  <GATTGGAAAATAGGATGGGGTTTGTTTATATTTGGGTTGAGGGATGTCCCACCTTGTCGTCCTTACGTTTCCGGAAGGGAGTGGTTGCTCGAAGCCCA
  >CCGTAGAGTTAGACACTCGACCGTGGTGAATCCGCGACCACCGCTTTGACGGGCGCTCTACGGTATCCCGCGATTTGTGTACGTGAAGCAGTGATTAAAC
  <CCTAGAGTTAGACACTCGACCGTGGTGAATCCGCGATCTACCGCTTTGACGGGCGCTCTACGGTATCCCGCGATTTGTGTACGTGAAGCGAGTGATTAAAC

Usage:
  1. Align two sequences from the positional arguments.

        wfa-go [options] <query seq> <target seq>

  2. Align sequence pairs from the input file (described above).

        wfa-go [options] -i input.txt

Options/Flags:
  -N    do not output alignment (for benchmark)
  -a    do not use adaptive reduction
  -g    do not use global alignment
  -h    print help message
  -i string
        input file.
  -m    mem pprof. go tool pprof -http=:8080 mem.pprof
  -p    cpu pprof. go tool pprof -http=:8080 cpu.pprof
```
</details>

## Benchmark

Generate datasets with WFA2-lib (v2.3.5):

    ./bin/generate_dataset -n 100000 -l 1000 -e 0.05 -o l1000-e0.05.seq
    ./bin/generate_dataset -n 100000 -l 1000 -e 0.10 -o l1000-e0.10.seq
    ./bin/generate_dataset -n 100000 -l 1000 -e 0.20 -o l1000-e0.20.seq

Commands:

    # memusg: https://github.com/shenwei356/memusg

    # WFA2-lib
    memusg -t -s "./bin/align_benchmark -i l1000-e0.05.seq  -a gap-affine-wfa"

    # WFA-go (this package)
    # global alignment && do not output results
    memusg -t -s "go run ./benchmark/wfa-go.go -N -i /home/shenwei/Downloads/WFA2-lib/l1000-e0.05.seq"

Results:

|Seq-len|Seq-num|Error-rate|Package |Time   |Memory   |Time-ratio|
|:-----:|:-----:|:--------:|:-------|:-----:|--------:|---------:|
|1000   |100000 |0.05      |WFA2-lib|5.762s |4.79 MB  |1x        |
|       |       |          |WFA-go  |1m:19s |51.3 MB  |14x       |
|1000   |100000 |0.10      |WFA2-lib|14.762s|8.04 MB  |1x        |
|       |       |          |WFA-go  |4m:22s |112.73 MB|18x       |
|1000   |100000 |0.20      |WFA2-lib|47.714s|9.27 MB  |1x        |
|       |       |          |WFA-go  |12m:26s|107.97 MB|16x       |


## Reference

- **Santiago Marco-Sola, Juan Carlos Moure, Miquel Moreto, Antonio Espinosa**. ["Fast gap-affine pairwise alignment using the wavefront algorithm."](https://doi.org/10.1093/bioinformatics/btaa777) Bioinformatics, 2020.
- https://github.com/smarco/WFA-paper/
- https://github.com/smarco/WFA2-lib

## Support

Please [open an issue](https://github.com/shenwei356/wfa/issues) to report bugs,
propose new functions or ask for help.

## License

[MIT License](https://github.com/shenwei356/wfa/blob/master/LICENSE)
