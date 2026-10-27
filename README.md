# Wavefront alignment algorithm (WFA) in Golang

[![Go Reference](https://pkg.go.dev/badge/github.com/shenwei356/wfa.svg)](https://pkg.go.dev/github.com/shenwei356/wfa)

This golang package implements [Wavefront alignment algorithm (WFA)](https://doi.org/10.1093/bioinformatics/btaa777),
not [BiWFA](https://doi.org/10.1093/bioinformatics/btad074) (maybe in the future).

Implemented features:

- Distance metrics: gap-affine.
- Alignment types: global (end-to-end), semi-global (free ends).
- Heuristics: wf-adaptive. 

Executable binaries for common operating systems are available to [download](https://github.com/shenwei356/wfa/releases),
for benchmarking with lots of sequences, or for fast aligning two sequences in the command line.

```
$ wfa-go -g "Bioinformatics helps Biology" "We learn bioinformatics to help biologists"
query   ---------Bioinformatics ---helps Biology---
                  ||||||||||||||   |||| | |||||
target  We learn bioinformatics to help- biologists
cigar   9I1X14M3I4M1D1M1X5M1X3I

align-score : 32
match-region: q[2, 27]/28 vs t[11, 38]/42
align-length: 29, matches: 24 (82.76%), gaps: 4, gap regions: 2
```

## Table of Contents

+ [History](#history)
+ [Details](#details)
+ [Examples](#examples)
+ [Usages](#usages)
+ [CLI](#cli)
+ [Benchmark](#benchmark)
+ [Reference](#reference)

## History

1. I need a fast DNA alignment package in Golang for my project.
1. WFA "looks easy" to implement as it does not heavily reply on SIMD intrinsics,
   though there are some other algorithms performing well in [a benchmark](https://github.com/rchikhi/rust-alignbench).
1. I'm not familiar with C++, and I found it difficult to understand the official code.
   I only found [one](https://github.com/cschin/wavefront-aln) (in Rust) 3rd party implementation.
1. After reading the WFA paper, I thought the algorithm is easy, so I implemented WFA from the scratch.
1. Later I found it not easy, there were so many details.
   After reading the `next` step in the rust version, I realised it's because I don't know gap-affine penalties.
1. The `backtrace` step is the most difficult part. In v0.1.0, I used 3 extra components to store the source offsets
   for each cell in I, D, M components. But it needs more memory. And the speed is also not ideal, 1/20 of official version.
    - Besides, I checked the bases again when tracing back matches. WFA did this too, but WFA2 did not.
1. Next, aftering reading thofficial implementation, I rewrote the whole project, using similar backtrace workfow with WFA2.
   The speed increased to 1/10 of the official version.
1. C++ is wild, it even support accessing list/array elements with negative indexes ([the diagonal k](https://github.com/smarco/WFA2-lib/issues/94)).

## Details

- A WFA component is saved as a list of WaveFront `[]*WaveFront`.
    - Score `s` is the index.
    - If the value is `nil`, it means the score does not exsit.
- A WaveFront is saved as a list of offsets `[]uint32`.
  We preset the list with a big length (2048) to avoid frequent `append` operations.
    - Diagonal `k` is the index.
      To support negative `k` values, we use [ZigZag encoding](https://gist.github.com/mfuerstenau/ba870a29e16536fdbaba):

          index: 0,  1,  2,  3,  4,  5,  6
          k:     0, -1,  1, -2,  2, -3,  3

      Value `0` means the `k` does not exist.

    - Offsets are saved with `uint32` integers, with the lower 3 bits for
      saving 6 possible paths which are used for backtrace.

          wfaInsertOpen uint32 = iota + 1
          wfaInsertExt
          wfaDeleteOpen
          wfaDeleteExt
          wfaMismatch
          wfaMatch

- Maximum sequence length: 512 Mb, `536,870,911` (1<<(32-3) - 1).
- All objects are saved in object pools for computation efficiency.
  Just don't forget to recycle them.

more details

- `Mismatch` and `GapExt` must be positive. `GapOpen` may be zero, and
  `GapOpen + GapExt` must fit in `uint32`. `Align` and `AlignPointers` return
  `ErrInvalidPenalties` for unsupported penalties.
- `AlignmentResult.Score` is the gap-affine cost of the returned CIGAR.
  Adaptive reduction is heuristic and may return a nonoptimal alignment.
- Semi-global alignment retains its established endpoint constraint: the shorter
  sequence is consumed, and the consumed offset in the other sequence must be
  at least that length. This is a constraint on endpoint positions, not a minimum
  length of the matched region. Only actual free flanks are excluded from the score;
  terminal gaps needed to reach an eligible endpoint remain charged.
- Coordinates and alignment statistics describe the region from the first to
  last exact match. With no exact matches, they are zero; `CIGAR(true)` and
  `AlignmentText(..., true)` return empty output. Untrimmed output remains available.


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

Flanking bases are excluded from the reported aligned region while the endpoint
constraint described above is preserved.

```
$ wfa-go -g ACGATCTCG CAGGCTCCTCGG
query   -ACGATC-TCG-
         | | || |||
target  CAGGCTCCTCGG
cigar   1I1M1X1M1X2M1I3M1I

align-score : 16
match-region: q[1, 9]/9 vs t[2, 11]/12
align-length: 10, matches: 7 (70.00%), gaps: 1, gap regions: 1
```

## Usages

```
import "github.com/shenwei356/wfa"

// ------------------[ initialization ]------------------

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

// ------------------[ align one pair of seqs ]------------------

q := []byte("ACCATACTCG")
t := []byte("AGGATGCTCG")

// align
result, err := algn.Align(q, t)
checkErr(err)

// score table of M
algn.Plot(&q, &t, os.Stdout, algn.M, true, -1)

if outputAlignment {
    fmt.Println()
    fmt.Printf("CIGAR:  %s\n", result.CIGAR(false))

    Q, A, T := result.AlignmentText(&q, &t, false)
    fmt.Printf("query   %s\n", *Q)
    fmt.Printf("        %s\n", *A)
    fmt.Printf("target  %s\n", *T)

    fmt.Println()
    fmt.Printf("align-score : %d\n", result.Score)
    fmt.Printf("align-region: q[%d, %d] vs t[%d, %d]\n",
        result.QBegin, result.QEnd, result.TBegin, result.TEnd)
    fmt.Printf("align-length: %d, matches: %d (%.2f%%), gaps: %d, gapRegions: %d\n",
        result.AlignLen, result.Matches, float64(result.Matches)/float64(result.AlignLen)*100,
        result.Gaps, result.GapRegions)
    fmt.Println()

    
    wfa.RecycleAlignmentText(Q, A, T) // !! important, recycle objects
}

wfa.RecycleAlignmentResult(result) // !! important, recycle objects

// ------------------[ clean ]------------------

wfa.RecycleAligner(algn) // !! important, recycle objects

```

## CLI

A [CLI](https://github.com/shenwei356/wfa/blob/main/wfa-go/wfa-go.go)
(download [binaries](https://github.com/shenwei356/wfa/releases)) is available to
align two sequences from either positional arguments or an input file
([format](https://github.com/smarco/WFA-paper?tab=readme-ov-file#41-introduction-to-benchmarking-wfa-simple-tests)).
Each pair requires a nonempty `>` query line followed by a nonempty `<` target
line. Blank lines, invalid prefixes, and incomplete pairs are errors. Read and
write failures, including a failed output flush, produce a nonzero exit status.

<details>
<summary>Example</summary>

Fast alignment.

```
$ wfa-go AGCTAGTGTCAATGGCTACTTTTCAGGTCCT AACTAAGTGTCGGTGGCTACTATATATCAGGTCCT
query   AGCTA-GTGTCAATGGCTACT---TTTCAGGTCCT
        | ||| |||||  ||||||||   | |||||||||
target  AACTAAGTGTCGGTGGCTACTATATATCAGGTCCT
cigar   1M1X3M1I5M2X8M3I1M1X9M

align-score : 36
match-region: q[1, 31]/31 vs t[1, 35]/35
align-length: 35, matches: 27 (77.14%), gaps: 4, gap regions: 2
```

From a input file, for benchmark.

```
$ wfa-go -i seqs.txt
query   -ATTGGAAAATAGGATTGGGGTTTGTTTATATTTGGGTTGAGGGATGTCCCACCTTCGTCGTCCTTACGTTTCCGGAAGGGAGTGGTTAGCTCGAAGCCCA
         ||||||||||||||| ||||||||||||||||||||||||||||||||||||||| ||||||||||||||||||||||||||||||| ||||||||||||
target  GATTGGAAAATAGGAT-GGGGTTTGTTTATATTTGGGTTGAGGGATGTCCCACCTT-GTCGTCCTTACGTTTCCGGAAGGGAGTGGTT-GCTCGAAGCCCA
cigar   1I15M1D39M1D31M1D12M

align-score : 32
match-region: q[1, 100]/100 vs t[2, 98]/98
align-length: 100, matches: 97 (97.00%), gaps: 3, gap regions: 3

query   CCGTAGAGTTAGACACTCGACCGTGGTGAATCCGCGA-CCACCGCTTTGACGGGCGCTCTACGGTATCCCGCGATTTGTGTACGTGAAGC-AGTGATTAAAC
        || |||||||||||||||||||||||||||||||||| | |||||||||||||||||||||||||||||||||||||||||||||||||| |||||||||||
target  CC-TAGAGTTAGACACTCGACCGTGGTGAATCCGCGATCTACCGCTTTGACGGGCGCTCTACGGTATCCCGCGATTTGTGTACGTGAAGCGAGTGATTAAAC
cigar   2M1D34M1I1M1X50M1I11M

align-score : 28
match-region: q[1, 100]/100 vs t[1, 101]/101
align-length: 102, matches: 98 (96.08%), gaps: 3, gap regions: 3
```
</details>

<details>
<summary>Usage</summary>

```
$ wfa-go -h
WFA alignment in Golang

 Author: Wei Shen <shenwei356@gmail.com>
   Code: https://github.com/shenwei356/wfa
Version: v0.5.1

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

        wfa-go [options] <query seq> <target seq>

  2. Align sequence pairs from the input file (described above).

        wfa-go [options] -i input.txt

Options/Flags:
  -N    do not output alignment (for benchmark)
  -a    do not use adaptive reduction
  -g    use semi-global alignment with free ends
  -h    print help message
  -i string
        input file.
  -m    mem pprof. go tool pprof -http=:8080 mem.pprof
  -p    cpu pprof. go tool pprof -http=:8080 cpu.pprof
  -t    only show the region from the first to last exact match (empty if none)
```
</details>

## Benchmark

Generate datasets with WFA (e634175) or WFA2-lib (v2.3.5):

    ./bin/generate_dataset -n 100000 -l 1000 -e 0.05 -o l1000-e0.05.seq
    ./bin/generate_dataset -n 100000 -l 1000 -e 0.10 -o l1000-e0.10.seq
    ./bin/generate_dataset -n 100000 -l 1000 -e 0.20 -o l1000-e0.20.seq

    ./bin/generate_dataset -n 500 -l 50000 -e 0.05 -o l50000-e0.05.seq
    ./bin/generate_dataset -n 500 -l 50000 -e 0.10 -o l50000-e0.10.seq
    ./bin/generate_dataset -n 500 -l 50000 -e 0.20 -o l50000-e0.20.seq

Commands (all tools use wfa-adaptive heuristic 10,50,1):

    # memusg: https://github.com/shenwei356/memusg

    # WFA
    memusg -t -s "./bin/align_benchmark -i l1000-e0.05.seq -a gap-affine-wfa-adaptive"

    # WFA2-lib
    memusg -t -s "./bin/align_benchmark -i l1000-e0.05.seq -a gap-affine-wfa --wfa-heuristic wfa-adaptive --wfa-heuristic-parameters 10,50,1"

    # WFA-go (this package, binary files are availabe in the release page)
    # global alignment && do not output results
    memusg -t -s "wfa-go -N -i l1000-e0.05.seq"

    csvtk csv2md -t benchmark.tsv -a c,c,c,l,r,r

Results (wfa-go v0.5.0):

|Seq-len|Seq-num|Error-rate|Package|Time   |Memory   |
|:-----:|:-----:|:--------:|:------|------:|--------:|
|1000   |100000 |0.05      |WFA1   |4.523s |3.21 MB  |
|       |       |          |WFA2   |3.597s |4.90 MB  |
|       |       |          |WFA-go |8.937s |13.22 MB |
|1000   |100000 |0.1       |WFA1   |8.031s |2.04 MB  |
|       |       |          |WFA2   |6.973s |6.55 MB  |
|       |       |          |WFA-go |23.561s|13.51 MB |
|1000   |100000 |0.2       |WFA1   |15.538s|4.04 MB  |
|       |       |          |WFA2   |13.450s|9.24 MB  |
|       |       |          |WFA-go |64s    |15.48 MB |
|50000  |500    |0.05      |WFA1   |2.180s |56.12 MB |
|       |       |          |WFA2   |1.481s |107.76 MB|
|       |       |          |WFA-go |3.395s |65.51 MB |
|50000  |500    |0.1       |WFA1   |4.144s |55.3 MB  |
|       |       |          |WFA2   |3.296s |190.32 MB|
|       |       |          |WFA-go |9.735s |126.41 MB|
|50000  |500    |0.2       |WFA1   |7.574s |81.15 MB |
|       |       |          |WFA2   |6.842s |314.08 MB|
|       |       |          |WFA-go |28.144s|260.57 MB|

Run in a laptop PC, with single-thread.

## Reference

- **Santiago Marco-Sola, Juan Carlos Moure, Miquel Moreto, Antonio Espinosa**. ["Fast gap-affine pairwise alignment using the wavefront algorithm."](https://doi.org/10.1093/bioinformatics/btaa777) Bioinformatics, 2020.
- **Santiago Marco-Sola, Jordan M Eizenga, Andrea Guarracino, Benedict Paten, Erik Garrison, Miquel Moreto**. ["Optimal gap-affine alignment in O(s) space"](https://doi.org/10.1093/bioinformatics/btad074). Bioinformatics, 2023.
- https://github.com/smarco/WFA-paper/
- https://github.com/smarco/WFA2-lib
- https://github.com/rchikhi/rust-alignbench

## Support

Please [open an issue](https://github.com/shenwei356/wfa/issues) to report bugs,
propose new functions or ask for help.

## License

[MIT License](https://github.com/shenwei356/wfa/blob/master/LICENSE)
