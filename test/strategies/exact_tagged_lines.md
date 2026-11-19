# Full-Range Tagged Line Orientation

The `tagged_orientation.hpp` kernel computes the Cartesian sign
of three points defined by pairs of binary64 supporting lines. Each point
retains a projective intersection `(X,Y,W)`, without division or GCD. It
accepts all finite binary64 endpoint values, including signed zero,
subnormals and the largest finite numbers. It also accepts unique line
intersections outside the source segments or the finite binary64 range.
Nonfinite endpoints and zero weights (parallel, coincident, or degenerate
supports) are rejected. This is not a general segment-intersection API.

## Representation and Proof

Use radix `R = 2^21` and a compile-time-sized array of binary64 digits.
The signed integer digit at tag `k` denotes `digit * R^k`. Two indices
describe the active span; tags inside the span are implicit. There is no
floating computation of `R^k`. Inputs are decoded directly from their
binary64 bit patterns into at most four digits. Only the active spans
are copied or initialized. All successful kernel paths use fixed-size
local arrays and no wide integer, GCD, FMA residual, or heap allocation.

An accumulator adds signed digit convolutions in binary64, then sweeps
from least to greatest tag. For a raw coefficient plus incoming carry `v`,
compute `q = trunc(v/R)` and `r = v - q*R`. Replace the current block by
`r` and add `q` to the next. Power-of-two scaling is exact, and the signed
remainder has `|r| < R`. No transform drops a tail. Digits need not all have
the same sign: the absolute value of the sum below the highest nonzero tag
`h` is at most `R^h - R^first`, strictly less than `R^h`. The highest digit
therefore gives the exact sign, and an empty span is exact zero.

### Exponent and Coefficient Bounds

For a stage of polynomial degree `d` and absolute monomial count `M`,
the weighted norm `sum_k |digit[k]| R^k` is strictly less than
`M * 2^(1024*d)`. For an input the norm equals its magnitude, and convolution
is submultiplicative while addition is subadditive. Carrying cannot increase
this norm: `r` and `q` have the sign of `v`, so their weighted magnitudes add
to the old contribution, and adding the carry to its neighbor satisfies
the triangle inequality. Thus the bound holds during cancellation and for
every partial sum, not just for the final value. No nonzero block can reach
tag `ceil((1024*d + ceil(log2(M)))/21)` or higher. The lowest raw tag is
`d * floor(-1074/21) = -52*d`; convolution adds tags and carry only raises
them. The table gives the complete inclusive/exclusive allocation ranges.

| Stage | d | M | Tags `[first,end)` | Digits | Max raw terms per digit |
| --- | ---: | ---: | --- | ---: | ---: |
| Binary64 input | 1 | 1 | `[-52,49)` | 101 | decode |
| Line A, B | 1 | 2 | `[-52,49)` | 101 | 2 |
| Line C | 2 | 2 | `[-104,98)` | 202 | 202 |
| Intersection X, Y | 3 | 8 | `[-156,147)` | 303 | 202 |
| Intersection W | 2 | 8 | `[-104,98)` | 202 | 202 |
| XY minors of two records | 5 | 128 | `[-260,245)` | 505 | 404 |
| W minor of two records | 6 | 128 | `[-312,293)` | 605 | 606 |
| Determinant of three records | 8 | 3072 | `[-416,391)` | 807 | 808 |

For a product, at most the smaller operand capacity contributes to one
coefficient; sum those bounds for all products. Treat a linear term as
one product with 1. For fan-in `F`, every raw coefficient and partial
accumulation is bounded by `F(R-1)^2`. Starting with zero, carry magnitude
is at most `F(R-1)` by induction, because the sum before each split is
bounded by `F R(R-1)`. With `F <= 808`, this is below `2^52`. Thus all the
integer arithmetic on binary64 digits is exact (even in directed rounding
modes); `trunc` explicitly supplies the desired carry. `sum()` checks
degree, monomial budget and the binary64 integer limit at compile time.

### Geometry

A line through two points has coefficients `(y0-y1, x1-x0, x0*y1-x1*y0)`
of degrees `(1,1,2)` and two monomials each. The cross product of two
lines gives `(X,Y,W)` of degrees `(3,3,2)` and at most eight monomials
per component. The cross product of the first two homogeneous rows gives
minors of degrees `(5,5,6)` and at most 128 monomials per component.
Dotting with the third row therefore uses degree eight and at most 3072
monomials. Its sign times the three weight signs equals the Cartesian
orientation sign; a line reversal or support-order swap cannot change it.

An intersection holds 808 digit slots and three spans. The determinant
phase uses 2,422 digit slots for the three minors and final accumulation
when return copies are elided, plus their span metadata; this excludes
the three input records. This is stack scratch, not a cache layout
proposal. Static capacity is not a claim of small per-record footprint.

## Verification and Benchmark

These are kernel measurements. The separate [opt-in dispatch trial](exact_tagged_dispatch.md)
provides traffic and end-to-end benchmark commands.

A rational parameter-intersection oracle followed by a Cartesian
determinant checks 466 input triples. Readback checks X/W and Y/W
separately on 40 triples. Hand cases include exact zero, negative zero,
all weight-sign combinations, nearest-even ties, collinearity across
exponent extremes, half-subnormal intersections, maximal inputs and
cancellation of maximal products leaving a small nonzero sign. Directed
rounding and invalid inputs have separate checks.

```sh
CASE_FILTER=tagged_lines BOOST_ROOT=/path/to/boost \
  sh test/strategies/run_exact_primitives.sh /tmp/bg-tagged-lines
clang++ -std=c++17 -O2 -DNDEBUG -ffp-contract=off -Iinclude \
  -I/path/to/boost test/strategies/exact_tagged_lines.cpp -o /tmp/bg-tagged-lines-bench
/tmp/bg-tagged-lines-bench --benchmark
```

The fixed baseline is the current bounded kernel: construct normalized
homogeneous records with `direct_native_homogeneous`, then use compact
128-bit orientation if all three records fit, otherwise the wide
`fixed_orientation`. Cold timing includes all three constructions
(including normalization); prepared timing reuses records from both
implementations and measures the sign only. No earlier filters, cache
lookup, cache bound refresh, record allocation, or geometry traversal are
timed. Random full-exponent line pairs are a stress workload, not a
measure of fallback frequency in ordinary overlay.

### Results

Apple Clang 21 and Homebrew GCC 16.2, arm64; C++17,
`-O2 -DNDEBUG -ffp-contract=off`. Medians of five runs, ns per sign;
each cell is fixed / tagged. All checksums agree within every path pair.

| Corpus | Clang prepared | Clang cold | GCC prepared | GCC cold |
| --- | ---: | ---: | ---: | ---: |
| Integer grid | 107 / 34 | 3350 / 268 | 109 / 27 | 3630 / 204 |
| Near-collinear | 552 / 103 | 7095 / 647 | 532 / 87 | 7118 / 598 |
| Sparse extremes | 1426 / 635 | 7644 / 1096 | 1228 / 770 | 6777 / 1196 |
| Full exponent range | 90851 / 10491 | 510935 / 15748 | 68408 / 24203 | 424760 / 30996 |

An intersection record is 6,488 bytes tagged vs. 2,320 fixed on both
builds. Apple Clang also passed the registered test with assertions and a
separate ASan/UBSan build. GCC passed the exact checks in the timed build.

### Integration Recommendation

The same-state kernel timings justify an opt-in trial on real three-record
fallback traffic, starting at cold misses with original line supports. They
do not justify replacing the default or storing full-width tagged arrays
in every cache entry. Reconstructing tagged records is slower than reusing
prepared fixed records on the integer-grid and near-collinear workloads.
Measure cold/hot frequency, active-span storage and peak cache footprint,
and paired end-to-end overlay and relate before any default dispatch change.
Do not replace compact or cached fixed values blindly. Mixed native and
constructed operands, recursive constructions, comparisons and coordinate
emission are outside this orientation prototype.
