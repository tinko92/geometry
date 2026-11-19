# Adaptive Exact Cartesian Filters

This follows the [pre-filter benchmark](exact_cartesian_performance.md), on the
same local Apple M5 / Apple Clang 21 machine, 2026-11-11. All work is local on
`codex/exact-cartesian-overlay-relate`; nothing was pushed.
The measured headers are committed as `e00fff01f`. They were uncommitted when
the benchmark manifest was captured; its full header hashes identify that same
content, rather than the manifest's earlier HEAD alone.

## Filtering Changes

The earlier implementation still eagerly evaluated rational bounds, measures,
collinear projections and output coordinates. This pass replaces those paths:

- Native scalars stay inline. Checked integer arithmetic and error-free binary
  residuals preserve exact native results without creating expression objects.
- Bounds, distances, measures and other scalar expressions carry certified
  intervals and retain owning references to their exact operands.
- Ambiguous constructed coordinates can refine their intervals through floating
  expansions of the original line equations, without materializing rationals.
  Comparisons then use exact polynomial expansions or certified producer
  incidence, before falling back to integer homogeneous coordinates/rationals.
- Output rounding verifies the candidate's rounding cell with exact expansion
  comparisons, including ties-to-even. It does not simply trust a quotient.
- Ring measures use an analytic sum-error bound, then interval refinement and
  finally rational evaluation. Orientation and closure traits remain respected.
- Collinear classification and segment parameters use filtered scalars, retaining
  the original endpoint nodes. Collinear spike direction uses coordinate ordering
  rather than unnecessarily evaluating a rational dot product.
- The geometry/box covered-by helper now respects the strategy-selected box type,
  avoiding eager coordinate access on implicit points. Defaults are unchanged.

Floating overflow and underflow are outside the filter contract. Roundoff,
cancellation, ties and integer overflow are not. Conservative exponent guards
remain as fast-filter admission checks; uncertainty is never treated as zero.
No epsilon, snapping or omitted exact fallback was introduced.

## Correctness

All 18 registered exact targets pass, including new scalar, construction, box,
measure, collinear-segment and filter-routing tests. Arithmetic/EFT comparisons,
mixed 64-bit integer/binary boundaries and predicate oracles cover all four
standard rounding modes. ASan/UBSan passes include the new scalar, construction,
box-lifetime and collinear tests and the combined routing/predicate test.

All 44 default Cartesian relate/set-operation source units were rebuilt: 36
pass; four runtime failures reproduce byte-for-byte against PR1491, and four
compile failures reproduce against its headers. No new default regression.

Stable replay of 164 native/retained issue records, 82 groups / 47 selected
issues, finds no semantic changes. Six #888 records have 16 changed WKT fields;
every change is only component order or ring starting vertex. Exact point-sequence
matching verifies those permutations independently of overlay; the frozen
rational evaluator additionally confirms valid, equal sets. Raw WKT is retained.

The [issue classifications](exact_issue_comparison.md) remain unchanged:
**20 selected public-symptom improvements versus develop for PR1491 plus exact,
15 incremental versus PR1491**. These are not whole-issue closure claims or a
numerical-only causal attribution. The five already fixed by PR1491 are #1490,
#1421, #1420, #1238 and #1178.

Reproduce the stable issue comparison, preserving the frozen runners:

```sh
node test/strategies/exact_issue_comparison_optimized.js /tmp/bg-adaptive-issue-comparison
node test/strategies/exact_issue_comparison_set_check.js /tmp/bg-adaptive-issue-comparison
```

The first command intentionally reports raw WKT differences as a nonzero exit;
the second examines them without rewriting the results. `semantic-summary.json`
records the source/header/result hashes and exact equivalence evidence; a copy
is [committed here](exact_adaptive_issue_results.json).

## Measurement

The benchmark uses the same five workloads, nine API calls, 64 outer vertices
and eight prepared input pairs as before. Double inputs are integer-valued;
intersection coordinates need not be integral. This is not a representative
sample of all fractional floating-point inputs. Generation/validation are outside
timing; output allocation, conversion and destruction are included.

Pinned baselines are develop `0b4851a77`, PR1491 `72b2d832f`, initial exact
`bc65c6d60`, and the previous optimization `7f2001390`. Builds use C++14,
`-O3 -DNDEBUG -ffp-contract=fast -pthread`, with frame pointers and line tables.
Three sequential timing rounds rotate the variants, at least 0.2 seconds per
case. Allocation measurements are separate. All our compile/test jobs are stopped
during timing; unrelated macOS storage services remain active. CPU time and
per-case min/max values are retained as well as wall-clock medians.

There are 1,755 timing measurements and 90 separate allocation measurements.
The [CSV](exact_adaptive_benchmark_results.csv) contains all 675 summaries;
the [manifest](exact_adaptive_benchmark_manifest.json) pins source/header hashes.
Hashes match before and after timing. Exact before/previous/adaptive output
digests agree for every measured case. These digests are not full set oracles.
The geometric mean per-case maximum/minimum timing ratio is 1.030; aggregate
CPU-time ratios agree with wall-time ratios to the shown precision.

## Results

Geometric means across the five workloads, not an application-mix prediction:

| Group | Double Speedup vs Previous | Int64 Speedup vs Previous | Double Slowdown vs Develop | Int64 Slowdown vs Develop |
| --- | ---: | ---: | ---: | ---: |
| Four polygon set operations | 2.93x | 2.41x | 13.7x | 12.6x |
| Three relation matrices | 1.84x | 1.61x | 14.6x | 16.5x |
| Linear/areal intersection | 2.55x | 2.22x | 11.8x | 13.0x |
| Point/polygon covered_by | 1.00x | 0.92x | 1.85x | 1.52x |
| All 45 cases | 2.19x | 1.87x | 11.0x | 11.0x |

Across all cases, adaptive exact is 7.66x / 8.05x faster than the initial exact
implementation (double / int64), but still 10.3x slower than PR1491 for both.
The tiny int64 point-query microbenchmark does not improve in this pass; it is
about 8% slower than the previous exact build. These ratios include all cases,
not only the ones that improved. This is not yet default-strategy performance.

Median microseconds per convex double call:

| Operation | Develop | PR1491 | Previous Exact | Adaptive Exact |
| --- | ---: | ---: | ---: | ---: |
| Union | 4.75 | 4.79 | 195.9 | 53.5 |
| Intersection | 4.60 | 4.57 | 181.1 | 48.4 |
| Difference | 5.48 | 5.48 | 201.8 | 56.5 |
| Symmetric difference | 15.5 | 15.6 | 556.1 | 217.7 |
| Areal/areal relate | 2.43 | 2.51 | 92.9 | 26.5 |
| Linear/areal relate | 9.64 | 10.1 | 185.0 | 99.7 |
| Linear/linear relate | 6.08 | 8.54 | 140.3 | 112.5 |
| Point/polygon covered_by | 0.150 | 0.154 | 0.300 | 0.299 |
| Linear/areal intersection | 11.5 | 12.9 | 357.0 | 137.2 |

Convex int64 union goes from 120.2 to 48.6 us; develop takes 4.12 us.

Disabling the two algorithm-side exactness optimizations stays within about 1%
overall, with mixed signs. They are not the cause of this pass's speedup. Their
previously verified cluster-comparison bound remains O(n log n) rather than
quadratic; no further topology/traversal decisions were removed here.

## Filters And Remaining Costs

The separately instrumented eight-pair convex corpus evaluates **zero** native
rational lifts, lazy rational expressions, point-coordinate caches, homogeneous
point caches or ratio caches for union, intersection, difference, all three
relation matrices, covered_by and linear/areal intersection. Symmetric difference
still evaluates exact values when it re-overlays coincident rational boundaries
from its two difference results. That final fallback is retained.

The new six-second convex-union sample has 2,368 main-thread stack samples.
Exclusive top-of-stack categories include about 7.4% multiprecision (small
rational cache construction/copy), 11.0% nextafter, 9.6% interval comparison,
8.8% copying and 6.7% allocation/freeing. In the pre-filter profile,
multiprecision alone was about 62.7%. Inlined work can be attributed to its
caller; these are sampled categories, not a precise cost decomposition.

Convex-union allocations fall from about 930 per call in the previous exact
build to 325 (initial exact: 4,028). However, allocated bytes rise from about
162 KB to 211 KB (initial: 560 KB): owning filtered records are larger.
Remaining opportunities include smaller scalar/cache storage, fewer temporary
measure-vector allocations and cheaper certified interval evaluation. Faster
predicates do not justify removing genuine ring-parent or traversal choices.

Raw timings, hashes, compiler commands, profiles and allocation runs live in
`/tmp/bg-adaptive-performance`; regression artifacts are in
`/tmp/bg-adaptive-default-regressions` and `/tmp/bg-adaptive-exact-ctest.log`.

```sh
BOOST_ROOT=/Users/tinko/dev/boost node test/strategies/run_exact_benchmark.js build /tmp/bg-adaptive-performance
BOOST_ROOT=/Users/tinko/dev/boost node test/strategies/run_exact_benchmark.js run /tmp/bg-adaptive-performance
BOOST_ROOT=/Users/tinko/dev/boost node test/strategies/run_exact_benchmark.js allocations /tmp/bg-adaptive-performance
BOOST_ROOT=/Users/tinko/dev/boost node test/strategies/run_exact_benchmark.js summary /tmp/bg-adaptive-performance
```

The standalone `exact_benchmark_filter_counts.cpp` runs each prepared pair once
with opt-in diagnostics. These count native rational lifts, expression/coordinate
evaluations and construction filters, not every trivial rational cache object.
Its [double](exact_adaptive_filter_counts_double.csv) and
[int64](exact_adaptive_filter_counts_int64.csv) results are separate from timings.
