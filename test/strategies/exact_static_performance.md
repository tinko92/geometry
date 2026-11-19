# Typed Exact Cartesian Strategies

Historical measurements. The shared-boundary scratch regression is addressed in
[exact_compact_performance.md](exact_compact_performance.md).

This replaces the separate segment parameter, producer pools and type-erased
scalar implementation at `b518175cf8d7c676bb9fbcea538f3c71fd49595a`.
The complete predicate formulas, data members and operation-specific pseudocode
are in [exact_static_kernel.md](exact_static_kernel.md).

## Changes

- Turns store their exact point, not an additional intersection fraction.
  Strategy-owned ordering chooses a nonconstant coordinate and reverses its
  comparison when the segment runs in the decreasing direction.
- Predicates and constructions are statically typed. Native intersections own
  their supporting endpoints; no producer pool or type-erased scalar callbacks
  remain. Transient predicates borrow existing exact points so refinements stay
  cached instead of being copied into short-lived temporaries.
- Numerical allocation is confined to indexed precision caches. Cache entries
  hold statically bounded expansions and checked fixed-capacity integers.
  Explicit rational I/O and the algorithms' ordinary containers are separate
  allocation boundaries; whole API calls are not allocation-free.
- Ring orientation, cancellation and parenting use boundary topology, not
  scalar area or area ordering. Ordinary rings, hole-free polygons, boxes and
  single positive output rings avoid the cancellation/parenting workspaces.
- Search boxes hold four outward-rounded doubles, only for conservative pruning.
  Exact geometry boxes and section working points remain separate. This avoids
  rounding implicit vertices together before turns are computed.
- Whole-endpoint adapters preserve construction provenance for clipping,
  segment predicates and boxes, without extracting allocating rational scalars.
  Default strategy paths and C++14 aggregate turn initialization are preserved.

## Correctness

All 23 registered exact-strategy tests pass. They include independent rational
controls for constructions, ordering, clipping and spatial results; 4000 random
clipping cases; close implicit points; integer extremes; cancellation; directed
rounding; cache growth; retained-output lifetime; and concurrent reads.
Warmed primitive-path allocation assertions pass, including implicit segment
adapters and simple topological ring analysis/parenting.

ASan/UBSan checks pass for combined relate/overlay, typed predicates, exact points,
borrowed point factories and the final ring fast paths. Separate precision-cache
concurrency/lifetime checks include TSan. These are local ARM checks, not a fresh
x86 run or an exhaustive proof of correctness.

The frozen issue replay has 164 records, 82 fixture groups and 47 selected issues.
There are no process failures, validator changes or non-WKT semantic changes.
Six #888 records change 16 WKT fields only by ring rotations/component order;
the independent rational evaluator confirms valid outputs and exact set equality.
[Evidence](exact_static_issue_comparison.json) retains hashes and the individual
checks, without rewriting the raw WKT.

The broader Cartesian suites were rerun against isolated baseline headers:

| Suite | Passing Sources | Unchanged Runtime Failures | Existing Compile Failures |
|---|---:|---:|---:|
| Exact-enabled, versus previous exact | 34 | 7 | 3 |
| Default, versus PR1491 | 37 | 4 | 3 |

There are no new compile/runtime failures or differing failure logs.
The [exact](exact_static_exact_regressions.tsv) and
[default](exact_static_default_regressions.tsv) source lists retain every result.
The seven exact failures are relate_areal_areal, difference, difference_multi,
intersection_integer, intersection_multi, union and union_multi. The default
failures are difference, difference_multi, union and union_multi. Both runs have
the existing intersection, intersection_segment and set_ops_ll compile failures.
These are not reported as passing tests or silently changed expectations.

## Benchmark Method

Apple M5, macOS Darwin 25.6.0, Apple Clang 21.0.0; C++14, `-O3 -DNDEBUG
-ffp-contract=fast -pthread`, identical line-table/frame-pointer flags.
The develop baseline is `aaf0560e520041a833c11ad5d61017b684b46f9a`.
Previous exact is `b518175cf8d7c676bb9fbcea538f3c71fd49595a`.

The unchanged benchmark uses five workloads, nine calls, eight prepared pairs
and 64 exterior vertices. Double inputs are integer-valued; their intersections
need not be integral. Linear/linear inputs repeat across shape labels, so this
is not 45 independent topologies or a fractional-input performance survey.
Generation/validation is outside timing; output allocation, conversion and
destruction are inside. Corpus digests observe output structure/matrices, not
full geometric equivalence; the issue replay supplies separate semantic checks.

Six serial rounds rotate all six type/baseline variants, at least 0.3 seconds
per case; summaries use the upper median of six samples. Our compilation and
correctness jobs were stopped for the recorded
run. A preceding pass overlapping suite compilation was retained only as warmup
and excluded. Allocation counters run separately after timings, not concurrently.
Other desktop applications and macOS services remain outside our control.

The [CSV](exact_static_benchmark_results.csv) contains case summaries, CPU times,
variation, allocation counts and digests; the
[manifest](exact_static_benchmark_manifest.json) pins flags and header hashes.
Raw rounds and profiles are in `/tmp/bg-static-performance`.
The manifest was captured before committing the redesign, so its dirty-header
list and complete header hashes, not its recorded HEAD alone, identify the code.

## Measurements

Geometric means over the indicated cases, not an application-mix prediction:

| Group | Double Speedup vs Previous Exact | Int64 Speedup vs Previous Exact | Double Slowdown vs Develop | Int64 Slowdown vs Develop |
|---|---:|---:|---:|---:|
| Four polygon set operations | 1.38x | 1.60x | 7.05x | 5.14x |
| Three relation matrices | 1.82x | 2.10x | 5.05x | 4.71x |
| Linear/areal intersection | 1.51x | 1.83x | 5.39x | 4.45x |
| Point/polygon covered_by | 0.99x | 0.98x | 1.86x | 1.36x |
| All 45 cases | 1.47x | 1.68x | 5.28x | 4.24x |

CPU-time aggregate speedups are also 1.47x and 1.68x. The geometric-mean
per-case maximum/minimum ratio for the new strategy is 1.022 for double and
1.030 for int64. Covered_by has no demonstrated improvement.
There are 1620 timing records, 135 separate allocation records and 405 summaries;
all previous/new exact corpus digests match. Develop differs on double star
symmetric difference, so its performance does not imply identical output.

Upper-median microseconds per convex double call:

| Operation | Develop | Previous Exact | New Exact |
|---|---:|---:|---:|
| Union | 4.83 | 31.18 | 18.55 |
| Intersection | 4.57 | 29.48 | 17.41 |
| Difference | 5.44 | 33.27 | 20.31 |
| Symmetric difference | 15.45 | 210.46 | 93.63 |
| Areal/areal relate | 2.41 | 17.18 | 8.42 |
| Linear/areal relate | 9.62 | 50.38 | 31.72 |
| Linear/linear relate | 6.06 | 71.31 | 36.89 |
| Point/polygon covered_by | 0.150 | 0.299 | 0.300 |
| Linear/areal intersection | 11.48 | 90.54 | 60.69 |

Separate whole-call allocation census, double coordinates:

| Call | Develop Allocations / Bytes | Previous Exact Allocations / Bytes | New Exact Allocations / Bytes |
|---|---:|---:|---:|
| Convex union | 79 / 31107 | 90 / 110027 | 89 / 96299 |
| Convex areal/areal relate | 12 / 7992 | 16 / 21440 | 16 / 8432 |
| Linear/areal intersection, convex workload | 88 / 63424 | 116 / 170728 | 116 / 72104 |
| Point/polygon covered_by | 0 / 0 | 0 / 0 | 0 / 0 |

These include algorithm/output containers and are not counts of numerical
allocations. Numerical scratch arithmetic is fixed-size, while indexed precision
caches may grow. The warmed primitive allocation assertions test that distinction.

## Remaining Costs

The shared-boundary rectangle union/intersection workloads regress, by 2.73x /
2.75x for double and 1.31x / 1.33x for int64 versus previous exact. Double union
is 644.5 us versus 236.2 us previously; intersection is 640.5 us versus 233.1 us.
They are included in the aggregate, not omitted as inconvenient outliers.

Five-second CPU samples of the final double binaries contain 3813 convex-union
samples and 3873 shared-boundary-union samples. In convex union, exclusive hot
frames include turn extraction (312), native orientation (262), input validation
(200), exact coordinate comparison (182) and expansion products (113).
Memory copying/clearing accounts for 325 samples; fegetround accounts for 99.
The ordinary nondegenerate path therefore still pays for larger exact points,
input checks, strategy calls and certified construction/output work.

The shared-boundary profile instead spends 2074 samples (53.6%) in memmove and
memset. Its call stacks lead through traversal's tied squared-distance/orientation
decisions into wide fixed-capacity fallback arithmetic. Proving an exact zero
does not pass a nonzero-sign filter. Removing dynamic scalar scratch allocation
did not remove the cost of clearing/copying conservative worst-case buffers.
The profile identifies a remaining degenerate-case optimization target; it does
not justify relaxing exact comparisons or omitting required traversal decisions.
Exclusive frame counts and inclusive call paths are not added together.

This is not yet a performance-compatible replacement for default strategies.
The kernel's numerical allocation goal and faster ordinary calls do not imply
uniform speedups, unbounded materialized precision, or allocation-free overlays.

## Reproduce

From this worktree with the pinned Git objects available:

```sh
export BOOST_ROOT=/Users/tinko/dev/boost
export BENCH_SECONDS=0.3 BENCH_REPETITIONS=6
node test/strategies/run_exact_benchmark.js build /tmp/bg-static-performance-new
node test/strategies/run_exact_benchmark.js run /tmp/bg-static-performance-new
node test/strategies/run_exact_benchmark.js allocations /tmp/bg-static-performance-new
node test/strategies/run_exact_benchmark.js summary /tmp/bg-static-performance-new
```

Do not run other builds/tests during timing. `BENCH_INCLUDE` can select frozen
headers; header hashes are verified before and after the run. `CXX` selects the
compiler. All runs are finite; failures abort the runner.
