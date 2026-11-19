# Forward-Error Exact Cartesian Strategies

Historical measurements. The subsequent typed-kernel redesign is described in
[exact_static_kernel.md](exact_static_kernel.md), with new measurements in
[exact_static_performance.md](exact_static_performance.md).

This replaces the normal interval-evaluation path in the experimental exact
Cartesian strategies with explicit forward-error formulas. The previous version
is `50c99d791`; the develop comparison is `0b4851a77`.
The measured implementation is committed as `769181f1c`. The manifest was
captured before that commit, so its recorded HEAD is the previous version;
its complete header hashes identify the new implementation, not that HEAD.

## Representation And Filters

Native points and intersections are values, not individually allocated nodes.
An intersection owns the eight native coordinates of its two supporting
segments, two approximate coordinates and their absolute error bounds, and an
epoch/index cache handle. Owning the small endpoint tuple, rather than borrowing
iterators, also supports retained exact output after the inputs are destroyed.
Intersections on inherited supporting lines retain this representation during
subsequent overlays. A genuinely new line between constructed points, or an
arbitrary rational literal, uses the explicitly separate compatibility path.

First filters evaluate the predicate and its nonnegative magnitude polynomial
with a static error coefficient. They reuse differences/products where possible;
there is no runtime expression tree. Rounded wide-integer input conversion is
included in the magnitude formulas. Uncertainty requests more precision, never
an epsilon-based answer. Directed rounding bypasses the nearest-rounding filters.
Floating overflow and underflow are outside the filter contract.

Ambiguous intersection coordinates use cached, statically sized numerator and
denominator expansions, refined approximation bounds, then bounded exact binary
integers. For native coordinates with B bits on the common dyadic lattice and
scale S, homogeneous intersection numerators need at most 3B+3 bits and the
denominator at most 2B+S+3 bits. Binary64 gives 6297 and 5273 bits respectively;
integer coordinate types have S=0. Checked, allocator-free fixed backends are
chosen at compile time. Native output rounding also uses bounded arithmetic,
including nearest/ties-to-even verification.

The precision cache is an operation-scoped, thread-local vector of fixed-size
entries, addressed by epoch/index. It grows only when fallback is needed and
retains capacity between calls. Copies and retained results cannot reuse a stale
entry from another operation or thread. Separate paged pools recycle fixed-size
scalar/parameter producer records; there is no allocation per native producer
after their high-water capacity is reached. Those pools retain their high-water
pages through process termination to support retained static objects safely.

There is no claim that arbitrary rational literals or unbounded global rational
sums fit a native primitive's fixed precision bound. Mixed coordinate types that
do not match the umbrella's native type also keep their lossless literal path.
These compatibility APIs
and exceptional fallback still use dynamic rationals. Native primitive
predicates and construction caches do not need that representation.

On this Apple M5 build, an implicit double or int64 point occupies 144 bytes,
and a filtered scalar 56 bytes. A segment parameter occupies 176 bytes: its
optional compatibility value is larger than the previous pair, but is no longer
constructed/copied on the native endpoint or deferred paths. Precision-cache
entries occupy 3216 bytes for double and 1040 for int64, only on fallback.

Fixed precision alone does not guarantee allocation-free evaluation:
Boost.Multiprecision's large multiplication can allocate Karatsuba scratch even
with an allocator-free fixed backend. Native formulas therefore use the bounded
digit-multiplication helper. Tests count allocations at full-width arithmetic
boundaries, not just small or sparse operands.

Common scalar expressions first use a 512-bit exact tier and compact inline
producer payloads; values that do not fit promote to the worst-case fixed tier.
This is a checked precision choice, not an approximate answer. Dyadic input
tuples also cancel a shared power of two before intersection products are
formed. Both changes avoid work on worst-case-size operands for ordinary input
values, without changing their exact homogeneous coordinates.

## Correctness

All 24 registered exact-strategy tests pass on the frozen final headers.
The stable issue replay covers 164 native/retained records, 82 fixture groups
and 47 selected issues, with no semantic changes or failing records. Six #888
records change 16 WKT fields only by component permutations or ring starting
points; independent exact point-sequence matching checks those changes, and the
frozen rational evaluator confirms set equality. Raw results are preserved.

The default-strategy Cartesian regression run has 36 passing sources, four
byte-identical existing runtime failures and four existing compile failures.
The run enabling the exact pass in the existing helpers, compared against the
previous exact version, has 33 passes, seven byte-identical existing assertion
failures and four existing compile failures. Neither run introduces a new
regression. Existing epsilon/rounded-output expectations have not been changed
to make the exact run pass.

Native arithmetic tests include int64/uint64 extremes, dense full-width
operands, cancellation, all standard rounding modes, output rounding ties,
subnormals, cache growth, retained-output lifetime and concurrent reads. Native
construction and warm precision-cache paths have allocation assertions. The
full-width scalar batch verifies zero allocations for addition, subtraction,
multiplication, division, comparisons and aliased outputs. Sanitizer checks
cover the new kernels, cache storage and combined predicate/overlay routing.

## Reproduction

```sh
export BENCH_SECONDS=0.3 BENCH_REPETITIONS=6
node test/strategies/run_exact_benchmark.js build /tmp/bg-forward-performance
node test/strategies/run_exact_benchmark.js run /tmp/bg-forward-performance
node test/strategies/run_exact_benchmark.js allocations /tmp/bg-forward-performance
node test/strategies/run_exact_benchmark.js summary /tmp/bg-forward-performance
```

The default selection is develop, previous exact, and the current exact headers,
for double and int64. `BENCH_VARIANTS` can select the other historical variants.
The unchanged benchmark has five workloads and nine API calls, with eight
prepared input pairs and 64 exterior vertices. Its double inputs are integer
valued; this is not a representative survey of fractional inputs. Generation
and validation are outside timing; output conversion/allocation/destruction are
inside. Six sequential rounds rotate all six variants; allocation runs are separate.
The manifest records compiler flags and source/header hashes.

## Measurements

Apple M5, Apple Clang 21, C++14, `-O3 -DNDEBUG -ffp-contract=fast -pthread`.
All our builds and correctness tests were idle during timings. Other desktop
applications and macOS services remained active. There are 1620 timing records,
135 separate allocation records and 405 case summaries. Header hashes agree
before and after timing, and all previous/new exact corpus digests agree.
Those digests observe result structure/matrices, not full geometric set equality;
the separate issue replay supplies the stronger semantic checks described above.

Geometric means over the five workloads, not an application-mix prediction:

| Group | Double Speedup vs Previous | Int64 Speedup vs Previous | Double Slowdown vs Develop | Int64 Slowdown vs Develop |
| --- | ---: | ---: | ---: | ---: |
| Four polygon set operations | 1.58x | 1.53x | 9.32x | 8.22x |
| Three relation matrices | 1.75x | 1.69x | 8.67x | 9.77x |
| Linear/areal intersection | 1.62x | 1.61x | 8.09x | 7.91x |
| Point/polygon covered_by | 1.10x | 1.17x | 1.78x | 1.32x |
| All 45 cases | 1.57x | 1.54x | 7.45x | 7.08x |

CPU-time aggregate speedups are 1.57x / 1.55x, agreeing with wall time.
However, the geometric-mean per-case maximum/minimum ratio is about 2.0 for
both types. These are noisy local estimates, not precise throughput guarantees.
The [full CSV](exact_forward_benchmark_results.csv) retains every median,
minimum, maximum, CPU time and allocation count; the
[manifest](exact_forward_benchmark_manifest.json) pins the experiment.

Median microseconds for the convex double workload:

| Operation | Develop | Previous Exact | New Exact |
| --- | ---: | ---: | ---: |
| Union | 11.55 | 120.01 | 70.17 |
| Intersection | 10.31 | 108.65 | 66.96 |
| Difference | 12.07 | 125.59 | 67.07 |
| Symmetric difference | 34.59 | 487.63 | 419.23 |
| Areal/areal relate | 5.34 | 59.18 | 33.16 |
| Linear/areal relate | 20.77 | 221.36 | 91.70 |
| Linear/linear relate | 12.58 | 246.98 | 131.38 |
| Point/polygon covered_by | 0.297 | 0.666 | 0.616 |
| Linear/areal intersection | 20.01 | 315.87 | 203.09 |

Convex double union goes from 325 allocations / 211235 bytes per call to
90 / 110027; develop uses 79 / 31107. These whole-operation counts include
output containers and algorithm storage. They do not contradict the zero
allocation assertions for warmed native primitive/cache paths.

## Remaining Costs

The concave-star symmetric-difference workload regresses: about 1.46x / 1.50x
slower than the previous exact implementation (double / int64). Eight separate
alternating double pairs confirm the direction, with a geometric-mean slowdown
of 1.41x, though individual ratios range from 0.89 to 2.37. The
[paired records](exact_forward_targeted_results.json) are retained, not omitted
from the overall results. This is not default-strategy performance readiness.

Profiling first exposed full-width scalar scratch clearing/copying and eager
rational endpoint caches. The compact tier, dyadic-factor cancellation and lazy
endpoint storage address those costs without removing exact fallback. In the
final star profile, 4243 main-thread samples include about 24% exclusive exact
orientation expansion work, 15% memmove, 5% memset and 4% compact scalar
arithmetic. Inlining can attribute work to callers; these are sampled
categories, not a complete cost decomposition. Coincident intermediate
boundaries and repeated producer/predicate work remain expensive.

No additional traversal or ring-parent decisions were removed in this pass.
The algorithm changes only establish the numeric cache lifetime. Stronger
predicates do not justify dropping genuine topological choices.

Raw measurements and profiles are in `/tmp/bg-forward-performance`; the final
CTest log is `/tmp/bg-forward-final-ctest.log`, complete existing-source exact
comparison `/tmp/bg-forward-existing-final`, and default comparison
`/tmp/bg-forward-default-regressions`. Issue evidence is
`/tmp/bg-forward-issue-final`, with its compact
[semantic report](exact_forward_issue_results.json) committed here. No public
issue-improvement classifications change: the earlier selected-symptom counts
remain 20 versus develop for PR1491 plus exact, 15 incremental versus PR1491;
these are not whole-issue closure claims.
