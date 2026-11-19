# Borrowed Exact Cartesian Implementation

Implemented locally on `codex/exact-cartesian-static-kernel`. Nothing committed
or pushed. This follows the pre-borrowing redesign checkpoint, not a new change
to the default Boost.Geometry strategies.

## Representation and Construction Closure

Native inputs are borrowed. A constructed point refers to its retained actual
turn record; its two existing operations identify its original native supports.
Sorting and filtering move active pointers, not the owning records. Intermediate
edges retain original support labels through reversal, closure, copying,
duplicate/spike removal, clipping and composed operations. Native source objects
and the operation context must outlive every borrowed value. Export finishes
inside that context; implicit results are not persistent between public calls.

There is no producer pool, type-erased numerical callback, per-point atomic or
global numerical store. Fixed-size precision entries live in typed,
operation-owned indexed caches. Algorithm work containers also allocate; that
is distinct from numerical scratch allocation.

Typical polygon/linestring Root layouts on Apple arm64:

| Object | Ordinary Approximation | Eager/Deferred Nearest |
| --- | ---: | ---: |
| Point reference | 16 bytes | 16 bytes |
| Numerical state | 40 bytes | 24 bytes |
| Construction record | 72 bytes | 56 bytes |

The lean states contain two binary64 approximations and one cache index, with
no stored coordinate errors. Other native accessor alternatives, notably indexed
box vertices, can enlarge the point-reference variant.

The [capacity proof](exact_borrowed_capacity.md) gives raw formula bounds,
predicate scratch bounds and the induction over composed stages. Retained points
are native vertices or original-native-support intersections. Midpoint/splice
queries are distinct one-level types and cannot become geometry points/supports.
For grid magnitude bits `b` and scale `2^s`, intersection capacities are
`(3b+3,2b+s+3)` and midpoint capacities `(n+d+1,2d+1)`. Scratch follows actual
operand formats, including floating search metadata; it is not narrowed to the
native Root's format. This proves bounded construction families under the
algorithm's support/lifetime invariants, not global algorithmic correctness.

Ordinary output emits the current approximation without output-only refinement.
The two optional nearest-output policies establish nearest binary64 immediately,
or one-ULP accuracy immediately and nearest rounding at floating emission.
Integral output always uses ordinary BG conversion of the stored approximation.
Filters use exact fallback; overflow/underflow, fast-math and non-default rounding
modes remain outside the agreed floating execution contract.

## Correctness Qualification

- All 56 focused sources pass, including the corrected segment-box rerun.
- The additional 16-case public/type/policy qualification passes, including
  ASan/UBSan. The full nested-operation sanitizer test also passes.
- Default C++14 suite: 37 passing sources, four baseline-identical runtime
  failures and three matching baseline compile failures. No new regression.
- Frozen grid runs: 99,000 calls, no input errors or mismatches. Set operations
  use 3x3 tiles with four/six directions; relate uses 5x5 with eight directions.
  These use int64 and eager-nearest construction so integral grid intersections
  survive ordinary integer output conversion.
- Pinned issue replay: 656 records, 82 fixtures, 47 issues. No execution/evaluator
  failures, input mismatches or current internal-exact observations. Native
  representation/feedback observations remain separate: 27 approximate-output
  fixtures and 23 nearest-output fixtures. They are not automatically certified
  harmless, nor counted as internal topology failures.
- Standalone CMake configure/build/ctest passes for the borrowed-record target.
  Exact test targets use C++17 and controlled contraction; ordinary paths stay
  C++14. New tests are registered with CMake and b2.

The segment-box test originally called the public integer-input API with an
internal floating search box. Its mixed arithmetic/oracle/allocation tests remain;
public checks now require two losslessly admitted native inputs. A final missing
predicate declaration was added after the main frozen runs; the standalone
primitive test then passed with `-pedantic-errors` and ASan/UBSan, without warnings.
The benchmark binaries include that declaration. The later output-policy comment
correction changes no declarations or executable code.

The replay pins previous exact `2886f2acd`, PR1491 `72b2d832` and develop
`aaf0560e5`. Develop moved during this work; the benchmarks below instead pin
the newer `47286e8ed` (including #1514/#1515). Replay issue counts must not be
attributed to this newer develop revision.

## Whole-Operation Benchmarks

Apple M5, Apple Clang 21, C++17, `-O3 -DNDEBUG -ffp-contract=off`. Five workloads,
nine operations, eight inputs per workload, three serial rotated repetitions,
0.2 seconds per case/variant. There are 2,430 uninstrumented timing rows and
540 separately instrumented allocation rows. Ratios below are geometric means
over the same 45 cases; **larger than one means slower**.

| Ordinary Output Inputs | vs. Older Exact `2886f2acd` | vs. Pre-Borrowing Snapshot | vs. Develop `47286e8ed` |
| --- | ---: | ---: | ---: |
| Whole-number double | 0.60 | 1.70 | 2.70 |
| int64 | 0.93 | 1.48 | 3.22 |
| Fractional double | 0.59 | 1.63 | 2.70 |

Both lean nearest policies cost about 25-28% over ordinary output overall.
Deferred/eager time ratios are 0.997, 0.992 and 0.998 respectively for the rows
above: effectively tied at this measurement resolution, not a clear winner.

The smaller records are **not a throughput improvement over the immediate
checkpoint**. Ordinary whole-call allocation counts average 8.9% higher, while
allocated bytes average 7.6% lower (double corpus). These are allocation traffic,
not peak-memory figures, and include retained algorithm owners/support vectors.
The largest normal-double regressions are shared-boundary intersection (10.20x)
and union (9.41x). Both lean policies allocate about 1.97x ordinary-output bytes
overall, reaching 14.45x on some relate cases despite their smaller records.

A short sample profile of new convex union attributes 34.7% to ring topology
(8.5% in the pre-borrowing profile); BGL is 4.3% inclusive. New star relate spends
31.9% in coordinate forward comparisons and 14.4% in the construction forward
bound. Some symbols are folded/unclassified, so these are leads, not a complete
causal attribution. Further common-path optimization is still needed.

All exact variants preserve the structural control digests. Current and
pre-borrowing spatial digests match on both corpora. Native integer output differs
from the older exact/develop construction-rounding convention, as agreed. Develop
also differs on the star symmetric-difference structural/spatial digest; this
benchmark digest is not an independent spatial oracle and does not establish
which result is correct.

## Known Separate Limitation and Evidence

Mixed collection relate still misses a boundary/interior contact at `(2,1)`:
`1F2F01102` instead of `1F2001102`. This is reproducible in the previous native
and rational baselines and remains in the routing-control test. It is an
algorithmic issue, not a numerical-capacity failure. Some collection Box tuple
intersection dispatches are also unsupported by the existing BG algorithms.

Raw builds, manifests, profiles and logs are in
`/tmp/bg-exact-borrowing-20261116`. Compact evidence is retained in the parent
checkout's `.geometry-plans/exact-borrowing-evidence-20261116`. The archive
`/tmp/bg-exact-before-borrowing-20261116.tar` preserves the immediate checkpoint.
Benchmark commands are in `run_exact_benchmark.js`; use
`BENCH_VARIANTS=develop,previous,stage1,optimized,eager,deferred` and
`BENCH_DEVELOP_REVISION=47286e8edccdeb1eb74a5182f7570b7feb67e0cf`, then run
`build`, `paired`, `allocations`, `summary`. Fractional controls additionally set
`BENCH_FRACTIONAL=1 BENCH_COORDINATES=double` and use a separate output directory.
