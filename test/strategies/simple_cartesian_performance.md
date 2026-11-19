# Compact Cartesian strategy benchmark

Historical GCC measurement. The subsequent double-only Clang/x86-64-v2 rerun,
including both correctly rounded modes and removal of numerical volatile, is in
`exact_no_volatile_performance.md`. The latest matched-configuration GCC 13 rerun
is in `exact_gcc13_performance.md`.

Local, uncommitted work on `codex/exact-cartesian-static-kernel`. The preceding
performance fixes and qualification are recorded in `exact_performance_status.md`.
No existing checkout or campaign on the server was changed.

## Corpus and execution

Nine literal WKT input pairs: three with shared edges, collinear vertices or
contacts; three ordinary general-position examples, including holes and disjoint
polygons; three decimal-coordinate examples, including a shared-boundary case.
Every coordinate in the decimal category has three digits after the decimal
point and is parsed directly as double. Parsing and input validation are outside
the timers. There are no star/comb workloads, area or length measurements.

The nine operations are polygon union, intersection, difference and symmetric
difference; areal/areal, linear/areal and linear/linear relation matrices; point
covered-by; and line/polygon intersection. Each timed batch visits all three
inputs in its category. Output containers are reused, but every library call
performs its ordinary operation setup and cleanup.

CPU-only on `10.26.12.47`, Intel Xeon Gold 6530, GCC 13.3.0. Processes run serially
on CPU 2, with three rotated-order repetitions and 0.1 seconds per group/variant.
Flags: C++17, `-O3 -DNDEBUG -ffp-contract=off -fno-omit-frame-pointer -pthread`.
No GPU code is used. The final timing matrix takes 25 seconds, excluding builds
and untimed checks: 243 samples, 27 groups per variant, 79,040,079 total calls.
Median within-group sample spreads are 0.94% exact, 0.67% develop and 0.66%
simple, with individual worst spreads up to 10.25%. Repetition-wide exact/develop
ratios range from 1.2458 to 1.2586; CPU-based and wall-time ratios agree.

Develop is pinned to `87f45bc3705b56e347e47f3f9809fe32b179fa7e`; its headers are
identical to the preceding benchmark's `47286e8ed` pin. Exact uses the current
borrowed implementation with ordinary stored-approximation output. Its 1,226
previously qualified headers remain byte-for-byte unchanged; only the separate
simple-strategy header was added.

## Results

Ratios are geometric means of matched per-group median wall times, with equal
weight per operation/category; greater than one means slower. Incorrect simple
outputs remain in the timing table and are flagged below, not silently excluded.

| Category | Exact / Develop | Simple / Develop | Exact / Simple |
| --- | ---: | ---: | ---: |
| Degenerate | 1.096 | 0.585 | 1.873 |
| General | 1.320 | 0.635 | 2.080 |
| Decimal | 1.359 | 0.597 | 2.276 |
| All | 1.253 | 0.605 | 2.070 |

Geometric mean microseconds/call across the three categories:

| Operation | Develop | Exact | Simple |
| --- | ---: | ---: | ---: |
| Union | 6.799 | 9.855 | 5.236 |
| Intersection | 5.482 | 6.624 | 3.928 |
| Difference | 6.803 | 10.491 | 5.259 |
| Symmetric difference | 23.440 | 48.571 | 17.104 |
| Relate areal/areal | 2.198 | 2.704 | 1.363 |
| Relate linear/areal | 3.468 | 4.601 | 1.952 |
| Relate linear/linear | 3.289 | 3.892 | 1.837 |
| Point covered-by | 0.076 | 0.042 | 0.023 |
| Line intersection | 3.381 | 4.297 | 1.998 |

Exact and develop agree on all 27 structural/spatial controls. Simple agrees on
25 groups but produces an invalid self-touching symmetric-difference ring for
the first general and first decimal input, merging two components into one.
These are not merely different ring starting points or output rounding hashes.

For the general input, segments `(2,8)-(7,7)` and `(1,3)-(4,9)` meet at the exact
point `(37/11,85/11)`. The ordinary intersection formula evaluated in opposite
segment orders gives ordinates `7.7272727272727275` and `7.7272727272727266`.
Symmetric difference unions its two directional differences; those rounded
intermediates no longer have consistent coincident boundaries under strict
double predicates. This demonstrates the intended non-robust baseline's limits,
not a reason to add epsilon or exact fallback to it. Timing comparisons involving
those two groups are not comparisons of equally correct outputs.

Digests retain component/hole/point counts and relation outcomes, plus sorted
native edges quantized to `2^-16`. They are diagnostics, not independent proofs
of spatial equivalence. Actual output inspection and validity checks establish
the two simple-strategy discrepancies described above.

## Strategy and qualification

`strategies::relate::simple_cartesian` is a stateless, C++14-compatible baseline
for native two-dimensional Cartesian doubles. Orientation, equality, coordinate
ordering, winding, clustering, distance signs and intersection construction use
ordinary double arithmetic and raw comparisons. Intersection ordering uses a
coordinate along the segment, without a stored ratio. There are no numerical
allocations, precision caches, exact headers or adaptive filters. The certified
predicate flag remains false; ordinary floating recovery paths remain enabled.

It uses stock signed-area ring parenting with raw comparisons, while exact uses
topological parenting. These are whole-operation comparisons of strategy
configurations, not isolated measurements of numerical-kernel overhead.
The public `equals()` algorithm has separate legacy tolerance shortcuts outside
this benchmark's relate/set-operation call graph; that API is explicitly not
covered by this experimental baseline's strictness claim.

The new regression covers analytic set results, relation matrices, line clipping,
segment endpoint/overlap/degeneracy permutations, strict adjacent-double
comparisons, box helpers, and open counterclockwise inputs. It deliberately
checks that a rounded determinant can cancel to zero, rather than claiming
robustness. Normal C++14, ASan/UBSan, GCC/x86 C++14, and CMake/CTest pass.
Dependency inspection finds no exact-strategy headers in the new test.
The existing exact predicate regression also passes on the server: 872 rational
triples, 160 measure comparisons, 35 allocation checks and zero allocations.

## Running

Provide the current and pinned develop include trees and the existing Boost
umbrella. The POSIX shell runner requires no Node on the server:

```sh
EXACT_INCLUDE=/path/to/exact/include DEVELOP_INCLUDE=/path/to/develop/include \
BOOST_ROOT=/path/to/boost CXX=g++ SIMPLE=1 TASKSET_CPU=2 \
sh test/strategies/run_simple_cartesian_benchmark.sh all /tmp/cartesian-run
```

Omit `TASKSET_CPU` on macOS and use `CXX=clang++`. Use `run` instead of `all` to
reuse checked binaries; choose a fresh output directory for a new build/run.
The runner refuses to overwrite raw samples. Analyze copied results locally:

```sh
node test/strategies/summarize_simple_cartesian_benchmark.js /tmp/cartesian-run
```

Server artifacts: `/home/tinko/dev/bg-exact-compact-20261116/{initial,final}`.
Local raw results: `/tmp/bg-compact-20261116`. Durable evidence is retained in
the umbrella checkout's `.geometry-plans/exact-compact-benchmark-20261116`.
