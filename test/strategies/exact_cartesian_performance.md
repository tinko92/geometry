# Exact Cartesian Performance

These are the pre-adaptive-filter measurements. The later filtering pass and
fresh comparisons are recorded in [exact_adaptive_performance.md](exact_adaptive_performance.md).

Measured locally on Apple M5 / macOS 26.6.2, Apple Clang 21.0.0, 2026-11-11.
The optimized engine is commit `7f2001390` on
`codex/exact-cartesian-overlay-relate`. Nothing was pushed.

## Baselines And Method

- Official develop: `0b4851a77d91633cf86c3a5e679e67d8c60200e9`.
- PR1491: `72b2d832fb9596a2a23c0c164a5426ae720d0cd6`.
- Initial exact implementation: `bc65c6d60fe2d891ed978d0c8cdf91bb678712b9`.
- Both public remote heads were checked against these IDs on the measurement day.
- `-std=c++14 -O3 -DNDEBUG -ffp-contract=fast -pthread`, frame pointers and line tables enabled identically.
- Five workloads: disjoint convex polygons, intersecting convex polygons, stars,
  polygons with holes, and subdivided orthogonal boundaries with shared edges.
- 64 outer-ring vertices, eight prepared input pairs per workload, nine operations.
  Line tests use eight two-point members per input. Hole rings have 32 vertices.
  Linear/linear inputs repeat across shape labels; disjoint/convex also share A
  for point/areal queries. The matrix is not 45 independent topologies.
- Integer-valued `double` and `int64_t` inputs are identical between all builds.
  Their intersection coordinates need not be integral. These measurements do not
  establish performance on arbitrary fractional floating-point input distributions.
- Three sequential rounds per binary, at least 0.2 seconds per case per round;
  variants rotated between rounds. No compilation or other test jobs during timing.
  macOS background storage services remained active. The geometric mean of the
  per-case maximum/minimum timing ratio was 1.030.
- Geometry generation, input validation and parsing are outside timing. Complete
  public calls, output allocation/conversion/destruction, and a small observable
  output checksum are included. Operation names are decoded outside timing.
- Single-thread throughput: 1,485 timing measurements, plus 90 separate allocation
  measurements. Allocation counters were disabled in the timed binaries.

The [CSV](exact_cartesian_benchmark_results.csv) contains medians, minima, maxima,
CPU time, allocation counts and output digests for all 585 variant/case summaries.
The [manifest](exact_cartesian_benchmark_manifest.json) records compiler, revisions,
flags and complete-header-tree hashes. Raw rounds and CPU profiles are in
`/tmp/bg-exact-performance`. Hashes were checked before and after timing.

## Aggregate Results

Each entry is a geometric mean over the five workloads and the indicated
operations, not an extrapolation to an application mix. A slowdown above one
means the optimized exact strategy is slower than the baseline.

| Operation Group | Double Speedup Over Initial Exact | Int64 Speedup Over Initial Exact | Double Slowdown vs Develop | Int64 Slowdown vs Develop |
| --- | ---: | ---: | ---: | ---: |
| Four polygon set operations | 2.44x | 3.07x | 39.6x | 30.1x |
| Three relation matrices | 2.46x | 2.67x | 26.5x | 26.0x |
| Linear/areal intersection | 4.44x | 5.72x | 29.9x | 28.2x |
| Point/polygon covered_by | 33.7x | 54.8x | 1.89x | 1.39x |
| All 45 cases | 3.50x | 4.33x | 23.9x | 20.2x |

Across all 45 cases, slowdown versus PR1491 is 22.5x for double and 18.9x for
int64. Against develop with its existing robust-side option enabled, the double
baseline remains in the same performance range as ordinary develop; enabling a
robust orientation predicate is not equivalent to exact constructions and exact
predicates on those constructions.

### Intersecting Convex Polygons

Median microseconds per complete API call, double coordinates:

| Operation | Develop | PR1491 | Initial Exact | Optimized Exact |
| --- | ---: | ---: | ---: | ---: |
| Union | 4.83 | 4.84 | 554 | 195 |
| Intersection | 4.60 | 4.60 | 526 | 182 |
| Difference | 5.47 | 5.47 | 552 | 200 |
| Symmetric difference | 15.5 | 15.7 | 1598 | 558 |
| Areal/areal relation | 2.41 | 2.52 | 365 | 92.6 |
| Linear/areal relation | 9.62 | 10.2 | 470 | 181 |
| Linear/linear relation | 6.00 | 8.54 | 247 | 133 |
| Point/polygon covered_by | 0.150 | 0.150 | 9.92 | 0.302 |
| Linear/areal intersection | 11.4 | 12.7 | 1567 | 350 |

For the same convex union with int64 coordinates: develop 4.19 us, PR1491
4.19 us, initial exact 441 us, optimized exact 120 us. Int64 covered_by improves
from 4.19 us to 0.0797 us, versus develop's 0.0546 us.

## Implemented Optimizations

1. Homogeneous native point equality, ordering, orientation and winding use
   stack-only views. Mixed types and implicit constructed points keep the owning
   exact representation. Finite-input checks remain in place.
2. Native float/double orientation first tries a certified `4u` determinant bound.
   Nearest rounding and exponent guards prevent invalid use; explicit rounding
   barriers prevent contraction from changing the proof. Uncertain cases retain
   interval, expansion and exact fallback stages. Integers retain checked-width
   arithmetic and arbitrary-precision fallback.
3. Native segment pairs rejected by exact side signs allocate no producer nodes.
   Potential intersections still own their producing endpoints.
4. Winding rejects segments strictly above/below the query using exact ordinate
   comparisons, and handles horizontal boundary membership without orientation.
5. Leaf coordinate caches read the exact stored scalar directly. Full-point
   conversions preserve input leaves and producer ownership rather than rebuilding
   literal nodes through coordinate setters. Partial-coordinate copies are unchanged.
   Matching-type input leaves return to native output without rational conversion;
   constructed and mixed-type outputs still use correctly rounded conversion.
6. Bounds store exact rational scalars directly, without constructing implicit
   point nodes for every bound update. They do not use approximate bounds.
7. Strategies may explicitly certify predicate consistency. Such strategies skip
   the numerical zero-side recovery determinant, and may provide an exact ordering
   whose equivalence classes equal cluster equality. Clusters then form contiguous
   equality runs after sorting. Default strategies retain their previous paths.

All these changes retain exact answers. No epsilon, heuristic sign, omitted
uncertain fallback, dangling producer reference or fast-math assumption was added.

### Algorithm-Only Contribution

`BOOST_GEOMETRY_EXACT_DISABLE_ALGORITHM_OPTIMIZATIONS` leaves the exact kernel
unchanged but restores the previous side recovery and cluster search. Those two
algorithm changes alone stay within about 1% in aggregate here, with mixed signs
and no convincing aggregate speedup above timing noise. The largest measured
shared-boundary double improvements were about 4%; no int64 case exceeded 2.3%.
They are not the main source of the overall speedup in this matrix.

Their asymptotic benefit is real: the regression with 256 distinct collinear turn
points needs 32,640 cluster-equality calls previously, versus 255 after sorting.
Total sorting/grouping work is O(n log n), rather than the previous quadratic
same-ordinate search. Shared endpoints, duplicate constructions, discarded turns,
node unification, and points closer than one output ULP are covered by tests.

### Allocation And CPU Profile

For convex union, allocations decrease from approximately 4,028 / 560 KB per call
to 930 / 162 KB. Covered_by decreases from 192 allocations / 44.5 KB to zero.
These counts come from separately instrumented runs and are not timing results.

Six-second CPU samples of final convex union contained 2,369 exact and 2,364
develop stack samples. Exclusive top-of-stack categories for exact were about
62.7% multiprecision arithmetic, 16.2% memory copying, and 7.6% allocation/freeing;
nextafter accounted for 0.4%. Develop used no multiprecision arithmetic and spent
about 6.6% in copying and 25.5% in allocation/freeing. Inclusive call-tree costs
overlap and must not be summed with these exclusive categories.

The remaining overhead is therefore not principally the graph search or removed
numerical recovery in this workload. Rational bound updates/comparisons, exact
ring measures, scalar conversion and their copies remain expensive. Exactness
does not justify removing genuine clustering, containment, degeneracy, traversal
or ring-parent decisions. Further filtered/lazy scalar and measure work is needed
before considering this strategy as a performance-compatible default.

## Correctness Checks And Qualifications

- Twelve exact test targets pass, including native filter/grid equivalence,
  algorithm optimizations ON/OFF, exact constructions, output rounding, mixed
  types, concurrency, collections, and issue replays.
- ASan/UBSan point tests and UBSan native-filter tests pass.
- All 164 optimized issue records across 82 groups / 47 issues exactly match the
  frozen exact implementation's rational WKT, predicates and other semantic fields.
- All 44 default Cartesian set-operation/relate compilation units were rerun:
  36 pass, four runtime failures match PR1491 byte-for-byte, and four compile
  failures have matching baseline diagnostics. No new default regression.
- Those pre-existing runtime failures are difference, difference_multi, union and
  union_multi. Compile failures are intersection, intersection_segment, set_ops_ll
  and union_other_types. These are not presented as clean-suite passes.
- No new x86 run was performed; these measurements and the fresh regressions are local ARM results.

Before/after exact checksums match for every timed case. One develop-versus-exact
checksum differs: double star symmetric difference. In the eight-pair diagnostic,
all exact outputs have 80 components and pass both native and rational validation;
all eight develop rounded outputs fail rational validation, despite seven passing
the tolerant native validator. One develop result has 79 components. The registered
`exact_benchmark_check` preserves this construction-consistency regression; compile
the same source without `BG_BENCH_EXACT` to inspect a native baseline, with `--wkt`
for full output. This is not an additional counted GitHub issue, and timing does
not imply that native and exact outputs are geometrically identical.

Final coordinate rounding remains capable of collapsing topology in other cases.
The [issue comparison](exact_issue_comparison.md) records **20 selected public
symptom improvements for combined PR1491-plus-exact versus develop**, of which
**15 are incremental versus PR1491**. Five already belong to PR1491. Five further
retained-only successes are not counted as public fixes. These are selected-symptom
counts, not claims that entire issue reports can be closed.

## Reproduce

From this worktree, with the pinned Git objects available:

```sh
node test/strategies/run_exact_benchmark.js build /tmp/bg-exact-performance
node test/strategies/run_exact_benchmark.js run /tmp/bg-exact-performance
node test/strategies/run_exact_benchmark.js allocations /tmp/bg-exact-performance
node test/strategies/run_exact_benchmark.js summary /tmp/bg-exact-performance
```

`BOOST_ROOT` selects the umbrella include directory and `CXX` the compiler.
Builds are serial and cached by source/header/compiler hashes. Keep other build
and test jobs stopped during `run`; all runs are finite. Compilation failures and
nonzero benchmark exits abort the runner; performance checksums are not complete
correctness oracles. The benchmark is deliberately not registered as a unit test.

## Projected Fixture Pair

These opt-in workloads are not part of the historical results above. They model
two overlapping, meandering corridor footprints a few kilometres long, with
bent crossing routes, roughly parallel routes, and inside/outside point queries.
They are synthetic planar metre coordinates, not data from a particular CRS or
survey. Each corpus member varies location and bend phase; `--vertices` controls
the footprint detail and `--corpus` controls the number of pairs.

- `projected` adds `(450000, 5200000)` metres. All easting and northing values
  lie in `[2^18, 2^19)` and `[2^22, 2^23)`, respectively, so Sterbenz guarantees
  exact subtraction of any two input coordinates on the same axis.
- `projected-local` is the identical input before that translation and straddles
  zero on both axes. It lacks the positive-box guarantee; that does not mean its
  differences are inexact. Both frames use the same vertices snapped once to the
  projected binary64 grids (`2^-34` and `2^-30` metres). The generator checks
  lossless translation for every vertex and query, not just for the envelope.

Use `--workload projected-pair` for both workloads, or either workload name for
one. They are double-only; `--workload all` still selects the original five
workloads. For runner use, start with a separate output directory:

```sh
export BENCH_PROJECTED=1 BENCH_COORDINATES=double
node test/strategies/run_exact_benchmark.js build /tmp/bg-exact-projected
node test/strategies/run_exact_benchmark.js paired /tmp/bg-exact-projected
node test/strategies/run_exact_benchmark.js summary /tmp/bg-exact-projected
```

`run` and `allocations` also honor this selection. Do not combine it with
`BENCH_FRACTIONAL`; that control applies to the historical workloads.

The exact `exact_benchmark_check` target (or the benchmark binary with
`--projected-self-test 1`) compares the frames on every corpus member and all
nine operations. It checks structural digests and predicates, including the
expected point-query answers, plus spatial digests after rebasing the projected
outputs to local coordinates. The normal CSV still hashes output in its native
frame for comparisons between strategies on the *same* workload. The normalized
check is an export regression signal, not a proof of exact constructed geometry
or a claim that ordinary Cartesian overlay is translation invariant.

In the local smoke run, `optimized` and `previous` agreed on both digests for
all 18 workload/operation combinations. The pinned `develop` baseline differed
from exact on local-frame symmetric difference in both timing and allocation
builds. Keep that discrepancy visible when comparing timings; strict
`BENCH_REQUIRE_GEOMETRY_MATCH=1` reports it as a failure.
