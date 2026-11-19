# Double-coordinate qualification and strict x86-64-v2 timings

Local, uncommitted work, 2026-11-16. This is a Clang/x86-64-v2 measurement, not
a direct continuation of the original GCC timings. The subsequent GCC 13 rerun
matching the original configuration is in `exact_gcc13_performance.md`.

## Change and floating-point contract

Removed ten numerical volatile temporaries from four exact-kernel headers:
the four basic arithmetic helpers, FMA residual/coordinate temporaries, and the
explicit output cast. No formulas, bounds, fallback stages, caches, borrowed
supports or rounding policies changed. Corrected the orientation helper's stale
comment about supporting implicit contraction.

The contract requires round-to-nearest binary64 evaluation without excess
precision, fast-math, reassociation or implicit contraction. Explicit `std::fma`
remains fused. FP overflow/underflow remain outside scope. No release rounding-mode
checks were added. Benchmark volatile sinks remain because they prevent dead-code
elimination; unrelated geographic/projection code was not changed.

Clang 20 assembly confirms that a basic product changes from multiply, stack
store, stack reload, return to multiply, return. Error-free sums retain their
individual binary64 operations, and the explicit product-residual FMA remains.
The v2 target does not require hardware FMA; the probe uses the native libm FMA.

## Double-only qualification

Only double geometry coordinates/outputs were tested in this rerun. Existing
integer and float support was not removed, but is not qualified by these runs.

- New `exact_double_qualification`: 2,952 orientation cases, 2,822 unique
  intersections under all three output policies, 130 parallel/collinear cases,
  and four nonzero FMA residual examples. Includes decimal inputs, determinant
  cancellation, near-parallel supports, midpoint ties, public set operations,
  relation matrices and open counterclockwise inputs. Enclosures and rounding
  cells are checked against independent exact rational equations, including
  ties-to-even, not a floating conversion used as an oracle.
- Clang 20/x86-64-v2: the new qualification passes against both frozen before
  and current headers. Existing double-only ULP filters, rounding policies,
  simple strategy, and public predicate tests under all three policies pass.
  The rounding regression verifies 3,795 intersections under both nearest modes.
- Apple Clang 21/ARM64: the new qualification passes normally and with ASan/UBSan.
- All eight benchmark variants pass their digest self-tests. Exact before/after
  under every policy agrees with develop on all 27 structural/spatial controls.
  These controls are diagnostics, not independent correctness proofs.

An initial isolated Clang 18 attempt encountered an inherited `relate()` overload
ambiguity in both the unchanged checkpoint and current source. No workaround was
added to the code; all reported timings and x86 qualification use Clang 20.
This is selected-case qualification, not exhaustive testing of the double domain
or a rerun of the full mixed-coordinate regression suite.

## Benchmark

CPU-only on the RTX server, Intel Xeon Gold 6530, Ubuntu Clang 20.1.2. Develop is
the local develop tip `87f45bc3705b56e347e47f3f9809fe32b179fa7e`.
Before is the frozen pre-removal header tree; its only numerical differences
from current are the ten volatile sites. No existing checkout was overwritten.

Same compact nine-pair corpus and executor as before: three degenerate/shared
boundary cases, three general-position cases, three three-decimal WKT cases
parsed as double. Nine operations, three inputs per group, parsing/validation
outside timers. All variants use:

```text
-std=c++17 -O3 -DNDEBUG -fno-fast-math -ffp-contract=off
-march=x86-64-v2 -fno-omit-frame-pointer -pthread
```

Five serial, case-paired, rotated-order repetitions, pinned to CPU 2. No concurrent
compilation or other timing jobs from this task. Each sample runs for at least 0.1 seconds
and visits complete three-input batches. 1,080 samples, 261,655,284 calls; timed
campaign from 13:04:24 to 13:06:16 UTC. Source and header checksums are verified
before and after running.

Every timing factor below is elapsed time divided by develop, which is 1.000.
Greater than one means slower. Each aggregate is the equal-weight geometric mean
of matched per-group median times, not a pooled-throughput average.

| Policy | Before | Without numerical volatile |
| --- | ---: | ---: |
| Stored approximation | 1.549 | 1.521 |
| Correctly rounded at construction | 3.514 | 3.256 |
| Correctly rounded at output | 3.609 | 3.315 |
| Simple Cartesian | - | 0.713 |

The two nearest policies also require one-ULP initial approximations, unlike the
ordinary error-bound layout. Their extra cost is therefore not solely final
coordinate export; it is present in relation calls too. Eager is slightly faster
overall here, while deferred has an advantage on some relation workloads. Neither
nearest mode was replaced by ordinary stored-approximation output.

| Category | Simple | Stored Approx. | Nearest Eager | Nearest Deferred |
| --- | ---: | ---: | ---: | ---: |
| Degenerate | 0.714 | 1.337 | 2.641 | 2.678 |
| General | 0.720 | 1.595 | 3.129 | 3.173 |
| Decimal | 0.706 | 1.649 | 4.177 | 4.287 |
| All | 0.713 | 1.521 | 3.256 | 3.315 |

| Operation | Simple | Stored Approx. | Nearest Eager | Nearest Deferred |
| --- | ---: | ---: | ---: | ---: |
| Union | 0.845 | 1.697 | 3.364 | 3.497 |
| Intersection | 0.805 | 1.476 | 3.637 | 3.792 |
| Difference | 0.856 | 1.794 | 3.445 | 3.587 |
| Symmetric difference | 0.799 | 2.311 | 4.497 | 4.514 |
| Relate areal/areal | 0.712 | 1.554 | 4.241 | 4.233 |
| Relate linear/areal | 0.612 | 1.564 | 4.473 | 4.385 |
| Relate linear/linear | 0.660 | 1.566 | 3.450 | 3.424 |
| Point covered-by | 0.584 | 0.733 | 0.727 | 0.729 |
| Line/polygon intersection | 0.610 | 1.501 | 4.563 | 4.856 |

Simple retains its two known invalid symmetric-difference results in the general
and decimal groups. These timings are flagged, not silently excluded. On the 25
matching-control groups alone, factors are simple 0.709, stored approximation
1.446, eager 3.111, deferred 3.171. Other guarantees and parenting algorithms also
differ between simple, develop and exact: these are whole-strategy comparisons.

Across the five repetition-wide aggregates, factors range from 1.515 to 1.530
for stored approximation, 3.235 to 3.266 for eager, and 3.288 to 3.325 for deferred.
The corresponding before ranges are 1.534 to 1.554, 3.505 to 3.533, and 3.556 to
3.617. Small per-operation differences should not be treated as reliable wins.

## Reproducing and evidence

The runner now defaults to develop, simple and all three exact policies;
`BEFORE_INCLUDE` adds the three checkpoint variants. Its analyzer reports only
factors of develop and checks the declared variants, double coordinate type,
sample counts, repetitions and within-variant controls.

```sh
EXACT_INCLUDE=/path/to/current/include DEVELOP_INCLUDE=/path/to/develop/include \
BEFORE_INCLUDE=/path/to/before/include BOOST_ROOT=/path/to/boost CXX=clang++ \
BENCH_ARCH=x86-64-v2 TASKSET_CPU=2 BENCH_REPETITIONS=5 \
sh test/strategies/run_simple_cartesian_benchmark.sh all /tmp/cartesian-double
node test/strategies/summarize_simple_cartesian_benchmark.js /tmp/cartesian-double
```

Omit `BENCH_ARCH` and `TASKSET_CPU` for portable local runs. Server artifacts are
in `/home/tinko/dev/bg-exact-compact-20261116/no-volatile`, under `benchmark20` and
`qualification20`. Compilers were extracted in isolation, not system-installed.
Raw CSVs, checksums, numerical diff, assembly probes, sanitizer logs and analysis
are retained locally in the umbrella checkout's
`.geometry-plans/exact-no-volatile-20261116`. No commits or pushes were made.
