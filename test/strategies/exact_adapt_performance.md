# Mechanical Adaptive-Fallback Split: GCC 13 / C++20

Experiment on `experiment/exact-and-fast-cartesian-strategies`, 2026-11-16.
Checkpoint `d71bf1fe1` preserves the previous implementation; `1a057b2b8` contains
only the mechanical split. **The split regresses the ordinary exact mode in this
benchmark. It is retained as a separate experiment, not recommended as a
performance improvement.** No follow-up tuning or selective rollback was made.

## Scope

Seven exact-kernel headers changed. Initial filter failure paths move into
`BOOST_NOINLINE` helpers with `_adapt` suffixes; first-stage successes use
`[[likely]]` and failure branches use `[[unlikely]]`. Caller operands, mutable
results and other required locals are explicit arguments; member helpers also
use the original `this`. Fixed two-stage filter loops are split after their
first iteration with the same refinement and comparison order.

The files cover borrowed and legacy point orientation, coordinate/equality and
direction predicates, measure comparisons, intersection approximation/refinement,
nearest output, borrowed midpoint/splice construction and box-coordinate output.
No formulas, bounds, capacities, policies, numerical representations, operation
ownership, source supports or algorithmic decisions changed. Existing separated
Shewchuk code and pure `try_*` certificates remain untouched. Static reviews
found no mechanical-equivalence or helper-visibility errors.

Only benchmark output checks were rerun, not the unit suites, issue reproducers,
sanitizers or random-grid campaigns. Historical issue status is summarized
separately in `exact_issue_status.md`; this experiment adds no issue-fix claims.

## Measurement

Intel Xeon Gold 6530 server, GCC 13.3.0, default target, CPU 2. All eight variants
are rebuilt in **C++20**, including the unsplit checkpoint and develop baseline.
Develop is `87f45bc3705b56e347e47f3f9809fe32b179fa7e`.

```text
-std=c++20 -O3 -DNDEBUG -fno-fast-math -ffp-contract=off
-fno-omit-frame-pointer -pthread
```

Same compact double-only corpus: three degenerate/shared-boundary, three general
and three three-decimal WKT pairs. Nine operations, three inputs per group.
Parsing/validation outside timers. Five serial case-paired rotated repetitions,
0.1 seconds per sample, no concurrent compilation from this task. The timed run
lasted 14:00:30-14:02:22 UTC: 1,080 samples and 314,067,684 complete calls. Source
and header checksums are checked before and after execution.

Every factor below is elapsed time divided by develop, which is 1.000. Lower is
faster. Aggregation is an equal-weight geometric mean of matched group medians.
Earlier C++17 numbers are not used as the unsplit reference.

| Policy | Unsplit C++20 | Split C++20 |
| --- | ---: | ---: |
| Stored approximation | 1.205 | 1.250 |
| Correctly rounded at construction | 2.854 | 2.870 |
| Correctly rounded at output | 2.891 | 2.910 |
| Simple Cartesian | - | 0.595 |

Both nearest policies also require one-ULP initial approximations. Their cost is
not solely output rounding. Small nearest-policy differences should not be
treated as reliably established improvements or regressions.

| Category | Approx. Before | Approx. After | Eager Before | Eager After | Deferred Before | Deferred After |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Degenerate | 1.071 | 1.115 | 2.465 | 2.489 | 2.492 | 2.510 |
| General | 1.252 | 1.296 | 2.651 | 2.652 | 2.674 | 2.686 |
| Decimal | 1.307 | 1.351 | 3.559 | 3.580 | 3.626 | 3.655 |

| Operation | Approx. Before | Approx. After | Eager Before | Eager After | Deferred Before | Deferred After |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Union | 1.447 | 1.503 | 3.342 | 3.321 | 3.428 | 3.419 |
| Intersection | 1.245 | 1.290 | 3.616 | 3.591 | 3.736 | 3.709 |
| Difference | 1.566 | 1.617 | 3.476 | 3.470 | 3.563 | 3.564 |
| Symmetric difference | 2.041 | 2.101 | 4.709 | 4.748 | 4.724 | 4.767 |
| Relate areal/areal | 1.065 | 1.109 | 3.378 | 3.349 | 3.394 | 3.371 |
| Relate linear/areal | 1.164 | 1.205 | 3.427 | 3.418 | 3.387 | 3.367 |
| Relate linear/linear | 1.085 | 1.125 | 2.497 | 2.537 | 2.495 | 2.544 |
| Point covered-by | 0.547 | 0.577 | 0.547 | 0.577 | 0.547 | 0.576 |
| Line/polygon intersection | 1.267 | 1.303 | 4.020 | 4.009 | 4.166 | 4.177 |

Repetition-wide approximate/develop aggregates range from about 1.202 to 1.206 before
and 1.246 to 1.251 after: this is a consistent negative result on this corpus.
The experiment does not isolate individual helpers, attributes, inlining or
code-layout effects; no cause is inferred without further profiling.

All six exact variants match develop's 27 structural/spatial controls, and all
eight variants pass digest self-tests. Controls are diagnostics, not numerical
or topology proofs. Simple retains its two known invalid symmetric-difference
cases in the general and decimal groups; they remain included and flagged.

## Reproducing and Evidence

The production benchmark runner remains unchanged. For this experiment a
temporary copy changed only its two `-std=c++17` occurrences to `-std=c++20`:

```sh
sed 's/c++17/c++20/g' test/strategies/run_simple_cartesian_benchmark.sh > /tmp/run_cxx20.sh
BENCH_SOURCE="$PWD/test/strategies/simple_cartesian_benchmark.cpp" \
EXACT_INCLUDE=/path/to/split/include BEFORE_INCLUDE=/path/to/unsplit/include \
DEVELOP_INCLUDE=/path/to/develop/include BOOST_ROOT=/path/to/boost CXX=g++-13 \
BENCH_ARCH= TASKSET_CPU=2 BENCH_REPETITIONS=5 \
sh /tmp/run_cxx20.sh all /tmp/cartesian-adapt
node test/strategies/summarize_simple_cartesian_benchmark.js /tmp/cartesian-adapt
```

Server artifacts: `/home/tinko/dev/bg-exact-compact-20261116/split-adapt-gcc20`.
Raw CSVs, manifests, controls, numerical diff and `analysis.json` are retained
locally in the umbrella checkout's `.geometry-plans/exact-adapt-gcc20-20261116`.
