# Double-coordinate GCC 13 rerun

Local, uncommitted measurement, 2026-11-16. This is the continuation of the
original GCC benchmark in `simple_cartesian_performance.md`, not a comparison
against the intervening Clang/x86-64-v2 configuration.

## Setup and checks

Same Intel Xeon Gold 6530 server, CPU 2, GCC 13.3.0, default x86-64 target and
develop `87f45bc3705b56e347e47f3f9809fe32b179fa7e`. Same double-only nine-pair
corpus and executor. No library or benchmark source changes in this rerun.

```text
-std=c++17 -O3 -DNDEBUG -fno-fast-math -ffp-contract=off
-fno-omit-frame-pointer -pthread
```

No `-march` flag. Explicit `-fno-fast-math` restates the earlier default. `cmp`
confirms that newly built develop and pre-removal exact binaries are each
byte-for-byte identical to their original GCC-run counterparts. This removes
uncertainty about the minor benchmark-source additions for selecting policies.

Five serial, case-paired, rotated repetitions, 0.1 seconds per group, no competing
compilation from this task. 1,080 samples, 318,259,494 completed calls, all positive
complete-corpus batches. Timed run: 13:32:28-13:34:20 UTC. Source and header
checksums verified before and after. Eight variants include old/current headers
for every exact output policy, develop and simple.

The new double-only correctness qualification passes with both old and current
headers: 2,952 orientation cases, 2,822 intersections under all three policies,
130 parallel/collinear cases and four nonzero FMA residual examples. This is
selected-case qualification, not an exhaustive domain or full-suite rerun.
All benchmark digest self-tests pass. All six exact variants match develop's
27 structural/spatial controls; controls are diagnostics, not correctness proofs.

## Results

Every factor is elapsed time divided by develop, which is 1.000. Lower is faster.
Aggregates are equal-weight geometric means of per-group median times.

| Policy | Before | Without numerical volatile |
| --- | ---: | ---: |
| Stored approximation | 1.252 | 1.214 |
| Correctly rounded at construction | 3.101 | 2.864 |
| Correctly rounded at output | 3.154 | 2.890 |
| Simple Cartesian | - | 0.605 |

The earlier GCC stored-approximation factor was 1.253; the repeated checkpoint
factor of 1.252 reproduces it. The corresponding category factors before were
1.095, 1.317, 1.361, compared with the earlier 1.096, 1.320, 1.359.

| Category | Simple | Stored Approx. | Nearest Eager | Nearest Deferred |
| --- | ---: | ---: | ---: | ---: |
| Degenerate | 0.584 | 1.084 | 2.472 | 2.485 |
| General | 0.636 | 1.261 | 2.665 | 2.678 |
| Decimal | 0.597 | 1.310 | 3.567 | 3.628 |
| All | 0.605 | 1.214 | 2.864 | 2.890 |

| Operation | Simple | Stored Approx. | Nearest Eager | Nearest Deferred |
| --- | ---: | ---: | ---: | ---: |
| Union | 0.770 | 1.422 | 3.232 | 3.317 |
| Intersection | 0.718 | 1.185 | 3.431 | 3.535 |
| Difference | 0.781 | 1.518 | 3.348 | 3.426 |
| Symmetric difference | 0.726 | 1.940 | 4.534 | 4.530 |
| Relate areal/areal | 0.616 | 1.207 | 3.869 | 3.858 |
| Relate linear/areal | 0.562 | 1.245 | 3.678 | 3.582 |
| Relate linear/linear | 0.557 | 1.130 | 2.557 | 2.573 |
| Point covered-by | 0.305 | 0.555 | 0.552 | 0.552 |
| Line/polygon intersection | 0.592 | 1.226 | 3.840 | 3.940 |

Both nearest policies also require one-ULP initial approximations, so their extra
cost is not solely final export. Simple retains its two known invalid symmetric
difference cases in the general and decimal groups; their timings are included
and flagged. On the 25 matching-control groups alone, factors are simple 0.598,
stored approximation 1.151, eager 2.705, deferred 2.732. These are whole-strategy
comparisons with different guarantees and parenting, not isolated numeric costs.

The five repetition-wide current exact/develop aggregates range from 1.212 to
1.216 for stored approximation, 2.860 to 2.866 for eager, and 2.888 to 2.892 for
deferred. They do not establish broad-workload performance guarantees.

## Reproducing and evidence

```sh
EXACT_INCLUDE=/path/to/current/include DEVELOP_INCLUDE=/path/to/develop/include \
BEFORE_INCLUDE=/path/to/before/include BOOST_ROOT=/path/to/boost CXX=g++-13 \
BENCH_ARCH= TASKSET_CPU=2 BENCH_REPETITIONS=5 \
sh test/strategies/run_simple_cartesian_benchmark.sh all /tmp/cartesian-gcc13
node test/strategies/summarize_simple_cartesian_benchmark.js /tmp/cartesian-gcc13
```

Server artifacts: `/home/tinko/dev/bg-exact-compact-20261116/no-volatile`, under
`benchmark-gcc13` and `qualification-gcc13`. Raw samples, controls, checksums,
qualification logs and `analysis-gcc13.json` are retained in the umbrella
checkout's `.geometry-plans/exact-no-volatile-20261116`. No commits or pushes.
