# Bounded-coordinate filters and large-budget determinant retry

Two independent CPU experiments on Intel Xeon Gold 6530, GCC 13.3.0 and
Clang 22.1.8, binary64, C++17, -O3, x86-64-v2, nearest/ties-even,
no fast-math or contraction. Each cell is one public operation on one
workload, with three input pairs. There are nine operations, five workloads,
five interleaved repeats and a dedicated physical core per experiment.
Ratios are geometric means of 45 equally weighted ratios of cell medians;
less than one is faster.

## Predicate bounds

Compare `369510c05` to the new [filter and proof](exact_lean_n_ulp_filters.md)
with the same driver, changing no output budget or retry setting. These are
whole-operation timings, not the cost of a single predicate. The before/after
experiment used the driver from `369510c05`, with the old headers or the
filter at `bace68a89`, and .12-second windows; the retry matrix below used
the driver at `21133c8ec` and .15-second windows. Both use the same corpus.

| Policy | GCC new/old | Clang new/old |
|---|---:|---:|
| Lean 8, retry off | .9967 | .9968 |
| Lean 100, retry off | .9965 | .9991 |
| Eager correctly rounded | .9845 | 1.0000 |

The GCC eager improvement occurred in all five repetitions. The fractional
percent changes for bounded output are too small to infer a consistent
speedup. Inputs, relation results, component/hole/vertex counts and quantized
geometry digests matched before and after in every cell. Quantized checks are
not proofs of rounding: independent rational tests check absolute bounds,
signs, halfway cases, subnormal products and all native/constructed masks.

## Construction retry

On a miss with a certified denominator, the selective retry tightens the
error estimate without changing the original FMA coordinate candidate.
If the tighter bound still fails, the ordinary exact refinement runs. The
on/off comparison sets the separate residual-only option to false in both
arms and keeps denominator rescue enabled in both arms.

Replaying the **614 proper crossings of the original segments** gives:

| N | Axes certified off / 1,228 | Axes certified on / 1,228 | Crossings needing refinement off | On |
|---:|---:|---:|---:|---:|
| 8 | 945 | 1,181 | 251 | 47 |
| 100 | 1,184 | 1,216 | 44 | 12 |

These are source crossings, not dynamic fallback frequencies. At N=8,
158 of the 204 saved point refinements come from the projected-local data;
at N=100, 28 of 32 do. No axis needed retry at either budget on the
untranslated projected data. A separate rational-oracle regression includes
a crossing at x=3/7: retry certifies the original FMA, one step from nearest;
retry off refines to nearest. Both outputs satisfy the contract.

### Public-call timings

| N | Layout | GCC off/develop | GCC on/develop | GCC on/off | Clang off/develop | Clang on/develop | Clang on/off |
|---:|---|---:|---:|---:|---:|---:|---:|
| 8 | stored | 1.0640 | .9123 | .8574 | 1.4368 | 1.2055 | .8390 |
| 8 | lean | 1.0571 | .9076 | .8586 | 1.4161 | 1.1905 | .8407 |
| 100 | stored | .8907 | .8714 | .9783 | 1.1817 | 1.1485 | .9719 |
| 100 | lean | .8852 | .8667 | .9791 | 1.1630 | 1.1308 | .9723 |

All five aggregate repetitions favored retry for every row on both compilers.
At N=8 the decimal workload took .76 times the old time on GCC and .73-.74
on Clang. At N=100 the corresponding factors were .95 and .93-.94.
Untranslated projected operations stayed within 1% of parity on both.
CPU/wall time was .9999 for both compilers.

The original cutoff was a tuning choice, not a numerical precondition. The
default now retries after a miss at every budget; budgets above 100 were
not timed. Set the third template parameter to false to disable it, or use
`<N, true, false>` to keep just the residual retry. All variants preserve the
same N-step bound and exact predicate fallbacks.

Input fingerprints matched across compilers. All paired same-budget output
checks agreed (including DE-9IM and quantized geometry digests). The
projected-local symmetric difference still differs from develop, as in
the preceding experiment. Eleven distinct focused tests passed under both
compilers, including the rational oracles, borrowed operations and bounds,
rounding and output policies, and shared boundaries.

## Reproduction

Set `EXACT_INCLUDE`, `DEVELOP_INCLUDE`, `BOOST_ROOT`, `CXX` and `TASKSET_CPU`
as in the [original comparison](exact_lean_ulp_benchmark.md). For each
compiler use a fresh output directory:

```sh
BENCH_COMPARE_DET_RETRY=1 BENCH_REPETITIONS=5 BENCH_SECONDS=.15 \
BENCH_ARCH=x86-64-v2 DEVELOP_REVISION=082b99181 \
sh test/strategies/run_simple_cartesian_benchmark.sh all "$out"
```

This builds develop and eight bounded variants. `results.csv` holds the
2,025 samples per compiler; `comparison.csv` holds 180 matched cells.
