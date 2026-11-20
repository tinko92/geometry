# Lean versus stored bounded-ULP output

Implementation `8c729e7c8`, measurement driver `7169b4c2e`, upstream
develop `082b99181`. Native double inputs and outputs, Intel Xeon Gold 6530
CPU only, GCC 13.3.0 and Clang 22.1.8, C++17, -O3, x86-64-v2, no fast-math
or FP contraction. Lower ratios mean less time. Each compiler ran five
interleaved measurements per variant on one physical core, rotating order by
cell; its own develop binary was timed in the same rotation. Each ratio is the
geometric mean of the per-cell ratios of median wall nanoseconds per public
call. There are 45 equally weighted cells: nine operations, five workloads
(degenerate, general, decimal, projected and projected translated near the
origin), three geometries per workload. Each compiler produced 2,475 timed
rows. Process CPU / wall time was .9999 on both. No GPU was used.

| N | GCC stored/develop | GCC lean/develop | GCC lean/stored | Clang stored/develop | Clang lean/develop | Clang lean/stored |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 1.0492 | 1.0457 | .9967 | 1.4095 | 1.3979 | .9918 |
| 2 | 1.0009 | .9972 | .9963 | 1.3393 | 1.3270 | .9908 |
| 4 | .9457 | .9440 | .9982 | 1.2571 | 1.2464 | .9915 |
| 8 | 1.0639 | 1.0597 | .9961 | 1.4407 | 1.4219 | .9869 |
| 100 | .8895 | .8878 | .9981 | 1.1840 | 1.1658 | .9846 |

The layout saves 16 bytes per root-bound construction, not 16 bytes per
geometry vertex. The public crossing allocation probe on this x86 system
made 10 allocations in either layout (928 bytes lean, 960 stored).
The observed GCC aggregate benefit is only .2-.4%; three of the five GCC
budgets had individual repetition ratios on both sides of parity.
Clang improved .8-1.5%. At N=100, linear/linear relation took .957/.949
times stored and symmetric difference took 1.045/1.018 (GCC/Clang).
This is not evidence for unconditionally replacing the stored policy.

### Why eight is not faster than four

These timings predate the [independent retry on/off comparison](exact_large_budget_retry.md).
The following cutoff describes the implementation measured here.

The inherited default construction tries extra determinant corrections at
N <= 4. N=8 uses a looser threshold but does not get that retry. Outside the
timed loop, replaying the original proper segment crossings for the three
geometry pairings gives the following counts. Collinear overlaps, endpoint
contacts, subsequent constructions and repeated visits by the algorithms
are not included, so these are **not** dynamic fallback rates.

| Workload | Axes examined | Certified N=4 | Certified N=8 |
|---|---:|---:|---:|
| Degenerate | 20 | 17 | 12 |
| General | 66 | 61 | 30 |
| Decimal | 54 | 43 | 26 |
| Projected | 544 | 544 | 544 |
| Projected near origin | 544 | 475 | 333 |
| Total | 1,228 | 1,140 | 945 |

### Verification

The native input fingerprints agree across all binaries and both compilers.
Every same-N lean/stored pair agrees on components, holes, point counts,
relation results and quantized spatial checks. Develop differs on projected
near-origin symmetric difference, as before. These benchmark checks ignore
small differences in legal exported coordinates; they are not precision oracles.

Seven focused executables passed on each compiler: ULP filters, bounded-ULP
policy, borrowed predicates, construction rounding, borrowed bounds, borrowed
segment/box, and shared boundary. The numerical checks use independent
arbitrary-precision rationals for coordinate limits and sign, including
near-half-cell errors after moving the full N steps; they check signed
determinants and squared distances against their returned absolute bounds.
The public tests include a clipped triangle with distinct exact vertices
that coincide under nearest-double export, and check DE-9IM and the exact
rational output. The mixed owning-point test proves that an uncertified
coordinate head cannot acquire the enclosing policy's guarantee. See
[the filter proof](exact_lean_n_ulp_filters.md) for numerical assumptions.

### Reproduction

Set `EXACT_INCLUDE`, `DEVELOP_INCLUDE` and `BOOST_ROOT` to the respective
include trees. Run the following once per compiler with a fresh output
directory and `CXX` and `TASKSET_CPU` set appropriately:

```sh
BENCH_COMPARE_BOUNDED=1 BENCH_REPETITIONS=5 BENCH_SECONDS=.15 \
BENCH_ARCH=x86-64-v2 DEVELOP_REVISION=082b99181 \
sh test/strategies/run_simple_cartesian_benchmark.sh all "$out"
```

`comparison.csv` contains all 225 matched operation/workload/budget rows
for each compiler; `results.csv` contains the timing samples.
