# Operation-Owned Exact Caches

Historical checkpoint report. The in-progress static-kernel redesign uses one
compile-time coordinate type per context, preserves copied cache indices within
an operation, and does not support implicit values retained between independent
operations. The retained-value and mixed-format descriptions and benchmark
results below describe the earlier revision, not the current working tree.

This replaces the cache ownership in `ea8cf307877cbcc7b3bfc445f52b498502e1dc36`.
There is no internal parallelism in these algorithms. The former atomics and
locks supported extra sharing of internal exact-point objects, not an algorithm
requirement. They and the shared materialized-value store have been removed.

## Ownership And Bounds

The root operation owns typed vectors of precision entries; nested calls reuse
that context. The vectors and index-range bookkeeping are initialized only when
the operation first needs a cache entry. A successful native predicate filter
does not initialize them. Each native cache handle is one plain `size_t` index. Copies start
unbound; moves preserve the hint. Issued-index ranges prevent a retained native
point from accidentally using an unrelated entry in a later operation. Only
the current-context pointer, index counters, and optional diagnostic counters
are thread-local; none owns numerical storage.

For native coordinate bounds \(B,S\), the intersection record needs
\(N=3B+3\) numerator bits and \(D=2B+S+3\) denominator bits. A coordinate
splice needs \(N+D,2D\) bits; a midpoint needs \(N+D+1,2D+1\). These are
the auxiliary constructions used by the native-input algorithms.

| Coordinates | Native Numerator/Denominator Bits | Auxiliary Bits |
|---|---:|---:|
| int16 | 51 / 35 | 87 / 71 |
| int64 | 195 / 131 | 327 / 263 |
| binary32 | 834 / 706 | 1541 / 1413 |
| binary64 | 6297 / 5273 | 11571 / 10547 |

Predicate scratch is bounded separately from its polynomial expression, as
before. Integer arithmetic remains checked; narrowing checks active bits before
copying. There is no fixed 32768-bit compatibility record, reference count, free
list, or process-global numerical vector. Dynamic numerical storage is confined
to operation-owned vector growth. The warmed primitive allocation tests reserve
only those vectors before counting, not numerical temporaries.

Ordinary mixed native inputs retain their original bits inline. Native-support
intersection points remain self-contained and usable after an operation ends.
Auxiliary values and bounded rational compatibility literals require an enclosing
scope when used outside the public exact algorithm entry points. They cannot
outlive it or migrate between independent operations. Unbounded arbitrary
rational inputs and unlimited recursive auxiliary constructions are not part of
the native-coordinate contract.

Legacy setters preserve exactly representable native input coordinates inline.
This matters for box-corner construction: native box clipping results must not
accidentally acquire an auxiliary-cache lifetime. True constructed fractional
values are not rounded into that path.

## Correctness

All 23 exact test executables pass, including new capacity assertions, full-range
int64 splice/midpoint controls, native setters without a scope, nested scopes,
stale handles, retained native intersections, and independent concurrent callers.
ASan/UBSan point and predicate tests and ThreadSanitizer point/cache tests pass.
The public exact entry points for within, covered_by, envelope, validation and
correction now establish the same strategy-selected scope as overlay and relate;
ordinary strategies get a no-op scope.

The frozen replay has 164 records, 82 fixtures, 47 selected issues, and no failures.
Every record, including raw WKT and validator diagnostics, is identical to the
previous exact checkpoint ([hashes and comparison](exact_operation_previous_equivalence.json)).
Against the original evaluator, the only spatially
equivalent raw WKT changes remain the already documented #888 permutations.
The runner differs from its frozen source only by an enclosing cache scope for
legacy validation of internal exact points; removing that exact block reproduces
the original source hash. Fixtures and the independent evaluator are unchanged.

| Broader Suite | Passing Sources | Unchanged Runtime Failures | Existing Compile Failures |
|---|---:|---:|---:|
| Exact-enabled | 34 | 7 | 3 |
| Default | 37 | 4 | 3 |

All 88 source classifications match the preceding checkpoint. Existing failures
are not passes. The final exact-enabled run and issue replay use unchanged
headers throughout. See the [exact](exact_operation_exact_regressions.tsv),
[default](exact_operation_default_regressions.tsv), and
[issue](exact_operation_issue_comparison.json) evidence.
Raw logs remain under `/tmp/bg-operation-*`.

## Benchmarks

Apple M5, Apple Clang 21, C++14, `-O3 -DNDEBUG -ffp-contract=fast -pthread`.
Previous exact is `ea8cf3078`; develop is
`aaf0560e520041a833c11ad5d61017b684b46f9a`. The unchanged corpus has five
workloads and nine operations for double and int64. Six serial, case-interleaved
rounds use 0.2 seconds per case, with no concurrent builds or correctness tests.
The 1620 timed records and 135 separate allocation records produce
[405 summaries](exact_operation_benchmark_results.csv). All exact output digests
match the previous version. Inputs are integer-valued even for double; this is
not a fractional-input performance survey or a new x86 run.

Geometric-mean speedup against previous exact, using upper-median timings:

| Calls | Double | Int64 |
|---|---:|---:|
| Polygon set operations | 1.02x | 1.07x |
| Relation matrices | 1.04x | 1.10x |
| Linear/areal intersection | 1.03x | 1.11x |
| Point/polygon covered_by | 1.02x | 0.98x |
| All 45 cases | 1.02x | 1.07x |

Overall slowdowns against develop remain 4.73x for double and 3.71x for int64.
CPU-time ratios closely match wall-time ratios. Geometric-mean per-case
maximum/minimum ratios are 1.032 and 1.050; small differences, particularly the
2% int64 point-query difference, are not strong evidence of a speed change.
An initial eager-context version slowed those point queries by about 22%.
Lazy initialization removed most of that overhead without weakening predicates.

Whole-call allocation counts are unchanged in 42 of 45 double cases. Vector
growth is now paid inside each operation rather than amortized by retained
thread-local capacity. Convex symmetric difference adds 3 allocations, holes
symmetric difference adds 3, and star symmetric difference adds 16; the latter's
cumulative allocated bytes rise from 3,168,177 to 6,627,329 per call. These are
allocation-traffic figures, not peak live memory. Other cases generally allocate
fewer bytes with the smaller point objects, e.g. convex union 96,299 to 92,283,
with 89 allocations in both versions. No per-point heap allocation was added.
The [comparison](exact_operation_benchmark_comparison.json) retains every case,
including these increases. The [manifest](exact_operation_benchmark_manifest.json)
and [run record](exact_operation_benchmark_run.json) pin the source, headers and
timing procedure. This is a simpler ownership model with modest throughput gains,
not a claim that exact strategies have reached develop's performance.

```sh
export BOOST_ROOT=/Users/tinko/dev/boost BENCH_SECONDS=0.2 BENCH_REPETITIONS=6
node test/strategies/run_exact_benchmark.js build /tmp/bg-operation-bench-new
node test/strategies/run_exact_benchmark.js paired /tmp/bg-operation-bench-new
node test/strategies/run_exact_benchmark.js allocations /tmp/bg-operation-bench-new
node test/strategies/run_exact_benchmark.js summary /tmp/bg-operation-bench-new
```

Final raw timings are in `/tmp/bg-operation-lazy-performance`.
