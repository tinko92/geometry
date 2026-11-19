# Compact Exact Comparison Scratch

This fixes the shared-boundary regression in `ffa469495fd0642859329df0ad298d6a6a06c6fa`.
The pre-regression forward-filter checkpoint is
`b518175cf8d7c676bb9fbcea538f3c71fd49595a`.

## Change And Exactness

Tied distance comparisons cannot finish with a nonzero-sign filter. They were
forming worst-case-sized fixed integers, even for small integer coordinates.
Boost.Multiprecision's fixed backend initializes/copies its entire inline array
in several constructors, so allocator-free arithmetic was still expensive.

The new exact tier examines active limb counts and the highest used limb before
copying homogeneous coordinates into 128-bit components. Cached native and
materialized values are narrowed under the existing cache lock, avoiding an
intermediate worst-case point copy. Native inputs are built directly from their
binary significands/exponents; no 2^1074 grid is constructed and canceled first.
Point/cache layouts, forward error filters and the full fallback are unchanged.

Before refining an ambiguous squared-distance comparison, known native integer
coordinates in \([-2^{30},2^{30}]\) also have a direct exact uint64 path. Each
coordinate difference is at most \(2^{31}\), so the sum of two squared magnitudes
is at most \(2^{63}\). Both the integral-value and range checks must pass; no
fractional, large, or approximate implicit coordinate is accepted by rounding.
This removes the remaining small-integer tie overhead without replacing any
general predicate with a tolerance or a guessed answer.

For input components \(x=(-1)^s M2^e\), put
\(E=\min(0,e_x,e_y)\) and construct
\[
 (X:Y:W)=((-1)^{s_x}M_x2^{e_x-E}:
          (-1)^{s_y}M_y2^{e_y-E}:2^{-E}).
\]
This uses integer operations on the original stored bits, not rounded floating
arithmetic. Component bit lengths are checked before assignment or shifting.

For fitting points, \(|X|,|Y|,W<2^{128}\). Each affine difference numerator
has magnitude below \(2^{257}\). A sum of its two squares therefore fits 515
numerator bits, and the squared denominator product fits 512 bits. The existing
checked comparison gives their cross-product difference 1032 bits, including
spare carry bits. The signed-determinant comparison uses the same proven bound.
Thus the compact answer is exact, including zero; it is not another approximate
filter. If any component does not fit, the original full-capacity path runs.
The threshold is a scratch-size choice, not a new precision limit.

No numerical heap storage was added. Precision-cache growth remains allowed;
warmed primitive paths retain their zero-allocation assertions. This changes
only the comparison hot path, not every use of wide fixed integers in the kernel.

## Verification

All 23 registered exact tests pass. Added regressions cover signed and unsigned
limb boundaries, failed narrowing leaving output untouched, exact distance ties,
cached implicit and materialized points, full 128-bit components, and values just
outside the compact tier. Independent rational controls cover unequal results
that round together, mixed inputs and all standard rounding modes. Tests also
exercise both caches concurrently and assert zero warmed numerical allocations.
The native path includes tests at the \(2^{63}\) squared-distance boundary,
signed/unsigned mixtures, fractional rejection and just-outside-range rejection.

ASan/UBSan arithmetic and predicate tests and ThreadSanitizer predicate/cache
tests pass. The 164-record replay (82 fixtures, 47 selected issues) has no failures
or changed spatial results. The 16 raw WKT differences versus the original
frozen evaluator are the previously verified #888 rotations/component order.
The source headers remained unchanged through the replay.
Direct comparison against the hash-verified `ffa469495` replay also finds zero
differences across all 164 records, including raw WKT, not just set equality.
[Evidence](exact_compact_previous_equivalence.json) records both result hashes;
the [independent set check](exact_compact_issue_comparison.json) retains the
original evaluator's validation and permutation evidence.

The broader Cartesian suites also preserve all baseline results:

| Suite | Passing Sources | Unchanged Runtime Failures | Existing Compile Failures |
|---|---:|---:|---:|
| Exact-enabled | 34 | 7 | 3 |
| Default | 37 | 4 | 3 |

No new failures or changed failure logs. Every classification is retained in the
[exact](exact_compact_exact_regressions.tsv) and
[default](exact_compact_default_regressions.tsv) lists. Headers were unchanged
through both runs. Existing failures have not been relabeled as passes.

## Benchmark Method

The full matrix uses the unchanged five workloads and nine operations, for double
and int64, with six serial, case-interleaved rounds at 0.2 seconds per case. Previous
exact is `ffa469495`; develop is `aaf0560e520041a833c11ad5d61017b684b46f9a`.
Apple M5, Apple Clang 21, C++14, `-O3 -DNDEBUG -ffp-contract=fast -pthread`.
Input generation is excluded; complete output allocation/conversion/destruction
is included. Allocation measurements are separate. Builds and correctness jobs
are stopped during the recorded timing run. Other desktop activity is uncontrolled.

Whole-binary batches at 0.3 seconds per case were also run and are preserved in
`/tmp/bg-compact-performance/batch-run`. Those batches were too noisy for a reliable
aggregate conclusion: the new double build's geometric-mean maximum/minimum ratio
was 2.02, and even unchanged point-in-polygon calls appeared over 20% slower.
The final run alternates variants after each case, rather than after an entire
45-case binary, to reduce that timing drift. The C++ benchmark and inputs are
unchanged; the wrapper's `paired` mode only changes invocation order.

Summaries use the upper median of six samples. Aggregate ratios are geometric
means across cases, not estimates for an application mix. The unchanged benchmark
uses integer-valued double inputs; it is not a fractional-input performance survey.
Output digests are checked against previous exact but do not replace the separate
spatial correctness replay. Local ARM results do not substitute for a new x86 run.

## Shared-Boundary Comparison

Separate six-round measurements alternate the new binary, `ffa469495` and the
pre-regression `b518175cf` binaries, for both coordinate types and both operations,
at 0.5 seconds per case. The source, compiler and flags match. No build, test or
allocation-counting process ran concurrently. All 72 output digests match their
corresponding reference cases.

Upper-median microseconds per call:

| Coordinates / Operation | Before Regression | Last Version | Compact Scratch |
|---|---:|---:|---:|
| Double union | 330.29 | 866.92 | 196.64 |
| Double intersection | 330.44 | 862.67 | 186.24 |
| Int64 union | 132.05 | 171.26 | 60.85 |
| Int64 intersection | 129.04 | 171.68 | 60.63 |

The shared-boundary regression is removed for both types, including against the
older checkpoint, not only against the immediately preceding slower version.
The [targeted results](exact_compact_shared_results.json) retain minima and maxima;
desktop timing noise remains visible, so the last decimal is not a guarantee.

## Full Matrix And Profile

The [CSV](exact_compact_benchmark_results.csv) contains 405 case summaries:
1620 timed records and 135 separate allocation records. Exact output digests
match the previous version in every case. The
[manifest](exact_compact_benchmark_manifest.json) pins source/header/runner hashes,
compiler, flags and revisions.

Geometric means over the indicated cases:

| Group | Double Speedup vs Last Exact | Int64 Speedup vs Last Exact | Double Slowdown vs Develop | Int64 Slowdown vs Develop |
|---|---:|---:|---:|---:|
| Four polygon set operations | 1.17x | 1.11x | 6.17x | 4.82x |
| Three relation matrices | 1.03x | 1.04x | 4.77x | 4.47x |
| Linear/areal intersection | 1.01x | 0.95x | 5.32x | 4.73x |
| Point/polygon covered_by | 1.01x | 0.99x | 1.90x | 1.33x |
| All 45 cases | 1.09x | 1.05x | 4.89x | 4.07x |

CPU-time aggregate speedups are 1.09x and 1.06x. Case-level interleaving did not
eliminate desktop noise: geometric-mean per-case maximum/minimum ratios remain
2.02 for double and 2.03 for int64. These modest aggregate changes are estimates,
not evidence of uniform speedups. Apparent regressions, including in unchanged
point-in-polygon paths, remain in the full CSV rather than being discarded.
Small differences outside the shared-boundary cases cannot reliably be attributed
to this patch from these local measurements. Develop's double star symmetric
difference also has different output, as recorded in the prior report.

Whole-call allocation counts/bytes are unchanged in all 45 double cases. For
example, shared union uses 566 allocations / 130043 bytes in both exact builds;
convex union uses 89 / 96299. Those include algorithm and output containers.
The speedup comes from avoiding numerical scratch copying/clearing, not from
changing these heap counts. Warmed primitive allocation assertions remain zero.

A five-second sample of the final shared-boundary double union has 4119 samples.
Visible exclusive memmove/memset frames account for 18 samples, about 0.4%,
versus 53.6% before this patch. Frames with fewer than five samples are omitted
by the profiler's flat listing, so the meaningful conclusion is below 1%, not
the last decimal. Remaining hot frames are exact orientation (1056), nextafter
(533), interval subtraction (290), input validation (199) and native orientation
(175). These independent frame counts are not added to inclusive call-tree costs.
The raw profile is `/tmp/bg-compact-performance/profile-shared.txt`.

The shared-boundary regression is resolved, but the exact strategy remains
substantially slower than develop overall. This patch is not a claim of default
strategy performance readiness.

## Reproduce

```sh
export BOOST_ROOT=/Users/tinko/dev/boost
export BENCH_SECONDS=0.2 BENCH_REPETITIONS=6
node test/strategies/run_exact_benchmark.js build /tmp/bg-compact-performance-new
node test/strategies/run_exact_benchmark.js paired /tmp/bg-compact-performance-new
node test/strategies/run_exact_benchmark.js allocations /tmp/bg-compact-performance-new
node test/strategies/run_exact_benchmark.js summary /tmp/bg-compact-performance-new
```

Raw measurements and CPU profiles are in `/tmp/bg-compact-performance`.
Do not run compilation/tests during timing. Header and source hashes are checked
before and after measurement; the manifest identifies the pre-commit dirty headers
by content hash, not just its recorded HEAD.
