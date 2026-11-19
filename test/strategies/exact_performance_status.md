# Borrowed exact strategy: performance follow-up

Local, uncommitted work; default strategies are unchanged. This supersedes the
timings, not the contracts, in `exact_borrowed_status.md`.

## Changes

Profiles traced the shared-boundary regression to full-capacity homogeneous
arithmetic and scratch copying after inconclusive distance filters. Restore
checked small-integer distance evaluation and checked compact homogeneous
fallbacks before the unchanged full-capacity fallback. Coordinate filters read
only their requested axis. Identical original support pairs certify identical
intersections; this includes operation/stage provenance checks.

Ordinary construction initialization reuses its existing forward-error bound.
Private midpoint/splice queries retain exactly representable native results
without precision-cache entries and otherwise try checked compact arithmetic.
Borrowed inputs, original supports, point layouts, exact fallbacks, and the
stored-approximation output contract remain intact. No release rounding-mode
checks or numerical ownership pools were introduced.

## Qualification

All 57 focused sources pass, including the corrected `borrowed_operations`
rerun. Its original native-area checks accidentally required ordinary rounded
output to preserve exact area. The replacement checks retained internal results
against analytic rational geometries; nearest-output area checks remain. The
original failed log is preserved. Corrected operations also pass ASan/UBSan.

The new shared-boundary regression passes for double/int64 and all three output
policies, with 10,908 small-integer comparisons and no full measure comparisons
per combination. Predicate tests include conflicting local support IDs from
distinct stages. Sanitizer, allocation, and query-capacity tests pass.

Default C++14 qualification: 37 passes, four baseline-identical runtime failures,
three matching baseline compile failures; no new regression. Fresh grid tests:
99,000 calls and 18,000 input roundtrips, no errors. Issue replay: 656 records,
82 fixtures, 47 issues, no changed classifications or internal-exact observations.
Four ordinary native outputs differ because skipped refinements can change
rounding (#588, #852, #1235, #1236); retained and nearest results are unchanged.
The independent mixed-collection relate limitation remains.

## Timings

Apple M5, Clang 21, C++17, `-O3 -DNDEBUG -ffp-contract=off`. Five serial rotated
repetitions, 45 workload/operation combinations per coordinate corpus; every
timed batch traverses complete eight-input corpora. Ratios are geometric means
of matched median times; greater than one means slower.

| Inputs | vs. Regressed Borrowed | vs. Pre-Borrowing Checkpoint | vs. Develop |
| --- | ---: | ---: | ---: |
| Whole-number double | 0.687 | 1.210 | 1.890 |
| int64 | 0.785 | 1.144 | 2.542 |
| Fractional double | 0.699 | 1.138 | 1.886 |

Shared-boundary double union falls from 720.22 to 80.62 microseconds/call;
intersection from 652.10 to 70.80. The corresponding checkpoint times are
73.48 and 61.58 microseconds. There is still measurable overhead against the
checkpoint, especially some relate and artificial star cases.

All exact variants preserve structural and quantized spatial control digests.
The historical develop star symmetric-difference discrepancy and integer-output
conversion discrepancies remain diagnostics, not independent correctness oracles.
All 4,590 raw timing/allocation rows have positive complete-corpus call counts.

Develop is pinned to `47286e8edccdeb1eb74a5182f7570b7feb67e0cf`. The issue replay
still pins historical develop `aaf0560e5`; its issue counts must not be attributed
to this newer revision. Raw evidence and manifests: `/tmp/bg-exact-perf-20261116`.

The final shared-union profile reduces copy/clear samples from 46.1% to 3.4%;
full-width arithmetic leaves the hot-path groups. Allocation is now 18.7% and
turn generation 34.6%. These are sample shares, not additive explanations of
wall-time differences. Star symmetric difference still spends 16.5% copying and
13.8% in fixed arithmetic; it remains a useful stress case, not the intended
representative weighting for the next compact benchmark.
