# Cold tagged orientation

This ports the bounded kernel and cold dispatch from `357434d0a` to
`0f7ba7837`. `BOOST_GEOMETRY_EXACT_TAGGED_ORIENTATION=0` (the default)
keeps the fixed path. `=1` enables the tagged route. Use a consistent
setting across all translation units.

The route requires three borrowed proper-crossing records in a binary64
Root, after visitation, both sign filters, and the existing alias and
support exits. If any fixed homogeneous value is ready, keep the fixed
path. Otherwise compute the sign from the original supports with the
[bounded stack kernel](exact_tagged_lines.md). The readiness lookup is
read-only; neither it nor the tagged computation inserts cache entries.
Existing expansion/coordinate refinement may still populate its normal
cache before dispatch. No record layout or other predicate is changed.

## A/B and verification

Compile `test/strategies/exact_tagged_dispatch.cpp` at settings 0 and 1,
with identical C++17 optimization, `-ffp-contract=off -pthread`,
`-Iinclude -I$BOOST_ROOT`, and `-DBG_TAGGED_TIMINGS_ONLY` for timing.
`--benchmark` runs the original 30 cells (six fixtures times five
operations), including `near-parallel` intersection. `--check` prints
the same native WKT and DE-9IM snapshots without timing. Omit statistics
from timing builds. Add `-DBOOST_GEOMETRY_EXACT_FILTER_STATISTICS` in
separate builds for `--traffic`, which reports tagged calls per cell.

Without `BG_TAGGED_TIMINGS_ONLY`, running with no arguments compares
rational overlay sets, areas, validity, and DE-9IM on all fixtures. With
statistics, it also verifies the cold route without allocating or filling
per-point cache entries in the tagged stage, and fixed bypass for warm,
aliased, native and query cases. `exact_tagged_lines.cpp` checks the
kernel against independent rational line intersection and orientation.
