# Exact Cartesian Relate And Overlay

The opt-in umbrella is
`boost::geometry::strategies::relate::cartesian_exact<Coordinate>` in
`boost/geometry/strategies/relate/cartesian_exact.hpp`. Unconfigured defaults
are unchanged; see the [build-wide pilot](exact_default_pilot.md) for the
separate opt-in on a few public calls. The original prototype used PR 1491;
see the [current integration order](exact_integration_order.md) for the develop + PR 1517 series.

```cpp
#include <boost/geometry/strategies/relate/cartesian_exact.hpp>

boost::geometry::strategies::relate::cartesian_exact<double> strategy;
boost::geometry::intersection(a, b, output, strategy);
auto matrix = boost::geometry::relation(a, b, strategy);
```

## Mathematical Contract

Inputs are two-dimensional Cartesian geometries with finite integral, float or
double coordinates. Floating inputs mean their exact binary values, not an
intended decimal value or intended collinearity. Mixed native promotion must
preserve the input types; lossy type combinations and long double are outside
the interface. Decimal intent is a separate authoritative-scalar strategy.

Original supports and implicit intersections stay available until the operation
has finished, including intermediate set operations. The default emits the
stored approximation without further refinement at output. Alternative
policies specify a bounded number of representable binary64 steps or correct
rounding. Integral output uses ordinary BG conversion of the stored value.
Distinct vertices can coincide after export under any policy; do not feed a
rounded polygon back as though it represented the exact set.

For an eager binary64 step guarantee, choose
`strategies::relate::bounded_ulp_output<N>` (stored coordinate errors) or
`strategies::relate::lean_bounded_ulp_output<N>` (no coordinate error fields
in the root-bound construction). Pass the policy as the second template
argument to `cartesian_exact`. Both count steps from the correctly rounded
coordinate, use the same construction certificate and refine an uncertified
axis. Their predicate schedules can differ, so their admissible native output
coordinates need not be bit-identical.

The lean first stage compares coordinate ranks and uses an N-dependent
forward bound for side, direction, signed determinant and squared distance.
Equality and along-segment order use coordinate comparisons. A search envelope
still needs outward bounds, and a derived measure still needs its own error
bound; neither is an N-step coordinate. Midpoints are separately certified,
and a mixed owning point cannot borrow the enclosing operation's certificate.
See the [filter proof](exact_lean_n_ulp_filters.md).

Point-only contact is not a linear exact output. Legacy line/box calls retain
their established singleton convention; the two routes have separate tests.

## Algorithm And Strategy Boundary

Strategy hooks now own cluster equality and colocation, sweep ordering/windows,
direction/position predicates, segment-position equality, distance measures,
section predicates, clipping, box intersection, signed ring measures and their
comparisons, and envelope expansion. Legacy adapters preserve the old
predicates and tolerances when an umbrella does not supply an override.

Box/box touches and overlaps use a strategy-owned topology helper. Matching
native coordinate types keep their direct comparisons; mixed types are lifted
without narrowing by the exact override. Unary self-touches and self-intersects
expose only default-strategy APIs and are outside this opt-in scope. Optimized
`equals` and `is_valid` paths are also outside the relate/overlay umbrella scope.

Strategy services select intersection points, segment parameters, intersection
result containers, comparable-distance values, boxes and working/output geometry
types. Ring metadata reuses the strategy-selected signed-measure type; it does
not infer a new calculation type from an input point. Signed-measure zero is
not advertised as a general degeneracy test for arbitrary self-intersecting
rings.

Direct coordinate decisions have been removed from the overlay, relate, turn
sorting and sectionalization paths reviewed here. Coordinate access remains in
diagnostics and in shared representation utilities such as geometry conversion;
this is not a claim that every algorithm in Boost.Geometry is now an opaque
kernel client. Straight-segment turn policies still embody straight-segment
combinatorics. A strategy-selected result container alone does not make the
overlay engine support circular arcs.

## Arithmetic And Constructions

- Integral orientation uses a compile-time widening proof for small types,
  checked wider arithmetic where available, and bounded integer fallback.
- Native binary64 orientation uses filtered signs and exact fallback. In a
  certified same-sign box the Sterbenz dispatch knows the coordinate differences
  are exact; unequal rounded products then decide the sign. Ties do not.
- A constructed point borrows its original supports. Per-operation, indexed
  precision caches store reusable expansions or bounded integer homogeneous
  coordinates on demand. No per-point locks or process-wide precision store.
- Constructed predicates use forward-error filters with exact fallback on the
  implicit homogeneous values. Exact support incidence is distinct from
  rounding coincidence. Ring parenting uses topology, not numerical area.
- See the [borrowed representation](exact_borrowed_status.md) and
  [capacity proofs](exact_borrowed_capacity.md) for the numerical layouts.
  The certified algebraic buffer backend has a separate arithmetic contract.

Floating filters assume IEEE round-to-nearest/ties-to-even, gradual underflow,
no excess precision, no reassociation and no implicit contraction. Use strict
FP flags (e.g. `-fno-fast-math -ffp-contract=off`). An unsupported runtime
rounding mode is a precondition violation, not a request for dynamic dispatch.
Do not infer an all-exponent filter guarantee from the bounded integer fallback.

This is not yet a proposal to make the strategy the default. Allocation,
coordinate materialization and conservative filters still need performance
work and broader compiler/platform verification.

## Verification

The registered `exact_*` strategy tests cover arithmetic, implicit points,
segment ordering, exhaustive small-grid segment classification, set operations, relate,
geometry collections and issue fixtures. They include integer ranges above
2^53, mixed coordinates, exact one-third intersections, a genuine hole of area
2^-107, delayed collection rounding, partial box output on failure, and box
arithmetic before narrowing. The segment test checks 13,622 pairs against a
rational-coordinate reference. Sanitizer and concurrency checks supplement the
registered tests.

The latest representation, allocation boundaries and reproducible measurements
are described in [the forward-error report](exact_forward_performance.md).

With a configured Geometry CMake test build:

```sh
cmake --build "$build" --target \
  boost_geometry_strategies_exact_arithmetic \
  boost_geometry_strategies_exact_point \
  boost_geometry_strategies_exact_segment_position \
  boost_geometry_strategies_exact_segment \
  boost_geometry_strategies_exact_relate_overlay \
  boost_geometry_strategies_exact_gc \
  boost_geometry_strategies_exact_issues \
  boost_geometry_strategies_exact_additional_issues \
  boost_geometry_strategies_exact_additional_issue_api -j2
ctest --test-dir "$build" -R '^boost_geometry_strategies_exact_' --output-on-failure
```

`BOOST_GEOMETRY_TEST_EXACT_CARTESIAN` adds an exact pass to selected existing
relate/set-operation helpers without replacing their default pass or changing
their expected results. Some old expectations assume epsilon-based topology;
some validate a rounded result rather than the exact working result. Their raw
assertion counts must not be reported as a clean pass. Conversely, rational
reference agreement shares the generic overlay engine and is not an independent
proof of that engine. `run_exact_regressions.sh` checks default-strategy
compatibility against the original headers.

### Earlier Verification Snapshot

Before the adaptive scalar/construction filters, local Apple Clang, C++14,
2026-11-11 (later regression measurements supersede these counts):

- All nine registered exact tests pass after the final predicate changes;
  source hashes are unchanged during that build and test run.
- Forty-one existing source units execute the exact strategy. All eight public
  predicates in the custom/read-only-model source compile and execute without
  substitutions. This is not a claim that all original assertions pass under
  the new mathematical contract: the observed numerical-oracle/representation
  differences and invalid-input cases were individually diagnosed, with no
  unclassified observed failures remaining.
- Default-strategy compatibility: 42 of 44 Cartesian sources pass. Two legacy
  source-build problems reproduce against the original headers: a removed
  intersection header and helper redefinitions in the linear aggregate source.
  Seven non-Cartesian smokes and all 22 disjoint/touches/overlaps units pass.
- The last bounded release grid runs complete 3,708 set-operation calls on a
  3x3 six-direction grid and 3,159 relation calls on a 4x4 eight-direction grid,
  with zero input errors, native errors or mismatches. Buffer is not exercised.
- Arithmetic/point/segment sanitizer checks pass, including the checked integer
  path without native 128-bit arithmetic and contraction-enabled builds.

Detailed local verification artifacts are in
`/tmp/bg-exact-additional-cmake/exact-final-verification`,
`/tmp/bg-exact-final-post-box`, `/tmp/bg-exact-routing`, and
`/tmp/boost-relate-hooks.0aUF6K`. These are development records, not additional
committed test expectations. Jam targets are registered; this verification used
CMake and direct Clang builds, not a Jam test run.

The [issue audit](exact_cartesian_issues.md) records which complete or selected
report inputs were replayed, unavailable data, and distinctions between exact
topology and final output representation. No GitHub issue is declared fixed
solely because a selected fixture passes.
