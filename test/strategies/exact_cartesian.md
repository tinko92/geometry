# Exact Cartesian Relate And Overlay

The opt-in umbrella is
`boost::geometry::strategies::relate::cartesian_exact<Coordinate>` in
`boost/geometry/strategies/relate/cartesian_exact.hpp`. Unconfigured defaults
are unchanged; see the [build-wide pilot](exact_default_pilot.md) for the
separate opt-in on a few public calls. This work starts at PR1491 commit
`72b2d832fb9596a2a23c0c164a5426ae720d0cd6`.

```cpp
#include <boost/geometry/strategies/relate/cartesian_exact.hpp>

boost::geometry::strategies::relate::cartesian_exact<double> strategy;
boost::geometry::intersection(a, b, output, strategy);
auto matrix = boost::geometry::relation(a, b, strategy);
```

## Mathematical Contract

Inputs are two-dimensional Cartesian geometries with finite integral, float or
double coordinates. Floating inputs mean their exact binary values, not an
intended decimal value or intended collinearity. Mixed native input types are
not narrowed to the umbrella's `Coordinate` type. Exact rational literals can
also enter predicates and working geometries. Long double is not supported.
Floating-point overflow and underflow are outside the promised domain of the
filters; integer overflow is not. Conservative exponent guards may still select
the exact fallback rather than attempt an uncertified floating expansion.

Constructed intersections remain exact until the complete operation finishes,
including intermediate geometry-collection merges and subtractions. Final
integral and floating outputs use nearest, ties-to-even rounding. Distinct
vertices may collapse, or a rounded output may become invalid. Such an output
must not be fed back as if it were the unrounded exact set. Rational output
coordinates avoid that final representation loss.

This does not change existing dimensional output conventions. For example,
line/box clipping still emits a singleton linestring at a tangent contact,
as required by the existing box test. The resulting validity discrepancy in
issue #1448 is not a floating-point rounding problem.

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
- Native orientation starts with an analytic determinant error bound, then uses
  outward intervals and error-free expansions. Unresolved cases use exact
  arithmetic. Native integer comparisons and checked arithmetic stay inline.
- Native intersections own the coordinates of their two supporting segments,
  approximate coordinates and error bounds. Operation-scoped indexed caches
  hold fixed-size expansions and bounded integer homogeneous coordinates.
  Equality, ordering and orientation have exact homogeneous fallbacks.
- Constructed-point filters first use coordinate error bounds, then refine those
  bounds with certified expansion evaluation of the producing line equations.
  Exact polynomial expansion comparisons and producer-incidence certificates
  resolve further ambiguities without constructing rational coordinates. Only
  unresolved cases materialize homogeneous coordinates or rational values.
- Final coordinate conversion also verifies the nearest representable value
  using expansion comparisons against its rounding cell, including ties-to-even.
- Envelope coordinates, comparable distances, ring measures and segment
  parameters use explicit forward-error filters. Native scalar arithmetic
  retains exact error-free residuals inline. Closed primitive sources own their
  operands without a runtime expression tree; ambiguous native values use
  fixed-size expansions or checked bounded integers. Arbitrary rational literals
  and unbounded global rational sums retain a compatibility fallback.

Certification requires IEEE binary arithmetic and no fast-math/reassociation.
Floating overflow/underflow and flush-to-zero modes are outside this contract.
Basic-operation rounding barriers permit
contraction-enabled builds. Exact floating expansions require round-to-nearest;
other standard rounding modes take the exact fallback. Caches use
`std::call_once`; older Unix toolchains can require native thread link flags.

This is not yet a proposal to make the strategy the default. Allocation,
coordinate materialization and conservative filters still need performance
work and broader compiler/platform verification.

## Verification

Twenty-four registered `exact_*` strategy tests cover arithmetic, implicit points,
ratios, exhaustive small-grid segment classification, set operations, relate,
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
  boost_geometry_strategies_exact_ratio \
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
