# Decimal Polygon Overlay and Relate

The existing `prepare_rational_polygon` path is unchanged. Opt in to
certified predicate filters with `prepare_filtered_decimal_polygon`
and `strategies::relate::cartesian_decimal`. Both paths compute rational
geometry from the same exact input decimals.

```cpp
#include <boost/geometry/algorithms/intersection.hpp>
#include <boost/geometry/algorithms/prepare_filtered_decimal_polygon.hpp>
#include <boost/geometry/algorithms/prepare_rational_polygon.hpp>
#include <boost/geometry/algorithms/relation.hpp>
#include <boost/geometry/io/wkt/read.hpp>
#include <boost/geometry/strategies/relate/cartesian.hpp>
#include <boost/geometry/strategies/relate/cartesian_decimal.hpp>

namespace bg = boost::geometry;
using decimal = bg::exact_coordinate<std::string>;
using point = bg::model::point<decimal, 2, bg::cs::cartesian>;
using polygon = bg::model::polygon<point>;

polygon a, b; // Valid closed clockwise WKT in this model.
bg::read_wkt(wkt_a, a);
bg::read_wkt(wkt_b, b);
auto ra = bg::prepare_rational_polygon(a);
auto rb = bg::prepare_rational_polygon(b);
bg::strategies::relate::cartesian<boost::multiprecision::cpp_rational> strategy;
bg::model::multi_polygon<decltype(ra)> output;
bg::intersection(ra, rb, output, strategy);
auto matrix = bg::relation(ra, rb, strategy);

auto fa = bg::prepare_filtered_decimal_polygon(a);
auto fb = bg::prepare_filtered_decimal_polygon(b);
bg::strategies::relate::cartesian_decimal filtered;
decltype(output) filtered_output;  // Still the plain rational polygon type.
bg::intersection(fa, fb, filtered_output, filtered);
auto filtered_matrix = bg::relation(fa, fb, filtered);
```

Use the same prepared inputs with `union_`, `difference`, `sym_difference`,
and mask-based `relate`. The boundary copies a 2D Cartesian polygon or
multipolygon once, preserving holes, ring order, and closure. Keep overlay
output in rational polygons to retain intersections such as 2/3. Prepared
inputs and rational results can be reused without parsing or rounding again.
Supply valid input; preparation neither corrects nor validates it. WKT reading
is supported for closed decimal polygon models; construct open wrapper rings
programmatically. Wrapper WKT writing is not supported.

## Predicate Scope

Unfiltered preparation copies authoritative rationals and does not read the
views. The filtered preparation additionally copies the nearest-even view
and an absolute error bound from each original `exact_coordinate<std::string>`
coordinate. Exactly representable views have zero error; otherwise the
larger adjacent spacing bounds the error (the rational lies between those
neighbors). Any write to a prepared point discards its certificate. Copies
of an unchanged original point retain it. Constructed points and outputs
do not acquire a certificate from rounding.

Only exact zero and finite views of magnitude `[2^-128, 2^128]` are admitted.
In this band, views and stored whole-spacing errors are on the `2^-180`
grid. The degree-two terms, including products of difference errors and
the roundoff allowances, stay normal when nonzero: the smallest scaled
difference term is at least `2^-232`, and even two of those multiplied
and scaled by a roundoff factor stay above `2^-520`. Magnitudes stay
well below overflow. The existing Cartesian forward-error orientation
bound applies under round-to-nearest, no fast-math and no implicit
contraction. On a miss use the rational side predicate. Range and equality
checks use monotonicity of nearest-even: strict comparisons prove order,
while equal views require rational fallback unless exactness proves equality.

## Public Coverage

With both explicit preparation and the decimal strategy, the public
polygon/multipolygon tests exercise `intersection`, `union_`, `difference`,
`sym_difference`, `relation`, and mask-based `relate`. Segment
classification uses filtered range, equality and side where it can.
The winding and overlay side hooks also use certified original operands
when available; otherwise they use rationals. A certain rejection avoids
computing the rational edge differences needed by Cramer's rule.
Intersection construction, segment ratios, envelope and other ordering
work, signed measures and output remain rational.

This is **predicate filtering inside selected public operations**, not a
floating-point overlay engine. The standalone array predicates in
`coordinate_predicates.hpp` do not provide public dispatch on raw wrappers.
There is no automatic strategy selection, fixed-width admission, mixed
wrapper/native overload, wrapper overlay output, or non-areal qualification.
The default Cartesian and native borrowed strategies still do what they did
before; unadmitted exact values can use the same rational engine.

## Reproductions and Timing

See the [historical decimal-intent audit](exact_decimal_issues.md) for other
report fixtures and for unsupported geometry pairs.

- The local #663 triangles overlap in both interpretations. Independent
  slope/height oracles give areas about 7.186e-32 for decimal text and
  4.149e-33 for exact binary64 input. Neither interpretation is empty.
- A decimal gap of 1e-20 between two squares vanishes in binary64: the
  rational inputs are disjoint and their union has two components; the
  binary inputs share an edge and their union has one.
- Both public decimal paths are tested on #663 and the gap, on holes and
  multipolygons, and on rational output at 2/3, including crossing inputs
  beyond binary64 overflow and underflow. Moving and mutating a prepared
  point cannot leave a stale certificate. The tests assume valid input.

From the Geometry checkout:

```sh
CASE_FILTER=decimal_overlay BOOST_ROOT=/path/to/boost \
  sh test/strategies/run_exact_primitives.sh /tmp/bg-decimal-overlay
/tmp/bg-decimal-overlay/decimal_overlay --benchmark
```

The benchmark reports parsing, preparation, prepared intersection, and
prepared relation separately. All three timed paths use the same exact
decimal inputs and rational output. The annotated/unfiltered control runs
the old rational strategy on the *same* annotated points as the filtered
strategy. The final column in each group measures both fast predicates
and deferred rational edge differences, not an isolated side primitive.

Medians of five local Apple Clang 21 arm64 runs, C++17,
`-O1 -DNDEBUG -ffp-contract=off`, microseconds per pair (1,000 repeats per
small case, 60 for the 52-point corridors):

| Pair | Prepare plain / annotated | Intersection plain / annotated-off / filtered | Relation plain / annotated-off / filtered |
| --- | ---: | ---: | ---: |
| Crossing | 0.73 / 0.68 | 24.46 / 24.35 / 18.48 | 16.14 / 16.47 / 10.83 |
| #663 | 2.52 / 2.49 | 46.37 / 46.14 / 40.20 | 24.21 / 24.47 / 19.17 |
| Decimal gap | 1.08 / 0.99 | 8.01 / 8.24 / 8.23 | 6.40 / 6.84 / 6.57 |
| Corridors | 9.28 / 9.17 | 249.85 / 257.53 / 189.23 | 172.11 / 174.30 / 104.26 |

GCC 16 corridor intersection/relation likewise went from 420.97/284.50 to
313.03/163.45 us. The gap did not benefit: ordinary envelope pruning can
already reject it without segment classification. These are local fixtures,
not general performance guarantees.

For routing evidence, compile the same test with the existing
`BOOST_GEOMETRY_EXACT_FILTER_STATISTICS` define and run it **without**
`--benchmark`. The polygon-with-hole test requires side, range and
equality hits separately in each of the four set operations and in
`relation`. For example, public intersection observed 60 side hits,
7 range hits, 64 equality hits and 15 side fallbacks. Do not time the
instrumented build.
