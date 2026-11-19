# Decimal Polygon Overlay and Relate

```cpp
#include <boost/geometry/algorithms/intersection.hpp>
#include <boost/geometry/algorithms/prepare_rational_polygon.hpp>
#include <boost/geometry/algorithms/relation.hpp>
#include <boost/geometry/io/wkt/read.hpp>
#include <boost/geometry/strategies/relate/cartesian.hpp>

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
```

Use the same prepared inputs with `union_`, `difference`, `sym_difference`,
and mask-based `relate`. The boundary copies a 2D Cartesian polygon or
multipolygon once, preserving holes, ring order, and closure. Keep overlay
output in rational polygons to retain intersections such as 2/3. Prepared
inputs and rational results can be reused without parsing or rounding again.
Supply valid input; preparation neither corrects nor validates it. WKT reading
is supported for closed decimal polygon models; construct open wrapper rings
programmatically. Wrapper WKT writing is not supported.

## What Is Filtered

**No floating-point numerical filters run on this decimal polygon path.**
Preparation converts from the authoritative scalars and ignores their cached
binary64 views. The generic Cartesian engine performs rational comparisons,
predicates, and constructions on the prepared geometry. Rational segment
ratios keep their exact quotient for ordering and equality, instead of
dividing a separately rounded numerator and denominator. Spatial pruning in
the generic engine still applies; it is not a binary64 sign filter.

The standalone wrapper filters remain available for coordinate comparison
and 2D orientation, but are not called by these operations. The native
`cartesian_exact<double>` borrowed kernel has its own certified filters and
fixed-width fallbacks; it still interprets inputs as exact binary values.
There is no automatic wrapper dispatch, fixed-width admission, mixed
wrapper/native overload, wrapper overlay output, or non-areal support here.

## Reproductions and Timing

- The local #663 triangles overlap in both interpretations. Independent
  slope/height oracles give areas about 7.186e-32 for decimal text and
  4.149e-33 for exact binary64 input. Neither interpretation is empty.
- A decimal gap of 1e-20 between two squares vanishes in binary64: the
  rational inputs are disjoint and their union has two components; the
  binary inputs share an edge and their union has one.
- Tests exercise holes, multipolygons, rational output at 2/3, and crossing
  inputs scaled beyond binary64 overflow and underflow. They are not a
  qualification of arbitrary invalid input or of every BG algorithm.

From the Geometry checkout:

```sh
CASE_FILTER=decimal_overlay BOOST_ROOT=/path/to/boost \
  sh test/strategies/run_exact_primitives.sh /tmp/bg-decimal-overlay
/tmp/bg-decimal-overlay/decimal_overlay --benchmark
```

The benchmark reports parsing, preparation, prepared intersection, and
prepared relation separately, in microseconds per polygon pair. Both
intersection paths emit rational polygons. Example local Apple Clang 21
arm64 run, C++17, `-O1 -DNDEBUG -ffp-contract=off`, 1,000 repeats per row:

| Pair | Parse decimal / binary | Prepare rational | Intersection rational / borrowed binary | Relation rational / borrowed binary |
| --- | ---: | ---: | ---: | ---: |
| Crossing | 3.27 / 2.27 | 0.83 | 31.00 / 10.11 | 18.61 / 0.85 |
| #663 | 6.42 / 4.53 | 2.92 | 52.19 / 36.64 | 27.92 / 13.22 |

These are two small fixtures, not a general throughput comparison. In
particular the two paths can perform different work because their input
sets differ.
