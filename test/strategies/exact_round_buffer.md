# Exact round-buffer kernel

`strategies/cartesian/exact/round_buffer.hpp` is an internal Cartesian 2D
building block, not a `buffer` strategy. It accepts `exact_number` (unbounded
rational) inputs. Each call performs a bounded number of rational operations;
integer storage and arithmetic costs depend on operand size. There is no
subdivision loop, floating arithmetic, root approximation, or tolerance-based
sign decision.

## Constructions and certificates

- For a nondegenerate directed edge `a -> b`, let `n` be its raw left or
  right normal and `D = n.n > 0`. The offset endpoints are `a + r*n/sqrt(D)`
  and `b + r*n/sqrt(D)`; the support is `n.(p-a) = r*sqrt(D)`. Both endpoints
  and the line use the same symbolic radical, including when `D` is a rational
  square. `side_of_line` classifies a rational query against the infinite
  support, not against the finite segment. Endpoint coordinate comparisons
  similarly reduce to comparing `A` with `B*sqrt(D)`. Different signs decide
  immediately; for like signs compare squares and reverse the order if both
  are negative. This preserves exact equality without evaluating a root.
- `rational_arc_sample(q,t)` rotates `((1-t^2)/(1+t^2), 2*t/(1+t^2))` in
  quarter turns. The vector has squared norm one, and the four closed charts
  cover the unit circle in counterclockwise order with shared axis endpoints.
  A rational center and radius turn its coordinates into a rational point
  by translation and scaling. The charts do not locate or order joins.
- For unit directions `u,v` at a common center and nonnegative radius `r`,
  the chord's nearest point to the center is its midpoint: the midpoint
  vector is perpendicular to the chord and has squared length
  `r^2*(1+u.v)/2`. The maximum radial deficit is therefore bounded by `e`
  exactly when `r^2*(1+u.v)/2 >= (r-e)^2` for `0 <= e < r`; `e >= r` always
  suffices. The dot product needs only the square root of the product of
  the two squared lengths. This includes rational samples, one rational
  and one irrational direction, and two irrational directions, with no
  rounded endpoints. Equality is accepted. The predicate covers the chord
  of the minor arc, including either semicircle at antipodes; it does not
  certify that a caller chose the intended arc or avoided a major sweep.

## Public Point Proof of Concept

Opt in with `<boost/geometry/algorithms/buffer_point_certified.hpp>`:

```cpp
namespace bg = boost::geometry;
using point = bg::model::point<double, 2, bg::cs::cartesian>;
bg::model::polygon<point> output;
auto status = bg::buffer_point_certified(point(0, 0), output, 2, 0.01);
```

Success certifies **strictly less** than the requested Euclidean Hausdorff
distance, both between the output boundary and the exact circle and between
the filled polygon and the exact disk. It applies to a single 2D Cartesian
Point, one positive radius, and a hole-free Polygon. It does not change any
existing `buffer` overload. Integral, IEEE `float`/`double`, and `cpp_rational`
values are supported. Native floating inputs denote their exact binary values
and require gradual underflow and no fast-math.

Each quadrant is recursively bisected in its rational chart parameter until
the kernel accepts the exact chord with half the requested tolerance. This
gives room for rounding but is not itself an export certificate. The algorithm
rounds coordinates to nearest, ties to even (or keeps them exact with rational
output), builds the requested open/closed and clockwise/counterclockwise ring,
then reads its stored coordinates back exactly. It checks two properties:

1. Each directed edge, read counterclockwise, advances by an angle strictly in
   `(0, pi)` about the center, and the cycle wraps exactly once. Each edge
   intersects every ray in its angular interval exactly once; those intervals
   cover the full turn with disjoint interiors. This rules out crossings,
   duplicates and collapses, and proves the polygon is star-shaped about the
   center. Consecutive collinear edges are allowed.
2. Every stored vertex is strictly inside the outer radius `r + error`.
   Every stored edge stays strictly outside the inner radius
   `max(0, r - error)`, using the exact minimum distance on the segment
   (including an interior projection when applicable). The maximum distance
   on each edge is at a vertex. The boundary radius on every ray is then
   within the strict error bound. Pair the circle and polygon boundary along
   the same ray; for the filled sets, clip to the closer boundary on that ray.
   The maxima are attained on compact sets, so the bound is strict globally.

`invalid_input` covers nonfinite inputs, a nonpositive radius, and negative
tolerance. `cannot_certify` covers zero tolerance, exhausted subdivision
limits (depth 16 per quadrant, 4096 vertices before closure), coordinate
overflow and a failed final-ring certificate, including on coarse grids.
Either returned failure leaves the caller's polygon unchanged. Failure does
not prove no admissible polygon exists on the target grid: the algorithm only
tries one sampling and nearest-rounding policy.

## Remaining Work Beyond Point Buffers

- General buffering needs strategy and input-policy integration, degeneracy
  handling, join/cap selection and exact sweep ordering, and a tessellator
  that preserves the exact offset endpoints and certifies every chord in
  the chosen sweep. It needs a termination/resource policy, including what
  to do when a positive sweep is requested with zero tolerance.
- Forming the actual buffer requires exact intersection construction,
  ordering and noding for offsets and chords, followed by topology assembly
  and cleanup. Primitives from different edges can introduce expressions
  outside the kernel's one-radical comparisons. Neither exact topology
  nor an end-to-end distance or containment bound after assembly is proved
  by these local certificates.
- General native-coordinate export needs certified rounding of algebraic
  vertices, consistent rounding of shared vertices, and an explicit budget
  for export displacement in addition to tessellation and any assembly
  error. A per-vertex displacement bound bounds corresponding straight
  segments by interpolation, but does not prevent collapsed edges,
  crossings, or other topology changes. A public contract needs a topology
  check or proof, and a defined failure or higher-precision path when the
  target coordinate grid cannot meet the requested guarantee.

## Focused Tests

From the Geometry checkout, with `BOOST_ROOT` pointing at the Boost headers:

```sh
CASE_FILTER=round_buffer sh test/strategies/run_exact_primitives.sh
c++ -std=c++17 -Iinclude -I"$BOOST_ROOT" \
    test/algorithms/buffer/buffer_point_certified.cpp -o /tmp/buffer_point_certified
/tmp/buffer_point_certified
```
