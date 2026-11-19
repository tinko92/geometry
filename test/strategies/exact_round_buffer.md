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

## Public Segment Capsules

`<boost/geometry/algorithms/buffer_segment_certified.hpp>` adds
`buffer_segment_certified(segment, polygon, radius, max_error)` for a 2D
Cartesian Segment. The Point and Segment APIs share `certified_buffer_status`
(`certified_point_buffer_status` still names the same status type). Equal
endpoints delegate to the Point path. `buffer_collinear_linestring_certified`
in the same header accepts exactly collinear nonempty Linestrings: a connected
trace on a line covers the entire segment between its extreme vertices, even
with repeats and backtracking. An all-equal trace gives a disk. A bent trace,
including a nonzero bend below binary64 resolution, returns `cannot_certify`.

### Exact Sides and Caps

Let `d = b - a`, `n = (dy,-dx)`, and `D = d.d`. Rotate the unit-circle chart
by using the right-handed frame `(n,d)/sqrt(D)`. For the rational chart sample
`(u,v)`, the ideal cap vertex is `center + r*(u*n + v*d)/sqrt(D)`. The first
two quadrants run around `b`; the last two run around `a`. The cap endpoints
are exactly the kernel's left and right offset endpoints. Connecting matching
normals at different centers gives the two true straight offset sides, not
chords between perturbed tangent directions. The cap chords use the exact
Point subdivision predicate in the rotated frame, with unchanged angles,
radius and sagitta. All ideal coordinates use the same symbolic square root.

For export, write `D = N/M` and bracket `1/sqrt(D)` from the integer square
root of `(N*M) << (2*bits)`. A zero remainder gives an exact rational scale,
including nonintegral Pythagorean edge lengths. Otherwise native integral or
IEEE floating export requires the two interval endpoints to round to the
same nearest/ties-even value. It refines the shared scale up to 4096 bits
and fails conservatively on an unresolved rounding or overflow. Rational
output takes the interval midpoint for irrational coordinates. The stored
polygon is not claimed to have exactly parallel sides after conversion.

### Export Certificate

Let `S` be the exact source segment, and let `rho(p) = dist(p,S)`. Re-read
the stored ring and require a simple single turn about the source midpoint,
using the Point topology check. Require every vertex to satisfy
`max(0,r-error) < rho(p) < r+error`. The upper bound extends to the edges
and the filled polygon because the outer capsule is convex. To prove the
lower bound on each edge, check source endpoint distances to that edge and
check for proper crossings between the edge and `S`. The output endpoint
distances were already checked. A minimum distance between disjoint planar
segments is attained with at least one endpoint; touches and overlaps give
zero endpoint distance. A proper crossing also gives zero, even when all
four endpoints are strictly clear.

The lower test puts the whole inner capsule in the polygon: it contains the
source midpoint and is connected and disjoint from the polygon boundary.
The bounds are strict on the compact boundary, so these inclusions hold
with some smaller error `e < error` (using inner radius zero if `e >= r`).
Along a ray from a closest point of `S`, that point stays closest for all
positive radii. On such rays, pair a polygon boundary point with the point
at radius `r`, or pair a capsule boundary point with a polygon boundary
crossing between the inner and outer radii. The same nesting bounds distances
for the filled sets as well. This proves the **strict Euclidean Hausdorff
bound after output conversion** for both boundary and filled region.

Zero tolerance, exhausted depth or vertex limits, rounding failure, and
failure of the stored-ring checks return `cannot_certify` without modifying
the output. This is not a search over alternative polygons on the target
grid. The numeric and geometry requirements are the same as for Point.

## Exact Chord Unions for Rational-Length Linestrings

`<boost/geometry/algorithms/buffer_linestring_certified.hpp>` adds
`buffer_linestring_certified(linestring, polygon, radius, max_error)`.
Every nonzero source edge must have rational length: squared length must
be a perfect square in the rationals. This covers axis-aligned edges and
rational multiples of Pythagorean triples, but rejects, for example, a
unit diagonal. Repeated vertices are skipped; an all-equal trace yields a
disk approximation. The radius and samples are rational and the exact-scale
path of the capsule kernel is required, so both true supports and all
selected cap chord endpoints are rational **before** any union operation.

### Exact Chord Topology and Lossless Export

The cap-and-side polygons are strictly convex. For each oriented piece edge,
clip its parameter interval against all the convex pieces using rational
half-plane inequalities. Split at every interval endpoint, including
singleton contacts and co-directed overlaps. The piece supplying an open
atom occupies its left side. Another piece covers its right side exactly
when the atom lies in the interior of that piece or on its oppositely
directed facet. Remove such atoms and deduplicate identical oriented atoms.
Every remaining atom bounds the union, and any boundary contact between
noncoincident atoms was noded by the clipping. Require one directed simple
cycle consuming all atoms. Holes, point-touch branches, and additional
cycles are refused. Limits are 32 nonzero pieces, 512 total piece edges,
and 4096 boundary atoms (in addition to the capsule subdivision limits).

`cpp_rational` output preserves the result. Native output succeeds only if
the stored ring reads back as the same exact sequence in the requested order
and closure. A failed round trip, topology or resource check returns
`cannot_certify` without changing the caller's polygon.

The topology claim is about the **exact union of these chord polygons**, not
the smooth capsule union. We do not union rounded pieces or round the union
and infer topology from a metric tolerance. In `test_rounding_can_change_contact`,
two exact capsules (and hence their inscribed chord polygons) are separated
by `1/12`; their separately certified integer exports share a boundary point.
Hausdorff error of a finite filled union is bounded by the maximum per-piece
error, but that bound does not stop a new contact, crossing, or hole change.
The same distinction applies between the chord union and the smooth buffer:
the inclusion guarantee below bounds the filled-set Hausdorff distance, not
boundary Hausdorff distance or the topology of the smooth buffer.

### Sufficient Nested-Linestring Condition

Write `B(L,r)` for the closed radius-`r` neighborhood of the trace `L`, and
`P(L,r,e)` for a successful export. For `0 < e < r`, each chord capsule is
inscribed in the exact capsule and strictly contains that segment's radius
`r-e` capsule by the lower-distance check. Unions distribute over both
collections of capsules, so

```
B(L,r-e) subset interior(P(L,r,e))
P(L,r,e) subset B(L,r)
```

Finite pieces and strict clearance permit a slightly smaller error budget,
which gives strict filled-set Hausdorff error `< e`. For `e >= r` the
implementation uses `r/2` instead to keep a positive inner radius.

For [#578](https://github.com/boostorg/geometry/issues/578), suppose the inner
trace is contained in the outer trace and `r_inner <= r_outer - e_outer`.
Then the inner output is contained in `B(L_inner,r_inner)`, which is contained
in `B(L_outer,r_outer-e_outer)` and therefore in the **interior** of the outer
output. The condition is sufficient, not necessary. No inner approximation
allowance is needed here because its output is inscribed and exported
losslessly; it would be needed for an output with outward displacement.

`test_578_nested_traces` directly buffers two noncollinear linestrings:
`LINESTRING(-6 -8,0 0,9 0,12 4,20 4)` and
`LINESTRING(-3 -4,0 0,9 0,10.5 2)`. It uses outer radius `4`, outer error
`1/16`, inner radius `63/16` (equality in the sufficient condition), and a
different inner error `1/8`. An independent rational ray test checks all
inner output vertices for strict containment, and rational pairwise
segment-intersection tests verify the two output boundaries never meet.
This is a regression for the inclusion property, not the historical issue
input and not a fix to the existing `boost::geometry::buffer` path.

## Rational-Length Polygon Buffers

`<boost/geometry/algorithms/buffer_polygon_certified.hpp>` adds
`buffer_polygon_certified(polygon, multi_polygon, radius, max_error, mode)`.
Modes are `dilation` (default), `erosion`, and `boundary_band`. The input is
one 2D Cartesian Polygon with a nonempty simple shell, correctly oriented
simple holes strictly inside the shell, and pairwise disjoint ring
boundaries. Exact input tests reject crossings, overlaps, touching rings,
nested holes, zero-length edges, and missing closure on closed models.
Every source edge must have rational length; otherwise the operation returns
`cannot_certify`. In particular there is no rounded-support fallback.

### Exact Arrangement and Export

Build the exact inscribed chord capsule for each shell and hole edge as for
linestrings. Call their union `T`. Node the source and capsule edges together
using exact rational segment parameters for crossings and both endpoints of
collinear overlaps. For each distinct open atom, classify both incident sides
of the atom against the Polygon (shell minus holes) and the capsule union.
The midpoint is not a vertex or a crossing after noding, so if a ring contains
it on its boundary, the oriented supporting edge determines the local inside
side. This avoids offset probes, epsilon nudges, or rounded joins. Select
faces by union (`P or T`), difference (`P and not T`), or band (`T`). Keeping
the boundary of the selected faces implements regularized set operations;
isolated points and zero-area lines are not emitted.

Orient each retained atom with material to its left and require disjoint
directed cycles. Positive cycles are shells, and each negative cycle is
attached to the innermost shell containing it. This permits holes, nested
disjoint components, and erosion that splits into several polygons or
disappears entirely. Branching and point-touching output cycles are refused.
Export re-reads all stored rings and requires lossless coordinates, with
the requested orientation and closure. Failure leaves output unchanged.
Limits are 32 source edges, 512 capsule edges, 32768 pair cut parameters,
and 4096 distinct atoms, in addition to the capsule subdivision limits.

### What Is Certified

Write `N_s` for the closed radius-`s` neighborhood of the source boundary,
and use effective error `e` (`r/2` if the requested error is at least `r`).
The capsule certificate gives `N_(r-e) subset interior(T)` and
`T subset N_r`. Therefore `boundary_band` and `dilation` (`P union T`) have
strict filled-set Hausdorff error `< e`: the strictly contained inner
neighborhood can be made uniformly a little larger across the finite
capsule set. This is **not** a boundary Hausdorff or smooth-topology claim.

For erosion let `E_s` be the closure of the interior points of `P` whose
distance from its boundary is strictly greater than `s`. The output is the
closure of the selected open faces of `P minus T`, so it contains `E_r`
and is contained in the set of points of `P` at clearance `> r-e`.
No general Hausdorff or component-correspondence bound follows near a
critical width or inradius. A neck whose entire cross-section is within
`N_(r-e)` must be blocked; chamber points at clearance `> r` must survive.
The dumbbell regression exercises both conditions independently.

For an accepted **convex shell without holes**, the regularized erosion is
exact. At an interior point the distance to the complement is the smallest
distance to a supporting half-plane's exterior. The perpendicular foot on
a minimizing support lies on a boundary edge: the disk of that minimum
radius fits in every supporting half-plane and touches the minimizing one.
The chord capsule contains the entire exact side strip over that edge, so
it contains the point whenever its clearance is `<= r`. Thus `T` and `N_r`
agree inside the convex Polygon. The 3-4-5 triangle regression checks the
exact inward half-planes, area, and disappearance at the inradius.

## Remaining Work

- Independently algebraic supports and noding are not implemented. The
  linestring API still accepts one exterior ring only; the polygon API
  supports multiple disjoint cycles but not boundary contacts. MultiPolygon
  input and lossy but topology-preserving native export are not implemented.
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
c++ -std=c++17 -Iinclude -I"$BOOST_ROOT" \
    test/algorithms/buffer/buffer_segment_certified.cpp -o /tmp/buffer_segment_certified
/tmp/buffer_segment_certified
c++ -std=c++17 -Iinclude -I"$BOOST_ROOT" \
    test/algorithms/buffer/buffer_linestring_certified.cpp -o /tmp/buffer_linestring_certified
/tmp/buffer_linestring_certified
c++ -std=c++17 -Iinclude -I"$BOOST_ROOT" \
    test/algorithms/buffer/buffer_polygon_certified.cpp -o /tmp/buffer_polygon_certified
/tmp/buffer_polygon_certified
```
