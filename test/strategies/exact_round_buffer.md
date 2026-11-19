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

## Original #578 Fixture and Export Obstacle

`buffer_issue_578_certified` transcribes the two WKT strings from the
[original Wandbox program](https://wandbox.org/permlink/wvGM3SFlM5MXLx5Z)
and reads them as binary64, like the original. There are 16 inner and 101
outer vertices. The source computes `R = 4`, `n = 180`, and
`r = R - R*(1-cos(pi/n)) - 1e-6`. The 16-point trace is **not** a subset of
the 101-point trace: inner vertex 5 is more than `0.39` from every outer
segment. The sufficient subtrace theorem above cannot certify this fixture.

The new test nevertheless proves direct nesting. It checks the original
binary64 formula and also bounds the real formula by rational inequalities
(using a bracket for pi and alternating cosine bounds). The inner radius
is less than `3.9993898`. The test certifies containment of the entire
radius-`3.9993898` smooth inner buffer in the open radius-`3.99939075`
neighborhood of the original outer trace. In particular, **any** outer
support-exact radius-4 chord union that contains its radius-`4-e` inner
neighborhood for `e = 0.00060925` contains the original inner smooth buffer
strictly. An inscribed inner chord approximation also lies strictly inside.
No radius change is needed at the ideal-set level. This does not certify the
legacy 180-point output or assume its actual chords meet the nominal sagitta.

For this certificate all the outer capsules have a common interior point,
so their union is star-shaped. Each inner vertex's circle is covered by
rational tangent triangles: two exact chart directions bound an arc between
their chord and tangents; put the triangle's three vertices strictly in one
convex outer capsule. Split the chart interval when that test fails (depth
12 and 4096 total arc tests). Each of the 30 true inner offset sides also
lies in some convex outer capsule, as verified by algebraic endpoint to
rational segment distance tests. This covers each inner capsule boundary;
radial segments from the common outer point then cover the filled inner
capsule. The tangent triangles are enclosures for the proof, **not** a
replacement for any offset support or exported chord geometry.

Both original traces have an exposed source edge of irrational length
(zero-based edges 0 and 29, respectively). All source vertices are on or to
the right; the entire true left offset is exposed in its normal direction.
On that side `cross(d,p-a) = r*sqrt(d.d)` is irrational for either original
positive binary64 radius, whereas the left side is rational at all rational
points. No finite polygon with rational or native coordinates can realize
that exposed face losslessly: its support maximum must occur at a vertex.
The test therefore checks `cannot_certify` and unchanged output for both
original lines under the current public contract. The outer line also
exceeds the current 32-piece limit independently of the algebraic issue.

### Bounded Algebraic Predicates and Pair Noding

The isolated `bounded_radical_field` keeps up to four independent positive
rational square classes (degree at most 16). An incoming root is tested
against each existing subset product by an integer perfect-square test of
the quotient's numerator times denominator, so dependent roots and exact
cancellations are identified without factoring. The test is exhaustive:
if a root of a rational already belongs to the field, every generator
sign-flip sends it to itself or its negative. Each simultaneous eigenspace
of these sign-flips is spanned by one basis monomial, so the root must be
a rational multiple of a subset product. A fifth independent root is
refused without changing the field or output. No interval refinement or
rounded root is used in the sign fallback.

Write a sign query as `A+B*sqrt(d)` over the lower field. Equal signs or a
zero term decide immediately; opposite signs reduce to the lower-field sign
of `A*A-d*B*B`, with the sign of `A` accounted for. At rank four this visits
at most 121 recursive calls (including rational leaves) and at most 350
coefficient-pair products across all those norm reductions. One ordinary
field multiplication uses at most 256 pairs; registering a root tests at
most 16 subsets. These are operation-count bounds, not constant bit-cost
bounds for arbitrary precision rational coefficients.

`radical_buffer_noding` lifts true edge offsets and nodes two nondegenerate
algebraic segments with exact signs, using projective parameters for proper
crossings, coincident endpoints and overlaps. It compares cuts without
evaluating an algebraic quotient. Tests cover independent roots, a gap of
`2^-300`, and the exact intersection of two independently algebraic cap
chords from the first two **original outer** edges, each within the proved
radius-4 chord-error budget. Their node is not rational. This checks a
local arrangement operation, not a full union or a repair of `buffer`.

## Common-Witness Native Export

`<boost/geometry/algorithms/buffer_linestring_star_certified.hpp>` adds the
separate opt-in `buffer_linestring_star_certified(linestring, witness,
polygon, radius, sagitta, max_displacement)`. The witness must be strictly
inside every radius-`r-sagitta` source capsule, verified by exact squared
point-to-segment distances. This implies strict interior to each ideal
inscribed chord capsule. Neither the legacy `buffer` nor the existing
lossless `buffer_linestring_certified` contract changes.

### Exact Ideal Contour

Use the true straight sides and algebraic tangent endpoints of each
capsule. The rotated chart vertices have one root per capsule, and exact
chord selection certifies their sagitta. The whole radius-`r-sagitta`
capsule is covered: the side strip is exact, and the forward/backward
semicircles have that radial lower bound on every chord. The ideal union
therefore is inscribed in the radius-`r` neighborhood and contains the
radius-`r-sagitta` neighborhood. It is star-shaped about the witness.

Sweep in angular order about the witness and merge two contours at a time.
Every open angular interval has exactly one active producer edge from each
contour. In witness-relative coordinates a CCW edge `a -> a+d` has radial
parameter `cross(a,d)/cross(u,d)` on ray `u`. Both denominators are positive
on the interval, so the difference has the sign of a linear form in `u`.
The two knot tests use an actual point on one active edge and an exact side
test against the other. Equal signs mean no switch; opposite signs mean one
switch at the exact producer intersection. A tie at both ends means the
active supports coincide on the interval. Angular comparison and side tests
may use outward-rounded binary64 intervals, but an uncertain sign rebuilds
the exact local field from the producers. Knots retain either a source
vertex, an exact pair of producer edges, or the rational positive-x seam.
No rounded coordinates become producer geometry.

The path is bounded to 128 nonzero pieces, 65536 primitive vertices, 16384
marks per contour, and 32768 exact filter fallbacks; the depth and local
four-root caps still apply. Each sweep merge visits the two angular lists
once. A non-interior witness, a bounded-work or precision failure, or an
invalid stored ring results in `cannot_certify` rather than a different
construction. Point traces are not handled here. The binary64 filters
require the exact Cartesian floating environment (no fast-math, implicit
contraction, or flush-to-zero).

### Displacement and Original #578 Nesting

A rational evaluation using 128-bit dyadic root brackets only proposes
the native storage values. Re-evaluate the ideal projective vertex in its
exact field, and compare each stored coordinate with it against half the
displacement budget `d`. This proves Euclidean vertex displacement `<= d`
by the L1 bound, separately from the ideal chord sagitta. The same bound
holds along corresponding edges by interpolation. The actual stored ring
must also make one simple turn about the same witness. Matching the two
triangle fans about that witness gives filled-set Hausdorff displacement
`<= d`. These tests do not say the native vertices themselves are exact
support incidences or even necessarily the nearest native values.

On the original fixture use the common witness `(85,-62)`, the original
radii, sagitta `0.00060925`, and displacement `d = 2^-40` independently
for each result. Both binary64 exports succeed. The earlier cover
certificate places `B(L_inner, 3.9993898)` strictly inside the ideal outer
region. Since the actual inner radius plus `2*d` is smaller than
`3.9993898`, every point of the exported inner polygon is at distance
strictly greater than `d` from the ideal outer boundary: its closed
`d`-neighborhood lies inside that larger certified inner buffer. The
straight homotopy between corresponding outer boundary points cannot
cross it, so its winding remains inside after export. Independently,
the fixture sweeps the **stored binary64 rings** using exact rational
tests at both sets of angular knots and proves their radial order strict.

The same test compares a small integral star merge with the existing
lossless rational union at zero displacement and rejects a bad witness,
zero displacement on the irrational fixture, and a float grid that cannot
meet `d`. It also rejects a collapsed float ring even with a generous
displacement budget. Run the fixture binary with `--emit-wkt` to print the certified
inner and outer polygons in that order, with binary64 round-trip precision.

## Bounded Non-Star Linestring Export

`<boost/geometry/algorithms/buffer_linestring_arrangement_certified.hpp>`
adds `buffer_linestring_arrangement_certified(line, polygon, r, e, d)`.
This is a separate opt-in for a nonempty 2D Cartesian linestring with at
least one nonzero edge. It requires `0 < e < r` and `0 <= d < r-e` and does
not require a common witness. The ideal and stored outputs must have one
simple shell and zero or more simple, pairwise disjoint, unnested holes
strictly within it. This admits a connected chord union with holes, not an
arbitrary arrangement of touching or branching rings. No existing buffer
contract changes.

### Local Algebraic Arrangement

Use the same exact capsule producers as the star path, without interpreting
the coordinate anchor as an interior point. The true sides and inscribed cap
chords are never rounded before noding. Conservative interval boxes reject
disjoint candidates; intersect the remaining producer pairs in a local
radical field. A cut keeps an original vertex or the pair of producers that
defines an exact projective intersection. Collinear overlaps supply both
ends. Lexicographic comparison is monotone on a straight segment, so exact
coordinate comparisons can sort and deduplicate cuts without rounding the
parameters. Comparing two cut coordinates uses at most four producer roots.

For an open atom, its owner capsule occupies the left neighborhood. Test
the midpoint against every candidate convex capsule's exact half-planes.
Another capsule hides that atom exactly when it covers the right
neighborhood: either the midpoint is interior or it lies on an oppositely
directed facet. Noding excludes any transverse crossing or vertex in the
open atom. The midpoint of two cuts on a shared producer uses at most
three source roots, and testing another facet adds at most one. Identical
directed atoms are deduplicated; a branch or a point contact in the retained
graph is refused. The sign of the turn at a ring's lexicographic minimum
distinguishes a shell from a hole without forming a global radical area
sum. Exact ray tests verify that holes lie inside the shell and do not nest.

Limits are 16 nonzero source edges, 384 primitive edges, 8192 pair checks
after box rejection, 4096 open atoms and 16 cycles. The four-root, chord
depth and 32768 exact filter-fallback limits also apply. Coordinate bit
costs still depend on the inputs. This is a bounded pairwise arrangement,
not a scalable general sweep. Its interval filters require the same
binary64 environment as the star path.

### Native Displacement, Topology and Containment

Only after assembly, propose native vertices and certify each against its
exact algebraic value. The projective coordinates share a positive
denominator `D`. For a proposed stored point, subtract its exact rational
coordinates (relative to the anchor) times `D` from the ideal numerators.
The sum of the squared residuals must be at most `(d*D)^2` in the local
radical field. This uses the specified Euclidean budget rather than the
more restrictive rule of allowing only `d/2` on each axis. Verify the
stored rings with exact rational simplicity, orientation, inter-ring and
nesting checks, then read back the actual Polygon before publishing it.
A failed check returns `cannot_certify` and leaves the output unchanged.
The ideal shell and holes correspond one-for-one to native rings, not to
the possibly different smooth-buffer topology. Each matched boundary pair
is within Euclidean displacement `d`, including along edges by interpolation.

For `N(L,s)` the closed radius-`s` neighborhood of the trace, the support-
exact chord union `U` satisfies `N(L,r-e) subset U subset N(L,r)`. Winding
at any point farther than `d` from the ideal boundary is unchanged along
the straight homotopy between paired ring vertices and edges. Thus the
stored region `P` satisfies

```
N(L,r-e-d) subset P subset N(L,r+d)
```

The lower inclusion is strict for smaller radii, and equality follows by
closure at the positive radius `r-e-d`. This gives a filled-set Hausdorff
bound `<= e+d`; `e` and the export bound `d` remain separate. It does not
certify boundary Hausdorff distance to the smooth buffer. For nested traces,
`r_inner+d_inner < r_outer-e_outer-d_outer` suffices for strict nesting of
the two stored sets.

For [#1233](https://github.com/boostorg/geometry/issues/1233), the reduced trace
`LINESTRING(0 100,-4503599627370496 0,0 -100)` returns empty with the
ordinary Cartesian buffer on the local replay. A unit offset to the left
of the tip is representable. The arrangement constructs from the original
sloping supports after exact rational translation; it does not flatten the
arms or round before noding. With `r=1`, `e=1/8`, and `d=1/2`, it exports
a simple binary64 ring covering both arms and preserving the notch.
The same fixture is refused with `d=1/4` or `d=0`. An exposed cap sample has
x offset strictly between `-3/4` and `-1/2` relative to the tip, so neither
adjacent binary64 x coordinate is within `1/4`. The displacement contract
does not promise nearest rounding. The existing star and lossless-linestring
entry points are unchanged.

The published report's far x coordinate rounds to `-2^63` for binary64.
Under those same budgets the minimum x of a certified output would have
to be in `[-2^63-3/2, -2^63-3/8]` by the outer and inner inclusions.
There is no binary64 number in that interval (the next smaller number is
`-2^63-2048`). The offending linestring component is a refusal test.
This is an output-grid obstruction distinct from loss of the representable
reduced trace in the ordinary algorithm.

The open fixture `LINESTRING(0 0,8 8,0 12,-8 4,-0.5 0.25)` has source
lengths in two independent square classes. Its opposite parallel sides
are `sqrt(72)` apart, so radius-1 capsules have no common witness. The end
gap closes inside even the smaller certified neighborhood. Tests export
both radii `1` and `3/4` with `e = 1/8`, `d = 2^-42`, and verify one hole
in each. An independent rational oracle tests material membership and
every pair of shell/hole edges for strict nesting after export. Separate
integer tests check a square and a connected two-hole trace at `d = 0`,
and refusal tests cover an exact pinch at a point, a tight/zero displacement
budget and work limits.

## Algebraic Polygon Arrangement and Native Export

`<boost/geometry/algorithms/buffer_polygon_arrangement_certified.hpp>`
adds `buffer_polygon_arrangement_certified(input, multi_polygon, distance,
e, d)`. Positive distance dilates; negative distance applies regularized
erosion. Let `r = abs(distance)`; require `0 < e < r` and `0 <= d < r-e`.
As in the rational path the input is one finite 2D Cartesian Polygon with
correctly oriented simple rings and pairwise disjoint, unnested holes
strictly inside its shell. There is no square-class restriction on source
edge lengths, and no common-witness premise. This is a native opt-in;
legacy, star, and lossless buffer contracts are unchanged. In particular
native lossless output can still be impossible for exposed irrational edges.

### Exact Face Selection

Node the rational source rings and exact algebraic chord capsule producers
together, including collinear overlaps. The same bounded local fields and
exposed-atom classification as the non-star linestring case find the chord
band boundary. Let `P` be the filled source and `T` the union of inscribed
boundary capsules. Since `T` contains the positive radius-`r-e` neighborhood
of the entire source boundary, both sides of a source-boundary atom are in
`T`. Such an atom cannot be on either output boundary. Any exposed band atom
is away from the source boundary, so exact rational-ray comparisons at its
algebraic midpoint determine a single material state on both sides. For
`P union T`, keep the exposed atom precisely when it is outside `P`. For
`P minus T`, keep it when it is inside `P` and reverse it so material stays
on the left. Closed selected faces give the regularized result; a zero-area
remnant is not exported. No source, facet, cut or midpoint is rounded before
these decisions.

Assemble disjoint exact cycles. Local extreme-turn and ray predicates build
their containment forest, with alternating shell and hole orientations.
This handles multiple components, holes and empty erosion. Franklin's
exact Euclidean vertex-displacement check is applied with budget `d`;
the bound extends along each pair of matched edges by interpolation.
Exact rational checks reject native edge collapse, crossings, contacts,
orientation changes or any difference from the ideal cycle containment
matrix. Every actual stored ring is read back before publishing. This
preserves chord-result topology, not necessarily the topology of the
smooth round result. On refusal the output is unchanged.

The arrangement caps remain 16 nonzero source edges across shell and holes,
384 combined producer edges (including source edges), 8192 pair checks after
box rejection, 4096 open band atoms and 16 cycles. The four-root cap is per
predicate, not per Polygon. Chord depth and exact filter-fallback limits
still apply. Output contacts and branch nodes are refused. Pairwise noding
and export validation are quadratic in their edge counts; this is not a
scalable sweep. The same binary64 environment is required for the filters.

### Radius Certificates

Write `N(P,s)` for the closed `s`-neighborhood of filled `P`, and `E_s` for
the closure of its interior points with clearance strictly greater than
`s` from the source boundary. After native export:

```
dilation: N(P,r-e-d) subset output subset N(P,r+d)
erosion:  E_(r+d) subset output subset {x in int(P): clearance(x) > r-e-d}
```

The ideal band contains the source boundary's radius-`r-e` neighborhood
and lies inside its radius-`r` neighborhood. These give the ideal bounds
with `d = 0`. Winding is unchanged outside the closed `d`-neighborhood of
the ideal boundary during the paired-edge homotopy. For dilation this yields
the native neighborhood bounds (the inner equality follows by closure),
so filled-set Hausdorff error is at most `e+d`. For erosion, every point of
clearance `> r+d` has its `d`-ball in the ideal erosion interior. Conversely,
the stored erosion lies within distance `d` of the ideal erosion, whose
clearance is strictly `> r-e`; since `d < r-e`, that short displacement
cannot cross the source boundary. These are inclusion bounds, not a
general erosion Hausdorff or component-correspondence bound to the smooth
erosion. A topology change near a smooth critical radius is still possible.

The focused tests use a concave shell with two oblique holes and five
independent length square classes across the source. They check positive
and negative outputs against independent exact rational ray and distance
samples, including dilation that fills the holes. A sheared dumbbell has
no common witness and yields two erosion components: chamber samples meet
the inner inclusion bound, and a narrow cross-section meets the exclusion
bound. A sheared slab gives empty erosion. Tests also exercise invalid
closure and lossless/native-grid refusal. Run the fixture binary with
`--benchmark` for end-to-end timing, excluding the sample oracle.

One local Darwin arm64 run (Apple Clang 21, `-O2 -DNDEBUG -ffp-contract=off`,
20 calls per fixture, including construction, certification and export):

| Fixture | Mean per call | Output vertices before closure |
| --- | ---: | ---: |
| Concave shell, two holes, dilation | 6.5 ms | 34 |
| Concave shell, two holes, erosion | 9.1 ms | 85 |
| Sheared dumbbell, split erosion | 8.1 ms | 28 |

These are small-fixture timings, not a general scaling claim or a speed
comparison with the different lossless or legacy contracts.

## Public Round Buffer Dispatch

`<boost/geometry/algorithms/buffer_round_certified.hpp>` provides
`buffer_round_certified(input, output, distance, arc_sagitta, max_displacement)`.
For a single Polygon it uses the following certificates, without changing `buffer()`:

| Case | Certificate |
| --- | --- |
| Zero distance | Validate the source using exact rational coordinates, then copy the same filled polygon losslessly with the output orientation and closure. |
| Convex erosion | Exact straight offsets, with an independent export bound. A failed exact construction or export is not replaced by a chord approximation. |
| Other erosion or dilation, native coordinates | Algebraic chord arrangement, including irrational edge lengths. If necessary, tighten the budgets to keep `0 < e < r` and `0 <= d < r-e`. |
| Rational chord arrangement | Lossless export for rational-length sources, including as a fallback on a native grid. `cpp_rational` output does not use the native algebraic exporter. |

On a native grid, `d > 0` permits bounded displacement; `d = 0` requires
every vertex losslessly. Zero-distance identity is lossless even if `d > 0`.
The native allowance does not relax the arc sagitta or permit dropping rings.
Unsupported geometry, a failed certificate, or a split that cannot fit a
single `Polygon` returns `cannot_certify` without changing the output.

For chord erosion, let `C = reg(P \ T)` where `T` is the closed inscribed
boundary band, and `reg(S)` is the closure of the interior of `S`. Since
`T` is contained in the true closed radius-`r` boundary neighborhood, every
point with clearance strictly greater than `r` is in the interior of
`P \ T`. Taking closures gives `E_r subset C`. A successful export preserves
all positive-area ideal components, including at `d > 0`. Thus an empty
successful result proves `C` and `E_r` empty. Failed export proves neither.
The reverse implication is false: a chord can miss part of a round cap.
Nonempty chord results have the radius bounds above; neither their
nonemptiness nor their component count is a claim about smooth erosion.

## Multi-Input Algebraic Arrangement

`<boost/geometry/algorithms/buffer_multi_arrangement_certified.hpp>`
provides `buffer_multi_arrangement_certified(input, output, distance, e, d)`
for 2D Cartesian MultiPoint, MultiLineString and an admitted subset of
MultiPolygon. Output is a native MultiPolygon. Point and line inputs require
positive distance; polygon input admits both signs. As in the single
arrangement paths, `r = abs(distance)`, `0 < e < r` and `0 <= d < r-e`.
A wholly empty collection succeeds with empty output on admissible budgets.
Empty line members and default-empty polygon members are ignored; singleton
or all-equal lines contribute point traces. A Polygon with an empty shell
and present hole rings is not an empty member.

Points append exact rational inscribed disk chords. Line members append the
same exact algebraic segment capsules as the single-line path. The union is
noded **once** over the entire collection, so components and holes do not
depend on export of individual members. Polygon members each contribute
their exact source rings and boundary capsules. The source face classifier
computes the union of the members' filled sets, each evaluated as its shell
minus its own holes. Cross-member source boundaries must be disjoint, and
neither shell may lie in another member's filled material. Exact pair and
membership tests enforce those conditions. An island strictly inside another
member's hole is admitted. Valid OGC configurations with source boundary
contacts are outside this opt-in class. The boundary of the admitted filled
union is the union of its individual boundaries, so the polygon face action
and radius proofs apply without changing their meaning.

All existing work caps are shared across the collection: 16 point disks
and/or nonzero source edges, 384 combined producer edges, 8192 pair checks
after box rejection, 4096 open atoms and 16 cycles. Coincident input points
can consume separate disk slots; coincident directed boundary atoms are
deduplicated exactly. Successful export preserves the ideal chord-result
shell/hole containment forest and separately checks Euclidean displacement
`d`; there is no rounded-member fallback. Exact or native boundary contacts
are refused. Dilation has filled-set error at most `e+d`, not a smooth
boundary or topology guarantee. Polygon erosion has the same radius
inclusions as the single-polygon path. On failure the output is unchanged.

The regression suite checks overlapping and disjoint disks and traces, a
hole formed by more than one line member, touching and native-grid refusal,
empty inputs and members, merged polygon dilations, and an island inside a
hole at both radius signs. The full [#1233 published MultiLineString]
(https://github.com/boostorg/geometry/issues/1233) is also tested. At radius
`1`, `e = 1/8`, `d = 1/2`, the far-left input parses as `x = -2^63`. The
certified inner and outer neighborhoods would force the minimum output
abscissa into `[x-3/2,x-3/8]`; the next smaller binary64 value is `x-2048`.
The original therefore must return `cannot_certify` for this native output
contract, independent of arrangement capacity. The reduced `2^52` hairpin
does export when split across two line members beside a separate component.

### Unified API Admission and Erosion

`buffer_round_certified` dispatches all three multi tags to the shared
arrangement, with the usual staged output. A `Polygon` destination is only
published when the entire result fits one component. The wrapper admits
empty multi inputs and members. Valid nonpositive point and linear inputs
regularize to empty, and positive buffering of a wholly empty collection
also succeeds with empty output. Nonempty, nonzero multi chord operations
require positive sagitta and native output coordinates; the wrapper can
tighten, but never relax, a supplied error budget. In particular there is
no new componentwise convex erosion path or nonzero `cpp_rational` multi
fallback. Caps apply to the collection, not individually to each member.

MultiPolygon uses the same exact source admission for all distance signs,
including zero: simple valid members with strictly disjoint filled sets and
disjoint boundaries. At zero, copy those source rings with the destination
orientation and closure, checking exact read-back even if the caller
provided a positive displacement allowance. Unrepresentable coordinates
or a destination with room for only one of several members are refused
transactionally. Nothing is tessellated at zero distance.

For the negative-distance domain, let `P = union P_i` be the compact
filled union. Pairwise disjoint compact members have positive separation,
so `boundary(P) = union boundary(P_i)`. For `x` in `int(P_i)`, a straight
segment to any other member's boundary exits `P_i` first. Thus its clearance
from `boundary(P)` equals its clearance from `boundary(P_i)`. Taking the
closure of the strict clearance sets gives `E_r(P) = union E_r(P_i)` for
the **smooth regularized** erosion on this admitted domain. This identity
does not justify componentwise chord construction and native rounding:
the backend nodes one `T` from all boundary capsules and selects
`reg(P minus T)` using the whole material union. Its certified native
inclusions remain `E_(r+d)(P) subset output` and `clearance > r-e-d`
inside `P` for every output point. Dilation similarly selects `P union T`,
so disjoint sources can merge into a single exported result.

For shared edges the premise fails. For example, adjacent `[0,2] x [0,2]`
and `[2,4] x [0,2]` have a shared seam; `(2,1)` survives radius-`1/2`
erosion of the filled union but survives neither member's erosion. The
exact source checks reject this case (and contained filled overlap) before
any output is published, rather than inserting a spurious seam in `T`.

The public dispatch tests cover all three dimensional tags and both output
forms. A fractional line member fully hidden by a closed line's buffer,
and a fractional island hidden by dilation that fills a hole, each refuse
lossless export alone but succeed inside the combined exact union at `d=0`.
Other cases check separated erosion with an island in a hole, shared-edge
refusal, exact zero-distance round trips, failed zero-distance native
conversion, and preservation of caller output on failure.

## Remaining Work

- The algebraic pair primitive is integrated into the bounded opt-in
  linestring and polygon paths, not the legacy buffer. General arrangement
  scaling and point-touching topology remain open. Native lossless output
  is still impossible on #578; its separate common-witness export is the
  supported path for that larger original fixture.
- The lossless linestring API still accepts one exterior ring only. The
  algebraic polygon APIs support multiple disjoint cycles but not boundary
  contacts, overlapping source MultiPolygon members or boundary-only mode.
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
c++ -std=c++17 -ffp-contract=off -Iinclude -I"$BOOST_ROOT" \
    test/algorithms/buffer/buffer_issue_578_certified.cpp -o /tmp/buffer_issue_578_certified
/tmp/buffer_issue_578_certified
c++ -std=c++17 -ffp-contract=off -Iinclude -I"$BOOST_ROOT" \
    test/algorithms/buffer/buffer_linestring_arrangement_certified.cpp -o /tmp/buffer_linestring_arrangement_certified
/tmp/buffer_linestring_arrangement_certified
c++ -std=c++17 -O2 -DNDEBUG -ffp-contract=off -Iinclude -I"$BOOST_ROOT" \
    test/algorithms/buffer/buffer_polygon_arrangement_certified.cpp -o /tmp/buffer_polygon_arrangement_certified
/tmp/buffer_polygon_arrangement_certified
/tmp/buffer_polygon_arrangement_certified --benchmark
c++ -std=c++17 -ffp-contract=off -Iinclude -I"$BOOST_ROOT" \
    test/algorithms/buffer/buffer_multi_arrangement_certified.cpp -o /tmp/buffer_multi_arrangement_certified
/tmp/buffer_multi_arrangement_certified
```
