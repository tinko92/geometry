# Certified convex polygon erosion

Explicitly include `boost/geometry/algorithms/erode_convex_polygon_certified.hpp`
and call `erode_convex_polygon_certified(polygon, result, radius, max_displacement)`.
This API does not participate in default `buffer` dispatch. It takes a
hole-free 2D Cartesian Polygon and writes an empty or singleton MultiPolygon.
The radius is positive, the displacement nonnegative, and failure leaves the
output unchanged. Input closure and orientation follow the geometry traits.
Consecutive duplicate positions and same-direction subdivisions are redundant;
backtracking, non-convexity and self-crossing orders are rejected. In particular,
a locally left-turning star is not equivalent to a convex cycle: the exact edge
directions must make one strictly increasing sweep.

## Exact topology

For an inward (not unit) normal `n` of an edge through `a`, a radius-`r` disk
centered at `p` fits in the supporting halfspace exactly when
`n dot (p - a) >= r sqrt(n dot n)`. Intersection over the edges gives the exact
Euclidean erosion. Sequential halfspace clipping starts with the original
convex cycle and retains the source side and shift state for each active edge.
Each vertex is a Cramer intersection of two such supports. In particular an
intersection position is not reused as a producer of another support.

A vertex is affine in at most two square roots. Its signed clearance from
another support is affine in at most three, independent of total side count.
Outward binary64 intervals answer separated signs; rational root enclosures
and then a local exact radical field resolve uncertain cases, including zero.
An interval miss does not imply empty. While the current intersection has
interior, a cut preserves interior iff some current vertex is strictly inside
the new halfspace. Any remaining set of dimension zero or one regularizes to
empty and subsequent cuts cannot restore interior. An on-line vertex reuses
the correct incident supports rather than creating a zero-length edge.

Clipping uses quadratic time and linear active storage. Limits are 512
input positions (plus declared closure) and 4096 bits in an input scalar's
numerator or denominator, including radius and displacement. Work in each
exact sign is local and bounded in degree and field rank.

## Numeric export

Compute a rational bracket for each algebraic coordinate and round a proposal.
Read back the stored points, prove their Euclidean displacement in the rational
enclosures or in the local two-root field, and independently certify a single
strictly convex cycle with the correct closure and point count. Vertex
correspondence then bounds Hausdorff distance for the boundary and filled
polygon by interpolation and convex combinations, respectively. This is not
an inward-rounding or correctly-rounded algebraic-coordinate guarantee. If
export loses topology or cannot meet the requested displacement, return
`cannot_certify`, not empty. A zero budget requires exact representability.

## #1281 replay

The regression stores the supplied original WKT, not a simplified or fitted
contour. It reads binary64 counterclockwise closed input and runs the legacy
`+9` buffer with 36 points per circle, then feeds that output directly to the
opt-in `-9` erosion. Both input cycles validate in exact rational coordinates
and are strictly convex. On the local replay the 103 source edges produced a
184-side intermediate and **one** valid 161-side erosion with area
`764.98706637614`. The native export passed a `1e-9` displacement certificate;
an independent rational check also verifies every stored vertex clears every
intermediate support by at least `9 - 1e-9`. The optimized erosion call alone
took approximately 6 to 12 ms on local Apple Clang arm64 (not the dilation or
independent audit). This confirms the opt-in negative stage on this fixture;
it does not change the legacy behavior.

`exact_convex_erosion` also covers irrational offsets, forced three-root
signs, dependent roots, binary64 filter overflow, duplicate and redundant
supports, concave and self-crossing input, both orientation and closure
conventions, true point/segment collapse, and an ideal positive-area result
that would collapse on native export. To replay from a separate WKT file,
pass its filename as the test executable's only argument.
