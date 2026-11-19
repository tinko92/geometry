# Exact Cartesian ring parenting

## What is certified

Stored native operand rings and normalized generated triangles selected in an
areal Boolean operation are candidates. A source ID establishes provenance,
not validity. The check is lazy, when boundary location or residual
normalization needs it. Buffer and dissolve do not use either certificate.

### Stored originals

After the usual removal of consecutive duplicates and closure normalization,
split the cycle at its lexicographic minimum and maximum. From minimum to
maximum each chain must increase strictly. Interior vertices on either chain
must lie strictly on one side of the endpoint chord, and the two chains must
use opposite sides (one chain may consist solely of the chord). Check the
resulting orientation against the stored exterior or interior ring role.

Every point of a segment lies within its endpoints' lexicographic interval,
also for a vertical edge. Intervals of nonadjacent segments of the same chain
are disjoint. Each open half-plane is convex, so the two chains do not meet
except at their endpoints. At least one nonempty chain is strictly off the
chord. This proves a simple, nondegenerate cycle with unique nonzero edges
*within that ring*. This is a sufficient test, not a general validity checker.

As a result, a point on its raw boundary is on its surviving boundary. If
classification of a touching inner ring needs normalized edges, reroute the
original edges into exactly the old lexicographic order without splitting or
adding their signed flows. Reordering the contact queries could change which
non-boundary witness wins, so keep the old order. Both containment directions,
vertex and edge midpoint contact tests, nearest-parent selection, and
coincident-ring multiplicity reduction still run.

### Generated triangles

Boundary construction removes consecutive duplicates, normalizes closure and
adds the cyclic edges of a single ring. When it has exactly three normalized
vertices and three closed edges, **one nonzero exact orientation** proves
three distinct vertices and a noncollinear triangle. Each pair of edges is
incident at one endpoint and has distinct supports. No vertex lies within
an edge, no geometric subedge can cancel, and the signed area is nonzero.
Winding is immaterial to this certificate: either direction has the same
boundary, and the direct residual path keeps the same signed flows and
lexicographic edge order as the general path.

The check is independent of traversal metadata, edge provenance and the
successor rule. It does not certify orientation roles, whole-output validity
or uniqueness **across** rings. Oppositely oriented copies of the same
triangle still go through the coincident-ring reduction and cancel.
Collinear triples and all other generated walks still follow the old path;
neither orientation analysis nor parent selection changes.

## What a terminal proof would need

For generated walks of arbitrary size, sufficient geometric invariants are:

1. **Full noding** of the final candidate boundary arcs, including endpoints
   of partial overlaps, contacts and proper crossings. Distinct open arcs
   meet only at nodes after coincident copies are identified.
2. **One geometric owner per surviving open arc**, globally across selected
   originals and generated rings, not one ID per turn. A constructive rule
   is to retain an arc exactly when the Boolean fill differs on its two
   sides, use the same filled-side convention for its direction, and emit
   it once. Shared operands with equal Boolean fill on both sides emit
   nothing. This rules out both reverse pairs and duplicates, even with
   different source subdivisions.
3. **Jordan cycles**: follow the same filled sector at each node. Split or
   reject a walk that revisits a geometric vertex other than its closure;
   require at least three distinct vertices and nonzero edges. With (1)
   and (2), a closed cycle satisfying this condition has no self-intersection
   or overlapping edges and encloses positive area. Separate cycles may
   touch at nodes; that is not a proof of parentage.

The actual traverser closes on the start node, rejects reusing a visited
turn-operation ID, uses biconnected components and local target selection,
removes local duplicates/spikes/collinearity, and checks the ring's point
count. Marking other operations that share a segment ID and target node is
not equivalent to marking a geometric arc across both operands. These
mechanical checks do not establish (2) or (3), and one noncollinear corner
does not establish nonzero signed area of a nonsimple walk.

## What is not certified

Two independently valid, coincident operand rings disprove global uniqueness
from input validity alone. A regression feeds two such selected originals into
parenting with opposite directions and verifies they cancel.

Closed walks in a biconnected graph with distinct visited arcs need not be
simple. For example, take four internally disjoint paths between U and V.
Start at an interior vertex of the first path, go to V, follow the next three
paths back and forth, and return through the unused portion of the first path.
Every arc can have a different turn-operation ID, yet U and V occur twice.
This is a counterexample to the proposed bookkeeping proof, not a claim the
overlay's angle policy emits this particular walk.

The tests add two concrete **synthetic terminal** witnesses, not reports of
valid-input traverser failures. Two equal squares with opposite windings
touch at a repeated vertex: the walk is fully noded, has distinct geometric
edges and a unique noncollinear lexicographic extremum, yet its signed area
is zero. Another walk connects two lobes by going along an edge and
returning much later along a differently subdivided copy. Local spike
cleanup cannot account for that cancellation: the midpoint is on the raw
boundary but not the reduced boundary. Both stay on the general path.

The public exact intersection test clips two rationally validated operand
triangles that touch at their tips. Neither output can be an untouched
source ring. The retained result has two valid triangles, a shared tip and
exact total area one, and agrees with the expected set by DE-9IM.

An instrumented audit observed the walks immediately after traversal and
before selection and parenting. Overlapping cell and triangle masks, a
checkerboard, tangent holes, and issue replays were checked after certifying
the operands with rational validity predicates. In 128,084 emitted-walk
observations, it found no non-simple ring (after correcting orientation) and
no positive-length collinear edge overlap within or between walks of the same
overlay invocation. This finite audit does not establish the invariants for
all generated walks. That audit was not repeated here. This change proves
only normalized generated triangles; no general repair is removed.

## Local cost

Apple Clang 21, arm64, C++17, `-O2 -DNDEBUG -ffp-contract=off`, five alternating
before/after rounds against `1b74f781c`. The synthetic input is a projected
rectangle with eight distinct tangent diamond holes; inset holes serve as the
no-contact control. Rational validity, hole counts and exact output area were
checked outside the timer.

| Case | 64 edges per source ring | 128 edges per source ring |
| --- | ---: | ---: |
| Parenting alone, tangent holes | 11.08x | 27.17x |
| Full union with a disjoint nonempty polygon | 6.77x | 14.83x |
| Full intersecting difference, generated walks | 0.98x | 1.00x |

These are median paired CPU speedups on the fixture, not expected gains for
arbitrary overlays. A separate probe certified nine originals in the first
two contact cases and none in the traversed difference. The no-contact
controls ranged from 0.98x to 1.00x. These are historical stored-original
measurements, not a performance claim for the triangle check.

## Paired public triangle check

The [standalone benchmark](exact_generated_triangle_benchmark.cpp) rotates
the touching-triangle fixture so the shared tip is the lexicographic minimum
of both clipped triangles. This makes the actual public
`intersection(polygon, multipolygon, output, cartesian_exact<double>{})`
query the shared point during parenting. The control clips an ordinary convex
parallelogram with a rectangle and produces one nontriangular ring.

Separate instrumented builds counted work inside the parent selection.
The off probe withheld only the triangle candidate marker; it still counted
the two shapes that would be candidates. The timing binaries had no counters.

| Per call | Candidate shapes | Certified | Boundary queries | Raw boundary kept | Direct residual | General split |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Touching triangles, off | 2 | 0 | 2 | 0 | 0 | 3 |
| Touching triangles, on | 2 | 2 | 2 | 2 | 1 | 0 |
| Convex clip, either | 0 | 0 | 0 | 0 | 0 | 0 |

Both uninstrumented builds produced byte-identical retained rational WKT
for both cases. The harness separately checks rational validity, expected
exact area (1 and 12), and DE-9IM set equality with independent coordinates.
The timed calls emit native output and reuse its outer container.

Apple Clang 21, arm64, `-std=c++17 -O2 -DNDEBUG -ffp-contract=off`.
Off used the headers at `cc7de2684`; on used `491f34e01`. Same source,
21 paired samples of 0.15 seconds per variant and case, alternating
variant/case order. Table times are sample medians in CPU us/call; the
ratio is the median of the within-pair time ratios.

| Public intersection | Off | On | Paired on/off (p10-p90) | On faster |
| --- | ---: | ---: | ---: | ---: |
| Touching triangles | 10.453 | 8.559 | 0.818 (0.775-0.831) | 21/21 |
| Convex clip | 3.725 | 3.735 | 0.997 (0.977-1.018) | 13/21 |

A separate 13-pair repeat measured 0.817x and 0.998x respectively. This
demonstrates a local gain when generated triangles reach the boundary path,
not a general overlay speedup or a hit-rate estimate.

To build either variant from this checkout, set `HEADERS` to the relevant
Geometry `include` directory and use:

```sh
clang++ -std=c++17 -O2 -DNDEBUG -ffp-contract=off \
  -I"$HEADERS" -isystem "$BOOST_ROOT" test/strategies/exact_generated_triangle_benchmark.cpp \
  -o /tmp/bg-triangle-variant
/tmp/bg-triangle-variant
/tmp/bg-triangle-variant touching-triangles 0.15
/tmp/bg-triangle-variant convex-clip 0.15
```

Local raw samples, the counter-only header override, and the pairing script
are in `/private/tmp/bg1517-ring-proof-bench`.
