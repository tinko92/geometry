# Exact Cartesian ring parenting

## What is certified

Only a stored native operand ring selected in an areal Boolean operation is a
candidate. A source ID establishes provenance, not validity. Generated rings,
exact working intermediates, temporaries, buffer, and dissolve do not carry
this certificate. The check is lazy, after a boundary result needs resolution.

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

An instrumented audit observed the walks immediately after traversal and
before selection and parenting. Overlapping cell and triangle masks, a
checkerboard, tangent holes, and issue replays were checked after certifying
the operands with rational validity predicates. In 128,084 emitted-walk
observations, it found no non-simple ring (after correcting orientation) and
no positive-length collinear edge overlap within or between walks of the same
overlay invocation. This finite audit does not establish the invariants for
all generated walks. The implementation does not skip their cancellation.

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
controls ranged from 0.98x to 1.00x.
