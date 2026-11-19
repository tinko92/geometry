# Exact Borrowed Cartesian Capacity

## Theorem and Domain

For native-input Root-bound overlay/relate, every retained geometry point is a
native vertex or an intersection of two original native edge supports. Under
the support, source-association and lifetime invariants below, any number of
composed stages preserves the same fixed capacities. A temporary query is one
midpoint or coordinate splice of geometry operands, never a new edge support.
This is numerical construction closure, not global algorithmic correctness.

Inputs are immutable finite Cartesian native integers up to 64 bits, IEEE
binary32 or binary64, losslessly admitted into the selected coordinate domain.
Arbitrary rational input/reimport, wider formats and recursively constructed
queries are excluded. Floating filters and output rounding retain their own
execution contracts; integer-capacity bounds do not prove those contracts.

## Construction Bounds

Write native coordinates as integer grid values divided by `Q=2^s`, with grid
magnitudes `<2^b`. All bit counts are magnitude bits; sign is stored separately.
A format `(N,D)` means `|X|,|Y|<2^N`, `0<W<2^D` for `(x,y)=(X/W,Y/W)`.
The existing `native_fixed_traits<C>` supplies `b,s`: integral types have
`s=0, b=digits+is_signed`; binary32 has `(b,s)=(277,149)` and binary64
`(2098,1074)`.

For an original grid edge, let

```text
U = ay-by, V = bx-ax, C = ax*by-ay*bx
|U|,|V| < 2^(b+1), |C| < 2^(2b+1)
line: Q*U*x + Q*V*y + C = 0

X = V1*C2-C1*V2
Y = C1*U2-U1*C2
W = Q*(U1*V2-V1*U2)
n = 3b+3, d = 2b+s+3
```

Each numerator subtracts two products bounded by `2^(3b+2)`; the denominator
subtracts two products bounded by `2^(2b+2)` and shifts by `s`. These bound the
actual raw formula in
[`direct_native_homogeneous`](../../include/boost/geometry/strategies/cartesian/exact/fixed_arithmetic.hpp),
not just its reduced result. Removing a common binary input scale only reduces
the bounds. Sign normalization and division by `gcd(|X|,|Y|,W)` cannot enlarge
them. `W=0` is rejected: parallel/collinear/degenerate segment cases select
existing endpoints instead of constructing a unique crossing.

For Root operands `P1,P2` in `(n,d)`, the actual query formulas are

```text
midpoint: (X1*W2+X2*W1, Y1*W2+Y2*W1, 2*W1*W2)
          capacity (n+d+1, 2d+1)
splice:   (selected(P1)*W2, selected(P2)*W1, W1*W2)
          capacity (n+d, 2d)
```

The carry bits bound raw intermediates; subsequent primitive normalization only
reduces them. For a native axis bound `p/Q`, substituting into an original line
and cancelling `Q` first gives vertical `(V*p,-(U*p+C),Q*V)`, or horizontal
`(-(V*p+C),U*p,Q*U)`, hence `(2b+2,b+s+1)`. This smaller bound applies to that
cancelled formula, not a generic four-point intersection's intermediate work.
Production [`borrowed_clipping`](../../include/boost/geometry/strategies/cartesian/exact/borrowed_clipping.hpp)
instead retains an original segment support and original native Box-face
support, using the ordinary `(n,d)` Root kernel. Synthetic metadata corners do
not acquire original-support provenance by fitting either budget.

| Family | Formula `(N,D)` | 64-bit integral | binary32 | binary64 |
| --- | --- | ---: | ---: | ---: |
| Native vertex | `(b,s+1)` | `(64,1)` | `(277,150)` | `(2098,1075)` |
| Original-support intersection | `(3b+3,2b+s+3)` | `(195,131)` | `(834,706)` | `(6297,5273)` |
| Midpoint query / auxiliary storage | `(n+d+1,2d+1)` | `(327,263)` | `(1541,1413)` | `(11571,10547)` |
| Splice query | `(n+d,2d)` | `(326,262)` | `(1540,1412)` | `(11570,10546)` |
| Cancelled native-axis query | `(2b+2,b+s+1)` | `(130,65)` | `(556,427)` | `(4198,3173)` |

These are the descriptors in
[`construction_capacity.hpp`](../../include/boost/geometry/strategies/cartesian/exact/construction_capacity.hpp).
Smaller integral widths use the same formulas with their own `b`.

## Exact-Native Query Storage

The 24-byte query still contains only a typed handle and two doubles. An empty
handle now means both doubles are certified exact, not merely approximate.
Scalar extraction accepts an original native coordinate or a producing
support's certified constant axis. Round-to-nearest and controlled contraction
are execution preconditions, checked only by debug assertions; release factories
do not inspect the floating environment. The native conversion check must report
zero error; it checks wide integer magnitudes before any backconversion, so
rounded int64/uint64 values cannot pass.

Splicing these scalars requires no arithmetic. A midpoint uses error-free
`two_sum` and requires a zero residual, then verifies exact halving by doubling
back. This rejects lost subnormal half-bits. If an operand exceeds half the
largest finite double, the potentially overflowing unscaled sum is avoided:
each operand must halve exactly and the scaled `two_sum` must have zero
residual. Finite extreme same-sign midpoints are therefore accepted without
overflow, while inexact binary64 half-subnormals use the unchanged fixed/cache
fallback. Binary32 half-subnormals are exactly representable in binary64.
The certificate also assumes the existing IEEE gradual-underflow/no-FTZ
contract. Public factory calls under other rounding modes are outside scope.

The certified scalars lie on the original grid. A native midpoint has a
representation `(ax+bx,ay+by,2Q)` in `(b+1,s+2)`; a native splice fits `(b,s+1)`.
Both fit even the smaller Root splice budget. Thus an exact-double pair can
reconstruct its primitive dyadic homogeneous value on demand within
`(b+1,s+2)`, checked by `try_fixed_input_homogeneous`, before promotion to the
unchanged auxiliary return type. No fixed arithmetic or cache allocation is
needed during certified construction. This storage choice neither narrows
predicate scratch nor admits query operands, arbitrary exact doubles, or
query-to-geometry conversion. Uncertified cases retain the normalized,
family-checked fixed construction and cached representation.
For the certified no-cache pair, `approximation` and `refine` return zero error
because both stored coordinates are exact. Cached queries retain their computed
one-ULP error. This distinction uses the existing handle, not extra error fields.

After the native certificate fails, both geometry operands are checked for an
exact `(128,128)` representation, using the borrowed-point compact accessor or
native dyadic decomposition. Success selects `fixed_midpoint`'s `(257,257)` or
`fixed_splice`'s `(256,256)` expression format. Failure retains the full-domain
geometry path. Both paths normalize their original expression, check X/Y/W
against the semantic family with `fixed_fits`, and round in that expression's
declared format before checked promotion into the single auxiliary cache entry.
No full family-sized temporary is needed merely to check capacity or round a
small rational query. This changes work/storage scheduling, not the formula,
family theorem, exact value, or one-level restriction.

`copy_compact_homogeneous<N,D>` copies cached active limbs directly through the
typed handle, or reconstructs an exact-native pair directly into the requested
format. It adds no fields and never constructs the full auxiliary value first.
A successful copy certifies the requested operand bounds; predicate scratch
must follow that actual format. An undersized format returns false without
changing the output, allowing the ordinary full-width fallback.

## Typed and Algorithm Invariants

1. [`borrowed_point<Root>`](../../include/boost/geometry/strategies/cartesian/exact/borrowed_records.hpp)
   holds original native accessors or retained record references, not an
   auxiliary handle or query. Each constructed record resolves two original
   supports to native endpoints before calling the direct native kernel.
2. Every output edge, including closing edges, retains its original support
   label through copying, reversal, traversal, simplification and composition.
   Boundary flattening includes empty ring slots; local reverse flags map via
   the model operation's source index. Endpoint association must match the
   actual input range, as checked by the four debug intersection assertions.
3. Endpoint/overlap cases copy stable references resolved through Stage IDs;
   proper crossings use one numerical seed. Turn insertion rebinds that seed
   to its retained owner, and finalization maps supports to the original domain
   before temporary stages are released. Predicate-only `at()` prvalues bound
   to local const references live through `apply`; their addresses are not
   retained by endpoint or crossing construction.
4. [`query_point<Root>`](../../include/boost/geometry/strategies/cartesian/exact/borrowed_query.hpp)
   has private construction. Factories accept only same-Root geometry references
   or eligible native points and reject query operands. The fallback normalizes
   and checks the family bound before caching; certified native queries use the
   proved no-cache representation above. Queries cannot be borrowed as geometry or used
   as output-edge supports. Both query kinds return the actual wider auxiliary
   homogeneous type. No unrestricted rational or recursive-midpoint closure is
   claimed; normalization alone would not establish it.
5. [`root_bounds`](../../include/boost/geometry/strategies/cartesian/exact/borrowed_bounds.hpp)
   selects four geometry coordinate origins; it constructs no geometry corner.
   Search/R-tree boxes are metadata, not original native sources. Working Box
   values are rejected by
   [`borrowed_working`](../../include/boost/geometry/strategies/cartesian/exact/borrowed_geometry.hpp);
   lazy standalone Box output and original native Box-face clipping remain
   separate supported paths.
6. Original sources outlive their accessors. The operation context outlives
   bound strategies, retained records and query handles; its precision scope
   remains active through nested callbacks and retained-owner destruction.
   Native/rational export completes before context exit. Public and nested
   calls preserve the bound strategy before helper/working-type selection;
   see [`borrowed_strategy.hpp`](../../include/boost/geometry/strategies/cartesian/exact/borrowed_strategy.hpp).

These invariants give induction over stages: copying preserves the two geometry
families, and a new crossing reconstructs from original native supports, not
lines inferred from arbitrary constructed points. Stage count therefore does
not increase `n,d`. Fixed integers and predicate temporaries have compile-time
storage, without arbitrary-bigint numerical fallback. Operation cache/container
growth may allocate; diagnostic exact-number getters and requested rational
export are separate boundaries.

## Predicate Scratch

For actual operand formats `(nA,dA)`, `(nB,dB)`, `(nC,dC)`, define

```text
KAB = max(nA+dB,nB+dA)+1
KBC = max(nB+dC,nC+dB)+1
coordinate comparison: KAB
orientation: max(nA+nB+dC,nA+dB+nC,dA+nB+nC)+3
direction: KAB+KBC+1
```

The orientation expansion has six signed triple products, requiring three
carry bits. Direction sums two products of coordinate-difference numerators.
Existing [`fixed_measures.hpp`](../../include/boost/geometry/strategies/cartesian/exact/fixed_measures.hpp)
supplies distance/determinant pair budgets; these are finite scratch expressions,
not new point families. Scratch must follow actual template types, including
auxiliary widening, rather than narrower semantic origins. In particular,
int64 Root/native-double metadata comparison requires **2230** bits, and
Root/Root/native-double orientation **2427**, not int64-only scratch. Metadata
is never cast to integral Root coordinates; infinite bounds are classified
before fixed conversion. Segment-box corner splices remain stack predicate
scratch and cannot become geometry supports.

## Executable Checks

The compact predicate stage first checks every homogeneous component against
`(128,128)`. Only successful exact copies use scratch derived from that format;
otherwise the existing full-format predicate runs. This is not a restriction
of the coordinate domain or a rounding of a large operand into 128 bits.
The same rule applies to compact midpoint/splice expressions before cache
promotion. Numerical storage and scratch remain fixed-capacity.

Squared-distance comparisons additionally admit exact native coordinates in
`[-2^30,2^30]`: differences fit signed 64 bits with magnitude at most `2^31`,
and their squared sum is at most `2^63`, fitting unsigned 64 bits. A coordinate
is accepted only if it is an original native coordinate, or constant on an
original producing support, and converts back exactly after the bounded integer
conversion. No rounded implicit coordinate supplies this certificate.

Distinct retained records can certify equal points from the same pair of
original supports, in either order. A constructed record represents a unique
finite crossing, so equality of those two lines determines equality of the
point. Context identity is checked first; local IDs are compared directly only
within one stage or after both have been finalized. Other cases resolve IDs
through their stages. Unequal IDs do not certify unequal points: arithmetic
fallback still handles coincident crossings with different producers.

Ordinary construction now uses the previously proved a-priori forward bound,
with the same initial approximation. A wider valid bound can cause more exact
fallbacks but cannot accept an incorrect predicate sign. Lean construction
retains its tighter residual calculation and ULP/nearest certificates.

- [`exact_construction_capacity.cpp`](exact_construction_capacity.cpp): native
  widths, formulas, actual return types, normalization, direct-axis equivalence,
  recursive-midpoint nonclosure, and allocation checks on numerical kernels.
- [`exact_borrowed_query.cpp`](exact_borrowed_query.cpp): three policies for
  int64/double, query values/layout, recursive-operand and geometry-conversion
  barriers, mixed-Root rejection, independent rational-oracle extreme/native
  cases (including float and wide integers), constant-axis origins, subnormal
  halving, exact non-dyadic compact construction, and zero allocations without
  cache reservation on certified paths (or after reservation for cached paths).
- [`exact_borrowed_capacity.cpp`](exact_borrowed_capacity.cpp): multiring source
  mapping, closing/open support transfer, helper/seed types, metadata/work-Box
  rejection modes, lazy native Box output and original Box clipping.
- [`exact_borrowed_bounds.cpp`](exact_borrowed_bounds.cpp),
  [`exact_borrowed_segment_box.cpp`](exact_borrowed_segment_box.cpp),
  [`exact_borrowed_root_ring.cpp`](exact_borrowed_root_ring.cpp) and
  [`exact_borrowed_reverse_dispatch.cpp`](exact_borrowed_reverse_dispatch.cpp):
  metadata/mixed-grid scratch, direct segment envelopes, exact separation,
  native/Root parenting and reversed source association with retained supports.

## Separate Algorithmic Limitation

The mixed-GC routing fixture in
[`exact_borrowed_gc.cpp`](exact_borrowed_gc.cpp) compares the current baseline,
not a mathematical topology oracle. Native, the previous rational baseline and
Root agree on `1F2F01102` where `1F2001102` is expected, missing B/I at `(2,1)`.
This known algorithmic matrix gap is separate from numerical capacity closure.
