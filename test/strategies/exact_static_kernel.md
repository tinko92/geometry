# Typed Exact Cartesian Kernel

Scope: Cartesian relate and straight-segment set operations on finite native
integral, float and double coordinates. Floating overflow and underflow are
excluded from the forward-filter contract. Long double is not supported.
No claim is made about the inherited metric algorithms.

## Predicates and Constructions

Write an exact point as \(P=(X_P:Y_P:W_P)\), with \(W_P>0\).
Mathematically a native input is \((x:y:1)\); its integer representation may
use a common binary scale to preserve floating-point inputs exactly.

| Operation | Mathematical decision or result |
|---|---|
| Coordinate order | \(\operatorname{sgn}(X_AW_B-X_BW_A)\), or the analogous expression for \(y\). |
| Point equality and clustering | Both coordinate-order predicates return zero. No epsilon. |
| Orientation | \(\operatorname{sgn}\det[H_A;H_B;H_P]\). |
| Continuation direction | \(\operatorname{sgn}((B-A)\cdot(P-B))\). |
| Squared-distance order | Compare \(\|B-A\|^2\) and \(\|D-C\|^2\), clearing positive denominators. |
| Segment order | Choose \(x\) when \(A_x\ne B_x\), otherwise \(y\); multiply coordinate order by \(\operatorname{sgn}(B_k-A_k)\). |
| Segment position | Endpoint equality, orientation zero, and segment-order comparisons against its two endpoints. A zero-length segment is classified separately. |
| Point on segment | Orientation zero and inclusive coordinate bounds in both dimensions. |
| Point in ring | Exact winding count; an upward crossing adds one when orientation is positive, a downward crossing subtracts one when negative. Boundary is separate. |
| Point/box and box/box relations | Coordinate-order conjunctions, with strict or inclusive bounds as required. |
| Ring orientation and degeneracy | Orientation at a lexicographic extreme. Repeated extrema and canceled excursions use directed boundary-edge cancellation, not signed area. This is not a full OGC-validity test. |
| Ring parent | Strict containment among selected boundaries; choose the innermost containing boundary. Coincident boundaries are reduced separately and cannot form parent cycles. |
| Intersection | \(H_{AB\cap CD}=(H_A\times H_B)\times(H_C\times H_D)\). Retain native supporting segments when available. |
| Coordinate splice | \((X_AW_B:Y_BW_A:W_AW_B)\). |
| Midpoint | \((X_AW_B+X_BW_A:Y_AW_B+Y_BW_A:2W_AW_B)\). Used only when containment samples lie on boundaries. |
| Segment/box clipping | Four exact half-plane tests; replace an outside endpoint by its intersection with the original supporting line and the selected box side. No clipping parameter. |
| Search enclosure | Outward-rounded coordinate bounds; false positives are permitted, false negatives are not. This is not an output geometry. |
| Output conversion | Round \(X/W,Y/W\) to the requested native type, nearest with ties to even. Rational output is an explicit compatibility boundary. |

In particular, a squared-distance value has the fixed rational form
\[
 \frac{(X_BW_A-X_AW_B)^2+(Y_BW_A-Y_AW_B)^2}{(W_AW_B)^2}.
\]
Its denominator is positive. There is no segment-ratio object.

## Filter Sequence

    if native coordinates permit a proved exact integer calculation:
        use that calculation
    evaluate the approximate polynomial and its forward error bound
    if the approximate sign exceeds the bound:
        return its sign
    refine ambiguous intersection coordinates using cached fixed-size expansions
    retry the forward filter
    evaluate a statically bounded exact expression if ambiguity remains

The common first stage does not evaluate an interval expression tree.
For native orientation, with \(u=2^{-53}\), the determinant filter uses
\(4u(|p|+|q|)\), where \(p,q\) are the two computed products.
The construction filter uses fixed final coefficients \(16u\) for exactly
converted inputs and \(32u\) for rounded integral inputs, multiplying the
proved magnitude expression. Constructed-point predicates also include the
stored coordinate errors. Conservative interval fallbacks remain secondary
stages, including non-default rounding modes.

Ambiguous distance/measure comparisons first try an exact compact tier: inspect
the active coordinate limbs and copy fitting homogeneous values into 128-bit
components. Native dyadics are assembled directly from their stored binary
significands/exponents, not through the worst-case binary64 grid. The resulting
distance numerators/denominators need 515/512 bits and their comparison 1032 bits.
Failure to fit selects the unchanged full-capacity exact path, never a rounded
answer. The compact tier adds neither stored point members nor heap allocation.
Squared-distance ties with known native integer coordinates in
\([-2^{30},2^{30}]\) are resolved even earlier, before refinement: coordinate
differences are at most \(2^{31}\), so the sum of their squares fits uint64
exactly (at most \(2^{63}\)). Integral-value/range checks are exact, not tolerances.

The operation-specific steps are:

| Family | Pseudocode |
|---|---|
| Coordinate order / equality | Compare native values directly when possible; otherwise filter the coordinate difference, refine cached coordinates if necessary, then compare the homogeneous cross-products exactly. Equality requires zero in both dimensions. |
| Orientation / continuation / squared distance | Evaluate the determinant, dot product or distance difference with its magnitude bound; refine ambiguous constructed points, then evaluate the cleared-denominator polynomial with fixed-capacity integers. |
| Along-segment order | Pick a nonconstant coordinate; compare the two points in it; reverse the sign when the segment decreases in that coordinate. A degenerate segment has one position. |
| Segment position / winding | Combine exact endpoint order and orientation tests; count signed ray crossings, returning boundary separately. |
| Intersection | Classify segment endpoints using exact orientation and order; copy an existing endpoint for touches; otherwise retain the two native supporting segments and approximate the intersection with a forward bound. Cache expanded or exact homogeneous coordinates only when needed. |
| Coordinate splice / midpoint | Form the homogeneous expression above, normalize it, and keep the checked fixed-size value in an indexed precision entry. |
| Clipping | Visit the four box half-planes; reject two outside endpoints, retain two inside endpoints, or replace the outside endpoint by a support-line/box-side intersection. Always use the original support line. |
| Search boxes | Enclose native or implicit coordinates outwards; only use these bounds to reject definitely disjoint candidates. Actual containment uses exact extrema instead. |
| Ring analysis | Try the unique-extreme orientation test; on ambiguity split coincident boundary edges at existing vertices, cancel reverse-directed excursions, and analyze the surviving boundary. No signed-area sum. |
| Ring parenting | Reduce coincident boundary multiplicities; test surviving boundaries for containment; select the innermost strict container; attach children according to effective orientation and the requested operation. |
| Output rounding | Compare the exact rational value against candidate rounding boundaries; choose the nearest representable coordinate, breaking ties to even. |

## Stored Data

Members below omit constructors, accessors and implementation helpers.
The fixed_integer types are checked, allocator-free fixed-capacity integers.

    struct coordinate_bound {
        double approximate;
        double error;
    };

    template <class C> struct exact_point {
        bool m_intersection;
        bool m_mixed, m_source_signed;
        operation_detail::format m_source_format;
        union {
            std::array<C, 8> m_values;
            std::array<std::uint64_t, 2> m_mixed_values;
        };
        std::array<interval, 2> m_bounds;
        point_precision_detail::point_handle m_precision;
        mutable point_cache_detail::native_cache_handle m_cache;
    };

    struct point_handle {
        std::size_t m_index;
        operation_detail::format m_format;
    };

    struct native_cache_handle {
        std::size_t m_index;
    };

    template <class C> struct native_precision_entry {
        typename native_fixed_traits<C>::intersection_type value;
        bool exact_ready;
        std::array<expanded_coordinate, 2> expanded;
    };

    struct expanded_coordinate {
        expansion<32> numerator;
        expansion<16> denominator;
        bool attempted;
        bool valid;
        bool bounded;
        coordinate_bound bounds;
    };

    template <class C> struct operation_point_storage {
        std::vector<native_precision_entry<C>> native;
        std::vector<auxiliary_point_type<C>> auxiliary;
        std::size_t native_base, auxiliary_base;
    };

    template <class A, class B> struct static_squared_distance_value {
        exact_point<A> a;
        exact_point<B> b;
        double approximate, error;
        bool valid, unset;
    };

    struct orientation_result { int side; };
    struct search_box_members { double min_x, min_y, max_x, max_y; };

    template <class P> struct segment_intersection_result {
        std::size_t count;
        P intersections[2];
    };

    template <class P> struct exact_ring_properties {
        bool valid;
        P point, boundary_point;
        bool clockwise;
        int direction;
        bool reversed, discarded;
        ring_identifier parent;
        std::vector<ring_identifier> children;
    };

An operation owns the cache vectors. Nested calls reuse the enclosing context;
the root scope initializes it on the first cache use and destroys every numerical
entry on return. Successful native filters do not initialize it. Only a non-owning
current-context pointer and issued-index counters are thread-local. There is no
global numerical store, lock, atomic, reference count, or free list. Point handles
are plain indices; disjoint issued-index ranges reject stale hints. Copies of
native-support points start unbound, while moves retain their cache hint.

Mixed native inputs are stored inline in their original integer/float/double
format, never rounded to the strategy's coordinate type. Auxiliary constructed
points and bounded rational compatibility literals require an explicit enclosing
operation scope if used outside the public algorithm entry points. They cannot
escape that scope or be shared across independent operations. Native-support
intersection points are self-contained and remain usable after it ends.

For a coordinate format with \(B\) integer-coordinate bits and binary scale
\(2^S\), a native intersection uses \(3B+3\) numerator bits and \(2B+S+3\)
denominator bits. Binary64 uses \(B=2098,S=1074\), hence 6297 and 5273 bits.
Writing these bounds as \(N,D\), auxiliary coordinate splices need \(N+D,2D\)
bits, and midpoint samples need \(N+D+1,2D+1\). Their operation-owned records
use the latter bounds: int64 327/263 bits, binary64 11571/10547 bits. Oversized
compatibility literals or recursive constructions are rejected, not truncated.
Predicate intermediate capacities are derived separately from their polynomial
degrees and term counts, rather than reusing a point's coordinate capacity.

Native supporting endpoint tuples are owned by the point, avoiding iterator
lifetime assumptions for retained outputs. Transient predicates borrow existing
exact points, retaining their refined cache entries instead of refining and
discarding temporary copies. Values retained beyond a predicate call, such as
squared-distance comparisons, own their point copies. No producer pool, general scalar
expression, type-erased numerical callback, or stored intersection fraction
remains. Dynamic numerical storage belongs only to these precision caches.
Algorithm containers and explicitly requested rational output can still allocate.

Search-box storage must not select the algorithm's working point type.
The umbrella supplies a separate section-point alias, preserving native or
implicit exact vertices through direction and duplicate classification.
Whole-geometry containment also requests exact extrema separately from search
enclosures. In particular, \(1\) and \(1+2^{-54}\) must remain distinct even
though both round to the same binary64 value.
