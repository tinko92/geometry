# Conventional binary multiprecision buffer: missing area

The repairs are opt-in. Include
`<boost/geometry/strategies/buffer/cartesian_mp.hpp>` and pass the
named strategy after the five usual buffer strategies:

```cpp
bg::buffer(input, output, distance, side, join, end, point,
           bg::strategies::buffer::cartesian_mp{});
```

This overload supports concrete Cartesian geometries. It does not
change the seven-argument overload or specialize its default strategy,
even in a translation unit that includes the optional header.

Legacy buffer retains two optional hooks. The turn hook returns a
certificate for the unique outgoing side; the buffer algorithm still
decides which turns can be repaired and applies the operation changes.
The round hook supplies the certified fan or declines the request. Both
hooks are absent from the default umbrella. The MP type checks, rational
conversion, ray predicates and cardinal sampling live in the optional
header, not in either legacy algorithm header.

On a concrete double multipolygon buffer smoke translation unit, Clang
21/C++14 sees the same 2,954 distinct headers as the parent, with an
identical normalized include list. The original repair added one header.
Preprocessed size is 13,419,327 bytes vs 13,416,063 at the parent and
13,425,567 with the original repair. Neither the optional strategy nor
the removed repair header appears in the refactored native include graph.
Seven alternating warm `-fsyntax-only` runs on local arm64 gave
median user CPU times of 1.74 s (parent) and 1.73 s (refactor), within
run variation.
The pre-existing transitive Boost.Multiprecision includes remain;
removing those would require separate changes to unrelated strategies
and core type promotion.

Base `d4f6bab74`. Both reductions reproduce with the independent
`known-buffer/mp-reduce.cpp` on this base and on local develop
`082b99181`. This uses conventional Cartesian `cpp_bin_float_100`
buffering, not the exact buffer pilot or the represented touch-placement
experiment. No epsilon, clustering, or graph rules changed.

## Miter: a required boundary switch was erased

There are three offset rings. The relevant operation pair is on ring 1
segment 5 and ring 2 segment 0, at their exact common endpoint
`(7+sqrt(2)/2, 4-sqrt(2)/2)` (as represented). Independently taking
rational determinants of the stored endpoints agreed with the native MP
side predicates: this is **not** a side-sign error. The touch handler
emits `intersection/blocked`; `get_turns()` erases that combination.
All three rings acquire intersection flags but traversal has no full
exterior cycle, so the final result is empty.

The repair is restricted to two different nondeflated, nonconcave
offset rings, an `intersection/blocked` end/end touch, and an exact
shared represented vertex. Four nonzero distinct incident rays are
required. The exact rational directions sort by half-plane and
determinant. For a clockwise ring its interior sector runs CCW from the
backward incoming ray to the forward outgoing ray. Another ring covers
a ray's adjacent angle just when it strictly covers the ray, because
coincident rays are excluded. Recover a union/intersection switch only
when there is exactly one exposed incoming and one exposed outgoing
ray, belonging to opposite rings. All other classifications are left
alone; subsequent original/piece coverage checks still run.

## Round: distinct intersections rounded to the same cluster point

Near `(0 6)`, the original round fan evaluates analytically cardinal
samples via trigonometry. Two generated vertices had x coordinates
about -9.86e-101 and -4.14e-101, at y=6. Different incidences with the
connecting ring round to the same turn point. Both turns survive piece
checks, but their angular cluster has `open_count == 0` and is marked
nontraversable. Ranking by the actual interior edge direction instead
of the displacement from the stored turn point still left this cluster
closed. The patch does not force clusters open or conflate nearby turns.

For built-in `join_round` with `side_straight` at positive equal
width, obtain the two exact left normals (up to positive scale) from the
original adjacent edges. Admit only normals with known integer octants,
a convex clockwise span of 1-3 octants, and a subdivision that actually
samples an interior cardinal direction. If the span is s octants and
the configured full fan has P points, use n=ceil(P*s/8). A cardinal
direction d octants after the start is sampled precisely if d*n is
divisible by s and 0<d<s. Construct that sample directly from the
center and radius, changing only one coordinate. The other samples
use trigonometry and the stored tangent endpoints are unchanged.
Noncertified joins stay on the original path. The two bounded
predicates use rational arithmetic on stored binary values; the band
[2^-512,2^512] in magnitude (or zero) limits conversion size, not
geometric discrimination within the band.

## Verification

The registered `buffer_mp_missing_area` regression checks valid public
results, component counts, area, and coverage of every original vertex
with each join, using the explicit strategy for MP and the default
for double. It checks the miter area
`11.5 + 5*sqrt(2)` and the round chord-fan area 24.44698493496707
in MP. Independent rational scanline integration of the offset rings
dumped with 110-digit coordinates gives the same areas: split at all vertices' x values and all
proper cross-ring edge intersections, union vertical intervals at
each slab midpoint, then multiply by the slab width. This arithmetic
does not use buffer turns, zoning, traversal, or parenting. Even
the pre-fix round fan's reference area rounds to the same value; the
lost ~10.93639 units were not a change in tessellation area. The
registered CTest passed, as did the standalone regression on GCC 16
(C++14) and Clang ASan/UBSan.

Apple Clang 21, arm64, C++17, release: both existing ordinary double
`buffer_polygon` and `buffer_multi_polygon` passed in both winding
orders. With the older suites' additional MP type enabled (clockwise)
and `BOOST_GEOMETRY_TEST_BUFFER_CARTESIAN_MP` defined,
the polygon suite passed and the multipolygon suite went from two
failures to one, with no changed expectations. `nores_wt_1` also
recovers all input vertices in the independent public reducer on the full
input (area 80.15797331784302). That suite uses the built-in asymmetric
distance with equal widths, while the original reduction uses the
symmetric strategy. Both are admitted; unequal widths are not.

The full historic `rt_w6` is still wrong in conventional MP: expected
79.799, observed 3.17157287525381 both before and after. The repaired
three-ring reduction is not evidence that every larger arrangement has
the same failure. After the patch the full case still reduces to a
different seven-component witness with area 19.48528137423857 and
20 omitted source-vertex entries. No general overlap/cluster repair
or change to the original area assertion is proposed here.

Local test logs and phase/turn/ray traces are in
`/private/tmp/bg1517-missing-buffer.MQgv9u`. The debug instrumentation
used to collect the traces is not in the patch.
The opt-in replay and compile manifests are in
`/private/tmp/bg1517-mp-strategy.FeC9m7`.

Recommendation: adopt the explicit strategy and the lightweight hooks,
not automatic MP dispatch. The repairs are useful on the proved cases,
but the remaining full miter failure precludes claiming complete
conventional MP buffering.
