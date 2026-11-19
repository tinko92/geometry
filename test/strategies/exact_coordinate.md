# Exact coordinate storage

`boost/geometry/util/exact_coordinate.hpp` provides `boost::geometry::exact_coordinate<Scalar>`.
`scalar()` exposes the owned source by const reference. `rational()` converts
that source without rounding. `binary64()` returns a cached nearest-even
view, and `binary64_is_exact()` says whether the view equals the source.
The binary64 view may be infinite or signed zero while the scalar remains
finite and nonzero. There is no implicit conversion to `double`.

The supplied `traits::exact_coordinate_conversion` specializations accept
`boost::multiprecision::cpp_rational` and `std::string`. Text is kept verbatim
and interpreted as an exact base-ten rational, not as the rational value of
a parsed `double`. For example, `"0.1"` and `cpp_rational(0.1)` have the same
binary64 view but different rational values. The decimal grammar and size
limits are documented on the specialization. Wrapper replacement validates
and rounds before changing the object. Copying from an rvalue leaves the
source usable, too.

## Predicates and limits

`strategies/cartesian/exact/coordinate_predicates.hpp` provides standalone
coordinate comparison and 2D orientation over arrays of wrapped coordinates.
The `compare` overload certifies order from distinct rounded values, and
equality only when both views are exact. `orientation_filter` uses outward
intervals under the existing strict floating-point contract (binary64,
gradual underflow, no fast-math or implicit contraction). `compare_coordinates`
and `orientation` fall back to rational arithmetic when the filter cannot
decide. The operands may use different wrapper specializations.

The wrapper is not a general arithmetic type. `model::point` storage and
`get`/`set` work, but the standalone predicates are not selected for raw
wrappers by BG algorithm or strategy dispatch. See
[decimal polygon overlay](exact_decimal_overlay.md) for the explicit
unfiltered and filtered rational boundaries and closed polygon WKT input.
Neither is direct wrapper support in the native borrowed kernel.
