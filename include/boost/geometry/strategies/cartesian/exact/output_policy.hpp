// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.

#ifndef BOOST_GEOMETRY_STRATEGIES_CARTESIAN_EXACT_OUTPUT_POLICY_HPP
#define BOOST_GEOMETRY_STRATEGIES_CARTESIAN_EXACT_OUTPUT_POLICY_HPP

#include <cstddef>
#include <type_traits>

namespace boost { namespace geometry { namespace strategies { namespace relate
{

// Emit the approximation already stored at the output boundary. Predicates
// remain exact, but emission neither refines nor certifies these coordinates.
struct approximate_output {};

// Round floating output to nearest at emission. Root-bound operations use the
// deferred one-ULP layout; legacy point compatibility keeps stored error bounds.
// Integral output follows BG numeric_cast conversion under every policy.
struct nearest_output {};

// Bound constructed binary64 output by N representable steps from nearest.
// The one-step policy selectively tightens determinant bounds on a miss.
// For larger budgets, Retry enables the exact FMA division residual stage.
template <std::size_t N, bool Retry = false>
struct bounded_ulp_output
{
    static_assert(N >= 1, "bounded_ulp_output requires at least one ULP");
    static constexpr std::size_t steps = N;
    static constexpr bool retry = Retry;
};

}}}} // namespace boost::geometry::strategies::relate

namespace boost { namespace geometry { namespace strategies { namespace detail
{ namespace exact_cartesian
{

template <typename Policy> struct is_bounded_ulp_output : std::false_type {};
template <std::size_t N, bool Retry>
struct is_bounded_ulp_output<relate::bounded_ulp_output<N, Retry>>
    : std::bool_constant<(relate::bounded_ulp_output<N, Retry>::steps > 0)> {};

}}}}}

#endif
