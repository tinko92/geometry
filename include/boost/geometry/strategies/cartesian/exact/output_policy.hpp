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
// DeterminantRetry selectively tightens the determinant bounds on a miss.
// Without it, Retry enables just the FMA division-residual retry.
template <std::size_t N, bool Retry = false, bool DeterminantRetry = (N <= 4)>
struct bounded_ulp_output
{
    static_assert(N >= 1, "bounded_ulp_output requires at least one ULP");
    static constexpr std::size_t steps = N;
    static constexpr bool retry = Retry;
    static constexpr bool determinant_retry = DeterminantRetry;
};

// Same construction and emission contract. Root-bound points retain no
// coordinate error fields; predicates use a bound specialized on N.
template <std::size_t N, bool Retry = false, bool DeterminantRetry = (N <= 4)>
struct lean_bounded_ulp_output : bounded_ulp_output<N, Retry, DeterminantRetry> {};

}}}} // namespace boost::geometry::strategies::relate

namespace boost { namespace geometry { namespace strategies { namespace detail
{ namespace exact_cartesian
{

template <typename Policy> struct is_bounded_ulp_output : std::false_type {};
template <std::size_t N, bool R, bool D>
struct is_bounded_ulp_output<relate::bounded_ulp_output<N, R, D>> : std::true_type {};
template <std::size_t N, bool R, bool D>
struct is_bounded_ulp_output<relate::lean_bounded_ulp_output<N, R, D>> : std::true_type {};

}}}}}

#endif
