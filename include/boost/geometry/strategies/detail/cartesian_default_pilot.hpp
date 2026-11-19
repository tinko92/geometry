// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.

#ifndef BOOST_GEOMETRY_STRATEGIES_DETAIL_CARTESIAN_DEFAULT_PILOT_HPP
#define BOOST_GEOMETRY_STRATEGIES_DETAIL_CARTESIAN_DEFAULT_PILOT_HPP

#if defined(BOOST_GEOMETRY_ENABLE_CARTESIAN_EXACT_DEFAULT_PILOT)
#include <boost/config.hpp>
#include <cfloat>

// The build macro attests to strict FP (including disabled implicit
// contraction). Reject configurations we can identify as incompatible
// before any legacy header can pull in the exact implementation.
#if BOOST_CXX_VERSION >= 201703L && !defined(BOOST_NO_CXX17_IF_CONSTEXPR) \
    && !defined(__FAST_MATH__) && !defined(_M_FP_FAST) \
    && !(defined(__FINITE_MATH_ONLY__) && __FINITE_MATH_ONLY__ > 0) \
    && defined(FLT_EVAL_METHOD) && FLT_EVAL_METHOD == 0 && FLT_RADIX == 2 \
    && DBL_MANT_DIG == 53 && DBL_MIN_EXP == -1021 && DBL_MAX_EXP == 1024
#define BOOST_GEOMETRY_DETAIL_CARTESIAN_EXACT_DEFAULT_PILOT

#include <cfenv>
#include <cstdint>
#include <cstring>
#include <limits>
#include <type_traits>

#include <boost/geometry/algorithms/for_each.hpp>
#include <boost/geometry/geometries/multi_polygon.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/strategies/relate/cartesian_exact.hpp>

namespace boost { namespace geometry { namespace detail
{
namespace cartesian_default_pilot
{

using polygon = model::polygon<model::d2::point_xy<double>>;
using polygons = model::multi_polygon<polygon>;
using strategy = strategies::relate::cartesian_exact<double, strategies::relate::approximate_output>;

// Match the original public types without instantiating geometry traits for
// other models, variant alternatives, or tuple output.
template <typename A, typename B>
constexpr bool admits_pair = std::is_same<A, polygon>::value && std::is_same<B, polygon>::value;

template <typename A, typename B, typename Out>
constexpr bool admits_output = admits_pair<A, B> && std::is_same<Out, polygons>::value;

inline bool admits_coordinate(double value)
{
    // IEEE binary64: |value| <= 2^16 and an integer multiple of 2^-64.
    // Inspect bits without changing exceptions, even on rejected values.
    std::uint64_t bits;
    std::memcpy(&bits, &value, sizeof(bits));
    constexpr int bias = 1023, fraction_bits = 52;
    bits &= (std::uint64_t{1} << 63) - 1;
    int const exponent = static_cast<int>(bits >> fraction_bits) - bias;
    if (exponent >= 16)
        return bits == (std::uint64_t{bias + 16} << fraction_bits);
    if (exponent < -64)
        return bits == 0;
    int const discarded = fraction_bits - 64 - exponent;
    return discarded <= 0 || (bits & ((std::uint64_t{1} << discarded) - 1)) == 0;
}

inline bool admits_values(polygon const& a, polygon const& b)
{
    if (!std::numeric_limits<double>::is_iec559 || std::fegetround() != FE_TONEAREST)
        return false;
    auto const admits_point = [](auto const& p)
    {
        return admits_coordinate(geometry::get<0>(p)) && admits_coordinate(geometry::get<1>(p));
    };
    return geometry::all_points_of(a, admits_point) && geometry::all_points_of(b, admits_point);
}

} // namespace cartesian_default_pilot
}}} // namespace boost::geometry::detail

#endif
#endif
#endif // BOOST_GEOMETRY_STRATEGIES_DETAIL_CARTESIAN_DEFAULT_PILOT_HPP
