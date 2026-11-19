// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.

#ifndef BOOST_GEOMETRY_STRATEGIES_CARTESIAN_EXACT_COORDINATE_PREDICATES_HPP
#define BOOST_GEOMETRY_STRATEGIES_CARTESIAN_EXACT_COORDINATE_PREDICATES_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

#include <boost/geometry/util/exact_coordinate.hpp>
#include <boost/geometry/strategies/cartesian/exact/forward_error.hpp>

namespace boost { namespace geometry { namespace strategies { namespace detail
{ namespace exact_cartesian
{

template <typename S, typename SC, typename T, typename TC>
inline filter_result compare(exact_coordinate<S, SC> const& a, exact_coordinate<T, TC> const& b)
{
    // Rounding is monotone, including at overflow and signed zero. Distinct
    // rounded values certify strict order, but equal views may hide a gap.
    if (a.binary64() < b.binary64()) { return {-1, true}; }
    if (a.binary64() > b.binary64()) { return {1, true}; }
    return {0, a.binary64_is_exact() && b.binary64_is_exact()};
}

template <typename S, typename SC, typename T, typename TC>
inline int compare_coordinates(exact_coordinate<S, SC> const& a, exact_coordinate<T, TC> const& b)
{
    auto const filtered = compare(a, b);
    if (filtered.certain) { return filtered.value; }
    auto const& x = a.rational();
    auto const& y = b.rational();
    return (x > y) - (x < y);
}

namespace coordinate_predicates_detail
{

template <typename S, typename C>
inline bool admit(exact_coordinate<S, C> const& coordinate)
{
    double const magnitude = std::abs(coordinate.binary64());
    // In this band the rounded views and their spacing bounds lie on the
    // 2^-180 grid. Nonzero degree-two terms and their error bounds are normal.
    return magnitude == 0 ? coordinate.binary64_is_exact()
        : magnitude >= 0x1p-128 && magnitude <= 0x1p128;
}

template <typename S, typename C>
inline interval enclosure(exact_coordinate<S, C> const& coordinate)
{
    double const value = coordinate.binary64();
    if (coordinate.binary64_is_exact()) { return {value, 0}; }
    if (!std::isfinite(value) || std::abs(value) == (std::numeric_limits<double>::max)())
    {
        return interval::uncertain();
    }
    // The exact value is between the two neighbours. Their distances from
    // value are exactly representable, including around zero. Use a whole
    // spacing so the bound itself cannot underflow at a subnormal tie.
    double const upper = arithmetic_detail::upward(value);
    double const lower = arithmetic_detail::downward(value);
    return {value, (std::max)(upper - value, value - lower)};
}

template <typename X, typename Y>
inline bool approximate(X const& x, Y const& y, forward_error::point_approximation& p)
{
    if (!admit(x) || !admit(y)) { return false; }
    auto const a = enclosure(x), b = enclosure(y);
    p = {{{a.approximate, b.approximate}}, {{a.error, b.error}}};
    return true;
}

} // namespace coordinate_predicates_detail

template <typename S, typename SC, typename T, typename TC, typename U, typename UC>
inline filter_result orientation_filter(std::array<exact_coordinate<S, SC>, 2> const& a,
                                        std::array<exact_coordinate<T, TC>, 2> const& b,
                                        std::array<exact_coordinate<U, UC>, 2> const& c)
{
    using coordinate_predicates_detail::approximate;
    forward_error::point_approximation x, y, z;
    int result;
    if (approximate(a[0], a[1], x) && approximate(b[0], b[1], y)
        && approximate(c[0], c[1], z) && forward_error::try_orientation(x, y, z, result))
    { return {result, true}; }
    return {0, false};
}

template <typename S, typename SC, typename T, typename TC, typename U, typename UC>
inline int orientation(std::array<exact_coordinate<S, SC>, 2> const& a,
                       std::array<exact_coordinate<T, TC>, 2> const& b,
                       std::array<exact_coordinate<U, UC>, 2> const& c)
{
    auto const filtered = orientation_filter(a, b, c);
    if (filtered.certain) { return filtered.value; }
    auto const& ax = a[0].rational();
    auto const& ay = a[1].rational();
    exact_number const dx1 = b[0].rational() - ax;
    exact_number const dy1 = b[1].rational() - ay;
    exact_number const dx2 = c[0].rational() - ax;
    exact_number const dy2 = c[1].rational() - ay;
    exact_number const determinant = dx1 * dy2 - dy1 * dx2;
    return (determinant > 0) - (determinant < 0);
}

}}}}} // namespace boost::geometry::strategies::detail::exact_cartesian

#endif
