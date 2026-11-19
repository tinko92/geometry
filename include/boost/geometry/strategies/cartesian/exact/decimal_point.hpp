// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.

#ifndef BOOST_GEOMETRY_STRATEGIES_CARTESIAN_EXACT_DECIMAL_POINT_HPP
#define BOOST_GEOMETRY_STRATEGIES_CARTESIAN_EXACT_DECIMAL_POINT_HPP

#include <array>
#include <cmath>
#include <cstddef>
#include <optional>
#include <string>
#include <type_traits>
#include <utility>

#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/strategies/cartesian/exact/coordinate_predicates.hpp>
#include <boost/geometry/strategies/cartesian/exact/forward_error.hpp>

namespace boost { namespace geometry
{
namespace detail { namespace filtered_decimal
{

using decimal = exact_coordinate<std::string>;
using rational = boost::multiprecision::cpp_rational;
using approximation = strategies::detail::exact_cartesian::forward_error::point_approximation;

// An owning rational point. Only preparation from the authoritative decimal
// wrapper creates a certificate; writes from the geometry engine discard it.
class point
{
public:
    point() = default;
    point(decimal const& x, decimal const& y) : m_value{{x.rational(), y.rational()}}
    {
        approximation p;
        if (strategies::detail::exact_cartesian::coordinate_predicates_detail::approximate(x, y, p))
        {
            m_approx = p;
        }
    }

    point(point const&) = default;
    point(point&& other) : point() { swap(other); }
    point& operator=(point other) noexcept
    {
        swap(other);
        return *this;
    }

    template <std::size_t D> rational const& get() const
    {
        static_assert(D < 2, "Decimal points are two-dimensional");
        return m_value[D];
    }
    template <std::size_t D> void set(rational const& value)
    {
        static_assert(D < 2, "Decimal points are two-dimensional");
        m_approx.reset();
        m_value[D] = value;
    }
    approximation const* certified_approximation() const { return m_approx ? &*m_approx : nullptr; }

private:
    void swap(point& other) noexcept
    {
        m_value[0].swap(other.m_value[0]);
        m_value[1].swap(other.m_value[1]);
        m_approx.swap(other.m_approx);
    }

    std::array<rational, 2> m_value{};
    std::optional<approximation> m_approx;
};

}} // namespace detail::filtered_decimal

namespace traits
{
template <> struct tag<detail::filtered_decimal::point> { using type = point_tag; };
template <> struct coordinate_type<detail::filtered_decimal::point>
{ using type = detail::filtered_decimal::rational; };
template <> struct coordinate_system<detail::filtered_decimal::point> { using type = cs::cartesian; };
template <> struct dimension<detail::filtered_decimal::point> : std::integral_constant<std::size_t, 2> {};
template <std::size_t D> struct access<detail::filtered_decimal::point, D>
{
    using point = detail::filtered_decimal::point;
    using rational = detail::filtered_decimal::rational;
    static rational const& get(point const& p) { return p.template get<D>(); }
    static void set(point& p, rational const& value) { p.template set<D>(value); }
};
} // namespace traits

}} // namespace boost::geometry

#endif
