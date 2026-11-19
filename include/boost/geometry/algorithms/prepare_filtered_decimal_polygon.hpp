// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.

#ifndef BOOST_GEOMETRY_ALGORITHMS_PREPARE_FILTERED_DECIMAL_POLYGON_HPP
#define BOOST_GEOMETRY_ALGORITHMS_PREPARE_FILTERED_DECIMAL_POLYGON_HPP

#include <type_traits>

#include <boost/range/size.hpp>
#include <boost/range/value_type.hpp>

#include <boost/geometry/core/closure.hpp>
#include <boost/geometry/core/exterior_ring.hpp>
#include <boost/geometry/core/interior_rings.hpp>
#include <boost/geometry/core/point_order.hpp>
#include <boost/geometry/geometries/multi_polygon.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/strategies/cartesian/exact/decimal_point.hpp>
#include <boost/geometry/util/type_traits.hpp>

namespace boost { namespace geometry
{
namespace detail { namespace filtered_decimal
{

template <typename Polygon>
using polygon = model::polygon<point,
    geometry::point_order<Polygon>::value == clockwise,
    geometry::closure<Polygon>::value == closed>;

template <typename Source, typename Destination>
void copy_ring(Source const& source, Destination& destination)
{
    destination.reserve(boost::size(source));
    for (auto const& p : source)
        destination.emplace_back(geometry::get<0>(p), geometry::get<1>(p));
}

}} // namespace detail::filtered_decimal

// Opt-in preparation for cartesian_decimal. Like prepare_rational_polygon,
// this makes an owning rational copy and does not validate or correct rings.
template <typename Polygon, std::enable_if_t<util::is_polygon<Polygon>::value, int> = 0>
auto prepare_filtered_decimal_polygon(Polygon const& source)
{
    namespace decimal = detail::filtered_decimal;
    static_assert(dimension<Polygon>::value == 2
                  && std::is_same<cs_tag_t<Polygon>, cartesian_tag>::value,
                  "Preparing filtered decimal polygons requires two-dimensional Cartesian input");
    static_assert(std::is_same<coordinate_type_t<Polygon>, decimal::decimal>::value,
                  "Preparing filtered decimal polygons requires exact_coordinate<std::string> input");
    decimal::polygon<Polygon> result;
    decimal::copy_ring(exterior_ring(source), exterior_ring(result));
    auto const& inner = interior_rings(source);
    auto& output = interior_rings(result);
    output.resize(boost::size(inner));
    auto dest = output.begin();
    for (auto const& ring : inner)
        decimal::copy_ring(ring, *dest++);
    return result;
}

template <typename Multi, std::enable_if_t<util::is_multi_polygon<Multi>::value, int> = 0>
auto prepare_filtered_decimal_polygon(Multi const& source)
{
    model::multi_polygon<detail::filtered_decimal::polygon<typename boost::range_value<Multi>::type>> result;
    result.reserve(boost::size(source));
    for (auto const& part : source)
        result.push_back(prepare_filtered_decimal_polygon(part));
    return result;
}

}} // namespace boost::geometry

#endif
