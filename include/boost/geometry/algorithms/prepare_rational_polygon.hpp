// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.

#ifndef BOOST_GEOMETRY_ALGORITHMS_PREPARE_RATIONAL_POLYGON_HPP
#define BOOST_GEOMETRY_ALGORITHMS_PREPARE_RATIONAL_POLYGON_HPP

#include <type_traits>

#include <boost/range/size.hpp>
#include <boost/range/value_type.hpp>

#include <boost/geometry/core/access.hpp>
#include <boost/geometry/core/closure.hpp>
#include <boost/geometry/core/coordinate_dimension.hpp>
#include <boost/geometry/core/coordinate_system.hpp>
#include <boost/geometry/core/coordinate_type.hpp>
#include <boost/geometry/core/cs.hpp>
#include <boost/geometry/core/exterior_ring.hpp>
#include <boost/geometry/core/interior_rings.hpp>
#include <boost/geometry/core/point_order.hpp>
#include <boost/geometry/geometries/multi_polygon.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/util/exact_coordinate.hpp>
#include <boost/geometry/util/type_traits.hpp>

namespace boost { namespace geometry
{

namespace detail { namespace rational_polygon
{

template <typename Coordinate>
struct is_wrapped : std::false_type {};

template <typename Scalar, typename Conversion>
struct is_wrapped<exact_coordinate<Scalar, Conversion>> : std::true_type {};

using point = model::point<boost::multiprecision::cpp_rational, 2, cs::cartesian>;

template <typename Polygon>
using polygon = model::polygon<point,
    geometry::point_order<Polygon>::value == clockwise,
    geometry::closure<Polygon>::value == closed>;

template <typename Source, typename Destination>
void copy_ring(Source const& source, Destination& destination)
{
    destination.reserve(boost::size(source));
    for (auto const& p : source)
    {
        destination.emplace_back(geometry::get<0>(p).rational(),
                                 geometry::get<1>(p).rational());
    }
}

}} // namespace detail::rational_polygon

// Explicit operation boundary for two-dimensional wrapped Cartesian polygons.
// The resulting owning geometry can be reused with public set operations and
// relate/relation and strategies::relate::cartesian<cpp_rational>. No binary64
// view is read. Ring order, closure, and holes are preserved, not corrected.
template <typename Polygon, std::enable_if_t<util::is_polygon<Polygon>::value, int> = 0>
auto prepare_rational_polygon(Polygon const& source)
{
    static_assert(dimension<Polygon>::value == 2
                  && std::is_same<cs_tag_t<Polygon>, cartesian_tag>::value,
                  "Preparing rational polygons requires two-dimensional Cartesian input");
    static_assert(detail::rational_polygon::is_wrapped<coordinate_type_t<Polygon>>::value,
                  "Preparing rational polygons requires exact_coordinate input");

    detail::rational_polygon::polygon<Polygon> output;
    using detail::rational_polygon::copy_ring;
    copy_ring(exterior_ring(source), exterior_ring(output));
    auto const& source_inners = interior_rings(source);
    auto& output_inners = interior_rings(output);
    output_inners.resize(boost::size(source_inners));
    auto output_ring = output_inners.begin();
    for (auto const& ring : source_inners)
    {
        copy_ring(ring, *output_ring++);
    }
    return output;
}

template <typename MultiPolygon, std::enable_if_t<util::is_multi_polygon<MultiPolygon>::value, int> = 0>
auto prepare_rational_polygon(MultiPolygon const& source)
{
    using polygon = typename boost::range_value<MultiPolygon>::type;
    model::multi_polygon<detail::rational_polygon::polygon<polygon>> output;
    output.reserve(boost::size(source));
    for (auto const& part : source)
    {
        output.push_back(prepare_rational_polygon(part));
    }
    return output;
}

}} // namespace boost::geometry

#endif
