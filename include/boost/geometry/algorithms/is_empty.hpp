// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2015-2023, Oracle and/or its affiliates.

// Contributed and/or modified by Vissarion Fysikopoulos, on behalf of Oracle
// Contributed and/or modified by Menelaos Karavelas, on behalf of Oracle
// Contributed and/or modified by Adam Wulkiewicz, on behalf of Oracle

// Licensed under the Boost Software License version 1.0.
// http://www.boost.org/users/license.html

#ifndef BOOST_GEOMETRY_ALGORITHMS_IS_EMPTY_HPP
#define BOOST_GEOMETRY_ALGORITHMS_IS_EMPTY_HPP

#include <algorithm>

#include <boost/range/begin.hpp>
#include <boost/range/empty.hpp>
#include <boost/range/end.hpp>

#include <boost/geometry/algorithms/not_implemented.hpp>
#include <boost/geometry/algorithms/detail/visit.hpp>

#include <boost/geometry/core/exterior_ring.hpp>
#include <boost/geometry/core/interior_rings.hpp>
#include <boost/geometry/core/tag.hpp>
#include <boost/geometry/core/tags.hpp>
#include <boost/geometry/core/visit.hpp>

#include <boost/geometry/geometries/adapted/std_variant.hpp>
#include <boost/geometry/geometries/concepts/check.hpp>

#include <boost/geometry/util/type_traits_std.hpp>

namespace boost { namespace geometry
{


#ifndef DOXYGEN_NO_DETAIL
namespace detail
{

template <concepts::ConstDynamicGeometry DynamicGeometry>
inline bool is_empty_impl(DynamicGeometry const& dynamic);

template <concepts::ConstGeometryCollection GeometryCollection>
inline bool is_empty_impl(GeometryCollection const& collection);

template <typename Geometry>
    requires concepts::ConstPoint<Geometry>
          || concepts::ConstBox<Geometry>
          || concepts::ConstSegment<Geometry>
inline bool is_empty_impl(Geometry const&)
{
    return false;
}

template <typename Geometry>
    requires concepts::ConstLinestring<Geometry>
          || concepts::ConstRing<Geometry>
          || concepts::ConstMultiPoint<Geometry>
inline bool is_empty_impl(Geometry const& geometry)
{
    return boost::empty(geometry);
}

template <concepts::ConstPolygon Polygon>
inline bool is_empty_impl(Polygon const& polygon)
{
    auto const& rings = interior_rings(polygon);
    return boost::empty(exterior_ring(polygon))
        && std::all_of(boost::begin(rings), boost::end(rings),
            [](auto const& ring) { return boost::empty(ring); });
}

template <concepts::ConstMultiPolygon MultiPolygon>
inline bool is_empty_impl(MultiPolygon const& multi)
{
    return std::all_of(boost::begin(multi), boost::end(multi),
        [](auto const& polygon) { return is_empty_impl(polygon); });
}

template <concepts::ConstMultiLinestring MultiLinestring>
inline bool is_empty_impl(MultiLinestring const& multi)
{
    return std::all_of(boost::begin(multi), boost::end(multi),
        [](auto const& linestring) { return boost::empty(linestring); });
}

template <concepts::ConstPolyhedralSurface PolyhedralSurface>
inline bool is_empty_impl(PolyhedralSurface const& surface)
{
    return std::all_of(boost::begin(surface), boost::end(surface),
        [](auto const& face) { return boost::empty(face); });
}

template <concepts::ConstDynamicGeometry DynamicGeometry>
inline bool is_empty_impl(DynamicGeometry const& dynamic)
{
    bool result = true;
    traits::visit<DynamicGeometry>::apply([&](auto const& geometry)
    {
        result = is_empty_impl(geometry);
    }, dynamic);
    return result;
}

template <concepts::ConstGeometryCollection GeometryCollection>
inline bool is_empty_impl(GeometryCollection const& collection)
{
    bool result = true;
    detail::visit_breadth_first([&](auto const& geometry)
    {
        result = is_empty_impl(geometry);
        return result;
    }, collection);
    return result;
}

} // namespace detail
#endif // DOXYGEN_NO_DETAIL


/*!
\brief \brief_check{is the empty set}
\ingroup is_empty
\tparam Geometry \tparam_geometry
\param geometry \param_geometry
\return \return_check{is the empty set}

\qbk{[include reference/algorithms/is_empty.qbk]}
*/
template <concepts::ConstGeometry Geometry>
inline bool is_empty(Geometry const& geometry)
{
    return detail::is_empty_impl(geometry);
}


}} // namespace boost::geometry


#endif // BOOST_GEOMETRY_ALGORITHMS_IS_EMPTY_HPP
