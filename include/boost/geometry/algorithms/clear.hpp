// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2007-2012 Barend Gehrels, Amsterdam, the Netherlands.
// Copyright (c) 2008-2012 Bruno Lalande, Paris, France.
// Copyright (c) 2009-2012 Mateusz Loskot, London, UK.
// Copyright (c) 2024 Adam Wulkiewicz, Lodz, Poland.

// This file was modified by Oracle on 2020-2023.
// Modifications copyright (c) 2020-2023, Oracle and/or its affiliates.
// Contributed and/or modified by Vissarion Fysikopoulos, on behalf of Oracle
// Contributed and/or modified by Adam Wulkiewicz, on behalf of Oracle

// Parts of Boost.Geometry are redesigned from Geodan's Geographic Library
// (geolib/GGL), copyright (c) 1995-2010 Geodan, Amsterdam, the Netherlands.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_ALGORITHMS_CLEAR_HPP
#define BOOST_GEOMETRY_ALGORITHMS_CLEAR_HPP


#include <concepts>
#include <type_traits>

#include <boost/geometry/algorithms/not_implemented.hpp>
#include <boost/geometry/core/exterior_ring.hpp>
#include <boost/geometry/core/interior_rings.hpp>
#include <boost/geometry/core/mutable_range.hpp>
#include <boost/geometry/core/tag_cast.hpp>
#include <boost/geometry/core/tags.hpp>
#include <boost/geometry/core/visit.hpp>
#include <boost/geometry/geometries/adapted/std_variant.hpp>
#include <boost/geometry/geometries/concepts/check.hpp>


namespace boost { namespace geometry
{

template <concepts::MutableGeometry Geometry>
inline void clear(Geometry& geometry);

#ifndef DOXYGEN_NO_DETAIL
namespace detail
{

template <concepts::MutableGeometry Polygon>
    requires std::same_as<tag_t<Polygon>, polygon_tag>
inline void clear_impl(Polygon& polygon)
{
    traits::clear
        <
            std::remove_reference_t
                <
                    typename traits::interior_mutable_type<Polygon>::type
                >
        >::apply(interior_rings(polygon));
    traits::clear
        <
            std::remove_reference_t
                <
                    typename traits::ring_mutable_type<Polygon>::type
                >
        >::apply(exterior_ring(polygon));
}

template <concepts::MutableGeometry Geometry>
    requires std::same_as<tag_t<Geometry>, linestring_tag>
          || std::same_as<tag_t<Geometry>, ring_tag>
          || std::same_as<tag_t<Geometry>, multi_point_tag>
          || std::same_as<tag_t<Geometry>, multi_linestring_tag>
          || std::same_as<tag_t<Geometry>, multi_polygon_tag>
          || std::same_as<tag_t<Geometry>, polyhedral_surface_tag>
          || std::same_as<tag_t<Geometry>, geometry_collection_tag>
inline void clear_impl(Geometry& geometry)
{
    traits::clear<Geometry>::apply(geometry);
}

template <concepts::MutableGeometry Geometry>
    requires std::same_as<tag_t<Geometry>, dynamic_geometry_tag>
inline void clear_impl(Geometry& geometry)
{
    traits::visit<Geometry>::apply([](auto& g)
    {
        geometry::clear(g);
    }, geometry);
}

template <concepts::MutableGeometry Geometry>
inline void clear_impl(Geometry&)
{}

} // namespace detail
#endif // DOXYGEN_NO_DETAIL

/*!
\brief Clears a linestring, ring or polygon (exterior+interiors) or multi*
\details Generic function to clear a geometry. All points will be removed from the collection or collections
    making up the geometry. In most cases this is equivalent to the .clear() method of a std::vector<...>. In
    the case of a polygon, this clear functionality is automatically called for the exterior ring, and for the
    interior ring collection. In the case of a point, boxes and segments, nothing will happen.
\ingroup clear
\tparam Geometry \tparam_geometry
\param geometry \param_geometry which will be cleared
\note points and boxes cannot be cleared, instead they can be set to zero by "assign_zero"

\qbk{[include reference/algorithms/clear.qbk]}
*/
template <concepts::MutableGeometry Geometry>
inline void clear(Geometry& geometry)
{
    detail::clear_impl(geometry);
}


}} // namespace boost::geometry


#endif // BOOST_GEOMETRY_ALGORITHMS_CLEAR_HPP
