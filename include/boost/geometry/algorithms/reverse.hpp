// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2007-2012 Barend Gehrels, Amsterdam, the Netherlands.
// Copyright (c) 2008-2012 Bruno Lalande, Paris, France.
// Copyright (c) 2009-2012 Mateusz Loskot, London, UK.
// Copyright (c) 2014 Adam Wulkiewicz, Lodz, Poland.

// This file was modified by Oracle on 2020-2023.
// Modifications copyright (c) 2020-2023 Oracle and/or its affiliates.
// Contributed and/or modified by Vissarion Fysikopoulos, on behalf of Oracle
// Contributed and/or modified by Adam Wulkiewicz, on behalf of Oracle

// Parts of Boost.Geometry are redesigned from Geodan's Geographic Library
// (geolib/GGL), copyright (c) 1995-2010 Geodan, Amsterdam, the Netherlands.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_ALGORITHMS_REVERSE_HPP
#define BOOST_GEOMETRY_ALGORITHMS_REVERSE_HPP

#include <algorithm>

#include <boost/range/begin.hpp>
#include <boost/range/end.hpp>

#include <boost/geometry/algorithms/detail/multi_modify.hpp>
#include <boost/geometry/algorithms/detail/visit.hpp>
#include <boost/geometry/core/interior_rings.hpp>
#include <boost/geometry/core/tags.hpp>
#include <boost/geometry/core/visit.hpp>
#include <boost/geometry/geometries/adapted/std_variant.hpp>
#include <boost/geometry/geometries/concepts/check.hpp>


namespace boost { namespace geometry
{

template <concepts::MutableGeometry Geometry>
inline void reverse(Geometry& geometry);

#ifndef DOXYGEN_NO_DETAIL
namespace detail
{

template <concepts::MutableGeometry Geometry>
    requires concepts::Ring<Geometry> || concepts::Linestring<Geometry>
inline void reverse_impl(Geometry& geometry)
{
    std::reverse(boost::begin(geometry), boost::end(geometry));
}

template <concepts::MutableGeometry Polygon>
    requires concepts::Polygon<Polygon>
inline void reverse_impl(Polygon& polygon)
{
    reverse_impl(exterior_ring(polygon));

    auto&& rings = interior_rings(polygon);
    auto const end = boost::end(rings);
    for (auto it = boost::begin(rings); it != end; ++it)
    {
        reverse_impl(*it);
    }
}

template <concepts::MutableGeometry MultiGeometry>
    requires concepts::MultiLinestring<MultiGeometry>
          || concepts::MultiPolygon<MultiGeometry>
inline void reverse_impl(MultiGeometry& geometry)
{
    for (auto it = boost::begin(geometry); it != boost::end(geometry); ++it)
    {
        reverse_impl(*it);
    }
}

template <concepts::MutableGeometry DynamicGeometry>
    requires concepts::DynamicGeometry<DynamicGeometry>
inline void reverse_impl(DynamicGeometry& geometry)
{
    traits::visit<DynamicGeometry>::apply([](auto& g)
    {
        geometry::reverse(g);
    }, geometry);
}

template <concepts::MutableGeometry GeometryCollection>
    requires concepts::GeometryCollection<GeometryCollection>
inline void reverse_impl(GeometryCollection& geometry)
{
    detail::visit_breadth_first([](auto& g)
    {
        geometry::reverse(g);
        return true;
    }, geometry);
}

template <concepts::MutableGeometry Geometry>
inline void reverse_impl(Geometry&)
{
}

} // namespace detail
#endif // DOXYGEN_NO_DETAIL


/*!
\brief Reverses the points within a geometry
\details Generic function to reverse a geometry. It resembles the std::reverse
   functionality, but it takes the geometry type into account. Only for a ring
   or for a linestring it is the same as the std::reverse.
\ingroup reverse
\tparam Geometry \tparam_geometry
\param geometry \param_geometry which will be reversed

\qbk{[include reference/algorithms/reverse.qbk]}
*/
template <concepts::MutableGeometry Geometry>
inline void reverse(Geometry& geometry)
{
    detail::reverse_impl(geometry);
}

}} // namespace boost::geometry


#endif // BOOST_GEOMETRY_ALGORITHMS_REVERSE_HPP
