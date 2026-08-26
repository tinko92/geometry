// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2008-2012 Bruno Lalande, Paris, France.
// Copyright (c) 2008-2012 Barend Gehrels, Amsterdam, the Netherlands.
// Copyright (c) 2009-2012 Mateusz Loskot, London, UK.

// This file was modified by Oracle on 2020-2021.
// Modifications copyright (c) 2020-2021, Oracle and/or its affiliates.
// Contributed and/or modified by Adam Wulkiewicz, on behalf of Oracle

// Parts of Boost.Geometry are redesigned from Geodan's Geographic Library
// (geolib/GGL), copyright (c) 1995-2010 Geodan, Amsterdam, the Netherlands.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_GEOMETRIES_CONCEPTS_POLYGON_CONCEPT_HPP
#define BOOST_GEOMETRY_GEOMETRIES_CONCEPTS_POLYGON_CONCEPT_HPP

#include <concepts>
#include <type_traits>

#include <boost/geometry/core/access.hpp>
#include <boost/geometry/core/exterior_ring.hpp>
#include <boost/geometry/core/interior_rings.hpp>
#include <boost/geometry/core/point_type.hpp>
#include <boost/geometry/core/ring_type.hpp>

#include <boost/geometry/geometries/concepts/concept_type.hpp>
#include <boost/geometry/geometries/concepts/point_concept.hpp>
#include <boost/geometry/geometries/concepts/ring_concept.hpp>


namespace boost { namespace geometry { namespace concepts
{

/*!
\brief Checks polygon concept
\ingroup concepts
*/
template <typename Geometry>
concept ConstPolygon =
    std::same_as<tag_t<geometry_type_t<Geometry>>, polygon_tag>
    && ConstPoint<point_type_t<geometry_type_t<Geometry>>>
    && ConstRing<ring_type_t<geometry_type_t<Geometry>>>
    && requires(geometry_type_t<Geometry> const& polygon)
    {
        { traits::exterior_ring<geometry_type_t<Geometry>>::get(polygon) }
            -> std::convertible_to<typename traits::ring_const_type<geometry_type_t<Geometry>>::type>;
        { traits::interior_rings<geometry_type_t<Geometry>>::get(polygon) }
            -> std::convertible_to<typename traits::interior_const_type<geometry_type_t<Geometry>>::type>;
    };


template <typename Geometry>
concept Polygon =
    ! std::is_const_v<std::remove_reference_t<Geometry>>
    && ConstPolygon<Geometry>
    && Point<point_type_t<geometry_type_t<Geometry>>>
    && Ring<ring_type_t<geometry_type_t<Geometry>>>
    && requires(geometry_type_t<Geometry>& polygon)
    {
        { traits::exterior_ring<geometry_type_t<Geometry>>::get(polygon) }
            -> std::convertible_to<typename traits::ring_mutable_type<geometry_type_t<Geometry>>::type>;
        { traits::interior_rings<geometry_type_t<Geometry>>::get(polygon) }
            -> std::convertible_to<typename traits::interior_mutable_type<geometry_type_t<Geometry>>::type>;
    };


template <typename Geometry>
struct concept_type<Geometry, polygon_tag>
    : std::bool_constant<Polygon<Geometry>>
{};

template <typename Geometry>
struct concept_type<Geometry const, polygon_tag>
    : std::bool_constant<ConstPolygon<Geometry>>
{};


}}} // namespace boost::geometry::concepts

#endif // BOOST_GEOMETRY_GEOMETRIES_CONCEPTS_POLYGON_CONCEPT_HPP
