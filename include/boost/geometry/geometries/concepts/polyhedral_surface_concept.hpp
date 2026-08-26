// Boost.Geometry

// Copyright (c) 2025 Siddharth Kumar, Roorkee, India.
// Copyright (c) 2025 Oracle and/or its affiliates.
// Contributed and/or modified by Vissarion Fysikopoulos, on behalf of Oracle

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_GEOMETRIES_CONCEPTS_POLYHEDRAL_SURFACE_CONCEPT_HPP
#define BOOST_GEOMETRY_GEOMETRIES_CONCEPTS_POLYHEDRAL_SURFACE_CONCEPT_HPP

#include <concepts>
#include <type_traits>

#include <boost/range/value_type.hpp>

#include <boost/geometry/core/access.hpp>
#include <boost/geometry/core/cs.hpp>
#include <boost/geometry/core/ring_type.hpp>
#include <boost/geometry/core/tags.hpp>
#include <boost/geometry/geometries/concepts/detail/mutable_range.hpp>
#include <boost/geometry/geometries/concepts/polygon_concept.hpp>

namespace boost { namespace geometry { namespace concepts
{

template <typename Geometry>
concept ConstPolyhedralSurface =
    std::same_as<tag_t<geometry_type_t<Geometry>>, polyhedral_surface_tag>
    && detail::ConstRandomAccessRange<geometry_type_t<Geometry>>
    && ConstPolygon<typename boost::range_value<geometry_type_t<Geometry>>::type>
    && dimension<geometry_type_t<Geometry>>::value == 3
    && std::same_as<cs_tag_t<geometry_type_t<Geometry>>, cartesian_tag>;


template <typename Geometry>
concept PolyhedralSurface =
    ! std::is_const_v<std::remove_reference_t<Geometry>>
    && ConstPolyhedralSurface<Geometry>
    && Polygon<typename boost::range_value<geometry_type_t<Geometry>>::type>
    && detail::MutableRange
        <geometry_type_t<Geometry>,
         typename boost::range_value<geometry_type_t<Geometry>>::type>;

template <typename Geometry>
struct concept_type<Geometry, polyhedral_surface_tag>
    : std::bool_constant<PolyhedralSurface<Geometry>>
{};

template <typename Geometry>
struct concept_type<Geometry const, polyhedral_surface_tag>
    : std::bool_constant<ConstPolyhedralSurface<Geometry>>
{};

}}} // namespace boost::geometry::concepts
#endif // BOOST_GEOMETRY_GEOMETRIES_CONCEPTS_POLYHEDRAL_SURFACE_CONCEPT_HPP
