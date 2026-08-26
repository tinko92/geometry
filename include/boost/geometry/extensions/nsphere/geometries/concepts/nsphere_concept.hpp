// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2007-2012 Barend Gehrels, Amsterdam, the Netherlands.
// Copyright (c) 2008-2012 Bruno Lalande, Paris, France.
// Copyright (c) 2009-2012 Mateusz Loskot, London, UK.

// This file was modified by Oracle on 2021.
// Modifications copyright (c) 2021, Oracle and/or its affiliates.
// Contributed and/or modified by Adam Wulkiewicz, on behalf of Oracle

// Parts of Boost.Geometry are redesigned from Geodan's Geographic Library
// (geolib/GGL), copyright (c) 1995-2010 Geodan, Amsterdam, the Netherlands.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_EXTENSIONS_NSPHERE_GEOMETRIES_CONCEPTS_NSPHERE_CONCEPT_HPP
#define BOOST_GEOMETRY_EXTENSIONS_NSPHERE_GEOMETRIES_CONCEPTS_NSPHERE_CONCEPT_HPP

#include <concepts>
#include <utility>

#include <boost/geometry/core/coordinate_dimension.hpp>
#include <boost/geometry/core/access.hpp>
#include <boost/geometry/core/point_type.hpp>
#include <boost/geometry/extensions/nsphere/core/radius.hpp>

namespace boost { namespace geometry { namespace concepts {

/*!
    \brief Checks Nsphere concept (const version)
    \ingroup concepts
    \details The ConstNsphere concept check the same as the Nsphere concept,
    but does not check write access.
*/
template <typename Geometry, std::size_t Dimension>
concept ConstNsphereCoordinate = requires(Geometry const& geometry)
{
    { geometry::get<Dimension>(geometry) }
        -> std::convertible_to<coordinate_type_t<Geometry>>;
};

template <typename Geometry, std::size_t... Dimensions>
constexpr bool const_nsphere_coordinates(std::index_sequence<Dimensions...>)
{
    return (ConstNsphereCoordinate<Geometry, Dimensions> && ...);
}

template <typename Geometry>
concept ConstNsphere =
    std::same_as<tag_t<Geometry>, nsphere_tag>
    && concepts::ConstPoint<point_type_t<Geometry>>
    && const_nsphere_coordinates<Geometry>(
        std::make_index_sequence<dimension<Geometry>::value>{})
    && requires(Geometry const& geometry)
    {
        { geometry::get_radius<0>(geometry) }
            -> std::convertible_to<radius_type_t<Geometry>>;
    };


/*!
    \brief Checks nsphere concept
    \ingroup concepts
*/
template <typename Geometry, std::size_t Dimension>
concept MutableNsphereCoordinate = requires(
    Geometry& geometry, coordinate_type_t<Geometry> value)
{
    geometry::set<Dimension>(geometry, value);
};

template <typename Geometry, std::size_t... Dimensions>
constexpr bool mutable_nsphere_coordinates(std::index_sequence<Dimensions...>)
{
    return (MutableNsphereCoordinate<Geometry, Dimensions> && ...);
}

template <typename Geometry>
concept Nsphere =
    ConstNsphere<Geometry>
    && concepts::Point<point_type_t<Geometry>>
    && mutable_nsphere_coordinates<Geometry>(
        std::make_index_sequence<dimension<Geometry>::value>{})
    && requires(Geometry& geometry, radius_type_t<Geometry> radius)
    {
        geometry::set_radius<0>(geometry, radius);
    };


template <typename Geometry>
struct concept_type<Geometry, nsphere_tag>
    : std::bool_constant<Nsphere<Geometry>>
{};

template <typename Geometry>
struct concept_type<Geometry const, nsphere_tag>
    : std::bool_constant<ConstNsphere<Geometry>>
{};


}}} // namespace boost::geometry::concepts

#endif // BOOST_GEOMETRY_EXTENSIONS_NSPHERE_GEOMETRIES_CONCEPTS_NSPHERE_CONCEPT_HPP
