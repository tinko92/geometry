// Boost.Geometry (aka GGL, Generic Geometry Library)
//
// Copyright (c) 2008-2014 Bruno Lalande, Paris, France.
// Copyright (c) 2008-2014 Barend Gehrels, Amsterdam, the Netherlands.
// Copyright (c) 2009-2014 Mateusz Loskot, London, UK.

// This file was modified by Oracle on 2014-2021.
// Modifications copyright (c) 2014-2021, Oracle and/or its affiliates.

// Contributed and/or modified by Menelaos Karavelas, on behalf of Oracle
// Contributed and/or modified by Adam Wulkiewicz, on behalf of Oracle

// Parts of Boost.Geometry are redesigned from Geodan's Geographic Library
// (geolib/GGL), copyright (c) 1995-2010 Geodan, Amsterdam, the Netherlands.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_GEOMETRIES_CONCEPTS_POINT_CONCEPT_HPP
#define BOOST_GEOMETRY_GEOMETRIES_CONCEPTS_POINT_CONCEPT_HPP

#include <concepts>
#include <cstddef>
#include <type_traits>
#include <utility>

#include <boost/geometry/core/access.hpp>
#include <boost/geometry/core/coordinate_dimension.hpp>
#include <boost/geometry/core/coordinate_system.hpp>

#include <boost/geometry/geometries/concepts/concept_type.hpp>


namespace boost { namespace geometry { namespace concepts
{

namespace detail
{

template <typename Geometry, std::size_t Dimension>
concept ConstPointCoordinate = requires(Geometry const& geometry)
{
    { geometry::get<Dimension>(geometry) }
        -> std::convertible_to<coordinate_type_t<Geometry>>;
};

template <typename Geometry, std::size_t Dimension>
concept MutablePointCoordinate =
    ConstPointCoordinate<Geometry, Dimension>
    && requires(Geometry& geometry, coordinate_type_t<Geometry> value)
    {
        geometry::set<Dimension>(geometry, value);
    };

template <typename Geometry, std::size_t... Dimensions>
constexpr bool const_point_coordinates(std::index_sequence<Dimensions...>)
{
    return (ConstPointCoordinate<Geometry, Dimensions> && ...);
}

template <typename Geometry, std::size_t... Dimensions>
constexpr bool mutable_point_coordinates(std::index_sequence<Dimensions...>)
{
    return (MutablePointCoordinate<Geometry, Dimensions> && ...);
}

template <typename Geometry, std::size_t Index, std::size_t Dimension>
concept ConstIndexedCoordinate = requires(Geometry const& geometry)
{
    { geometry::get<Index, Dimension>(geometry) }
        -> std::convertible_to<coordinate_type_t<Geometry>>;
};

template <typename Geometry, std::size_t Index, std::size_t Dimension>
concept MutableIndexedCoordinate =
    ConstIndexedCoordinate<Geometry, Index, Dimension>
    && requires(Geometry& geometry, coordinate_type_t<Geometry> value)
    {
        geometry::set<Index, Dimension>(geometry, value);
    };

template <typename Geometry, std::size_t Index, std::size_t... Dimensions>
constexpr bool const_indexed_coordinates(std::index_sequence<Dimensions...>)
{
    return (ConstIndexedCoordinate<Geometry, Index, Dimensions> && ...);
}

template <typename Geometry, std::size_t Index, std::size_t... Dimensions>
constexpr bool mutable_indexed_coordinates(std::index_sequence<Dimensions...>)
{
    return (MutableIndexedCoordinate<Geometry, Index, Dimensions> && ...);
}

} // namespace detail


/*!
\brief point concept (const version).

\ingroup const_concepts

\details The ConstPoint concept apply the same as the Point concept,
but does not apply write access.

*/
template <typename Geometry>
concept ConstPoint =
    std::same_as<tag_t<geometry_type_t<Geometry>>, point_tag>
    && requires
    {
        typename coordinate_type_t<geometry_type_t<Geometry>>;
        typename coordinate_system_t<geometry_type_t<Geometry>>;
        sizeof(coordinate_system_t<geometry_type_t<Geometry>>);
    }
    && detail::const_point_coordinates<geometry_type_t<Geometry>>(
        std::make_index_sequence<dimension<geometry_type_t<Geometry>>::value>{});


template <typename Geometry>
concept Point =
    ! std::is_const_v<std::remove_reference_t<Geometry>>
    && ConstPoint<Geometry>
    && detail::mutable_point_coordinates<geometry_type_t<Geometry>>(
        std::make_index_sequence<dimension<geometry_type_t<Geometry>>::value>{});


template <typename Geometry>
struct concept_type<Geometry, point_tag>
    : std::bool_constant<Point<Geometry>>
{};

template <typename Geometry>
struct concept_type<Geometry const, point_tag>
    : std::bool_constant<ConstPoint<Geometry>>
{};


}}} // namespace boost::geometry::concepts

#endif // BOOST_GEOMETRY_GEOMETRIES_CONCEPTS_POINT_CONCEPT_HPP
