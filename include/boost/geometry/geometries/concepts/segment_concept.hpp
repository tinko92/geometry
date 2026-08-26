// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2008-2012 Bruno Lalande, Paris, France.
// Copyright (c) 2008-2012 Barend Gehrels, Amsterdam, the Netherlands.
// Copyright (c) 2009-2012 Mateusz Loskot, London, UK.

// This file was modified by Oracle on 2021.
// Modifications copyright (c) 2021 Oracle and/or its affiliates.
// Contributed and/or modified by Adam Wulkiewicz, on behalf of Oracle

// Parts of Boost.Geometry are redesigned from Geodan's Geographic Library
// (geolib/GGL), copyright (c) 1995-2010 Geodan, Amsterdam, the Netherlands.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_GEOMETRIES_CONCEPTS_SEGMENT_CONCEPT_HPP
#define BOOST_GEOMETRY_GEOMETRIES_CONCEPTS_SEGMENT_CONCEPT_HPP


#include <concepts>
#include <type_traits>
#include <utility>

#include <boost/geometry/core/access.hpp>
#include <boost/geometry/core/point_type.hpp>

#include <boost/geometry/geometries/concepts/concept_type.hpp>
#include <boost/geometry/geometries/concepts/point_concept.hpp>


namespace boost { namespace geometry { namespace concepts
{

template <typename Geometry>
concept ConstSegment =
    std::same_as<tag_t<geometry_type_t<Geometry>>, segment_tag>
    && ConstPoint<point_type_t<geometry_type_t<Geometry>>>
    && detail::const_indexed_coordinates<geometry_type_t<Geometry>, 0>(
        std::make_index_sequence<dimension<point_type_t<geometry_type_t<Geometry>>>::value>{})
    && detail::const_indexed_coordinates<geometry_type_t<Geometry>, 1>(
        std::make_index_sequence<dimension<point_type_t<geometry_type_t<Geometry>>>::value>{});


template <typename Geometry>
concept Segment =
    ! std::is_const_v<std::remove_reference_t<Geometry>>
    && ConstSegment<Geometry>
    && Point<point_type_t<geometry_type_t<Geometry>>>
    && detail::mutable_indexed_coordinates<geometry_type_t<Geometry>, 0>(
        std::make_index_sequence<dimension<point_type_t<geometry_type_t<Geometry>>>::value>{})
    && detail::mutable_indexed_coordinates<geometry_type_t<Geometry>, 1>(
        std::make_index_sequence<dimension<point_type_t<geometry_type_t<Geometry>>>::value>{});


template <typename Geometry>
struct concept_type<Geometry, segment_tag>
    : std::bool_constant<Segment<Geometry>>
{};

template <typename Geometry>
struct concept_type<Geometry const, segment_tag>
    : std::bool_constant<ConstSegment<Geometry>>
{};


}}} // namespace boost::geometry::concepts


#endif // BOOST_GEOMETRY_GEOMETRIES_CONCEPTS_SEGMENT_CONCEPT_HPP
