// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2007-2012 Barend Gehrels, Amsterdam, the Netherlands.
// Copyright (c) 2008-2012 Bruno Lalande, Paris, France.
// Copyright (c) 2009-2012 Mateusz Loskot, London, UK.

// This file was modified by Oracle on 2020-2021.
// Modifications copyright (c) 2020-2021 Oracle and/or its affiliates.
// Contributed and/or modified by Adam Wulkiewicz, on behalf of Oracle

// Parts of Boost.Geometry are redesigned from Geodan's Geographic Library
// (geolib/GGL), copyright (c) 1995-2010 Geodan, Amsterdam, the Netherlands.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)


#ifndef BOOST_GEOMETRY_GEOMETRIES_CONCEPTS_MULTI_POINT_CONCEPT_HPP
#define BOOST_GEOMETRY_GEOMETRIES_CONCEPTS_MULTI_POINT_CONCEPT_HPP


#include <ranges>

#include <boost/geometry/core/mutable_range.hpp>

#include <boost/geometry/geometries/concepts/concept_type.hpp>
#include <boost/geometry/geometries/concepts/detail/mutable_range.hpp>
#include <boost/geometry/geometries/concepts/point_concept.hpp>


namespace boost { namespace geometry { namespace concepts
{

template <typename Geometry>
concept ConstMultiPoint =
    std::same_as<tag_t<geometry_type_t<Geometry>>, multi_point_tag>
    && detail::ConstRandomAccessRange<geometry_type_t<Geometry>>
    && ConstPoint<std::ranges::range_value_t<geometry_type_t<Geometry>>>;


template <typename Geometry>
concept MultiPoint =
    ! std::is_const_v<std::remove_reference_t<Geometry>>
    && ConstMultiPoint<Geometry>
    && Point<std::ranges::range_value_t<geometry_type_t<Geometry>>>
    && detail::MutableRange
        <geometry_type_t<Geometry>,
         std::ranges::range_value_t<geometry_type_t<Geometry>>>;


template <typename Geometry>
struct concept_type<Geometry, multi_point_tag>
    : std::bool_constant<MultiPoint<Geometry>>
{};

template <typename Geometry>
struct concept_type<Geometry const, multi_point_tag>
    : std::bool_constant<ConstMultiPoint<Geometry>>
{};


}}} // namespace boost::geometry::concepts


#endif // BOOST_GEOMETRY_GEOMETRIES_CONCEPTS_MULTI_POINT_CONCEPT_HPP
