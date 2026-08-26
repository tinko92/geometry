// Boost.Geometry

// Copyright (c) 2021, Oracle and/or its affiliates.

// Contributed and/or modified by Adam Wulkiewicz, on behalf of Oracle

// Licensed under the Boost Software License version 1.0.
// http://www.boost.org/users/license.html

#ifndef BOOST_GEOMETRY_GEOMETRIES_CONCEPTS_CONCEPT_TYPE_HPP
#define BOOST_GEOMETRY_GEOMETRIES_CONCEPTS_CONCEPT_TYPE_HPP


#include <concepts>
#include <type_traits>

#include <boost/geometry/core/tag.hpp>
#include <boost/geometry/core/tags.hpp>


namespace boost { namespace geometry { namespace concepts
{

template <typename Geometry, typename Tag = tag_t<Geometry>>
struct concept_type : std::false_type
{};

template <typename Geometry>
using geometry_type_t = std::remove_cvref_t<Geometry>;


template <typename Geometry>
concept GeometryType = concept_type<std::remove_reference_t<Geometry>>::value;

template <typename Geometry>
concept ConstGeometry = concept_type
    <geometry_type_t<Geometry> const>::value;

template <typename Geometry>
concept MutableGeometry =
    ! std::is_const_v<std::remove_reference_t<Geometry>>
    && concept_type<geometry_type_t<Geometry>>::value;

template <typename Geometry, typename Category>
concept GeometryCategory =
    ConstGeometry<Geometry>
    && std::derived_from<tag_t<geometry_type_t<Geometry>>, Category>;

template <typename Geometry>
concept ArealGeometry = GeometryCategory<Geometry, areal_tag>;


}}} // namespace boost::geometry::concepts


#endif // BOOST_GEOMETRY_GEOMETRIES_CONCEPTS_CONCEPT_TYPE_HPP
