// Boost.Geometry

// Copyright (c) 2021, Oracle and/or its affiliates.

// Contributed and/or modified by Adam Wulkiewicz, on behalf of Oracle

// Licensed under the Boost Software License version 1.0.
// http://www.boost.org/users/license.html

#ifndef BOOST_GEOMETRY_GEOMETRIES_CONCEPTS_DYNAMIC_GEOMETRY_CONCEPT_HPP
#define BOOST_GEOMETRY_GEOMETRIES_CONCEPTS_DYNAMIC_GEOMETRY_CONCEPT_HPP


#include <concepts>
#include <type_traits>
#include <utility>

#include <boost/geometry/core/geometry_types.hpp>
#include <boost/geometry/core/tags.hpp>
#include <boost/geometry/core/visit.hpp>

#include <boost/geometry/geometries/concepts/box_concept.hpp>
#include <boost/geometry/geometries/concepts/concept_type.hpp>
#include <boost/geometry/geometries/concepts/geometry_collection_concept.hpp>
#include <boost/geometry/geometries/concepts/linestring_concept.hpp>
#include <boost/geometry/geometries/concepts/multi_point_concept.hpp>
#include <boost/geometry/geometries/concepts/multi_linestring_concept.hpp>
#include <boost/geometry/geometries/concepts/multi_polygon_concept.hpp>
#include <boost/geometry/geometries/concepts/point_concept.hpp>
#include <boost/geometry/geometries/concepts/polygon_concept.hpp>
#include <boost/geometry/geometries/concepts/ring_concept.hpp>
#include <boost/geometry/geometries/concepts/segment_concept.hpp>


namespace boost { namespace geometry { namespace concepts
{

namespace detail
{

template <typename DynamicGeometry, typename SubGeometry>
concept MutableDynamicAlternative =
    is_recursive_geometry_v<SubGeometry>
    || (concepts::GeometryType<SubGeometry>
        && requires(DynamicGeometry& dynamic, SubGeometry&& geometry)
        {
            dynamic = std::move(geometry);
        });

template <typename DynamicGeometry, typename Sequence>
struct mutable_dynamic_alternatives : std::false_type
{};

template <typename DynamicGeometry, typename... SubGeometries>
struct mutable_dynamic_alternatives
    <DynamicGeometry, util::type_sequence<SubGeometries...>>
    : std::bool_constant
        <(MutableDynamicAlternative<DynamicGeometry, SubGeometries> && ...)>
{};

} // namespace detail


template <typename Geometry>
concept ConstDynamicGeometry =
    std::same_as<tag_t<geometry_type_t<Geometry>>, dynamic_geometry_tag>
    && requires
    {
        typename traits::geometry_types<geometry_type_t<Geometry>>::type;
        requires detail::const_collection_alternatives
            <typename traits::geometry_types<geometry_type_t<Geometry>>::type>::value;
    }
    && requires(geometry_type_t<Geometry> const& dynamic)
    {
        traits::visit<geometry_type_t<Geometry>>::apply([](auto&&) {}, dynamic);
    };


template <typename Geometry>
concept DynamicGeometry =
    ! std::is_const_v<std::remove_reference_t<Geometry>>
    && ConstDynamicGeometry<Geometry>
    && requires
    {
        requires detail::mutable_dynamic_alternatives
            <geometry_type_t<Geometry>,
             typename traits::geometry_types<geometry_type_t<Geometry>>::type>::value;
    }
    && requires(geometry_type_t<Geometry>& dynamic)
    {
        traits::visit<geometry_type_t<Geometry>>::apply([](auto&&) {}, dynamic);
    };


template <typename Geometry>
struct concept_type<Geometry, dynamic_geometry_tag>
    : std::bool_constant<DynamicGeometry<Geometry>>
{};

template <typename Geometry>
struct concept_type<Geometry const, dynamic_geometry_tag>
    : std::bool_constant<ConstDynamicGeometry<Geometry>>
{};


}}} // namespace boost::geometry::concepts


#endif // BOOST_GEOMETRY_GEOMETRIES_CONCEPTS_DYNAMIC_GEOMETRY_CONCEPT_HPP
