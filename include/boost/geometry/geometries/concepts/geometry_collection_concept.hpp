// Boost.Geometry

// Copyright (c) 2021, Oracle and/or its affiliates.

// Contributed and/or modified by Adam Wulkiewicz, on behalf of Oracle

// Licensed under the Boost Software License version 1.0.
// http://www.boost.org/users/license.html

#ifndef BOOST_GEOMETRY_GEOMETRIES_CONCEPTS_GEOMETRY_COLLECTION_CONCEPT_HPP
#define BOOST_GEOMETRY_GEOMETRIES_CONCEPTS_GEOMETRY_COLLECTION_CONCEPT_HPP


#include <concepts>
#include <ranges>
#include <type_traits>
#include <utility>

#include <boost/geometry/core/geometry_types.hpp>
#include <boost/geometry/core/mutable_range.hpp>
#include <boost/geometry/core/tag.hpp>
#include <boost/geometry/core/tags.hpp>
#include <boost/geometry/core/visit.hpp>

#include <boost/geometry/geometries/concepts/box_concept.hpp>
#include <boost/geometry/geometries/concepts/concept_type.hpp>
#include <boost/geometry/geometries/concepts/linestring_concept.hpp>
#include <boost/geometry/geometries/concepts/multi_point_concept.hpp>
#include <boost/geometry/geometries/concepts/multi_linestring_concept.hpp>
#include <boost/geometry/geometries/concepts/multi_polygon_concept.hpp>
#include <boost/geometry/geometries/concepts/point_concept.hpp>
#include <boost/geometry/geometries/concepts/polygon_concept.hpp>
#include <boost/geometry/geometries/concepts/ring_concept.hpp>
#include <boost/geometry/geometries/concepts/segment_concept.hpp>

#include <boost/geometry/util/sequence.hpp>
#include <boost/geometry/util/type_traits.hpp>


namespace boost { namespace geometry { namespace concepts
{

namespace detail
{

template <typename Geometry>
inline constexpr bool is_recursive_geometry_v =
    util::is_dynamic_geometry<Geometry>::value
    || util::is_geometry_collection<Geometry>::value;

template <typename Collection, typename SubGeometry>
concept MutableCollectionAlternative =
    is_recursive_geometry_v<SubGeometry>
    || (concepts::GeometryType<SubGeometry>
        && requires(Collection& collection, SubGeometry&& geometry)
        {
            traits::emplace_back<Collection>::apply(
                collection, std::move(geometry));
        });

template <typename SubGeometry>
concept ConstCollectionAlternative =
    is_recursive_geometry_v<SubGeometry>
    || concepts::GeometryType<SubGeometry const>;

template <typename Collection, typename Sequence>
struct mutable_collection_alternatives : std::false_type
{};

template <typename Collection, typename... SubGeometries>
struct mutable_collection_alternatives
    <Collection, util::type_sequence<SubGeometries...>>
    : std::bool_constant
        <(MutableCollectionAlternative<Collection, SubGeometries> && ...)>
{};

template <typename Sequence>
struct const_collection_alternatives : std::false_type
{};

template <typename... SubGeometries>
struct const_collection_alternatives<util::type_sequence<SubGeometries...>>
    : std::bool_constant<(ConstCollectionAlternative<SubGeometries> && ...)>
{};


} // namespace detail


template <typename Geometry>
concept ConstGeometryCollection =
    std::same_as<tag_t<geometry_type_t<Geometry>>, geometry_collection_tag>
    && detail::ConstForwardRange<geometry_type_t<Geometry>>
    && requires
    {
        typename traits::geometry_types<geometry_type_t<Geometry>>::type;
        requires detail::const_collection_alternatives
            <typename traits::geometry_types<geometry_type_t<Geometry>>::type>::value;
    }
    && requires(geometry_type_t<Geometry> const& collection)
    {
        traits::iter_visit<geometry_type_t<Geometry>>::apply(
            [](auto&&) {}, std::ranges::begin(collection));
    };


template <typename Geometry>
concept GeometryCollection =
    ! std::is_const_v<std::remove_reference_t<Geometry>>
    && ConstGeometryCollection<Geometry>
    && requires
    {
        requires detail::mutable_collection_alternatives
            <geometry_type_t<Geometry>,
             typename traits::geometry_types<geometry_type_t<Geometry>>::type>::value;
    }
    && requires(geometry_type_t<Geometry>& collection)
    {
        traits::clear<geometry_type_t<Geometry>>::apply(collection);
        traits::iter_visit<geometry_type_t<Geometry>>::apply(
            [](auto&&) {}, std::ranges::begin(collection));
    };


template <typename Geometry>
struct concept_type<Geometry, geometry_collection_tag>
    : std::bool_constant<GeometryCollection<Geometry>>
{};

template <typename Geometry>
struct concept_type<Geometry const, geometry_collection_tag>
    : std::bool_constant<ConstGeometryCollection<Geometry>>
{};


}}} // namespace boost::geometry::concepts


#endif // BOOST_GEOMETRY_GEOMETRIES_CONCEPTS_GEOMETRY_COLLECTION_CONCEPT_HPP
