// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2007-2012 Barend Gehrels, Amsterdam, the Netherlands.
// Copyright (c) 2008-2014 Bruno Lalande, Paris, France.
// Copyright (c) 2009-2014 Mateusz Loskot, London, UK.
// Copyright (c) 2013 Adam Wulkiewicz, Lodz, Poland.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_EXTENSIONS_ALGEBRA_GEOMETRIES_CONCEPTS_ROTATION_MATRIX_CONCEPT_HPP
#define BOOST_GEOMETRY_EXTENSIONS_ALGEBRA_GEOMETRIES_CONCEPTS_ROTATION_MATRIX_CONCEPT_HPP

#include <concepts>
#include <type_traits>
#include <utility>

#include <boost/geometry/core/cs.hpp>
#include <boost/geometry/extensions/algebra/core/access.hpp>
#include <boost/geometry/extensions/algebra/core/coordinate_dimension.hpp>
#include <boost/geometry/extensions/algebra/core/coordinate_system.hpp>
#include <boost/geometry/extensions/algebra/core/coordinate_type.hpp>
#include <boost/geometry/extensions/algebra/geometries/concepts/detail/coordinate_concepts.hpp>
#include <boost/geometry/geometries/concepts/concept_type.hpp>

namespace boost { namespace geometry { namespace concepts
{

namespace detail
{

template <typename Geometry>
consteval bool const_rotation_matrix_coordinates()
{
    constexpr std::size_t size = dimension<Geometry>::value;
    return const_algebra_indexed_coordinates<Geometry, size>(
        std::make_index_sequence<size * size>{});
}

template <typename Geometry>
consteval bool mutable_rotation_matrix_coordinates()
{
    constexpr std::size_t size = dimension<Geometry>::value;
    return mutable_algebra_indexed_coordinates<Geometry, size>(
        std::make_index_sequence<size * size>{});
}

} // namespace detail

template <typename Geometry>
concept ConstRotationMatrix =
    std::same_as<tag_t<geometry_type_t<Geometry>>, rotation_matrix_tag>
    && std::same_as<coordinate_system_t<geometry_type_t<Geometry>>, cs::cartesian>
    && detail::const_rotation_matrix_coordinates<geometry_type_t<Geometry>>();

template <typename Geometry>
concept RotationMatrix =
    ! std::is_const_v<std::remove_reference_t<Geometry>>
    && ConstRotationMatrix<Geometry>
    && detail::mutable_rotation_matrix_coordinates<geometry_type_t<Geometry>>();

template <typename Geometry>
struct concept_type<Geometry, rotation_matrix_tag>
    : std::bool_constant<RotationMatrix<Geometry>>
{};

template <typename Geometry>
struct concept_type<Geometry const, rotation_matrix_tag>
    : std::bool_constant<ConstRotationMatrix<Geometry>>
{};

}}} // namespace boost::geometry::concepts

#endif
