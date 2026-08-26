// Boost.Geometry

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_EXTENSIONS_ALGEBRA_GEOMETRIES_CONCEPTS_DETAIL_COORDINATE_CONCEPTS_HPP
#define BOOST_GEOMETRY_EXTENSIONS_ALGEBRA_GEOMETRIES_CONCEPTS_DETAIL_COORDINATE_CONCEPTS_HPP

#include <concepts>
#include <cstddef>
#include <utility>

#include <boost/geometry/core/access.hpp>
#include <boost/geometry/core/coordinate_type.hpp>
#include <boost/geometry/extensions/algebra/core/coordinate_dimension.hpp>

namespace boost { namespace geometry { namespace concepts { namespace detail
{

template <typename Geometry, std::size_t Dimension>
concept ConstAlgebraCoordinate = requires(Geometry const& geometry)
{
    { geometry::get<Dimension>(geometry) }
        -> std::convertible_to<coordinate_type_t<Geometry>>;
};

template <typename Geometry, std::size_t Dimension>
concept MutableAlgebraCoordinate =
    ConstAlgebraCoordinate<Geometry, Dimension>
    && requires(Geometry& geometry, coordinate_type_t<Geometry> value)
    {
        geometry::set<Dimension>(geometry, value);
    };

template <typename Geometry, std::size_t... Dimensions>
consteval bool const_algebra_coordinates(std::index_sequence<Dimensions...>)
{
    return (ConstAlgebraCoordinate<Geometry, Dimensions> && ...);
}

template <typename Geometry, std::size_t... Dimensions>
consteval bool mutable_algebra_coordinates(std::index_sequence<Dimensions...>)
{
    return (MutableAlgebraCoordinate<Geometry, Dimensions> && ...);
}

template <typename Geometry, std::size_t Row, std::size_t Column>
concept ConstAlgebraIndexedCoordinate = requires(Geometry const& geometry)
{
    { geometry::get<Row, Column>(geometry) }
        -> std::convertible_to<coordinate_type_t<Geometry>>;
};

template <typename Geometry, std::size_t Row, std::size_t Column>
concept MutableAlgebraIndexedCoordinate =
    ConstAlgebraIndexedCoordinate<Geometry, Row, Column>
    && requires(Geometry& geometry, coordinate_type_t<Geometry> value)
    {
        geometry::set<Row, Column>(geometry, value);
    };

template <typename Geometry, std::size_t Columns, std::size_t... Indices>
consteval bool const_algebra_indexed_coordinates(
    std::index_sequence<Indices...>)
{
    return (ConstAlgebraIndexedCoordinate
        <Geometry, Indices / Columns, Indices % Columns> && ...);
}

template <typename Geometry, std::size_t Columns, std::size_t... Indices>
consteval bool mutable_algebra_indexed_coordinates(
    std::index_sequence<Indices...>)
{
    return (MutableAlgebraIndexedCoordinate
        <Geometry, Indices / Columns, Indices % Columns> && ...);
}

}}}} // namespace boost::geometry::concepts::detail

#endif
