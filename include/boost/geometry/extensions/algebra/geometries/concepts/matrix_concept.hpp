// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2007-2012 Barend Gehrels, Amsterdam, the Netherlands.
// Copyright (c) 2008-2012 Bruno Lalande, Paris, France.
// Copyright (c) 2009-2012 Mateusz Loskot, London, UK.
// Copyright (c) 2013 Adam Wulkiewicz, Lodz, Poland.

// This file was modified by Oracle on 2018-2020.
// Modifications copyright (c) 2018-2020 Oracle and/or its affiliates.
// Contributed and/or modified by Adam Wulkiewicz, on behalf of Oracle

// Parts of Boost.Geometry are redesigned from Geodan's Geographic Library
// (geolib/GGL), copyright (c) 1995-2010 Geodan, Amsterdam, the Netherlands.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_EXTENSIONS_ALGEBRA_GEOMETRIES_CONCEPTS_MATRIX_CONCEPT_HPP
#define BOOST_GEOMETRY_EXTENSIONS_ALGEBRA_GEOMETRIES_CONCEPTS_MATRIX_CONCEPT_HPP

#include <concepts>
#include <type_traits>
#include <utility>

#include <boost/geometry/extensions/algebra/core/access.hpp>
#include <boost/geometry/extensions/algebra/core/coordinate_dimension.hpp>
#include <boost/geometry/extensions/algebra/core/coordinate_type.hpp>
#include <boost/geometry/extensions/algebra/geometries/concepts/detail/coordinate_concepts.hpp>
#include <boost/geometry/geometries/concepts/concept_type.hpp>

namespace boost { namespace geometry { namespace concepts {

namespace detail
{

template <typename Geometry>
consteval bool const_matrix_coordinates()
{
    constexpr std::size_t rows = traits::indexed_dimension<Geometry, 0>::value;
    constexpr std::size_t columns = traits::indexed_dimension<Geometry, 1>::value;
    return const_algebra_indexed_coordinates<Geometry, columns>(
        std::make_index_sequence<rows * columns>{});
}

template <typename Geometry>
consteval bool mutable_matrix_coordinates()
{
    constexpr std::size_t rows = traits::indexed_dimension<Geometry, 0>::value;
    constexpr std::size_t columns = traits::indexed_dimension<Geometry, 1>::value;
    return mutable_algebra_indexed_coordinates<Geometry, columns>(
        std::make_index_sequence<rows * columns>{});
}

} // namespace detail

template <typename Geometry>
concept ConstMatrix =
    std::same_as<tag_t<geometry_type_t<Geometry>>, matrix_tag>
    && requires
    {
        typename coordinate_type_t<geometry_type_t<Geometry>>;
        traits::indexed_dimension<geometry_type_t<Geometry>, 0>::value;
        traits::indexed_dimension<geometry_type_t<Geometry>, 1>::value;
    }
    && detail::const_matrix_coordinates<geometry_type_t<Geometry>>();

template <typename Geometry>
concept Matrix =
    ! std::is_const_v<std::remove_reference_t<Geometry>>
    && ConstMatrix<Geometry>
    && detail::mutable_matrix_coordinates<geometry_type_t<Geometry>>();

template <typename Geometry>
struct concept_type<Geometry, matrix_tag>
    : std::bool_constant<Matrix<Geometry>>
{};

template <typename Geometry>
struct concept_type<Geometry const, matrix_tag>
    : std::bool_constant<ConstMatrix<Geometry>>
{};

}}} // namespace boost::geometry::concepts

#endif // BOOST_GEOMETRY_EXTENSIONS_ALGEBRA_GEOMETRIES_CONCEPTS_MATRIX_CONCEPT_HPP
