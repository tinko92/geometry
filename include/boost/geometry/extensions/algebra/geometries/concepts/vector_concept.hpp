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

#ifndef BOOST_GEOMETRY_EXTENSIONS_ALGEBRA_GEOMETRIES_CONCEPTS_VECTOR_CONCEPT_HPP
#define BOOST_GEOMETRY_EXTENSIONS_ALGEBRA_GEOMETRIES_CONCEPTS_VECTOR_CONCEPT_HPP

#include <concepts>
#include <type_traits>
#include <utility>

#include <boost/geometry/core/cs.hpp>
#include <boost/geometry/extensions/algebra/core/access.hpp>
#include <boost/geometry/extensions/algebra/core/coordinate_dimension.hpp>
#include <boost/geometry/extensions/algebra/core/coordinate_type.hpp>
#include <boost/geometry/extensions/algebra/core/coordinate_system.hpp>
#include <boost/geometry/extensions/algebra/geometries/concepts/detail/coordinate_concepts.hpp>
#include <boost/geometry/geometries/concepts/concept_type.hpp>

namespace boost { namespace geometry { namespace concepts {

template <typename Geometry>
concept ConstVector =
    std::same_as<tag_t<geometry_type_t<Geometry>>, vector_tag>
    && std::same_as<coordinate_system_t<geometry_type_t<Geometry>>, cs::cartesian>
    && detail::const_algebra_coordinates<geometry_type_t<Geometry>>(
        std::make_index_sequence<dimension<geometry_type_t<Geometry>>::value>{});

template <typename Geometry>
concept Vector =
    ! std::is_const_v<std::remove_reference_t<Geometry>>
    && ConstVector<Geometry>
    && detail::mutable_algebra_coordinates<geometry_type_t<Geometry>>(
        std::make_index_sequence<dimension<geometry_type_t<Geometry>>::value>{});

template <typename Geometry>
struct concept_type<Geometry, vector_tag>
    : std::bool_constant<Vector<Geometry>>
{};

template <typename Geometry>
struct concept_type<Geometry const, vector_tag>
    : std::bool_constant<ConstVector<Geometry>>
{};

}}} // namespace boost::geometry::concepts

#endif // BOOST_GEOMETRY_EXTENSIONS_ALGEBRA_GEOMETRIES_CONCEPTS_VECTOR_CONCEPT_HPP
