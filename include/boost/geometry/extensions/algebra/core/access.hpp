// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2007-2012 Barend Gehrels, Amsterdam, the Netherlands.
// Copyright (c) 2008-2012 Bruno Lalande, Paris, France.
// Copyright (c) 2009-2012 Mateusz Loskot, London, UK.
// Copyright (c) 2013 Adam Wulkiewicz, Lodz, Poland.

// This file was modified by Oracle on 2020.
// Modifications copyright (c) 2020, Oracle and/or its affiliates.
// Contributed and/or modified by Adam Wulkiewicz, on behalf of Oracle

// Parts of Boost.Geometry are redesigned from Geodan's Geographic Library
// (geolib/GGL), copyright (c) 1995-2010 Geodan, Amsterdam, the Netherlands.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_EXTENSIONS_ALGEBRA_CORE_ACCESS_HPP
#define BOOST_GEOMETRY_EXTENSIONS_ALGEBRA_CORE_ACCESS_HPP

#include <boost/geometry/core/access.hpp>

#include <boost/geometry/extensions/algebra/core/tags.hpp>

namespace boost { namespace geometry
{

#ifndef DOXYGEN_NO_DISPATCH
namespace core_dispatch
{

template <typename Vector, std::size_t Dimension>
struct access<vector_tag, Vector, Dimension, false>
{
    static inline coordinate_type_t<Vector> get(Vector const& v)
    {
        return traits::access<Vector, Dimension>::get(v);
    }
    static inline void set(Vector& v, coordinate_type_t<Vector> const& value)
    {
        traits::access<Vector, Dimension>::set(v, value);
    }
};

template <typename Vector, std::size_t Dimension>
struct access<vector_tag, Vector, Dimension, true>
{
    static inline coordinate_type_t<Vector> get(Vector const* v)
    {
        return traits::access<std::remove_pointer_t<Vector>, Dimension>::get(*v);
    }
    static inline void set(Vector* v, coordinate_type_t<Vector> const& value)
    {
        traits::access<std::remove_pointer_t<Vector>, Dimension>::set(*v, value);
    }
};

template <typename Q, std::size_t Dimension>
struct access<quaternion_tag, Q, Dimension, false>
{
    static inline coordinate_type_t<Q> get(Q const& v)
    {
        return traits::access<Q, Dimension>::get(v);
    }
    static inline void set(Q& v, coordinate_type_t<Q> const& value)
    {
        traits::access<Q, Dimension>::set(v, value);
    }
};

template <typename Q, std::size_t Dimension>
struct access<quaternion_tag, Q, Dimension, true>
{
    static inline coordinate_type_t<Q> get(Q const* v)
    {
        return traits::access<std::remove_pointer_t<Q>, Dimension>::get(*v);
    }
    static inline void set(Q* v, coordinate_type_t<Q> const& value)
    {
        traits::access<std::remove_pointer_t<Q>, Dimension>::set(*v, value);
    }
};

template <typename M, std::size_t I, std::size_t J>
struct indexed_access<matrix_tag, M, I, J, false>
    : detail::indexed_access_non_pointer<M, I, J>
{};

template <typename M, std::size_t I, std::size_t J>
struct indexed_access<matrix_tag, M, I, J, true>
    : detail::indexed_access_pointer<M, I, J>
{};


template <typename Q, std::size_t Dimension>
struct access<rotation_quaternion_tag, Q, Dimension, false>
{
    static inline coordinate_type_t<Q> get(Q const& v)
    {
        return traits::access<Q, Dimension>::get(v);
    }
    static inline void set(Q& v, coordinate_type_t<Q> const& value)
    {
        traits::access<Q, Dimension>::set(v, value);
    }
};

template <typename Q, std::size_t Dimension>
struct access<rotation_quaternion_tag, Q, Dimension, true>
{
    static inline coordinate_type_t<Q> get(Q const* v)
    {
        return traits::access<std::remove_pointer_t<Q>, Dimension>::get(*v);
    }
    static inline void set(Q* v, coordinate_type_t<Q> const& value)
    {
        traits::access<std::remove_pointer_t<Q>, Dimension>::set(*v, value);
    }
};

template <typename RM, std::size_t I, std::size_t J>
struct indexed_access<rotation_matrix_tag, RM, I, J, false>
    : detail::indexed_access_non_pointer<RM, I, J>
{};

template <typename RM, std::size_t I, std::size_t J>
struct indexed_access<rotation_matrix_tag, RM, I, J, true>
    : detail::indexed_access_pointer<RM, I, J>
{};

} // namespace core_dispatch
#endif // DOXYGEN_NO_DISPATCH


}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_EXTENSIONS_ALGEBRA_CORE_ACCESS_HPP
