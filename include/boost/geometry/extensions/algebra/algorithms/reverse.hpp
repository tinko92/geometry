// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2007-2012 Barend Gehrels, Amsterdam, the Netherlands.
// Copyright (c) 2008-2012 Bruno Lalande, Paris, France.
// Copyright (c) 2009-2012 Mateusz Loskot, London, UK.
// Copyright (c) 2013 Adam Wulkiewicz, Lodz, Poland.

// Parts of Boost.Geometry are redesigned from Geodan's Geographic Library
// (geolib/GGL), copyright (c) 1995-2010 Geodan, Amsterdam, the Netherlands.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_EXTENSIONS_ALGEBRA_ALGORITHMS_REVERSE_HPP
#define BOOST_GEOMETRY_EXTENSIONS_ALGEBRA_ALGORITHMS_REVERSE_HPP

#include <boost/geometry/algorithms/reverse.hpp>

#include <boost/geometry/extensions/algebra/algorithms/detail.hpp>

namespace boost { namespace geometry
{

// This is experimental implementation of reverse() which negates vectors
// and inverses rotations. It doesn't work for them as for Geometries.

template <typename Vector>
    requires concepts::MutableGeometry<Vector> && concepts::Vector<Vector>
inline void reverse(Vector& vector)
{
    detail::algebra::neg<0, dimension<Vector>::value>(vector);
}

template <typename RotationQuaternion>
    requires concepts::MutableGeometry<RotationQuaternion>
          && concepts::RotationQuaternion<RotationQuaternion>
inline void reverse(RotationQuaternion& rotation)
{
    detail::algebra::neg<1, 4>(rotation);
}

template <typename RotationMatrix>
    requires concepts::MutableGeometry<RotationMatrix>
          && concepts::RotationMatrix<RotationMatrix>
inline void reverse(RotationMatrix& rotation)
{
    detail::algebra::matrix_transpose
        <RotationMatrix, 0, 0, dimension<RotationMatrix>::value>::apply(rotation);
}

}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_EXTENSIONS_ALGEBRA_ALGORITHMS_CLEAR_HPP
