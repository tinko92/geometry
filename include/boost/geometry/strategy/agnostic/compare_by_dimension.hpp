// Boost.Geometry (aka GGL, Generic Geometry Library)

// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_STRATEGY_AGNOSTIC_COMPARE_BY_DIMENSION_HPP
#define BOOST_GEOMETRY_STRATEGY_AGNOSTIC_COMPARE_BY_DIMENSION_HPP

#include <boost/geometry/core/access.hpp>

#include <boost/geometry/util/condition.hpp>

namespace boost { namespace geometry
{

namespace strategy { namespace compare_by_dimension
{

template <bool Less>
struct agnostic
{
    template <std::size_t Dimension, typename Point>
    static bool apply(Point const& p1, Point const& p2)
    {
        if ( BOOST_GEOMETRY_CONDITION(Less) )
        {
            return get<Dimension>(p1) < get<Dimension>(p2);
        }
        else
        {
            return get<Dimension>(p1) > get<Dimension>(p2);
        }
    }
};

}} // namespace strategy::less_by_dimension

}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_STRATEGY_AGNOSTIC_COMPARE_BY_DIMENSION_HPP
