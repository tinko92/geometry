// Boost.Geometry (aka GGL, Generic Geometry Library)

// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_STRATEGY_SPHERICAL_SIDE_VALUE_ZERO_HPP
#define BOOST_GEOMETRY_STRATEGY_SPHERICAL_SIDE_VALUE_ZERO_HPP

#include <boost/geometry/algorithms/detail/make/make.hpp>

#include <boost/geometry/arithmetic/infinite_line_functions.hpp>

#include <boost/geometry/util/select_coordinate_type.hpp>

namespace boost { namespace geometry
{

namespace strategy { namespace side_value_line
{

struct zero
{
    template <typename Point1, typename Point2>
    static inline auto apply(Point1 const& p1, Point1 const& p2, Point2 const& p)
    {
        typename select_coordinate_type<Point1, Point2>::type result(0);
        return result;
    }
};

}} // namespace strategy::direction

}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_STRATEGY_SPHERICAL_SIDE_VALUE_ZERO_HPP
