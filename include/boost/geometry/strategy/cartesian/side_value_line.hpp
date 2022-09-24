// Boost.Geometry (aka GGL, Generic Geometry Library)

// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_STRATEGY_CARTESIAN_SIDE_VALUE_LINE_HPP
#define BOOST_GEOMETRY_STRATEGY_CARTESIAN_SIDE_VALUE_LINE_HPP

#include <boost/geometry/algorithms/detail/make/make.hpp>

#include <boost/geometry/arithmetic/infinite_line_functions.hpp>

#include <boost/geometry/util/select_coordinate_type.hpp>

namespace boost { namespace geometry
{

namespace strategy { namespace side_value_line
{

struct cartesian
{
    template <typename Point1, typename Point2>
    static inline auto apply(Point1 const& p1, Point1 const& p2, Point2 const& p)
    {
        using coord_t = typename select_coordinate_type<Point1, Point2>::type;
        auto const line = geometry::detail::make::make_infinite_line<coord_t>(p1, p2);
        return arithmetic::side_value(line, p);
    }
};

}} // namespace strategy::direction

}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_STRATEGY_CARTESIAN_SIDE_VALUE_LINE_HPP
