// Boost.Geometry (aka GGL, Generic Geometry Library)

// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_STRATEGY_CARTESIAN_DIRECTION_HPP
#define BOOST_GEOMETRY_STRATEGY_CARTESIAN_DIRECTION_HPP

#include <boost/geometry/algorithms/detail/make/make.hpp>

#include <boost/geometry/arithmetic/infinite_line_functions.hpp>

#include <boost/geometry/util/select_coordinate_type.hpp>

namespace boost { namespace geometry
{

namespace strategy { namespace direction
{

struct cartesian
{
    template <typename Point1, typename Point2>
    static inline int apply(Point1 const& segment_a, Point1 const& segment_b, Point2 const& point)
    {
        using calc_t = typename geometry::select_coordinate_type<Point1, Point2>::type;
        
        using line_type = model::infinite_line<calc_t>;
                    
        // Situation and construction of perpendicular line
        //
        //     P1     a--------------->b   P2    
        //                             |
        //                             |
        //                             v
        //
        // P1 is located right of the (directional) perpendicular line
        // and therefore gets a negative side_value, and returns -1.
        // P2 is to the left of the perpendicular line and returns 1.
        // If the specified point is located on top of b, it returns 0.
         
        line_type const line
            = detail::make::make_perpendicular_line<calc_t>(segment_a, segment_b, segment_b);
    
        if (arithmetic::is_degenerate(line))
        {
            return 0;
        }
    
        calc_t const sv = arithmetic::side_value(line, point);
        static calc_t const zero = 0;
        return sv == zero ? 0 : sv > zero ? 1 : -1;  
    }
};

}} // namespace strategy::direction

}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_STRATEGY_CARTESIAN_DIRECTION_HPP
