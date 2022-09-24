// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2007-2015 Barend Gehrels, Amsterdam, the Netherlands.

// This file was modified by Oracle on 2015-2020.
// Modifications copyright (c) 2015-2020 Oracle and/or its affiliates.

// Contributed and/or modified by Menelaos Karavelas, on behalf of Oracle
// Contributed and/or modified by Adam Wulkiewicz, on behalf of Oracle

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_STRATEGY_CARTESIAN_INTERSECTION_HPP
#define BOOST_GEOMETRY_STRATEGY_CARTESIAN_INTERSECTION_HPP

#include <boost/geometry/core/access.hpp>

#include <boost/geometry/util/select_coordinate_type.hpp>
#include <boost/geometry/util/select_most_precise.hpp>

namespace boost { namespace geometry
{

namespace strategy { namespace intersection
{

/*!
    \brief Strategy: line clipping algorithm after Liang Barsky
    \ingroup overlay
    \details The Liang-Barsky line clipping algorithm clips a line with a clipping box.
    It is slightly adapted in the sense that it returns which points are clipped
    \tparam B input box type of clipping box
    \tparam P input/output point-type of segments to be clipped
    \note The algorithm is currently only implemented for 2D Cartesian points
    \author Barend Gehrels, and the following recourses
    - A tutorial: http://www.skytopia.com/project/articles/compsci/clipping.html
    - a German applet (link broken): http://ls7-www.cs.uni-dortmund.de/students/projectgroups/acit/lineclip.shtml
*/
class liang_barsky
{
private:
    template <typename CoordinateType, typename CalcType>
    static inline bool check_edge(CoordinateType const& p, CoordinateType const& q, CalcType& t1, CalcType& t2)
    {
        bool visible = true;

        if(p < 0)
        {
            CalcType const r = static_cast<CalcType>(q) / p;
            if (r > t2)
                visible = false;
            else if (r > t1)
                t1 = r;
        }
        else if(p > 0)
        {
            CalcType const r = static_cast<CalcType>(q) / p;
            if (r < t1)
                visible = false;
            else if (r < t2)
                t2 = r;
        }
        else
        {
            if (q < 0)
                visible = false;
        }

        return visible;
    }

public:
    template <typename Box, typename Segment>
    static inline bool apply(Box const& b, Segment& s, bool& sp1_clipped, bool& sp2_clipped)
    {
        using coordinate_type = typename select_coordinate_type<Box, Segment>::type;
        using calc_type = typename select_most_precise<coordinate_type, double>::type;

        calc_type t1 = 0;
        calc_type t2 = 1;

        coordinate_type const dx = get<1, 0>(s) - get<0, 0>(s);
        coordinate_type const dy = get<1, 1>(s) - get<0, 1>(s);

        coordinate_type const p1 = -dx;
        coordinate_type const p2 = dx;
        coordinate_type const p3 = -dy;
        coordinate_type const p4 = dy;

        coordinate_type const q1 = get<0, 0>(s) - get<min_corner, 0>(b);
        coordinate_type const q2 = get<max_corner, 0>(b) - get<0, 0>(s);
        coordinate_type const q3 = get<0, 1>(s) - get<min_corner, 1>(b);
        coordinate_type const q4 = get<max_corner, 1>(b) - get<0, 1>(s);

        if (check_edge(p1, q1, t1, t2)      // left
            && check_edge(p2, q2, t1, t2)   // right
            && check_edge(p3, q3, t1, t2)   // bottom
            && check_edge(p4, q4, t1, t2))   // top
        {
            sp1_clipped = t1 > 0;
            sp2_clipped = t2 < 1;

            if (sp2_clipped)
            {
                set<1, 0>(s, get<0, 0>(s) + t2 * dx);
                set<1, 1>(s, get<0, 1>(s) + t2 * dy);
            }

            if(sp1_clipped)
            {
                set<0, 0>(s, get<0, 0>(s) + t1 * dx);
                set<0, 1>(s, get<0, 1>(s) + t1 * dy);
            }

            return true;
        }

        return false;
    }
};


}} // namespace strategy::intersection

}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_STRATEGY_CARTESIAN_INTERSECTION_HPP
