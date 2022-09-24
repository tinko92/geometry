// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2007-2015 Barend Gehrels, Amsterdam, the Netherlands.

// This file was modified by Oracle on 2015-2020.
// Modifications copyright (c) 2015-2020 Oracle and/or its affiliates.

// Contributed and/or modified by Menelaos Karavelas, on behalf of Oracle
// Contributed and/or modified by Adam Wulkiewicz, on behalf of Oracle

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_ALGORITHMS_DETAIL_OVERLAY_CLIP_LINESTRING_HPP
#define BOOST_GEOMETRY_ALGORITHMS_DETAIL_OVERLAY_CLIP_LINESTRING_HPP

#include <boost/range/begin.hpp>
#include <boost/range/empty.hpp>
#include <boost/range/end.hpp>

#include <boost/geometry/algorithms/clear.hpp>
#include <boost/geometry/algorithms/convert.hpp>

#include <boost/geometry/algorithms/detail/overlay/append_no_duplicates.hpp>

#include <boost/geometry/geometries/segment.hpp>

namespace boost { namespace geometry
{

#ifndef DOXYGEN_NO_DETAIL
namespace detail { namespace intersection
{

/*!
    \brief Clips a linestring with a box
    \details A linestring is intersected (clipped) by the specified box
    and the resulting linestring, or pieces of linestrings, are sent to the specified output operator.
    \tparam OutputLinestring type of the output linestrings
    \tparam OutputIterator an output iterator which outputs linestrings
    \tparam Linestring linestring-type, for example a vector of points, matching the output-iterator type,
         the points should also match the input-iterator type
    \tparam Box box type
    \tparam Strategy strategy, a clipping strategy which should implement the methods "clip_segment" and "apply"
*/
template
<
    typename OutputLinestring,
    typename OutputIterator,
    typename Range,
    typename RobustPolicy,
    typename Box,
    typename Strategy
>
OutputIterator clip_range_with_box(Box const& b, Range const& range,
            RobustPolicy const&,
            OutputIterator out, Strategy const& strategy)
{
    if (boost::begin(range) == boost::end(range))
    {
        return out;
    }

    using point_type = typename point_type<OutputLinestring>::type;

    OutputLinestring line_out;

    auto vertex = boost::begin(range);
    for (auto previous = vertex++;
              vertex != boost::end(range);
              ++previous, ++vertex)
    {
        point_type p1, p2;
        geometry::convert(*previous, p1);
        geometry::convert(*vertex, p2);

        // Clip the segment. Five situations:
        // 1. Segment is invisible, finish line if any (shouldn't occur)
        // 2. Segment is completely visible. Add (p1)-p2 to line
        // 3. Point 1 is invisible (clipped), point 2 is visible. Start new line from p1-p2...
        // 4. Point 1 is visible, point 2 is invisible (clipped). End the line with ...p2
        // 5. Point 1 and point 2 are both invisible (clipped). Start/finish an independant line p1-p2
        //
        // This results in:
        // a. if p1 is clipped, start new line
        // b. if segment is partly or completely visible, add the segment
        // c. if p2 is clipped, end the line

        bool c1 = false;
        bool c2 = false;
        model::referring_segment<point_type> s(p1, p2);

        if (!strategy.intersection(b, s).apply(b, s, c1, c2))
        {
            if (!boost::empty(line_out))
            {
                *out++ = line_out;
                geometry::clear(line_out);
            }
        }
        else
        {
            // a. If necessary, finish the line and add a start a new one
            if (c1)
            {
                if (!boost::empty(line_out))
                {
                    *out++ = line_out;
                    geometry::clear(line_out);
                }
            }

            // b. Add p1 only if it is the first point, then add p2
            if (boost::empty(line_out))
            {
                detail::overlay::append_with_duplicates(line_out, p1);
            }
            detail::overlay::append_no_duplicates(line_out, p2, strategy);

            // c. If c2 is clipped, finish the line
            if (c2)
            {
                if (!boost::empty(line_out))
                {
                    *out++ = line_out;
                    geometry::clear(line_out);
                }
            }
        }

    }

    // Add last part
    if (!boost::empty(line_out))
    {
        *out++ = line_out;
        geometry::clear(line_out);
    }
    return out;
}

}} // namespace detail::intersection
#endif // DOXYGEN_NO_DETAIL

}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_ALGORITHMS_DETAIL_OVERLAY_CLIP_LINESTRING_HPP
