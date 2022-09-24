// Boost.Geometry (aka GGL, Generic Geometry Library)

// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_STRATEGY_AGNOSTIC_SUCCESSOR_APPROX_RANGE_POINT_HPP
#define BOOST_GEOMETRY_STRATEGY_AGNOSTIC_SUCCESSOR_APPROX_RANGE_POINT_HPP

#include <boost/geometry/core/tag.hpp>
#include <boost/geometry/core/tags.hpp>

#include <boost/geometry/strategy/detail/approximately_equals.hpp>

#include <boost/geometry/algorithms/detail/overlay/copy_segment_point.hpp>
#include <boost/geometry/algorithms/detail/overlay/segment_identifier.hpp>

#include <boost/geometry/util/condition.hpp>

namespace boost { namespace geometry
{

namespace strategy { namespace successor
{

template <typename CalculationType, int ToleranceMultiple = 1000000000, int Limit = 10>
struct range_point_tolerance
{

    template
    <
        bool Reverse,
        bool Forward = true,
        typename PointOut,
        typename Geometry,
        typename Point
    >
    static void apply(segment_identifier const& seg_id,
                      Geometry const& range,
                      Point const& point,
                      PointOut& out,
                      int offset = 0)
    {
        static CalculationType const tolerance = static_cast<CalculationType>(ToleranceMultiple);
        using csp = dispatch::copy_segment_point
            <
                typename tag<Geometry>::type,
                Geometry,
                Reverse,
                segment_identifier,
                PointOut
            >;
        if ( BOOST_GEOMETRY_CONDITION(Forward) )
        {
            csp::apply(range, seg_id, offset, out);
            while (detail::approximately_equals(out, point, tolerance) && offset < Limit)
            {
                csp::apply(range, seg_id, ++offset, out);
            }
        }
        else
        {
            csp::apply(range, seg_id, offset, out);
            while (detail::approximately_equals(out, point, tolerance) && offset > -Limit)
            {
                csp::apply(range, seg_id, --offset, out);
            }
        }
    }
};

}} // namespace strategy::successor

}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_STRATEGY_AGNOSTIC_SUCCESSOR_APPROX_RANGE_POINT_HPP
