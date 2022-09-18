// Boost.Geometry (aka GGL, Generic Geometry Library)

// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_STRATEGY_AGNOSTIC_DEGENERATE_SEGMENT_HPP
#define BOOST_GEOMETRY_STRATEGY_AGNOSTIC_DEGENERATE_SEGMENT_HPP

#include <boost/geometry/core/access.hpp>

#include <boost/geometry/util/math.hpp>
#include <boost/geometry/util/algorithm.hpp>

namespace boost { namespace geometry
{

namespace strategy { namespace degenerate_segment
{

struct math_equals
{
    template <typename Segment>
    static bool apply(Segment const& seg)
    {
        return detail::all_dimensions_of<Segment>([&](auto index)
        {
            return geometry::math::equals
                (
                    geometry::get<0, index>(seg),
                    geometry::get<1, index>(seg)
                );

        });
    }
};

}} // namespace strategy::degenerate_segment

}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_STRATEGY_AGNOSTIC_DEGENERATE_SEGMENT_HPP
