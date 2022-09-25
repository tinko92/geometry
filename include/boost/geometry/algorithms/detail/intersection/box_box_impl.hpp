// Boost.Geometry (aka GGL, Generic Geometry Library)

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_ALGORITHMS_DETAIL_INTERSECTION_BOX_BOX_IMPL_HPP
#define BOOST_GEOMETRY_ALGORITHMS_DETAIL_INTERSECTION_BOX_BOX_IMPL_HPP


#include <boost/geometry/core/access.hpp>

#include <boost/geometry/util/algorithm.hpp>

namespace boost { namespace geometry
{


#ifndef DOXYGEN_NO_DETAIL
namespace detail { namespace intersection
{

struct intersection_box_box
{
    template
    <
        typename Box1,
        typename Box2,
        typename BoxOut,
        typename Strategy
    >
    static inline bool apply(Box1 const& box1,
                             Box2 const& box2,
                             BoxOut& box_out,
                             Strategy const&)
    {
        return detail::all_dimensions_of<Box1>([&](auto index)
        {
            auto max1 = get<max_corner, index>(box1);
            auto min2 = get<min_corner, index>(box2);

            if (max1 < min2)
            {
                return false;
            }

            auto max2 = get<max_corner, index>(box2);
            auto min1 = get<min_corner, index>(box1);

            if (max2 < min1)
            {
                return false;
            }

            // Set dimensions of output coordinate
            set<min_corner, index>(box_out, min1 < min2 ? min2 : min1);
            set<max_corner, index>(box_out, max1 > max2 ? max2 : max1);

            return true;
        });
    }
};


}} // namespace detail::intersection
#endif // DOXYGEN_NO_DETAIL


}} // namespace boost::geometry


#endif // BOOST_GEOMETRY_ALGORITHMS_DETAIL_INTERSECTION_BOX_BOX_IMPL_HPP
