// Boost.Geometry (aka GGL, Generic Geometry Library)

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_STRATEGY_AGNOSTIC_CLUSTER_COLOCATE_FIRST_HPP
#define BOOST_GEOMETRY_STRATEGY_AGNOSTIC_CLUSTER_COLOCATE_FIRST_HPP


namespace boost { namespace geometry
{

namespace strategy { namespace cluster_colocate
{

/*!
\brief Assigns to a range of points that form a cluster in overlay the first point.
\ingroup strategies
*/

class first
{
public :
    template <typename PointIt>
    static inline void apply(PointIt begin, PointIt end)
    {
        auto it = begin;
        auto const& first_point = *it++;
        for (; it != end; ++it)
            *it = first_point;
    }
};

}} // namespace strategy::cluster_colocate


}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_STRATEGY_AGNOSTIC_CLUSTER_COLOCATE_FIRST_HPP
