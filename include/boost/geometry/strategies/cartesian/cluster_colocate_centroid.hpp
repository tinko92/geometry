// Boost.Geometry (aka GGL, Generic Geometry Library)

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_STRATEGY_CARTESIAN_CLUSTER_COLOCATE_CENTROID_HPP
#define BOOST_GEOMETRY_STRATEGY_CARTESIAN_CLUSTER_COLOCATE_CENTROID_HPP

#include <iterator>

#include <boost/geometry/core/access.hpp>

namespace boost { namespace geometry
{

namespace strategy { namespace cluster_colocate
{

/*!
\brief Assigns all points in a cluster to their centroid.
\ingroup strategies
*/

class centroid
{
public :
    template <typename PointIt>
    static inline void apply(PointIt begin, PointIt end)
    {
        auto it = begin;
        auto centroid_0 = geometry::get<0>(*it);
        auto centroid_1 = geometry::get<1>(*it);
        for (++it; it != end; ++it)
        {
            centroid_0 += geometry::get<0>(*it);
            centroid_1 += geometry::get<1>(*it);
        }
        centroid_0 /= std::distance(begin, end);
        centroid_1 /= std::distance(begin, end);
        for (auto it = begin; it != end; ++it)
        {
            geometry::set<0>(*it, centroid_0);
            geometry::set<1>(*it, centroid_1);
        }
    }
};

}} // namespace strategy::cluster_colocate


}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_STRATEGY_CARTESIAN_CLUSTER_COLOCATE_CENTROID_HPP
