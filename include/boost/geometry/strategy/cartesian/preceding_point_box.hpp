// Boost.Geometry (aka GGL, Generic Geometry Library)

// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_STRATEGY_CARTESIAN_PRECEDING_POINT_BOX_HPP
#define BOOST_GEOMETRY_STRATEGY_CARTESIAN_PRECEDING_POINT_BOX_HPP

#include <boost/geometry/core/access.hpp>
#include <boost/geometry/core/tags.hpp>

#include <boost/geometry/strategy/preceding.hpp>

namespace boost { namespace geometry
{

namespace strategy { namespace preceding
{

struct cartesian_point_box
{
    template <std::size_t Dimension, typename Point, typename Box>
    static bool apply(int dir, Point const& point, Box const& /*point_box*/, Box const& other_box)
    {
        return (dir == 1  && get<Dimension>(point) < get<min_corner, Dimension>(other_box))
            || (dir == -1 && get<Dimension>(point) > get<max_corner, Dimension>(other_box));
    }
};


#ifndef DOXYGEN_NO_STRATEGY_SPECIALIZATIONS

namespace services
{

template <typename CalculationType>
struct default_strategy<cartesian_tag, CalculationType>
{
    typedef cartesian_point_box type;
};


} // namespace services

#endif // DOXYGEN_NO_STRATEGY_SPECIALIZATIONS


}} // namespace strategy::preceding

}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_STRATEGY_CARTESIAN_PRECEDING_POINT_BOX_HPP
