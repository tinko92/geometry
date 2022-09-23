// Boost.Geometry (aka GGL, Generic Geometry Library)

// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_STRATEGY_AGNOSTIC_CLUSTER_HPP
#define BOOST_GEOMETRY_STRATEGY_AGNOSTIC_CLUSTER_HPP

#include <boost/geometry/core/access.hpp>

#include <boost/geometry/algorithms/detail/overlay/approximately_equals.hpp>

#include <boost/geometry/util/math.hpp>
#include <boost/geometry/util/algorithm.hpp>

namespace boost { namespace geometry
{

namespace strategy { namespace cluster
{

struct exact
{
    template <typename Point>
    static bool equals(Point const& p1, Point const& p2)
    {
        return detail::all_dimensions_of<Point>([&](auto index)
        {
            return get<index>(p1) == get<index>(p2);
        });
    }

    template <std::size_t Dimension, typename Point>
    static bool exceeds(Point const& p1, Point const& p2)
    {
        return get<Dimension>(p1) > get<Dimension>(p2);
    }
};

struct integral_neighbours
{
    template <typename Point>
    static bool equals(Point const& p1, Point const& p2)
    {
        return detail::all_dimensions_of<Point>([&](auto index)
        {
            return math::abs(get<index>(p1) - get<index>(p2)) <= 1;
        });
    }

    template <std::size_t Dimension, typename Point>
    static bool exceeds(Point const& p1, Point const& p2)
    {
        return get<Dimension>(p1) - get<Dimension>(p2) > 1;
    }
};

struct approx_equal
{
    template <typename Point>
    static bool equals(Point const& p1, Point const& p2)
    {
        using ct = typename coordinate_type<Point>::type;
        return detail::overlay::approximately_equals(p1, p2, ct(1000));
    }

    template <std::size_t Dimension, typename Point>
    static bool exceeds(Point const& p1, Point const& p2)
    {
        using ct = typename coordinate_type<Point>::type;
        return get<Dimension>(p1) - get<Dimension>(p2) > (ct(1)/ct(1000));
    }
};

}} // namespace strategy::cluster

}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_STRATEGY_AGNOSTIC_CLUSTER_HPP
