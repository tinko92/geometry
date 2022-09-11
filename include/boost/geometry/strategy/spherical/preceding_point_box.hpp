// Boost.Geometry (aka GGL, Generic Geometry Library)

// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_STRATEGY_SPHERICAL_PRECEDING_POINT_BOX_HPP
#define BOOST_GEOMETRY_STRATEGY_SPHERICAL_PRECEDING_POINT_BOX_HPP

#include <boost/geometry/core/access.hpp>
#include <boost/geometry/core/coordinate_type.hpp>
#include <boost/geometry/core/cs.hpp>
#include <boost/geometry/core/tags.hpp>

#include <boost/geometry/strategies/cartesian/point_in_box.hpp>

#include <boost/geometry/strategy/preceding.hpp>

#include <boost/geometry/util/normalize_spheroidal_coordinates.hpp>
#include <boost/geometry/util/select_coordinate_type.hpp>

namespace boost { namespace geometry
{

namespace strategy { namespace preceding
{

struct spherical_point_box
{
    template <std::size_t Dimension, typename Point, typename Box>
    static bool apply(int dir, Point const& point, Box const& point_box, Box const& other_box)
    {
        if ( BOOST_GEOMETRY_CONDITION(Dimension == 0) )
        {
            using calc_t = typename select_coordinate_type
                <    
                    Point, Box
                >::type;
            using units_t = typename detail::cs_angular_units<Point>::type;
                          
            calc_t const c0 = 0;
                          
            calc_t const value = get<0>(point);
            calc_t const other_min = get<min_corner, 0>(other_box);
            calc_t const other_max = get<max_corner, 0>(other_box);
     
            bool const pt_covered = strategy::within::detail::covered_by_range
                                        <
                                            Point, 0, spherical_tag
                                        >::apply(value,
                                                 other_min,
                                                 other_max);

            if (pt_covered)
            {
                return false;
            }

            if (dir == 1)
            {
                calc_t const diff_min = math::longitude_distance_signed
                                            <
                                                units_t, calc_t
                                            >(other_min, value);

                calc_t const diff_min_min = math::longitude_distance_signed
                                            <
                                                units_t, calc_t
                                            >(other_min, get<min_corner, 0>(point_box));

                return diff_min < c0 && diff_min_min <= c0 && diff_min_min <= diff_min;
            }
            else if (dir == -1)
            {
                calc_t const diff_max = math::longitude_distance_signed
                                            <
                                                units_t, calc_t
                                            >(other_max, value);

                calc_t const diff_max_max = math::longitude_distance_signed
                                            <
                                                units_t, calc_t
                                            >(other_max, get<max_corner, 0>(point_box));

                return diff_max > c0 && diff_max_max >= c0 && diff_max <= diff_max_max;
            }

            return false;
        }
        else
        {
            return (dir == 1  && get<Dimension>(point) < get<min_corner, Dimension>(other_box))
                || (dir == -1 && get<Dimension>(point) > get<max_corner, Dimension>(other_box));
        }
    }
};


#ifndef DOXYGEN_NO_STRATEGY_SPECIALIZATIONS

namespace services
{

template <typename CalculationType>
struct default_strategy<spherical_equatorial_tag, CalculationType>
{
    typedef spherical_point_box type;
};

template <typename CalculationType>
struct default_strategy<spherical_polar_tag, CalculationType>
{
    typedef spherical_point_box type;
};

template <typename CalculationType>
struct default_strategy<geographic_tag, CalculationType>
{
    typedef spherical_point_box type;
};

} // namespace services

#endif // DOXYGEN_NO_STRATEGY_SPECIALIZATIONS


}} // namespace strategy::preceding

}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_STRATEGY_SPHERICAL_PRECEDING_POINT_BOX_HPP
