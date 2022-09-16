// Boost.Geometry (aka GGL, Generic Geometry Library)

// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_STRATEGY_SPHERICAL_DIRECTION_HPP
#define BOOST_GEOMETRY_STRATEGY_SPHERICAL_DIRECTION_HPP

#include <boost/geometry/core/access.hpp>

#include <boost/geometry/util/math.hpp>
#include <boost/geometry/util/normalize_spheroidal_coordinates.hpp>
#include <boost/geometry/util/select_coordinate_type.hpp>

namespace boost { namespace geometry
{

namespace strategy { namespace direction
{

struct spherical_equatorial
{
    template <typename Point1, typename Point2>
    static inline int apply(Point1 const& segment_a, Point1 const& segment_b, Point2 const& p)
    {
        using coord1_t = typename coordinate_type<Point1>::type;
        using coord2_t = typename coordinate_type<Point2>::type;
        using units_t = typename geometry::detail::cs_angular_units<Point1>::type;
        using units2_t = typename geometry::detail::cs_angular_units<Point2>::type;
        BOOST_GEOMETRY_STATIC_ASSERT(
            (std::is_same<units_t, units2_t>::value),
            "Not implemented for different units.",
            units_t, units2_t);

        using calc_t = typename geometry::select_coordinate_type <Point1, Point2>::type;
        using constants1 = math::detail::constants_on_spheroid<coord1_t, units_t>;
        using constants2 = math::detail::constants_on_spheroid<coord2_t, units_t>;
        static coord1_t const pi_half1 = constants1::max_latitude();
        static coord2_t const pi_half2 = constants2::max_latitude();
        static calc_t const c0 = 0;
        
        auto const a0 = geometry::get<0>(segment_a);
        auto const a1 = geometry::get<1>(segment_a);
        auto const b0 = geometry::get<0>(segment_b);
        auto const b1 = geometry::get<1>(segment_b);
        auto const p0 = geometry::get<0>(p);   
        auto const p1 = geometry::get<1>(p);
                
        if ( (math::equals(b0, a0) && math::equals(b1, a1))
          || (math::equals(b0, p0) && math::equals(b1, p1)) )
        {
            return 0;
        }

        bool const is_a_pole = math::equals(pi_half1, math::abs(a1));
        bool const is_b_pole = math::equals(pi_half1, math::abs(b1));
        bool const is_p_pole = math::equals(pi_half2, math::abs(p1));

        if ( is_b_pole && ((is_a_pole && math::sign(b1) == math::sign(a1))
                        || (is_p_pole && math::sign(b1) == math::sign(p1))) )
        {
            return 0;
        }

        // NOTE: as opposed to the implementation for cartesian CS
        // here point b is the origin

        calc_t const dlon1 = math::longitude_distance_signed<units_t, calc_t>(b0, a0);
        calc_t const dlon2 = math::longitude_distance_signed<units_t, calc_t>(b0, p0);

        bool is_antilon1 = false, is_antilon2 = false;
        calc_t const dlat1 = latitude_distance_signed<units_t, calc_t>(b1, a1, dlon1, is_antilon1);
        calc_t const dlat2 = latitude_distance_signed<units_t, calc_t>(b1, p1, dlon2, is_antilon2);

        calc_t mx = is_a_pole || is_b_pole || is_p_pole ?
                    c0 :
                    (std::min)(is_antilon1 ? c0 : math::abs(dlon1),
                               is_antilon2 ? c0 : math::abs(dlon2));
        calc_t my = (std::min)(math::abs(dlat1),
                               math::abs(dlat2));

        int s1 = 0, s2 = 0;
        if (mx >= my)
        {
            s1 = dlon1 > 0 ? 1 : -1;
            s2 = dlon2 > 0 ? 1 : -1;
        }
        else
        {
            s1 = dlat1 > 0 ? 1 : -1;
            s2 = dlat2 > 0 ? 1 : -1;
        }

        return s1 == s2 ? -1 : 1;
    }

private:
    template <typename Units, typename T>
    static inline T latitude_distance_signed(T const& lat1, T const& lat2, T const& lon_ds,
                                             bool & is_antilon)
    {
        using constants = math::detail::constants_on_spheroid<T, Units>;
        static T const pi = constants::half_period();
        static T const c0 = 0;
        T res = lat2 - lat1;

        is_antilon = math::equals(math::abs(lon_ds), pi);
        if (is_antilon)
        {
            res = lat2 + lat1;
            if (res >= c0)
                res = pi - res;
            else
                res = -pi - res;
        }
        return res;
    }
};

struct spherical_polar
{
    template <typename Point1, typename Point2>
    static inline int apply(Point1 segment_a, Point1 segment_b,
                            Point2 p)
    {
        using constants1 = math::detail::constants_on_spheroid
            <
                typename coordinate_type<Point1>::type,
                typename geometry::detail::cs_angular_units<Point1>::type
            >;
        using constants2 = math::detail::constants_on_spheroid
            <
                typename coordinate_type<Point2>::type,
                typename geometry::detail::cs_angular_units<Point2>::type
            >;

        geometry::set<1>(segment_a, constants1::max_latitude() - geometry::get<1>(segment_a));
        geometry::set<1>(segment_b, constants1::max_latitude() - geometry::get<1>(segment_b));
        geometry::set<1>(p, constants2::max_latitude() - geometry::get<1>(p));

        return spherical_equatorial::apply(segment_a, segment_b, p);
    }
};

}} // namespace strategy::direction

}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_STRATEGY_SPHERICAL_DIRECTION_HPP
