// Boost.Geometry

// Copyright (c) 2007-2012 Barend Gehrels, Amsterdam, the Netherlands.

// This file was modified by Oracle on 2014-2023.
// Modifications copyright (c) 2014-2023, Oracle and/or its affiliates.
// Contributed and/or modified by Vissarion Fysikopoulos, on behalf of Oracle
// Contributed and/or modified by Adam Wulkiewicz, on behalf of Oracle

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_ALGORITHMS_AZIMUTH_HPP
#define BOOST_GEOMETRY_ALGORITHMS_AZIMUTH_HPP


#include <boost/geometry/algorithms/not_implemented.hpp>

#include <boost/geometry/core/radian_access.hpp>
#include <boost/geometry/core/tag.hpp>

#include <boost/geometry/geometries/concepts/check.hpp>

#include <boost/geometry/strategies/default_strategy.hpp>
#include <boost/geometry/strategies/azimuth/cartesian.hpp>
#include <boost/geometry/strategies/azimuth/geographic.hpp>
#include <boost/geometry/strategies/azimuth/spherical.hpp>

namespace boost { namespace geometry
{


#ifndef DOXYGEN_NO_DETAIL
namespace detail
{

} // namespace detail
#endif // DOXYGEN_NO_DETAIL


#ifndef DOXYGEN_NO_DISPATCH
namespace dispatch
{

template <concepts::ConstPoint Point1, concepts::ConstPoint Point2,
          typename Strategy>
inline auto azimuth(Point1 const& p1, Point2 const& p2,
                    Strategy const& strategy)
{
    auto azimuth_strategy = strategy.azimuth();
    using calc_t = typename decltype(azimuth_strategy)::template result_type
        <
            coordinate_type_t<Point1>,
            coordinate_type_t<Point2>
        >::type;

    calc_t result = 0;
    calc_t const x1 = geometry::get_as_radian<0>(p1);
    calc_t const y1 = geometry::get_as_radian<1>(p1);
    calc_t const x2 = geometry::get_as_radian<0>(p2);
    calc_t const y2 = geometry::get_as_radian<1>(p2);

    azimuth_strategy.apply(x1, y1, x2, y2, result);

        // NOTE: It is not clear which units we should use for the result.
        //   For now radians are always returned but a user could expect
        //   e.g. something like this:
        /*
        bool const both_degree = std::is_same
                <
                    typename detail::cs_angular_units<Point1>::type,
                    geometry::degree
                >::value
            && std::is_same
                <
                    typename detail::cs_angular_units<Point2>::type,
                    geometry::degree
                >::value;
        if (both_degree)
        {
            result *= math::r2d<calc_t>();
        }
        */

    return result;
}

} // namespace dispatch
#endif // DOXYGEN_NO_DISPATCH


namespace resolve_strategy
{

template <concepts::ConstPoint Point1, concepts::ConstPoint Point2,
          typename Strategy>
inline auto azimuth(Point1 const& point1, Point2 const& point2,
                    Strategy const& strategy)
{
    if constexpr (std::same_as<Strategy, default_strategy>)
    {
        using strategy_type = typename strategies::azimuth::services
            ::default_strategy<Point1, Point2>::type;
        return dispatch::azimuth(point1, point2, strategy_type());
    }
    else if constexpr (strategies::detail::is_umbrella_strategy<Strategy>::value)
    {
        return dispatch::azimuth(point1, point2, strategy);
    }
    else
    {
        using strategies::azimuth::services::strategy_converter;
        return dispatch::azimuth(point1, point2,
            strategy_converter<Strategy>::get(strategy));
    }
}


} // namespace resolve_strategy


namespace resolve_variant
{
} // namespace resolve_variant


/*!
\brief Calculate azimuth of a segment defined by a pair of points.
\ingroup azimuth
\tparam Point1 Type of the first point of a segment.
\tparam Point2 Type of the second point of a segment.
\param point1 First point of a segment.
\param point2 Second point of a segment.
\return Azimuth in radians.

\qbk{[include reference/algorithms/azimuth.qbk]}

\qbk{
[heading Example]
[azimuth]
[azimuth_output]
}
*/
template <concepts::ConstPoint Point1, concepts::ConstPoint Point2>
inline auto azimuth(Point1 const& point1, Point2 const& point2)
{
    return resolve_strategy::azimuth(point1, point2, default_strategy());
}


/*!
\brief Calculate azimuth of a segment defined by a pair of points.
\ingroup azimuth
\tparam Point1 Type of the first point of a segment.
\tparam Point2 Type of the second point of a segment.
\tparam Strategy Type of an umbrella strategy defining azimuth strategy.
\param point1 First point of a segment.
\param point2 Second point of a segment.
\param strategy Umbrella strategy defining azimuth strategy.
\return Azimuth in radians.

\qbk{distinguish,with strategy}
\qbk{[include reference/algorithms/azimuth.qbk]}

\qbk{
[heading Example]
[azimuth_strategy]
[azimuth_strategy_output]
}
*/
template <concepts::ConstPoint Point1, concepts::ConstPoint Point2,
          typename Strategy>
inline auto azimuth(Point1 const& point1, Point2 const& point2, Strategy const& strategy)
{
    return resolve_strategy::azimuth(point1, point2, strategy);
}


}} // namespace boost::geometry


#endif // BOOST_GEOMETRY_ALGORITHMS_AZIMUTH_HPP
