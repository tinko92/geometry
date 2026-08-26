// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2023-2024 Adam Wulkiewicz, Lodz, Poland.

// Copyright (c) 2018-2023 Oracle and/or its affiliates.
// Contributed and/or modified by Vissarion Fysikopoulos, on behalf of Oracle
// Contributed and/or modified by Adam Wulkiewicz, on behalf of Oracle

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_ALGORITHMS_LINE_INTERPOLATE_HPP
#define BOOST_GEOMETRY_ALGORITHMS_LINE_INTERPOLATE_HPP

#include <type_traits>

#include <boost/range/begin.hpp>
#include <boost/range/end.hpp>
#include <boost/range/value_type.hpp>
#include <boost/geometry/algorithms/detail/convert_point_to_point.hpp>
#include <boost/geometry/algorithms/detail/dummy_geometries.hpp>

#include <boost/geometry/core/exception.hpp>
#include <boost/geometry/core/static_assert.hpp>
#include <boost/geometry/core/tags.hpp>
#include <boost/geometry/core/visit.hpp>

#include <boost/geometry/geometries/concepts/check.hpp>

#include <boost/geometry/strategies/default_strategy.hpp>
#include <boost/geometry/strategies/detail.hpp>
#include <boost/geometry/strategies/line_interpolate/cartesian.hpp>
#include <boost/geometry/strategies/line_interpolate/geographic.hpp>
#include <boost/geometry/strategies/line_interpolate/spherical.hpp>

#include <boost/geometry/util/range.hpp>
#include <boost/geometry/util/type_traits.hpp>

#include <boost/geometry/views/segment_view.hpp>

namespace boost { namespace geometry
{


#ifndef DOXYGEN_NO_DETAIL
namespace detail { namespace line_interpolate
{

struct convert_and_push_back
{
    template <typename Range, typename Point>
    static inline void apply(Point const& p, Range& range)
    {
        typename boost::range_value<Range>::type p2;
        geometry::detail::conversion::convert_point_to_point(p, p2);
        range::push_back(range, p2);
    }
};

struct convert_and_assign
{
    template <typename Point1, typename Point2>
    static inline void apply(Point1 const& p1, Point2& p2)
    {
        geometry::detail::conversion::convert_point_to_point(p1, p2);
    }

};


/*!
\brief Internal, calculates interpolation point of a linestring using iterator pairs and
    specified strategy
*/
template <typename Policy>
struct interpolate_range
{
    template
    <
        typename Range,
        typename Distance,
        typename PointLike,
        typename Strategies
    >
    static inline void apply(Range const& range,
                             Distance const& max_distance,
                             PointLike & pointlike,
                             Strategies const& strategies)
    {
        typedef typename boost::range_value<Range const>::type point_t;

        auto it = boost::begin(range);
        auto const end = boost::end(range);

        if (it == end) // empty(range)
        {
            BOOST_THROW_EXCEPTION(empty_input_exception());
            return;
        }
        if (max_distance <= 0) //non positive distance
        {
            Policy::apply(*it, pointlike);
            return;
        }

        auto const pp_strategy = strategies.distance(dummy_point(), dummy_point());
        auto const strategy = strategies.line_interpolate(range);

        typedef decltype(pp_strategy.apply(
                    std::declval<point_t>(), std::declval<point_t>())) distance_type;

        auto prev = it++;
        distance_type repeated_distance = max_distance;
        distance_type prev_distance = 0;
        distance_type current_distance = 0;
        point_t start_p = *prev;

        for ( ; it != end ; ++it)
        {
            distance_type dist = pp_strategy.apply(*prev, *it);
            current_distance = prev_distance + dist;

            while (current_distance >= repeated_distance)
            {
                point_t p;
                distance_type diff_distance = current_distance - prev_distance;
                BOOST_ASSERT(diff_distance != distance_type(0));
                strategy.apply(start_p, *it,
                               (repeated_distance - prev_distance)/diff_distance,
                               p,
                               diff_distance);
                Policy::apply(p, pointlike);
                if constexpr (util::is_point<PointLike>::value)
                {
                    return;
                }
                else // else prevents unreachable code warning
                {
                    start_p = p;
                    prev_distance = repeated_distance;
                    repeated_distance += max_distance;
                }
            }
            prev_distance = current_distance;
            prev = it;
            start_p = *prev;
        }

        // case when max_distance is larger than linestring's length
        // return the last point in range (range is not empty)
        if (repeated_distance == max_distance)
        {
            Policy::apply(*(end-1), pointlike);
        }
    }
};

template <typename Policy>
struct interpolate_segment
{
    template <typename Segment, typename Distance, typename Pointlike, typename Strategy>
    static inline void apply(Segment const& segment,
                             Distance const& max_distance,
                             Pointlike & point,
                             Strategy const& strategy)
    {
        interpolate_range<Policy>().apply(segment_view<Segment>(segment),
                                          max_distance, point, strategy);
    }
};

}} // namespace detail::line_interpolate
#endif // DOXYGEN_NO_DETAIL


#ifndef DOXYGEN_NO_DISPATCH
namespace dispatch
{


template <concepts::ConstGeometry Geometry,
          concepts::MutableGeometry Pointlike,
          typename Distance,
          typename Strategies>
    requires (concepts::ConstLinestring<Geometry>
              || concepts::ConstSegment<Geometry>)
          && (concepts::Point<Pointlike>
              || concepts::MultiPoint<Pointlike>)
inline void line_interpolate(Geometry const& geometry,
                             Distance const& max_distance,
                             Pointlike& pointlike,
                             Strategies const& strategies)
{
    using policy = std::conditional_t
        <
            concepts::Point<Pointlike>,
            detail::line_interpolate::convert_and_assign,
            detail::line_interpolate::convert_and_push_back
        >;

    if constexpr (concepts::ConstLinestring<Geometry>)
    {
        detail::line_interpolate::interpolate_range<policy>::apply(
            geometry, max_distance, pointlike, strategies);
    }
    else
    {
        detail::line_interpolate::interpolate_segment<policy>::apply(
            geometry, max_distance, pointlike, strategies);
    }
}

} // namespace dispatch
#endif // DOXYGEN_NO_DISPATCH


namespace resolve_strategy {

template <concepts::ConstGeometry Geometry,
          typename Distance,
          concepts::MutableGeometry Pointlike,
          typename Strategy>
inline void line_interpolate(Geometry const& geometry,
                             Distance const& max_distance,
                             Pointlike& pointlike,
                             Strategy const& strategy)
{
    if constexpr (std::same_as<Strategy, default_strategy>)
    {
        using strategy_type = typename strategies::line_interpolate::services
            ::default_strategy<Geometry>::type;
        dispatch::line_interpolate(
            geometry, max_distance, pointlike, strategy_type());
    }
    else if constexpr (strategies::detail::is_umbrella_strategy<Strategy>::value)
    {
        dispatch::line_interpolate(
            geometry, max_distance, pointlike, strategy);
    }
    else
    {
        using strategies::line_interpolate::services::strategy_converter;
        dispatch::line_interpolate(
            geometry, max_distance, pointlike,
            strategy_converter<Strategy>::get(strategy));
    }
}

} // namespace resolve_strategy


namespace resolve_dynamic {

template <concepts::ConstGeometry Geometry,
          typename Distance,
          concepts::MutableGeometry Pointlike,
          typename Strategy>
inline void line_interpolate(Geometry const& geometry,
                             Distance const& max_distance,
                             Pointlike& pointlike,
                             Strategy const& strategy)
{
    if constexpr (concepts::ConstDynamicGeometry<Geometry>)
    {
        traits::visit<Geometry>::apply([&](auto const& g)
        {
            resolve_dynamic::line_interpolate(
                g, max_distance, pointlike, strategy);
        }, geometry);
    }
    else
    {
        resolve_strategy::line_interpolate(
            geometry, max_distance, pointlike, strategy);
    }
}

} // namespace resolve_dynamic

/*!
\brief     Returns one or more points interpolated along a LineString \brief_strategy
\ingroup line_interpolate
\tparam Geometry Any type fulfilling a LineString concept
\tparam Distance A numerical distance measure
\tparam Pointlike Any type fulfilling Point or Multipoint concept
\tparam Strategy A type fulfilling a LineInterpolatePointStrategy concept
\param geometry Input geometry
\param max_distance Distance threshold (in units depending on coordinate system)
representing the spacing between the points
\param pointlike Output: either a Point (exactly one point will be constructed) or
a MultiPoint (depending on the max_distance one or more points will be constructed)
\param strategy line_interpolate strategy to be used for interpolation of
points

\qbk{[include reference/algorithms/line_interpolate.qbk]}

\qbk{distinguish,with strategy}

\qbk{
[heading Available Strategies]
\* [link geometry.reference.strategies.strategy_line_interpolate_cartesian Cartesian]
\* [link geometry.reference.strategies.strategy_line_interpolate_spherical Spherical]
\* [link geometry.reference.strategies.strategy_line_interpolate_geographic Geographic]

[heading Example]
[line_interpolate_strategy]
[line_interpolate_strategy_output]

[heading See also]
\* [link geometry.reference.algorithms.densify densify]
}
 */
template <concepts::ConstGeometry Geometry,
          typename Distance,
          concepts::MutableGeometry Pointlike,
          typename Strategy>
    requires (concepts::Point<Pointlike>
              || concepts::MultiPoint<Pointlike>)
inline void line_interpolate(Geometry const& geometry,
                             Distance const& max_distance,
                             Pointlike & pointlike,
                             Strategy const& strategy)
{
    // detail::throw_on_empty_input(geometry);

    return resolve_dynamic::line_interpolate(
        geometry, max_distance, pointlike, strategy);
}


/*!
\brief     Returns one or more points interpolated along a LineString.
\ingroup line_interpolate
\tparam Geometry Any type fulfilling a LineString concept
\tparam Distance A numerical distance measure
\tparam Pointlike Any type fulfilling Point or Multipoint concept
\param geometry Input geometry
\param max_distance Distance threshold (in units depending on coordinate system)
representing the spacing between the points
\param pointlike Output: either a Point (exactly one point will be constructed) or
a MultiPoint (depending on the max_distance one or more points will be constructed)

\qbk{[include reference/algorithms/line_interpolate.qbk]

[heading Example]
[line_interpolate]
[line_interpolate_output]

[heading See also]
\* [link geometry.reference.algorithms.densify densify]
}
 */
template <concepts::ConstGeometry Geometry,
          typename Distance,
          concepts::MutableGeometry Pointlike>
    requires (concepts::Point<Pointlike>
              || concepts::MultiPoint<Pointlike>)
inline void line_interpolate(Geometry const& geometry,
                             Distance const& max_distance,
                             Pointlike & pointlike)
{
    // detail::throw_on_empty_input(geometry);

    return resolve_dynamic::line_interpolate(
        geometry, max_distance, pointlike, default_strategy());
}

}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_ALGORITHMS_LINE_INTERPOLATE_HPP
