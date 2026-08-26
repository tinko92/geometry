// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2007-2014 Barend Gehrels, Amsterdam, the Netherlands.
// Copyright (c) 2008-2014 Bruno Lalande, Paris, France.
// Copyright (c) 2009-2014 Mateusz Loskot, London, UK.

// This file was modified by Oracle on 2014-2020.
// Modifications copyright (c) 2014-2020, Oracle and/or its affiliates.

// Contributed and/or modified by Menelaos Karavelas, on behalf of Oracle
// Contributed and/or modified by Adam Wulkiewicz, on behalf of Oracle

// Parts of Boost.Geometry are redesigned from Geodan's Geographic Library
// (geolib/GGL), copyright (c) 1995-2010 Geodan, Amsterdam, the Netherlands.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_STRATEGIES_CONCEPTS_DISTANCE_CONCEPT_HPP
#define BOOST_GEOMETRY_STRATEGIES_CONCEPTS_DISTANCE_CONCEPT_HPP

#include <concepts>
#include <type_traits>

#include <boost/geometry/strategies/distance.hpp>
#include <boost/geometry/strategies/tags.hpp>

namespace boost { namespace geometry { namespace concepts
{


/*!
    \brief Checks strategy for point-point or point-box or box-box distance
    \ingroup distance
*/
template <typename Strategy, typename Point1, typename Point2>
concept PointDistanceStrategy =
    requires
    {
        typename strategy::distance::services
            ::return_type<Strategy, Point1, Point2>::type;
        typename strategy::distance::services
            ::comparable_type<Strategy>::type;
        typename strategy::distance::services::tag<Strategy>::type;
    }
    && (std::same_as
            <typename strategy::distance::services::tag<Strategy>::type,
             strategy_tag_distance_point_point>
        || std::same_as
            <typename strategy::distance::services::tag<Strategy>::type,
             strategy_tag_distance_point_box>
        || std::same_as
            <typename strategy::distance::services::tag<Strategy>::type,
             strategy_tag_distance_box_box>)
    && requires(Strategy const& instance,
                Point1 const& point1,
                Point2 const& point2)
    {
        { instance.apply(point1, point2) }
            -> std::convertible_to<typename strategy::distance::services
                ::return_type<Strategy, Point1, Point2>::type>;
        { strategy::distance::services::get_comparable<Strategy>::apply(instance) }
            -> std::convertible_to<typename strategy::distance::services
                ::comparable_type<Strategy>::type>;
        { strategy::distance::services::result_from_distance
            <Strategy, Point1, Point2>::apply(instance, 1.0) }
            -> std::convertible_to<typename strategy::distance::services
                ::return_type<Strategy, Point1, Point2>::type>;
    };


/*!
    \brief Checks strategy for point-segment distance
    \ingroup strategy_concepts
*/
template <typename Strategy, typename Point, typename PointOfSegment>
concept PointSegmentDistanceStrategy =
    requires
    {
        typename strategy::distance::services::tag<Strategy>::type;
        typename strategy::distance::services
            ::return_type<Strategy, Point, PointOfSegment>::type;
        typename strategy::distance::services
            ::comparable_type<Strategy>::type;
    }
    && std::same_as
        <typename strategy::distance::services::tag<Strategy>::type,
         strategy_tag_distance_point_segment>
    && requires(Strategy const& instance,
                Point const& point,
                PointOfSegment const& segment_point)
    {
        { instance.apply(point, segment_point, segment_point) }
            -> std::convertible_to<typename strategy::distance::services
                ::return_type<Strategy, Point, PointOfSegment>::type>;
        { strategy::distance::services::get_comparable<Strategy>::apply(instance) }
            -> std::convertible_to<typename strategy::distance::services
                ::comparable_type<Strategy>::type>;
        { strategy::distance::services::result_from_distance
            <Strategy, Point, PointOfSegment>::apply(instance, 1.0) }
            -> std::convertible_to<typename strategy::distance::services
                ::return_type<Strategy, Point, PointOfSegment>::type>;
    };


}}} // namespace boost::geometry::concepts


#endif // BOOST_GEOMETRY_STRATEGIES_CONCEPTS_DISTANCE_CONCEPT_HPP
