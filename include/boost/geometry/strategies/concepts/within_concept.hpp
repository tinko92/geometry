// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2007-2012 Barend Gehrels, Amsterdam, the Netherlands.
// Copyright (c) 2008-2012 Bruno Lalande, Paris, France.
// Copyright (c) 2009-2012 Mateusz Loskot, London, UK.

// This file was modified by Oracle on 2018-2020.
// Modifications copyright (c) 2018-2020 Oracle and/or its affiliates.
// Contributed and/or modified by Adam Wulkiewicz, on behalf of Oracle

// Parts of Boost.Geometry are redesigned from Geodan's Geographic Library
// (geolib/GGL), copyright (c) 1995-2010 Geodan, Amsterdam, the Netherlands.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_STRATEGIES_CONCEPTS_WITHIN_CONCEPT_HPP
#define BOOST_GEOMETRY_STRATEGIES_CONCEPTS_WITHIN_CONCEPT_HPP


#include <concepts>
#include <type_traits>
#include <utility>

#include <boost/geometry/core/point_type.hpp>

#include <boost/geometry/geometries/concepts/box_concept.hpp>
#include <boost/geometry/geometries/concepts/point_concept.hpp>

#include <boost/geometry/strategies/detail.hpp>

namespace boost { namespace geometry { namespace concepts
{


namespace detail
{


template <typename Point, typename Geometry, typename Strategy>
consteval auto relate_strategy_identity()
{
    if constexpr (strategies::detail::is_umbrella_strategy<Strategy>::value)
    {
        using type = decltype(std::declval<Strategy>().relate(
            std::declval<Point>(), std::declval<Geometry>()));
        return std::type_identity<type>{};
    }
    else
    {
        return std::type_identity<Strategy>{};
    }
}

template <typename Point, typename Geometry, typename Strategy>
using relate_strategy_t = typename decltype(
    relate_strategy_identity<Point, Geometry, Strategy>())::type;

template <typename Point, typename Geometry, typename Strategy>
consteval auto within_strategy_identity()
{
    if constexpr (strategies::detail::is_umbrella_strategy<Strategy>::value)
    {
        using type = decltype(std::declval<Strategy>().within(
            std::declval<Point>(), std::declval<Geometry>()));
        return std::type_identity<type>{};
    }
    else
    {
        return std::type_identity<Strategy>{};
    }
}

template <typename Point, typename Geometry, typename Strategy>
using within_strategy_t = typename decltype(
    within_strategy_identity<Point, Geometry, Strategy>())::type;


} // namespace detail


/*!
\brief Checks strategy for within (point-in-polygon)
\ingroup within
*/
template <typename Point, typename Polygonal, typename Strategy>
concept WithinStrategyPolygonal =
    concepts::ConstPoint<Point>
    && requires(detail::relate_strategy_t<Point, Polygonal, Strategy> const& strategy,
                Point const& point,
                point_type_t<Polygonal> const& segment_point,
                typename detail::relate_strategy_t
                    <Point, Polygonal, Strategy>::state_type& state)
    {
        typename detail::relate_strategy_t
            <Point, Polygonal, Strategy>::state_type;
        { strategy.apply(point, segment_point, segment_point, state) }
            -> std::same_as<bool>;
        { strategy.result(state) } -> std::same_as<int>;
    };

template <typename Point, typename Box, typename Strategy>
concept WithinStrategyPointBox =
    concepts::ConstPoint<Point>
    && concepts::ConstBox<Box>
    && requires(detail::within_strategy_t<Point, Box, Strategy> const& strategy,
                Point const& point,
                Box const& box)
    {
        { strategy.apply(point, box) } -> std::same_as<bool>;
    };

template <typename Box1, typename Box2, typename Strategy>
concept WithinStrategyBoxBox =
    concepts::ConstBox<Box1>
    && concepts::ConstBox<Box2>
    && requires(detail::within_strategy_t<Box1, Box2, Strategy> const& strategy,
                Box1 const& box1,
                Box2 const& box2)
    {
        { strategy.apply(box1, box2) } -> std::same_as<bool>;
    };

// So now: boost::geometry::concepts::within
namespace within
{

/*!
\brief Checks, in compile-time, the concept of any within-strategy
\ingroup concepts
*/
template <concepts::ConstGeometry Geometry1,
          concepts::ConstGeometry Geometry2,
          typename Strategy>
constexpr void check()
{
    if constexpr (concepts::ConstPoint<Geometry1>
                  && concepts::ConstBox<Geometry2>)
    {
        static_assert(WithinStrategyPointBox
            <Geometry1, Geometry2, Strategy>);
    }
    else if constexpr (concepts::ConstBox<Geometry1>
                       && concepts::ConstBox<Geometry2>)
    {
        static_assert(WithinStrategyBoxBox
            <Geometry1, Geometry2, Strategy>);
    }
    else if constexpr (concepts::ConstPoint<Geometry1>)
    {
        static_assert(WithinStrategyPolygonal
            <Geometry1, Geometry2, Strategy>);
    }
}


}}}} // namespace boost::geometry::concepts::within


#endif // BOOST_GEOMETRY_STRATEGIES_CONCEPTS_WITHIN_CONCEPT_HPP
