// Boost.Geometry

// Copyright (c) 2021, Oracle and/or its affiliates.

// Contributed and/or modified by Adam Wulkiewicz, on behalf of Oracle

// Licensed under the Boost Software License version 1.0.
// http://www.boost.org/users/license.html

#ifndef BOOST_GEOMETRY_STRATEGIES_DISTANCE_CARTESIAN_HPP
#define BOOST_GEOMETRY_STRATEGIES_DISTANCE_CARTESIAN_HPP


//#include <boost/geometry/strategies/cartesian/azimuth.hpp>

#include <boost/geometry/strategies/cartesian/distance_projected_point.hpp>
#include <boost/geometry/strategies/cartesian/distance_pythagoras.hpp>
#include <boost/geometry/strategies/cartesian/distance_pythagoras_box_box.hpp>
#include <boost/geometry/strategies/cartesian/distance_pythagoras_point_box.hpp>
#include <boost/geometry/strategies/cartesian/distance_segment_box.hpp>

#include <boost/geometry/strategies/detail.hpp>
#include <boost/geometry/strategies/distance/comparable.hpp>
#include <boost/geometry/strategies/distance/detail.hpp>
#include <boost/geometry/strategies/distance/services.hpp>

//#include <boost/geometry/strategies/normalize.hpp>
#include <boost/geometry/strategies/relate/cartesian.hpp>

#include <boost/geometry/util/type_traits.hpp>


namespace boost { namespace geometry
{

namespace strategies { namespace distance
{

template <typename CalculationType = void>
struct cartesian
    : public strategies::relate::cartesian<CalculationType>
{
    template <typename Geometry1, typename Geometry2>
        requires detail::geometry_pair<Geometry1, Geometry2>
    static auto distance(Geometry1 const&, Geometry2 const&)
    {
        if constexpr (detail::is_pp_v<Geometry1, Geometry2>)
            return strategy::distance::pythagoras<CalculationType>();
        else if constexpr (detail::is_ps_v<Geometry1, Geometry2>)
            return strategy::distance::projected_point
                <CalculationType, strategy::distance::pythagoras<CalculationType>>();
        else if constexpr (detail::is_pb_v<Geometry1, Geometry2>)
            return strategy::distance::pythagoras_point_box<CalculationType>();
        else if constexpr (detail::is_sb_v<Geometry1, Geometry2>)
            return strategy::distance::cartesian_segment_box
                <CalculationType, strategy::distance::pythagoras<CalculationType>>();
        else
            return strategy::distance::pythagoras_box_box<CalculationType>();
    }
};


namespace services
{

template <typename Geometry1, typename Geometry2>
struct default_strategy<Geometry1, Geometry2, cartesian_tag, cartesian_tag>
{
    using type = strategies::distance::cartesian<>;
};


template <typename CT>
struct strategy_converter<strategy::distance::pythagoras<CT> >
{
    template <typename S>
    static auto get(S const&)
    {
        return strategies::distance::cartesian<CT>();
    }
};

template <typename CT, typename PPS>
struct strategy_converter<strategy::distance::projected_point<CT, PPS> >
    : strategy_converter<PPS>
{};

template <typename CT>
struct strategy_converter<strategy::distance::pythagoras_point_box<CT> >
{
    static auto get(strategy::distance::pythagoras_point_box<CT> const&)
    {
        return strategies::distance::cartesian<CT>();
    }
};

template <typename CT, typename PPS>
struct strategy_converter<strategy::distance::cartesian_segment_box<CT, PPS> >
    : strategy_converter<PPS>
{};

template <typename CT>
struct strategy_converter<strategy::distance::pythagoras_box_box<CT> >
{
    static auto get(strategy::distance::pythagoras_box_box<CT> const&)
    {
        return strategies::distance::cartesian<CT>();
    }
};


template <typename CT>
struct strategy_converter<strategy::distance::comparable::pythagoras<CT> >
{
    template <typename S>
    static auto get(S const&)
    {
        return strategies::distance::detail::make_comparable(
                    strategies::distance::cartesian<CT>());
    }
};

template <typename CT>
struct strategy_converter<strategy::distance::comparable::pythagoras_point_box<CT> >
{
    static auto get(strategy::distance::comparable::pythagoras_point_box<CT> const&)
    {
        return strategies::distance::detail::make_comparable(
                    strategies::distance::cartesian<CT>());
    }
};

template <typename CT>
struct strategy_converter<strategy::distance::comparable::pythagoras_box_box<CT> >
{
    static auto get(strategy::distance::comparable::pythagoras_box_box<CT> const&)
    {
        return strategies::distance::detail::make_comparable(
                    strategies::distance::cartesian<CT>());
    }
};


} // namespace services

}} // namespace strategies::distance

}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_STRATEGIES_DISTANCE_CARTESIAN_HPP
