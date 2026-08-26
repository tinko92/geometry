// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2007-2012 Barend Gehrels, Amsterdam, the Netherlands.
// Copyright (c) 2008-2012 Bruno Lalande, Paris, France.
// Copyright (c) 2009-2012 Mateusz Loskot, London, UK.
// Copyright (c) 2017-2024 Adam Wulkiewicz, Lodz, Poland.

// This file was modified by Oracle on 2017-2023.
// Modifications copyright (c) 2017-2023 Oracle and/or its affiliates.
// Contributed and/or modified by Vissarion Fysikopoulos, on behalf of Oracle
// Contributed and/or modified by Adam Wulkiewicz, on behalf of Oracle

// Parts of Boost.Geometry are redesigned from Geodan's Geographic Library
// (geolib/GGL), copyright (c) 1995-2010 Geodan, Amsterdam, the Netherlands.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_ALGORITHMS_AREA_HPP
#define BOOST_GEOMETRY_ALGORITHMS_AREA_HPP

#include <boost/core/ignore_unused.hpp>
#include <boost/range/begin.hpp>
#include <boost/range/end.hpp>
#include <boost/range/size.hpp>
#include <boost/range/value_type.hpp>

#include <boost/geometry/core/closure.hpp>
#include <boost/geometry/core/exterior_ring.hpp>
#include <boost/geometry/core/interior_rings.hpp>
#include <boost/geometry/core/point_order.hpp>
#include <boost/geometry/core/point_type.hpp>
#include <boost/geometry/core/ring_type.hpp>
#include <boost/geometry/core/tags.hpp>
#include <boost/geometry/core/visit.hpp>

#include <boost/geometry/algorithms/detail/calculate_null.hpp>
#include <boost/geometry/algorithms/detail/calculate_sum.hpp>
// #include <boost/geometry/algorithms/detail/throw_on_empty_input.hpp>
#include <boost/geometry/algorithms/detail/multi_sum.hpp>
#include <boost/geometry/algorithms/detail/visit.hpp>

#include <boost/geometry/algorithms/area_result.hpp>
#include <boost/geometry/algorithms/default_area_result.hpp>

#include <boost/geometry/geometries/adapted/std_variant.hpp>
#include <boost/geometry/geometries/concepts/check.hpp>

#include <boost/geometry/strategies/area/services.hpp>
#include <boost/geometry/strategies/area/cartesian.hpp>
#include <boost/geometry/strategies/area/geographic.hpp>
#include <boost/geometry/strategies/area/spherical.hpp>
#include <boost/geometry/strategies/concepts/area_concept.hpp>
#include <boost/geometry/strategies/default_strategy.hpp>

#include <boost/geometry/views/detail/closed_clockwise_view.hpp>


namespace boost { namespace geometry
{


#ifndef DOXYGEN_NO_DETAIL
namespace detail { namespace area
{

struct box_area
{
    template <typename Box, typename Strategies>
    static inline auto
    apply(Box const& box, Strategies const& strategies)
    {
        // Currently only works for 2D Cartesian boxes
        assert_dimension<Box, 2>();

        return strategies.area(box).apply(box);
    }
};


struct ring_area
{
    template <typename Ring, typename Strategies>
    static inline typename area_result<Ring, Strategies>::type
    apply(Ring const& ring, Strategies const& strategies)
    {
        using strategy_type = decltype(strategies.area(ring));

        BOOST_CONCEPT_ASSERT( (geometry::concepts::AreaStrategy<Ring, strategy_type>) );
        assert_dimension<Ring, 2>();

        // Ignore warning (because using static method sometimes) on strategy
        boost::ignore_unused(strategies);

        // An open ring has at least three points,
        // A closed ring has at least four points,
        // if not, there is no (zero) area
        if (boost::size(ring) < detail::minimum_ring_size<Ring>::value)
        {
            return typename area_result<Ring, Strategies>::type();
        }

        detail::closed_clockwise_view<Ring const> const view(ring);
        auto it = boost::begin(view);
        auto const end = boost::end(view);

        strategy_type const strategy = strategies.area(ring);
        typename strategy_type::template state<Ring> state;

        for (auto previous = it++; it != end; ++previous, ++it)
        {
            strategy.apply(*previous, *it, state);
        }

        return strategy.result(state);
    }
};


}} // namespace detail::area


#endif // DOXYGEN_NO_DETAIL


#ifndef DOXYGEN_NO_DISPATCH
namespace dispatch
{

template <concepts::ConstGeometry Geometry, typename Strategy>
inline auto area(Geometry const& geometry, Strategy const& strategy)
{
    if constexpr (concepts::ConstBox<Geometry>)
    {
        return detail::area::box_area::apply(geometry, strategy);
    }
    else if constexpr (concepts::ConstRing<Geometry>)
    {
        return detail::area::ring_area::apply(geometry, strategy);
    }
    else if constexpr (concepts::ConstPolygon<Geometry>)
    {
        return detail::calculate_polygon_sum::apply
            <
                typename area_result<Geometry, Strategy>::type,
                detail::area::ring_area
            >(geometry, strategy);
    }
    else if constexpr (concepts::ConstMultiPolygon<Geometry>)
    {
        typename area_result<Geometry, Strategy>::type result = 0;
        for (auto it = boost::begin(geometry); it != boost::end(geometry); ++it)
        {
            result += dispatch::area(*it, strategy);
        }
        return result;
    }
    else
    {
        return detail::calculate_null::apply
            <typename area_result<Geometry, Strategy>::type>(geometry, strategy);
    }
}


} // namespace dispatch
#endif // DOXYGEN_NO_DISPATCH


namespace resolve_strategy
{

template <concepts::ConstGeometry Geometry, typename Strategy>
inline auto area(Geometry const& geometry, Strategy const& strategy)
{
    if constexpr (std::same_as<Strategy, default_strategy>)
    {
        using strategy_type = typename strategies::area::services::default_strategy
            <Geometry>::type;
        return dispatch::area(geometry, strategy_type());
    }
    else if constexpr (strategies::detail::is_umbrella_strategy<Strategy>::value)
    {
        return dispatch::area(geometry, strategy);
    }
    else
    {
        using strategies::area::services::strategy_converter;
        return dispatch::area(geometry,
            strategy_converter<Strategy>::get(strategy));
    }
}


} // namespace resolve_strategy


namespace resolve_dynamic
{

template <concepts::ConstGeometry Geometry, typename Strategy>
inline auto area(Geometry const& geometry, Strategy const& strategy)
{
    if constexpr (concepts::ConstDynamicGeometry<Geometry>)
    {
        typename area_result<Geometry, Strategy>::type result = 0;
        traits::visit<Geometry>::apply([&](auto const& g)
        {
            result = resolve_dynamic::area(g, strategy);
        }, geometry);
        return result;
    }
    else if constexpr (concepts::ConstGeometryCollection<Geometry>)
    {
        typename area_result<Geometry, Strategy>::type result = 0;
        detail::visit_breadth_first([&](auto const& g)
        {
            result += resolve_dynamic::area(g, strategy);
            return true;
        }, geometry);
        return result;
    }
    else
    {
        return resolve_strategy::area(geometry, strategy);
    }
}

} // namespace resolve_dynamic


/*!
\brief \brief_calc{area}
\ingroup area
\details \details_calc{area}. \details_default_strategy

The area algorithm calculates the surface area of all geometries having a surface, namely
box, polygon, ring, multipolygon. The units are the square of the units used for the points
defining the surface. If subject geometry is defined in meters, then area is calculated
in square meters.

The area calculation can be done in all three common coordinate systems, Cartesian, Spherical
and Geographic as well.

\tparam Geometry \tparam_geometry
\param geometry \param_geometry
\return \return_calc{area}

\qbk{[include reference/algorithms/area.qbk]}
\qbk{[heading Examples]}
\qbk{[area] [area_output]}
*/
template <concepts::ConstGeometry Geometry>
inline auto area(Geometry const& geometry)
{
    // detail::throw_on_empty_input(geometry);

    return resolve_dynamic::area(geometry, default_strategy());
}

/*!
\brief \brief_calc{area} \brief_strategy
\ingroup area
\details \details_calc{area} \brief_strategy. \details_strategy_reasons
\tparam Geometry \tparam_geometry
\tparam Strategy \tparam_strategy{Area}
\param geometry \param_geometry
\param strategy \param_strategy{area}
\return \return_calc{area}

\qbk{distinguish,with strategy}

\qbk{
[include reference/algorithms/area.qbk]

[heading Available Strategies]
\* [link geometry.reference.strategies.strategy_area_cartesian Cartesian]
\* [link geometry.reference.strategies.strategy_area_spherical Spherical]
\* [link geometry.reference.strategies.strategy_area_geographic Geographic]

[heading Example]
[area_with_strategy]
[area_with_strategy_output]
}
 */
template <concepts::ConstGeometry Geometry, typename Strategy>
inline auto area(Geometry const& geometry, Strategy const& strategy)
{
    // detail::throw_on_empty_input(geometry);

    return resolve_dynamic::area(geometry, strategy);
}


}} // namespace boost::geometry


#endif // BOOST_GEOMETRY_ALGORITHMS_AREA_HPP
