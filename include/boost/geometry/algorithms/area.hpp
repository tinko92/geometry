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

#include <type_traits>

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

// #include <boost/geometry/algorithms/detail/throw_on_empty_input.hpp>
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
namespace detail
{

template <typename Geometry, typename Strategy>
using area_result_t = typename area_result<Geometry, Strategy>::type;

template <concepts::ConstGeometry Geometry>
inline auto resolve_area_strategy(Geometry const&, default_strategy)
{
    using strategy_type = typename strategies::area::services::default_strategy
        <Geometry>::type;
    return strategy_type();
}

template <concepts::ConstGeometry Geometry, typename Strategy>
    requires strategies::detail::is_umbrella_strategy<Strategy>::value
inline Strategy const& resolve_area_strategy(Geometry const&, Strategy const& strategy)
{
    return strategy;
}

template <concepts::ConstGeometry Geometry, typename Strategy>
    requires (! std::same_as<Strategy, default_strategy>)
          && (! strategies::detail::is_umbrella_strategy<Strategy>::value)
inline auto resolve_area_strategy(Geometry const&, Strategy const& strategy)
{
    using strategies::area::services::strategy_converter;
    return strategy_converter<Strategy>::get(strategy);
}

template <concepts::ConstBox Box, typename Strategies>
inline area_result_t<Box, Strategies>
area_impl(Box const& box, Strategies const& strategies)
{
    assert_dimension<Box, 2>();
    return strategies.area(box).apply(box);
}

template <concepts::ConstRing Ring, typename Strategies>
inline area_result_t<Ring, Strategies>
area_impl(Ring const& ring, Strategies const& strategies)
{
    using strategy_type = decltype(strategies.area(ring));

    static_assert(concepts::AreaStrategy<Ring, strategy_type>);
    assert_dimension<Ring, 2>();

    boost::ignore_unused(strategies);

    if (boost::size(ring) < detail::minimum_ring_size<Ring>::value)
    {
        return area_result_t<Ring, Strategies>();
    }

    detail::closed_clockwise_view<Ring const> const view(ring);
    auto it = boost::begin(view);
    auto const end = boost::end(view);

    strategy_type const ring_strategy = strategies.area(ring);
    typename strategy_type::template state<Ring> state;

    for (auto previous = it++; it != end; ++previous, ++it)
    {
        ring_strategy.apply(*previous, *it, state);
    }

    return ring_strategy.result(state);
}

template <concepts::ConstPolygon Polygon, typename Strategies>
inline area_result_t<Polygon, Strategies>
area_impl(Polygon const& polygon, Strategies const& strategies)
{
    area_result_t<Polygon, Strategies> result
        = area_impl(exterior_ring(polygon), strategies);
    auto const& rings = interior_rings(polygon);
    for (auto it = boost::begin(rings); it != boost::end(rings); ++it)
    {
        result += area_impl(*it, strategies);
    }
    return result;
}

template <concepts::ConstMultiPolygon MultiPolygon, typename Strategies>
inline area_result_t<MultiPolygon, Strategies>
area_impl(MultiPolygon const& multi, Strategies const& strategies)
{
    area_result_t<MultiPolygon, Strategies> result = 0;
    for (auto it = boost::begin(multi); it != boost::end(multi); ++it)
    {
        result += area_impl(*it, strategies);
    }
    return result;
}

template <concepts::ConstGeometry Geometry, typename Strategies>
    requires (! concepts::ArealGeometry<Geometry>)
          && (! concepts::ConstDynamicGeometry<Geometry>)
          && (! concepts::ConstGeometryCollection<Geometry>)
inline auto area_impl(Geometry const& geometry, Strategies const& strategies)
{
    boost::ignore_unused(geometry, strategies);
    return area_result_t<Geometry, Strategies>();
}

template <concepts::ConstDynamicGeometry DynamicGeometry, typename Strategy>
inline area_result_t<DynamicGeometry, Strategy>
area_resolved(DynamicGeometry const& dynamic, Strategy const& strategy);

template <concepts::ConstGeometryCollection GeometryCollection, typename Strategy>
inline area_result_t<GeometryCollection, Strategy>
area_resolved(GeometryCollection const& collection, Strategy const& strategy);

template <concepts::ConstGeometry Geometry, typename Strategy>
    requires (! concepts::ConstDynamicGeometry<Geometry>)
          && (! concepts::ConstGeometryCollection<Geometry>)
inline auto area_resolved(Geometry const& geometry, Strategy const& strategy)
{
    auto&& strategies = resolve_area_strategy(geometry, strategy);
    return area_impl(geometry, strategies);
}

template <concepts::ConstDynamicGeometry DynamicGeometry, typename Strategy>
inline area_result_t<DynamicGeometry, Strategy>
area_resolved(DynamicGeometry const& dynamic, Strategy const& strategy)
{
    area_result_t<DynamicGeometry, Strategy> result = 0;
    traits::visit<DynamicGeometry>::apply([&](auto const& geometry)
    {
        result = area_resolved(geometry, strategy);
    }, dynamic);
    return result;
}

template <concepts::ConstGeometryCollection GeometryCollection, typename Strategy>
inline area_result_t<GeometryCollection, Strategy>
area_resolved(GeometryCollection const& collection, Strategy const& strategy)
{
    area_result_t<GeometryCollection, Strategy> result = 0;
    detail::visit_breadth_first([&](auto const& geometry)
    {
        result += area_resolved(geometry, strategy);
        return true;
    }, collection);
    return result;
}

} // namespace detail
#endif // DOXYGEN_NO_DETAIL


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

    return detail::area_resolved(geometry, default_strategy());
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

    return detail::area_resolved(geometry, strategy);
}


}} // namespace boost::geometry


#endif // BOOST_GEOMETRY_ALGORITHMS_AREA_HPP
