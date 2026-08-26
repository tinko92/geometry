// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2007-2014 Barend Gehrels, Amsterdam, the Netherlands.
// Copyright (c) 2008-2014 Bruno Lalande, Paris, France.
// Copyright (c) 2009-2014 Mateusz Loskot, London, UK.

// This file was modified by Oracle on 2014-2023.
// Modifications copyright (c) 2014-2023, Oracle and/or its affiliates.
// Contributed and/or modified by Vissarion Fysikopoulos, on behalf of Oracle
// Contributed and/or modified by Menelaos Karavelas, on behalf of Oracle
// Contributed and/or modified by Adam Wulkiewicz, on behalf of Oracle

// Parts of Boost.Geometry are redesigned from Geodan's Geographic Library
// (geolib/GGL), copyright (c) 1995-2010 Geodan, Amsterdam, the Netherlands.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_ALGORITHMS_LENGTH_HPP
#define BOOST_GEOMETRY_ALGORITHMS_LENGTH_HPP

#include <concepts>

#include <boost/range/begin.hpp>
#include <boost/range/end.hpp>
#include <boost/range/value_type.hpp>

#include "boost/geometry/algorithms/detail/assign_indexed_point.hpp"
#include <boost/geometry/algorithms/detail/dummy_geometries.hpp>
// #include <boost/geometry/algorithms/detail/throw_on_empty_input.hpp>
#include <boost/geometry/algorithms/detail/visit.hpp>

#include <boost/geometry/core/closure.hpp>
#include <boost/geometry/core/tag.hpp>
#include <boost/geometry/core/tags.hpp>
#include <boost/geometry/core/visit.hpp>

#include <boost/geometry/geometries/adapted/std_variant.hpp>
#include <boost/geometry/geometries/concepts/check.hpp>

#include <boost/geometry/strategies/default_strategy.hpp>
#include <boost/geometry/strategies/default_length_result.hpp> // TODO: Move to algorithms
#include <boost/geometry/strategies/detail.hpp>
#include <boost/geometry/strategies/length/cartesian.hpp>
#include <boost/geometry/strategies/length/geographic.hpp>
#include <boost/geometry/strategies/length/spherical.hpp>

#include <boost/geometry/views/closeable_view.hpp>

namespace boost { namespace geometry
{


#ifndef DOXYGEN_NO_DETAIL
namespace detail
{

template <closure_selector Closure, typename Range, typename Strategies>
inline typename default_length_result<Range>::type
range_length(Range const& range, Strategies const& strategies)
{
    typename default_length_result<Range>::type sum = 0;
    detail::closed_view<Range const> const view(range);
    auto it = boost::begin(view);
    auto const end = boost::end(view);
    if (it != end)
    {
        auto const strategy = strategies.distance(dummy_point(), dummy_point());
        for (auto previous = it++; it != end; ++previous, ++it)
        {
            sum += strategy.apply(*previous, *it);
        }
    }
    return sum;
}

template <concepts::ConstGeometry Geometry>
inline auto resolve_length_strategy(Geometry const&, default_strategy)
{
    using strategies_type = typename strategies::length::services::default_strategy
        <Geometry>::type;
    return strategies_type();
}

template <concepts::ConstGeometry Geometry, typename Strategy>
    requires strategies::detail::is_umbrella_strategy<Strategy>::value
inline Strategy const& resolve_length_strategy(
    Geometry const&, Strategy const& strategy)
{
    return strategy;
}

template <concepts::ConstGeometry Geometry, typename Strategy>
    requires (! std::same_as<Strategy, default_strategy>)
          && (! strategies::detail::is_umbrella_strategy<Strategy>::value)
inline auto resolve_length_strategy(Geometry const&, Strategy const& strategy)
{
    using strategies::length::services::strategy_converter;
    return strategy_converter<Strategy>::get(strategy);
}

template <concepts::ConstDynamicGeometry DynamicGeometry, typename Strategy>
inline typename default_length_result<DynamicGeometry>::type
length_impl(DynamicGeometry const& dynamic, Strategy const& strategy);

template <concepts::ConstGeometryCollection GeometryCollection, typename Strategy>
inline typename default_length_result<GeometryCollection>::type
length_impl(GeometryCollection const& collection, Strategy const& strategy);

template <concepts::ConstLinestring Linestring, typename Strategy>
inline typename default_length_result<Linestring>::type
length_impl(Linestring const& linestring, Strategy const& strategy)
{
    auto&& strategies = resolve_length_strategy(linestring, strategy);
    return range_length<closed>(linestring, strategies);
}

template <concepts::ConstSegment Segment, typename Strategy>
inline typename default_length_result<Segment>::type
length_impl(Segment const& segment, Strategy const& strategy)
{
    auto&& strategies = resolve_length_strategy(segment, strategy);
    point_type_t<Segment> p1, p2;
    geometry::detail::assign_point_from_index<0>(segment, p1);
    geometry::detail::assign_point_from_index<1>(segment, p2);
    return strategies.distance(p1, p2).apply(p1, p2);
}

template <concepts::ConstMultiLinestring MultiLinestring, typename Strategy>
inline typename default_length_result<MultiLinestring>::type
length_impl(MultiLinestring const& multi, Strategy const& strategy)
{
    auto&& strategies = resolve_length_strategy(multi, strategy);
    typename default_length_result<MultiLinestring>::type result = 0;
    for (auto it = boost::begin(multi); it != boost::end(multi); ++it)
    {
        result += length_impl(*it, strategies);
    }
    return result;
}

template <concepts::ConstGeometry Geometry, typename Strategy>
    requires (! concepts::GeometryCategory<Geometry, linear_tag>)
          && (! concepts::ConstDynamicGeometry<Geometry>)
          && (! concepts::ConstGeometryCollection<Geometry>)
inline typename default_length_result<Geometry>::type
length_impl(Geometry const&, Strategy const&)
{
    return 0;
}

template <concepts::ConstDynamicGeometry DynamicGeometry, typename Strategy>
inline typename default_length_result<DynamicGeometry>::type
length_impl(DynamicGeometry const& dynamic, Strategy const& strategy)
{
    typename default_length_result<DynamicGeometry>::type result = 0;
    traits::visit<DynamicGeometry>::apply([&](auto const& geometry)
    {
        result = length_impl(geometry, strategy);
    }, dynamic);
    return result;
}

template <concepts::ConstGeometryCollection GeometryCollection, typename Strategy>
inline typename default_length_result<GeometryCollection>::type
length_impl(GeometryCollection const& collection, Strategy const& strategy)
{
    typename default_length_result<GeometryCollection>::type result = 0;
    detail::visit_breadth_first([&](auto const& geometry)
    {
        result += length_impl(geometry, strategy);
        return true;
    }, collection);
    return result;
}

} // namespace detail
#endif // DOXYGEN_NO_DETAIL


/*!
\brief \brief_calc{length}
\ingroup length
\details \details_calc{length, length (the sum of distances between consecutive points)}. \details_default_strategy
\tparam Geometry \tparam_geometry
\param geometry \param_geometry
\return \return_calc{length}

\qbk{[include reference/algorithms/length.qbk]}
\qbk{[length] [length_output]}
 */
template<concepts::ConstGeometry Geometry>
inline typename default_length_result<Geometry>::type
length(Geometry const& geometry)
{
    // detail::throw_on_empty_input(geometry);

    return detail::length_impl(geometry, default_strategy());
}


/*!
\brief \brief_calc{length} \brief_strategy
\ingroup length
\details \details_calc{length, length (the sum of distances between consecutive points)} \brief_strategy. \details_strategy_reasons
\tparam Geometry \tparam_geometry
\tparam Strategy \tparam_strategy{distance}
\param geometry \param_geometry
\param strategy \param_strategy{distance}
\return \return_calc{length}

\qbk{distinguish,with strategy}
\qbk{[include reference/algorithms/length.qbk]}
\qbk{[length_with_strategy] [length_with_strategy_output]}
 */
template<concepts::ConstGeometry Geometry, typename Strategy>
inline typename default_length_result<Geometry>::type
length(Geometry const& geometry, Strategy const& strategy)
{
    // detail::throw_on_empty_input(geometry);

    return detail::length_impl(geometry, strategy);
}


}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_ALGORITHMS_LENGTH_HPP
