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

#ifndef BOOST_GEOMETRY_ALGORITHMS_PERIMETER_HPP
#define BOOST_GEOMETRY_ALGORITHMS_PERIMETER_HPP

#include <boost/range/value_type.hpp>

#include <boost/geometry/algorithms/length.hpp>
#include <boost/geometry/algorithms/detail/calculate_null.hpp>
#include <boost/geometry/algorithms/detail/calculate_sum.hpp>
#include <boost/geometry/algorithms/detail/multi_sum.hpp>
// #include <boost/geometry/algorithms/detail/throw_on_empty_input.hpp>
#include <boost/geometry/algorithms/detail/visit.hpp>

#include <boost/geometry/core/closure.hpp>
#include <boost/geometry/core/tags.hpp>
#include <boost/geometry/core/visit.hpp>

#include <boost/geometry/geometries/adapted/std_variant.hpp>
#include <boost/geometry/geometries/concepts/check.hpp>

#include <boost/geometry/strategies/default_length_result.hpp>
#include <boost/geometry/strategies/default_strategy.hpp>
#include <boost/geometry/strategies/detail.hpp>
#include <boost/geometry/strategies/length/cartesian.hpp>
#include <boost/geometry/strategies/length/geographic.hpp>
#include <boost/geometry/strategies/length/spherical.hpp>

namespace boost { namespace geometry
{

#ifndef DOXYGEN_NO_DISPATCH
namespace dispatch
{

template <concepts::ConstGeometry Geometry, typename Strategy>
inline typename default_length_result<Geometry>::type
perimeter(Geometry const& geometry, Strategy const& strategy)
{
    if constexpr (concepts::ConstRing<Geometry>)
    {
        return detail::length::range_length
            <Geometry, closure<Geometry>::value>::apply(geometry, strategy);
    }
    else if constexpr (concepts::ConstPolygon<Geometry>)
    {
        using policy = detail::length::range_length
            <ring_type_t<Geometry>, closure<Geometry>::value>;
        return detail::calculate_polygon_sum::apply
            <typename default_length_result<Geometry>::type, policy>(
                geometry, strategy);
    }
    else if constexpr (concepts::ConstMultiPolygon<Geometry>)
    {
        typename default_length_result<Geometry>::type result = 0;
        for (auto it = boost::begin(geometry); it != boost::end(geometry); ++it)
        {
            result += dispatch::perimeter(*it, strategy);
        }
        return result;
    }
    else
    {
        return detail::calculate_null::apply
            <typename default_length_result<Geometry>::type>(geometry, strategy);
    }
}

} // namespace dispatch
#endif // DOXYGEN_NO_DISPATCH


namespace resolve_strategy {

template <concepts::ConstGeometry Geometry, typename Strategy>
inline typename default_length_result<Geometry>::type
perimeter(Geometry const& geometry, Strategy const& strategy)
{
    if constexpr (std::same_as<Strategy, default_strategy>)
    {
        using strategies_type = typename strategies::length::services::default_strategy
            <Geometry>::type;
        return dispatch::perimeter(geometry, strategies_type());
    }
    else if constexpr (strategies::detail::is_umbrella_strategy<Strategy>::value)
    {
        return dispatch::perimeter(geometry, strategy);
    }
    else
    {
        using strategies::length::services::strategy_converter;
        return dispatch::perimeter(
            geometry, strategy_converter<Strategy>::get(strategy));
    }
}

} // namespace resolve_strategy


namespace resolve_dynamic {

template <concepts::ConstGeometry Geometry, typename Strategy>
inline typename default_length_result<Geometry>::type
perimeter(Geometry const& geometry, Strategy const& strategy)
{
    if constexpr (concepts::ConstDynamicGeometry<Geometry>)
    {
        typename default_length_result<Geometry>::type result = 0;
        traits::visit<Geometry>::apply([&](auto const& g)
        {
            result = resolve_dynamic::perimeter(g, strategy);
        }, geometry);
        return result;
    }
    else if constexpr (concepts::ConstGeometryCollection<Geometry>)
    {
        typename default_length_result<Geometry>::type result = 0;
        detail::visit_breadth_first([&](auto const& g)
        {
            result += resolve_dynamic::perimeter(g, strategy);
            return true;
        }, geometry);
        return result;
    }
    else
    {
        return resolve_strategy::perimeter(geometry, strategy);
    }
}

} // namespace resolve_dynamic


/*!
\brief \brief_calc{perimeter}
\ingroup perimeter
\details The function perimeter returns the perimeter of a geometry,
    using the default distance-calculation-strategy
\tparam Geometry \tparam_geometry
\param geometry \param_geometry
\return \return_calc{perimeter}

\qbk{[include reference/algorithms/perimeter.qbk]}
\qbk{
[heading Example]
[perimeter]
[perimeter_output]
}
 */
template<concepts::ConstGeometry Geometry>
inline typename default_length_result<Geometry>::type perimeter(
        Geometry const& geometry)
{
    // detail::throw_on_empty_input(geometry);
    return resolve_dynamic::perimeter(geometry, default_strategy());
}

/*!
\brief \brief_calc{perimeter} \brief_strategy
\ingroup perimeter
\details The function perimeter returns the perimeter of a geometry,
    using specified strategy
\tparam Geometry \tparam_geometry
\tparam Strategy \tparam_strategy{distance}
\param geometry \param_geometry
\param strategy strategy to be used for distance calculations.
\return \return_calc{perimeter}

\qbk{distinguish,with strategy}
\qbk{[include reference/algorithms/perimeter.qbk]}
 */
template<concepts::ConstGeometry Geometry, typename Strategy>
inline typename default_length_result<Geometry>::type perimeter(
        Geometry const& geometry, Strategy const& strategy)
{
    // detail::throw_on_empty_input(geometry);
    return resolve_dynamic::perimeter(geometry, strategy);
}

}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_ALGORITHMS_PERIMETER_HPP
