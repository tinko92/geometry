// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2014, Oracle and/or its affiliates.

// Contributed and/or modified by Menelaos Karavelas, on behalf of Oracle

// Licensed under the Boost Software License version 1.0.
// http://www.boost.org/users/license.html

#ifndef BOOST_GEOMETRY_ALGORITHMS_IS_SIMPLE_HPP
#define BOOST_GEOMETRY_ALGORITHMS_IS_SIMPLE_HPP

#include <boost/geometry/algorithms/detail/is_simple/implementation.hpp>

#include <boost/geometry/core/visit.hpp>
#include <boost/geometry/geometries/adapted/std_variant.hpp>
#include <boost/geometry/geometries/concepts/check.hpp>
#include <boost/geometry/strategies/default_strategy.hpp>
#include <boost/geometry/strategies/detail.hpp>
#include <boost/geometry/strategies/relate/services.hpp>

namespace boost { namespace geometry
{

namespace resolve_strategy
{

template <concepts::ConstGeometry Geometry, typename Strategy>
inline bool is_simple(Geometry const& geometry, Strategy const& strategy)
{
    if constexpr (std::same_as<Strategy, default_strategy>)
    {
        using strategy_type = typename strategies::relate::services
            ::default_strategy<Geometry, Geometry>::type;
        return detail::is_simple::apply(geometry, strategy_type());
    }
    else if constexpr (strategies::detail::is_umbrella_strategy<Strategy>::value)
    {
        return detail::is_simple::apply(geometry, strategy);
    }
    else
    {
        using strategies::relate::services::strategy_converter;
        return detail::is_simple::apply(
            geometry, strategy_converter<Strategy>::get(strategy));
    }
}

} // namespace resolve_strategy

namespace resolve_dynamic
{

template <concepts::ConstGeometry Geometry, typename Strategy>
inline bool is_simple(Geometry const& geometry, Strategy const& strategy)
{
    if constexpr (concepts::ConstDynamicGeometry<Geometry>)
    {
        bool result = false;
        traits::visit<Geometry>::apply([&](auto const& concrete)
        {
            result = resolve_strategy::is_simple(concrete, strategy);
        }, geometry);
        return result;
    }
    else
    {
        return resolve_strategy::is_simple(geometry, strategy);
    }
}

} // namespace resolve_dynamic

/*!
\brief \brief_check{is simple}
\ingroup is_simple
\tparam Geometry \tparam_geometry
\tparam Strategy \tparam_strategy{Is_simple}
\param geometry \param_geometry
\param strategy \param_strategy{is_simple}
\return \return_check{is simple}

\qbk{distinguish,with strategy}
\qbk{[include reference/algorithms/is_simple.qbk]}
*/
template <concepts::ConstGeometry Geometry, typename Strategy>
inline bool is_simple(Geometry const& geometry, Strategy const& strategy)
{
    return resolve_dynamic::is_simple(geometry, strategy);
}

/*!
\brief \brief_check{is simple}
\ingroup is_simple
\tparam Geometry \tparam_geometry
\param geometry \param_geometry
\return \return_check{is simple}

\qbk{[include reference/algorithms/is_simple.qbk]}
*/
template <concepts::ConstGeometry Geometry>
inline bool is_simple(Geometry const& geometry)
{
    return resolve_dynamic::is_simple(geometry, default_strategy());
}

}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_ALGORITHMS_IS_SIMPLE_HPP
