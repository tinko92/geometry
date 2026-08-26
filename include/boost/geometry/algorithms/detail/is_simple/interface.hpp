// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2014-2020, Oracle and/or its affiliates.

// Contributed and/or modified by Menelaos Karavelas, on behalf of Oracle
// Contributed and/or modified by Adam Wulkiewicz, on behalf of Oracle

// Licensed under the Boost Software License version 1.0.
// http://www.boost.org/users/license.html

#ifndef BOOST_GEOMETRY_ALGORITHMS_DETAIL_IS_SIMPLE_INTERFACE_HPP
#define BOOST_GEOMETRY_ALGORITHMS_DETAIL_IS_SIMPLE_INTERFACE_HPP

#include <type_traits>
#include <variant>

#include <boost/geometry/geometries/concepts/check.hpp>

#include <boost/geometry/algorithms/dispatch/is_simple.hpp>
#include <boost/geometry/core/cs.hpp>
#include <boost/geometry/strategies/default_strategy.hpp>
#include <boost/geometry/strategies/detail.hpp>
#include <boost/geometry/strategies/relate/services.hpp>


namespace boost { namespace geometry
{

namespace resolve_strategy
{

template
<
    typename Strategy,
    bool IsUmbrella = strategies::detail::is_umbrella_strategy<Strategy>::value
>
struct is_simple
{
    template <typename Geometry>
    static inline bool apply(Geometry const& geometry,
                             Strategy const& strategy)
    {
        return dispatch::is_simple<Geometry>::apply(geometry, strategy);
    }
};

template <typename Strategy>
struct is_simple<Strategy, false>
{
    template <typename Geometry>
    static inline bool apply(Geometry const& geometry,
                             Strategy const& strategy)
    {
        using strategies::relate::services::strategy_converter;
        return dispatch::is_simple
            <
                Geometry
            >::apply(geometry, strategy_converter<Strategy>::get(strategy));
    }
};

template <>
struct is_simple<default_strategy, false>
{
    template <typename Geometry>
    static inline bool apply(Geometry const& geometry,
                             default_strategy)
    {
        // NOTE: Currently the strategy is only used for Linear geometries
        typedef typename strategies::relate::services::default_strategy
            <
                Geometry, Geometry
            >::type strategy_type;

        return dispatch::is_simple<Geometry>::apply(geometry, strategy_type());
    }
};

} // namespace resolve_strategy

namespace resolve_variant
{

template <typename Geometry>
struct is_simple
{
    template <typename Strategy>
    static inline bool apply(Geometry const& geometry, Strategy const& strategy)
    {
        concepts::check<Geometry const>();

        return resolve_strategy::is_simple<Strategy>::apply(geometry, strategy);
    }
};

template <typename ...Ts>
struct is_simple<std::variant<Ts...>>
{
    template <typename Strategy>
    static inline bool
    apply(std::variant<Ts...> const& geometry,
          Strategy const& strategy)
    {
        return std::visit([&strategy](auto const& concrete)
        {
            return is_simple<std::decay_t<decltype(concrete)>>::apply(
                concrete, strategy);
        }, geometry);
    }
};

} // namespace resolve_variant


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
template <typename Geometry, typename Strategy>
inline bool is_simple(Geometry const& geometry, Strategy const& strategy)
{
    return resolve_variant::is_simple<Geometry>::apply(geometry, strategy);
}


/*!
\brief \brief_check{is simple}
\ingroup is_simple
\tparam Geometry \tparam_geometry
\param geometry \param_geometry
\return \return_check{is simple}

\qbk{[include reference/algorithms/is_simple.qbk]}
*/
template <typename Geometry>
inline bool is_simple(Geometry const& geometry)
{
    return resolve_variant::is_simple<Geometry>::apply(geometry, default_strategy());
}


}} // namespace boost::geometry


#endif // BOOST_GEOMETRY_ALGORITHMS_DETAIL_IS_SIMPLE_INTERFACE_HPP
