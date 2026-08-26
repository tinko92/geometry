// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2007-2012 Barend Gehrels, Amsterdam, the Netherlands.
// Copyright (c) 2008-2012 Bruno Lalande, Paris, France.
// Copyright (c) 2009-2012 Mateusz Loskot, London, UK.
// Copyright (c) 2014 Samuel Debionne, Grenoble, France.

// This file was modified by Oracle on 2014-2022.
// Modifications copyright (c) 2014-2022 Oracle and/or its affiliates.

// Contributed and/or modified by Adam Wulkiewicz, on behalf of Oracle

// Parts of Boost.Geometry are redesigned from Geodan's Geographic Library
// (geolib/GGL), copyright (c) 1995-2010 Geodan, Amsterdam, the Netherlands.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_ALGORITHMS_CROSSES_HPP
#define BOOST_GEOMETRY_ALGORITHMS_CROSSES_HPP

#include <cstddef>

#include <boost/geometry/algorithms/detail/gc_topological_dimension.hpp>
#include <boost/geometry/algorithms/detail/relate/relate_impl.hpp>
#include <boost/geometry/algorithms/relate.hpp>
#include <boost/geometry/core/access.hpp>
#include <boost/geometry/geometries/adapted/std_variant.hpp>
#include <boost/geometry/geometries/concepts/check.hpp>
#include <boost/geometry/strategies/default_strategy.hpp>
#include <boost/geometry/strategies/detail.hpp>
#include <boost/geometry/strategies/relate/cartesian.hpp>
#include <boost/geometry/strategies/relate/geographic.hpp>
#include <boost/geometry/strategies/relate/spherical.hpp>
#include <boost/geometry/views/detail/geometry_collection_view.hpp>


namespace boost { namespace geometry
{

#ifndef DOXYGEN_NO_DISPATCH
namespace dispatch
{


template <concepts::ConstGeometry Geometry1,
          concepts::ConstGeometry Geometry2,
          typename Strategy>
inline bool crosses(Geometry1 const& geometry1, Geometry2 const& geometry2,
                    Strategy const& strategy)
{
    if constexpr (concepts::ConstGeometryCollection<Geometry1>
                  && concepts::ConstGeometryCollection<Geometry2>)
    {
        int const dimension1 = detail::gc_topological_dimension(geometry1);
        int const dimension2 = detail::gc_topological_dimension(geometry2);

        if (dimension1 >= 0 && dimension2 >= 0)
        {
            if (dimension1 < dimension2)
            {
                return detail::relate::relate_impl
                    <
                        detail::de9im::static_mask_crosses_d1_le_d2_type,
                        Geometry1,
                        Geometry2
                    >::apply(geometry1, geometry2, strategy);
            }
            else if (dimension1 > dimension2)
            {
                return detail::relate::relate_impl
                    <
                        detail::de9im::static_mask_crosses_d2_le_d1_type,
                        Geometry1,
                        Geometry2
                    >::apply(geometry1, geometry2, strategy);
            }
            else if (dimension1 == 1 && dimension2 == 1)
            {
                return detail::relate::relate_impl
                    <
                        detail::de9im::static_mask_crosses_d1_1_d2_1_type,
                        Geometry1,
                        Geometry2
                    >::apply(geometry1, geometry2, strategy);
            }
        }

        return false;
    }
    else if constexpr (concepts::ConstGeometryCollection<Geometry2>)
    {
        using view_type = detail::geometry_collection_view<Geometry1>;
        return dispatch::crosses(
            view_type(geometry1), geometry2, strategy);
    }
    else if constexpr (concepts::ConstGeometryCollection<Geometry1>)
    {
        using view_type = detail::geometry_collection_view<Geometry2>;
        return dispatch::crosses(
            geometry1, view_type(geometry2), strategy);
    }
    else
    {
        return detail::relate::relate_impl
            <
                detail::de9im::static_mask_crosses_type,
                Geometry1,
                Geometry2
            >::apply(geometry1, geometry2, strategy);
    }
}


} // namespace dispatch
#endif // DOXYGEN_NO_DISPATCH


namespace resolve_strategy
{

template <concepts::ConstGeometry Geometry1,
          concepts::ConstGeometry Geometry2,
          typename Strategy>
inline bool crosses(Geometry1 const& geometry1,
                    Geometry2 const& geometry2,
                    Strategy const& strategy)
{
    if constexpr (std::same_as<Strategy, default_strategy>)
    {
        using strategy_type = typename strategies::relate::services
            ::default_strategy<Geometry1, Geometry2>::type;
        return dispatch::crosses(
            geometry1, geometry2, strategy_type());
    }
    else if constexpr (strategies::detail::is_umbrella_strategy<Strategy>::value)
    {
        return dispatch::crosses(geometry1, geometry2, strategy);
    }
    else
    {
        using strategies::relate::services::strategy_converter;
        auto const converted = strategy_converter<Strategy>::get(strategy);
        return dispatch::crosses(geometry1, geometry2, converted);
    }
}

} // namespace resolve_strategy


namespace resolve_dynamic
{

template <concepts::ConstGeometry Geometry1,
          concepts::ConstGeometry Geometry2,
          typename Strategy>
inline bool crosses(Geometry1 const& geometry1,
                    Geometry2 const& geometry2,
                    Strategy const& strategy)
{
    if constexpr (concepts::ConstDynamicGeometry<Geometry1>
                  && concepts::ConstDynamicGeometry<Geometry2>)
    {
        bool result = false;
        traits::visit<Geometry1, Geometry2>::apply(
            [&](auto const& g1, auto const& g2)
            {
                result = resolve_strategy::crosses(g1, g2, strategy);
            }, geometry1, geometry2);
        return result;
    }
    else if constexpr (concepts::ConstDynamicGeometry<Geometry1>)
    {
        bool result = false;
        traits::visit<Geometry1>::apply([&](auto const& g1)
        {
            result = resolve_strategy::crosses(g1, geometry2, strategy);
        }, geometry1);
        return result;
    }
    else if constexpr (concepts::ConstDynamicGeometry<Geometry2>)
    {
        bool result = false;
        traits::visit<Geometry2>::apply([&](auto const& g2)
        {
            result = resolve_strategy::crosses(geometry1, g2, strategy);
        }, geometry2);
        return result;
    }
    else
    {
        return resolve_strategy::crosses(geometry1, geometry2, strategy);
    }
}


} // namespace resolve_dynamic


/*!
\brief \brief_check2{crosses}
\ingroup crosses
\tparam Geometry1 \tparam_geometry
\tparam Geometry2 \tparam_geometry
\tparam Strategy \tparam_strategy{Crosses}
\param geometry1 \param_geometry
\param geometry2 \param_geometry
\param strategy \param_strategy{crosses}
\return \return_check2{crosses}

\qbk{distinguish,with strategy}
\qbk{[include reference/algorithms/crosses.qbk]}
*/
template <concepts::ConstGeometry Geometry1,
          concepts::ConstGeometry Geometry2,
          typename Strategy>
inline bool crosses(Geometry1 const& geometry1,
                    Geometry2 const& geometry2,
                    Strategy const& strategy)
{
    return resolve_dynamic::crosses(geometry1, geometry2, strategy);
}

/*!
\brief \brief_check2{crosses}
\ingroup crosses
\tparam Geometry1 \tparam_geometry
\tparam Geometry2 \tparam_geometry
\param geometry1 \param_geometry
\param geometry2 \param_geometry
\return \return_check2{crosses}

\qbk{[include reference/algorithms/crosses.qbk]}
\qbk{
[heading Examples]
[crosses]
[crosses_output]
}
*/
template <concepts::ConstGeometry Geometry1,
          concepts::ConstGeometry Geometry2>
inline bool crosses(Geometry1 const& geometry1, Geometry2 const& geometry2)
{
    return resolve_dynamic::crosses(
        geometry1, geometry2, default_strategy());
}

}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_ALGORITHMS_CROSSES_HPP
