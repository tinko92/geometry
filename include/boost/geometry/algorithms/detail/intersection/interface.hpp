// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2007-2012 Barend Gehrels, Amsterdam, the Netherlands.

// This file was modified by Oracle on 2014-2024.
// Modifications copyright (c) 2014-2024, Oracle and/or its affiliates.
// Contributed and/or modified by Vissarion Fysikopoulos, on behalf of Oracle
// Contributed and/or modified by Adam Wulkiewicz, on behalf of Oracle

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_ALGORITHMS_DETAIL_INTERSECTION_INTERFACE_HPP
#define BOOST_GEOMETRY_ALGORITHMS_DETAIL_INTERSECTION_INTERFACE_HPP


#include <boost/geometry/algorithms/detail/overlay/intersection_insert.hpp>
#include <boost/geometry/algorithms/detail/tupled_output.hpp>
#include <boost/geometry/geometries/adapted/std_variant.hpp>
#include <boost/geometry/strategies/default_strategy.hpp>
#include <boost/geometry/strategies/detail.hpp>
#include <boost/geometry/strategies/relate/services.hpp>
#include <boost/geometry/util/range.hpp>
#include <boost/geometry/util/type_traits_std.hpp>


namespace boost { namespace geometry
{


#ifndef DOXYGEN_NO_DISPATCH
namespace dispatch
{

template
<
    typename Geometry1, typename Geometry2,
    typename Tag1 = geometry::tag_t<Geometry1>,
    typename Tag2 = geometry::tag_t<Geometry2>,
    bool Reverse = reverse_dispatch<Geometry1, Geometry2>::type::value
>
struct intersection
{
    template <typename GeometryOut, typename Strategy>
    static inline bool apply(Geometry1 const& geometry1,
                             Geometry2 const& geometry2,
                             GeometryOut& geometry_out,
                             Strategy const& strategy)
    {
        using single_out = typename geometry::detail::output_geometry_value
            <
                GeometryOut
            >::type;

        intersection_insert
            <
                Geometry1, Geometry2, single_out,
                overlay_intersection
            >::apply(geometry1, geometry2,
                     geometry::detail::output_geometry_back_inserter(geometry_out),
                     strategy);

        return true;
    }
};

template <typename Geometry1, typename Geometry2,
          typename Tag1, typename Tag2>
struct intersection<Geometry1, Geometry2, Tag1, Tag2, true>
    : intersection<Geometry2, Geometry1, Tag2, Tag1, false>
{
    template <typename GeometryOut, typename Strategy>
    static inline bool apply(Geometry1 const& geometry1,
                             Geometry2 const& geometry2,
                             GeometryOut& geometry_out,
                             Strategy const& strategy)
    {
        return intersection<Geometry2, Geometry1, Tag2, Tag1, false>::apply(
            geometry2, geometry1, geometry_out, strategy);
    }
};

} // namespace dispatch
#endif // DOXYGEN_NO_DISPATCH


namespace resolve_collection
{

template
<
    typename Geometry1, typename Geometry2, typename GeometryOut,
    typename Tag1 = geometry::tag_t<Geometry1>,
    typename Tag2 = geometry::tag_t<Geometry2>,
    typename TagOut = geometry::tag_t<GeometryOut>
>
struct intersection
{
    template <typename Strategy>
    static bool apply(Geometry1 const& geometry1,
                      Geometry2 const& geometry2,
                      GeometryOut& geometry_out,
                      Strategy const& strategy)
    {
        return dispatch::intersection<Geometry1, Geometry2>::apply(
            geometry1, geometry2, geometry_out, strategy);
    }
};

} // namespace resolve_collection


namespace resolve_strategy
{

template <typename Geometry1, typename Geometry2,
          typename GeometryOut, typename Strategy>
inline bool intersection(Geometry1 const& geometry1,
                         Geometry2 const& geometry2,
                         GeometryOut& geometry_out,
                         Strategy const& strategy)
{
    if constexpr (std::same_as<Strategy, default_strategy>)
    {
        using strategy_type = typename strategies::relate::services::default_strategy
            <
                Geometry1, Geometry2
            >::type;
        return resolve_collection::intersection
            <
                Geometry1, Geometry2, GeometryOut
            >::apply(geometry1, geometry2, geometry_out, strategy_type());
    }
    else if constexpr (strategies::detail::is_umbrella_strategy<Strategy>::value)
    {
        return resolve_collection::intersection
            <
                Geometry1, Geometry2, GeometryOut
            >::apply(geometry1, geometry2, geometry_out, strategy);
    }
    else
    {
        using strategies::relate::services::strategy_converter;
        return resolve_collection::intersection
            <
                Geometry1, Geometry2, GeometryOut
            >::apply(geometry1, geometry2, geometry_out,
                     strategy_converter<Strategy>::get(strategy));
    }
}

} // resolve_strategy


namespace resolve_dynamic
{

template <concepts::ConstGeometry Geometry1,
          concepts::ConstGeometry Geometry2,
          typename GeometryOut, typename Strategy>
inline bool intersection(Geometry1 const& geometry1,
                         Geometry2 const& geometry2,
                         GeometryOut& geometry_out,
                         Strategy const& strategy)
{
    if constexpr (concepts::ConstDynamicGeometry<Geometry1>
                  && concepts::ConstDynamicGeometry<Geometry2>)
    {
        bool result = false;
        traits::visit<Geometry1, Geometry2>::apply(
            [&](auto const& g1, auto const& g2)
            {
                result = resolve_strategy::intersection(
                    g1, g2, geometry_out, strategy);
            }, geometry1, geometry2);
        return result;
    }
    else if constexpr (concepts::ConstDynamicGeometry<Geometry1>)
    {
        bool result = false;
        traits::visit<Geometry1>::apply([&](auto const& g1)
        {
            result = resolve_strategy::intersection(
                g1, geometry2, geometry_out, strategy);
        }, geometry1);
        return result;
    }
    else if constexpr (concepts::ConstDynamicGeometry<Geometry2>)
    {
        bool result = false;
        traits::visit<Geometry2>::apply([&](auto const& g2)
        {
            result = resolve_strategy::intersection(
                geometry1, g2, geometry_out, strategy);
        }, geometry2);
        return result;
    }
    else
    {
        return resolve_strategy::intersection(
            geometry1, geometry2, geometry_out, strategy);
    }
}

} // namespace resolve_dynamic


/*!
\brief \brief_calc2{intersection}
\ingroup intersection
\details \details_calc2{intersection, spatial set theoretic intersection}.
\tparam Geometry1 \tparam_geometry
\tparam Geometry2 \tparam_geometry
\tparam GeometryOut Collection of geometries (e.g. std::vector, std::deque, boost::geometry::multi*) of which
    the value_type fulfills a \p_l_or_c concept, or it is the output geometry (e.g. for a box)
\tparam Strategy \tparam_strategy{Intersection}
\param geometry1 \param_geometry
\param geometry2 \param_geometry
\param geometry_out The output geometry, either a multi_point, multi_polygon,
    multi_linestring, or a box (for intersection of two boxes)
\param strategy \param_strategy{intersection}

\qbk{distinguish,with strategy}
\qbk{[include reference/algorithms/intersection.qbk]}
*/
template
<
    typename Geometry1,
    typename Geometry2,
    typename GeometryOut,
    typename Strategy
>
inline bool intersection(Geometry1 const& geometry1,
                         Geometry2 const& geometry2,
                         GeometryOut& geometry_out,
                         Strategy const& strategy)
{
    return resolve_dynamic::intersection(
        geometry1, geometry2, geometry_out, strategy);
}


/*!
\brief \brief_calc2{intersection}
\ingroup intersection
\details \details_calc2{intersection, spatial set theoretic intersection}.
\tparam Geometry1 \tparam_geometry
\tparam Geometry2 \tparam_geometry
\tparam GeometryOut Collection of geometries (e.g. std::vector, std::deque, boost::geometry::multi*) of which
    the value_type fulfills a \p_l_or_c concept, or it is the output geometry (e.g. for a box)
\param geometry1 \param_geometry
\param geometry2 \param_geometry
\param geometry_out The output geometry, either a multi_point, multi_polygon,
    multi_linestring, or a box (for intersection of two boxes)

\qbk{[include reference/algorithms/intersection.qbk]}
*/
template
<
    typename Geometry1,
    typename Geometry2,
    typename GeometryOut
>
inline bool intersection(Geometry1 const& geometry1,
                         Geometry2 const& geometry2,
                         GeometryOut& geometry_out)
{
    return resolve_dynamic::intersection(
        geometry1, geometry2, geometry_out, default_strategy());
}


}} // namespace boost::geometry


#endif // BOOST_GEOMETRY_ALGORITHMS_DETAIL_INTERSECTION_INTERFACE_HPP
