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
// #include <boost/geometry/algorithms/detail/throw_on_empty_input.hpp>
#include <boost/geometry/algorithms/detail/visit.hpp>

#include <boost/geometry/core/closure.hpp>
#include <boost/geometry/core/exterior_ring.hpp>
#include <boost/geometry/core/interior_rings.hpp>
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

#ifndef DOXYGEN_NO_DETAIL
namespace detail
{

template <concepts::ConstDynamicGeometry DynamicGeometry, typename Strategy>
inline typename default_length_result<DynamicGeometry>::type
perimeter_impl(DynamicGeometry const& dynamic, Strategy const& strategy);

template <concepts::ConstGeometryCollection GeometryCollection, typename Strategy>
inline typename default_length_result<GeometryCollection>::type
perimeter_impl(GeometryCollection const& collection, Strategy const& strategy);

template <concepts::ConstRing Ring, typename Strategy>
inline typename default_length_result<Ring>::type
perimeter_impl(Ring const& ring, Strategy const& strategy)
{
    auto&& strategies = resolve_length_strategy(ring, strategy);
    return range_length<closure<Ring>::value>(ring, strategies);
}

template <concepts::ConstPolygon Polygon, typename Strategy>
inline typename default_length_result<Polygon>::type
perimeter_impl(Polygon const& polygon, Strategy const& strategy)
{
    auto&& strategies = resolve_length_strategy(polygon, strategy);
    typename default_length_result<Polygon>::type result
        = perimeter_impl(exterior_ring(polygon), strategies);
    auto const& rings = interior_rings(polygon);
    for (auto it = boost::begin(rings); it != boost::end(rings); ++it)
    {
        result += perimeter_impl(*it, strategies);
    }
    return result;
}

template <concepts::ConstMultiPolygon MultiPolygon, typename Strategy>
inline typename default_length_result<MultiPolygon>::type
perimeter_impl(MultiPolygon const& multi, Strategy const& strategy)
{
    auto&& strategies = resolve_length_strategy(multi, strategy);
    typename default_length_result<MultiPolygon>::type result = 0;
    for (auto it = boost::begin(multi); it != boost::end(multi); ++it)
    {
        result += perimeter_impl(*it, strategies);
    }
    return result;
}

template <concepts::ConstGeometry Geometry, typename Strategy>
    requires (! concepts::ConstRing<Geometry>)
          && (! concepts::ConstPolygon<Geometry>)
          && (! concepts::ConstMultiPolygon<Geometry>)
          && (! concepts::ConstDynamicGeometry<Geometry>)
          && (! concepts::ConstGeometryCollection<Geometry>)
inline typename default_length_result<Geometry>::type
perimeter_impl(Geometry const&, Strategy const&)
{
    return 0;
}

template <concepts::ConstDynamicGeometry DynamicGeometry, typename Strategy>
inline typename default_length_result<DynamicGeometry>::type
perimeter_impl(DynamicGeometry const& dynamic, Strategy const& strategy)
{
    typename default_length_result<DynamicGeometry>::type result = 0;
    traits::visit<DynamicGeometry>::apply([&](auto const& geometry)
    {
        result = perimeter_impl(geometry, strategy);
    }, dynamic);
    return result;
}

template <concepts::ConstGeometryCollection GeometryCollection, typename Strategy>
inline typename default_length_result<GeometryCollection>::type
perimeter_impl(GeometryCollection const& collection, Strategy const& strategy)
{
    typename default_length_result<GeometryCollection>::type result = 0;
    detail::visit_breadth_first([&](auto const& geometry)
    {
        result += perimeter_impl(geometry, strategy);
        return true;
    }, collection);
    return result;
}

} // namespace detail
#endif // DOXYGEN_NO_DETAIL


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
    return detail::perimeter_impl(geometry, default_strategy());
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
    return detail::perimeter_impl(geometry, strategy);
}

}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_ALGORITHMS_PERIMETER_HPP
