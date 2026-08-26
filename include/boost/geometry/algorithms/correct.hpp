// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2007-2012 Barend Gehrels, Amsterdam, the Netherlands.
// Copyright (c) 2008-2012 Bruno Lalande, Paris, France.
// Copyright (c) 2009-2012 Mateusz Loskot, London, UK.
// Copyright (c) 2014-2024 Adam Wulkiewicz, Lodz, Poland.

// This file was modified by Oracle on 2017-2023.
// Modifications copyright (c) 2017-2023 Oracle and/or its affiliates.
// Contributed and/or modified by Vissarion Fysikopoulos, on behalf of Oracle
// Contributed and/or modified by Adam Wulkiewicz, on behalf of Oracle

// Parts of Boost.Geometry are redesigned from Geodan's Geographic Library
// (geolib/GGL), copyright (c) 1995-2010 Geodan, Amsterdam, the Netherlands.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_ALGORITHMS_CORRECT_HPP
#define BOOST_GEOMETRY_ALGORITHMS_CORRECT_HPP


#include <algorithm>
#include <functional>

#include <boost/range/begin.hpp>
#include <boost/range/end.hpp>

#include <boost/geometry/algorithms/area.hpp>
#include <boost/geometry/algorithms/correct_closure.hpp>
#include <boost/geometry/algorithms/detail/multi_modify.hpp>
#include <boost/geometry/algorithms/detail/visit.hpp>

#include <boost/geometry/core/exterior_ring.hpp>
#include <boost/geometry/core/interior_rings.hpp>
#include <boost/geometry/core/tags.hpp>
#include <boost/geometry/core/visit.hpp>

#include <boost/geometry/geometries/adapted/std_variant.hpp>
#include <boost/geometry/geometries/concepts/check.hpp>

#include <boost/geometry/strategies/area/cartesian.hpp>
#include <boost/geometry/strategies/area/geographic.hpp>
#include <boost/geometry/strategies/area/spherical.hpp>
#include <boost/geometry/strategies/detail.hpp>

#include <boost/geometry/util/algorithm.hpp>


namespace boost { namespace geometry
{

// Silence warning C4127: conditional expression is constant
#if defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable : 4127)
#endif

#ifndef DOXYGEN_NO_DETAIL
namespace detail { namespace correct
{

struct correct_nop
{
    template <typename Geometry, typename Strategy>
    static inline void apply(Geometry& , Strategy const& )
    {}
};


// Correct a box: make min/max correct
struct correct_box
{
    template <typename Box, typename Strategy>
    static inline void apply(Box& box, Strategy const& )
    {
        // Currently only for Cartesian coordinates
        // (or spherical without crossing dateline)
        // Future version: adapt using strategies
        detail::for_each_dimension<Box>([&](auto dimension)
        {
            if (get<min_corner, dimension>(box) > get<max_corner, dimension>(box))
            {
                // Swap the coordinates
                auto max_value = get<min_corner, dimension>(box);
                auto min_value = get<max_corner, dimension>(box);
                set<min_corner, dimension>(box, min_value);
                set<max_corner, dimension>(box, max_value);
            }
        });
    }
};


// Close a ring, if not closed
template <typename Predicate = std::less<>>
struct correct_ring
{
    template <typename Ring, typename Strategy>
    static inline void apply(Ring& r, Strategy const& strategy)
    {
        // Correct closure if necessary
        detail::correct_closure::close_or_open_ring::apply(r);

        // NOTE: calculate_point_order should probably be used here instead.

        // Check area
        using area_t = typename area_result<Ring, Strategy>::type;
        area_t const zero = 0;
        if (Predicate()(detail::area_impl(r, strategy), zero))
        {
            std::reverse(boost::begin(r), boost::end(r));
        }
    }
};

// Correct a polygon: normalizes all rings, sets outer ring clockwise, sets all
// inner rings counter clockwise (or vice versa depending on orientation)
struct correct_polygon
{
    template <typename Polygon, typename Strategy>
    static inline void apply(Polygon& poly, Strategy const& strategy)
    {
        correct_ring<std::less<>>::apply(exterior_ring(poly), strategy);

        auto&& rings = interior_rings(poly);
        auto const end = boost::end(rings);
        for (auto it = boost::begin(rings); it != end; ++it)
        {
            correct_ring<std::greater<>>::apply(*it, strategy);
        }
    }
};


}} // namespace detail::correct
#endif // DOXYGEN_NO_DETAIL


namespace resolve_strategy
{

template <concepts::MutableGeometry Geometry, typename Strategy>
inline void correct(Geometry& geometry, Strategy const& strategy)
{
    if constexpr (concepts::DynamicGeometry<Geometry>)
    {
        traits::visit<Geometry>::apply([&](auto& g)
        {
            resolve_strategy::correct(g, strategy);
        }, geometry);
    }
    else if constexpr (concepts::GeometryCollection<Geometry>)
    {
        detail::visit_breadth_first([&](auto& g)
        {
            resolve_strategy::correct(g, strategy);
            return true;
        }, geometry);
    }
    else if constexpr (std::same_as<Strategy, default_strategy>)
    {
        using strategy_type = typename strategies::area::services::default_strategy
            <Geometry>::type;
        resolve_strategy::correct(geometry, strategy_type());
    }
    else if constexpr (! strategies::detail::is_umbrella_strategy<Strategy>::value)
    {
        using geometry::strategies::area::services::strategy_converter;
        resolve_strategy::correct(
            geometry, strategy_converter<Strategy>::get(strategy));
    }
    else if constexpr (concepts::Box<Geometry>)
    {
        detail::correct::correct_box::apply(geometry, strategy);
    }
    else if constexpr (concepts::Ring<Geometry>)
    {
        detail::correct::correct_ring<>::apply(geometry, strategy);
    }
    else if constexpr (concepts::Polygon<Geometry>)
    {
        detail::correct::correct_polygon::apply(geometry, strategy);
    }
    else if constexpr (concepts::MultiPolygon<Geometry>)
    {
        for (auto it = boost::begin(geometry); it != boost::end(geometry); ++it)
        {
            detail::correct::correct_polygon::apply(*it, strategy);
        }
    }
}

} // namespace resolve_strategy


/*!
\brief Corrects a geometry
\details Corrects a geometry: all rings which are wrongly oriented with respect
    to their expected orientation are reversed. To all rings which do not have a
    closing point and are typed as they should have one, the first point is
    appended. Also boxes can be corrected.
\ingroup correct
\tparam Geometry \tparam_geometry
\param geometry \param_geometry which will be corrected if necessary

\qbk{[include reference/algorithms/correct.qbk]}
*/
template <concepts::MutableGeometry Geometry>
inline void correct(Geometry& geometry)
{
    resolve_strategy::correct(geometry, default_strategy());
}

/*!
\brief Corrects a geometry
\details Corrects a geometry: all rings which are wrongly oriented with respect
    to their expected orientation are reversed. To all rings which do not have a
    closing point and are typed as they should have one, the first point is
    appended. Also boxes can be corrected.
\ingroup correct
\tparam Geometry \tparam_geometry
\tparam Strategy \tparam_strategy{Area}
\param geometry \param_geometry which will be corrected if necessary
\param strategy \param_strategy{area}

\qbk{distinguish,with strategy}

\qbk{[include reference/algorithms/correct.qbk]}
*/
template <concepts::MutableGeometry Geometry, typename Strategy>
inline void correct(Geometry& geometry, Strategy const& strategy)
{
    resolve_strategy::correct(geometry, strategy);
}

#if defined(_MSC_VER)
#pragma warning(pop)
#endif

}} // namespace boost::geometry


#endif // BOOST_GEOMETRY_ALGORITHMS_CORRECT_HPP
