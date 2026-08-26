// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2015 Barend Gehrels, Amsterdam, the Netherlands.

// This file was modified by Oracle on 2017-2023.
// Modifications copyright (c) 2017-2023 Oracle and/or its affiliates.
// Contributed and/or modified by Vissarion Fysikopoulos, on behalf of Oracle
// Contributed and/or modified by Adam Wulkiewicz, on behalf of Oracle

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_ALGORITHMS_IS_CONVEX_HPP
#define BOOST_GEOMETRY_ALGORITHMS_IS_CONVEX_HPP


#include <boost/range/empty.hpp>
#include <boost/range/size.hpp>

#include <boost/geometry/algorithms/detail/equals/point_point.hpp>
#include <boost/geometry/algorithms/detail/dummy_geometries.hpp>
#include <boost/geometry/algorithms/detail/visit.hpp>
#include <boost/geometry/core/closure.hpp>
#include <boost/geometry/core/exterior_ring.hpp>
#include <boost/geometry/core/interior_rings.hpp>
#include <boost/geometry/core/visit.hpp>
#include <boost/geometry/geometries/adapted/boost_variant.hpp> // For backward compatibility
#include <boost/geometry/geometries/concepts/check.hpp>
#include <boost/geometry/iterators/ever_circling_iterator.hpp>
#include <boost/geometry/strategies/default_strategy.hpp>
#include <boost/geometry/strategies/is_convex/cartesian.hpp>
#include <boost/geometry/strategies/is_convex/geographic.hpp>
#include <boost/geometry/strategies/is_convex/spherical.hpp>
#include <boost/geometry/views/detail/closed_clockwise_view.hpp>


namespace boost { namespace geometry
{


#ifndef DOXYGEN_NO_DETAIL
namespace detail { namespace is_convex
{

struct ring_is_convex
{
    template <typename Ring, typename Strategies>
    static inline bool apply(Ring const& ring, Strategies const& strategies)
    {
        std::size_t n = boost::size(ring);
        if (n < detail::minimum_ring_size<Ring>::value)
        {
            // (Too) small rings are considered as non-concave, is convex
            return true;
        }

        // Walk in clockwise direction, consider ring as closed
        // (though closure is not important in this algorithm - any dupped
        //  point is skipped)
        using view_type = detail::closed_clockwise_view<Ring const>;
        view_type const view(ring);

        using it_type = geometry::ever_circling_range_iterator<view_type const>;
        it_type previous(view);
        it_type current(view);
        current++;

        auto const equals_strategy = strategies.relate(dummy_point(), dummy_point());

        std::size_t index = 1;
        while (equals::equals_point_point(*current, *previous, equals_strategy)
            && index < n)
        {
            current++;
            index++;
        }

        if (index == n)
        {
            // All points are apparently equal
            return true;
        }

        it_type next = current;
        next++;
        while (equals::equals_point_point(*current, *next, equals_strategy))
        {
            next++;
        }

        auto const side_strategy = strategies.side();

        // We have now three different points on the ring
        // Walk through all points, use a counter because of the ever-circling
        // iterator
        for (std::size_t i = 0; i < n; i++)
        {
            int const side = side_strategy.apply(*previous, *current, *next);
            if (side == 1)
            {
                // Next is on the left side of clockwise ring:
                // the piece is not convex
                return false;
            }

            previous = current;
            current = next;

            // Advance next to next different point
            // (because there are non-equal points, this loop is not infinite)
            next++;
            while (equals::equals_point_point(*current, *next, equals_strategy))
            {
                next++;
            }
        }
        return true;
    }
};


struct polygon_is_convex
{
    template <typename Polygon, typename Strategies>
    static inline bool apply(Polygon const& polygon, Strategies const& strategies)
    {
        return boost::empty(interior_rings(polygon))
            && ring_is_convex::apply(exterior_ring(polygon), strategies);
    }
};

struct multi_polygon_is_convex
{
    template <typename MultiPolygon, typename Strategies>
    static inline bool apply(MultiPolygon const& multi_polygon, Strategies const& strategies)
    {
        auto const size = boost::size(multi_polygon);
        // TODO: this looks wrong, it should only return convex if all its rings are convex
        return size == 0 // For consistency with ring_is_convex
            || (size == 1 && polygon_is_convex::apply(range::front(multi_polygon), strategies));
    }
};


}} // namespace detail::is_convex
#endif // DOXYGEN_NO_DETAIL


namespace resolve_strategy {

template <concepts::ConstGeometry Geometry, typename Strategy>
inline bool is_convex(Geometry const& geometry, Strategy const& strategy)
{
    if constexpr (concepts::ConstDynamicGeometry<Geometry>)
    {
        bool result = false;
        traits::visit<Geometry>::apply([&](auto const& g)
        {
            result = resolve_strategy::is_convex(g, strategy);
        }, geometry);
        return result;
    }
    else if constexpr (concepts::ConstGeometryCollection<Geometry>)
    {
        bool result = false;
        bool is_first = true;
        detail::visit_breadth_first([&](auto const& g)
        {
            result = is_first && resolve_strategy::is_convex(g, strategy);
            is_first = false;
            return result;
        }, geometry);
        return result;
    }
    else if constexpr (std::same_as<Strategy, default_strategy>)
    {
        using strategy_type = typename strategies::is_convex::services::default_strategy
            <Geometry>::type;
        return resolve_strategy::is_convex(geometry, strategy_type());
    }
    else if constexpr (! strategies::detail::is_umbrella_strategy<Strategy>::value)
    {
        using strategies::is_convex::services::strategy_converter;
        return resolve_strategy::is_convex(
            geometry, strategy_converter<Strategy>::get(strategy));
    }
    else if constexpr (concepts::ConstBox<Geometry>)
    {
        return true;
    }
    else if constexpr (concepts::ConstRing<Geometry>)
    {
        return detail::is_convex::ring_is_convex::apply(geometry, strategy);
    }
    else if constexpr (concepts::ConstPolygon<Geometry>)
    {
        return detail::is_convex::polygon_is_convex::apply(geometry, strategy);
    }
    else if constexpr (concepts::ConstMultiPolygon<Geometry>)
    {
        return detail::is_convex::multi_polygon_is_convex::apply(
            geometry, strategy);
    }
    else
    {
        return false;
    }
}

} // namespace resolve_strategy

// TODO: documentation / qbk
template<concepts::ConstGeometry Geometry>
inline bool is_convex(Geometry const& geometry)
{
    return resolve_strategy::is_convex(geometry, geometry::default_strategy());
}

// TODO: documentation / qbk
template<concepts::ConstGeometry Geometry, typename Strategy>
inline bool is_convex(Geometry const& geometry, Strategy const& strategy)
{
    return resolve_strategy::is_convex(geometry, strategy);
}


}} // namespace boost::geometry


#endif // BOOST_GEOMETRY_ALGORITHMS_IS_CONVEX_HPP
