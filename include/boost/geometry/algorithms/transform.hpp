// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2007-2012 Barend Gehrels, Amsterdam, the Netherlands.
// Copyright (c) 2008-2012 Bruno Lalande, Paris, France.
// Copyright (c) 2009-2012 Mateusz Loskot, London, UK.
// Copyright (c) 2014 Adam Wulkiewicz, Lodz, Poland.

// This file was modified by Oracle on 2020-2023.
// Modifications copyright (c) 2020-2023 Oracle and/or its affiliates.
// Contributed and/or modified by Vissarion Fysikopoulos, on behalf of Oracle
// Contributed and/or modified by Adam Wulkiewicz, on behalf of Oracle

// Parts of Boost.Geometry are redesigned from Geodan's Geographic Library
// (geolib/GGL), copyright (c) 1995-2010 Geodan, Amsterdam, the Netherlands.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_ALGORITHMS_TRANSFORM_HPP
#define BOOST_GEOMETRY_ALGORITHMS_TRANSFORM_HPP

#include <type_traits>

#include <boost/range/begin.hpp>
#include <boost/range/end.hpp>
#include <boost/range/size.hpp>
#include <boost/range/value_type.hpp>

#include <boost/geometry/algorithms/clear.hpp>
#include "boost/geometry/algorithms/detail/assign_indexed_point.hpp"
#include "boost/geometry/algorithms/detail/assign_values.hpp"
#include <boost/geometry/algorithms/num_interior_rings.hpp>

#include <boost/geometry/core/cs.hpp>
#include <boost/geometry/core/exterior_ring.hpp>
#include <boost/geometry/core/interior_rings.hpp>
#include <boost/geometry/core/mutable_range.hpp>
#include <boost/geometry/core/tag_cast.hpp>
#include <boost/geometry/core/tags.hpp>
#include <boost/geometry/core/visit.hpp>
#include <boost/geometry/geometries/concepts/check.hpp>
#include <boost/geometry/strategies/default_strategy.hpp>
#include <boost/geometry/strategies/transform.hpp>
#include <boost/geometry/util/range.hpp>


namespace boost { namespace geometry
{

#ifndef DOXYGEN_NO_DETAIL
namespace detail { namespace transform
{

struct transform_point
{
    template <typename Point1, typename Point2, typename Strategy>
    static inline bool apply(Point1 const& p1, Point2& p2,
                Strategy const& strategy)
    {
        return strategy.apply(p1, p2);
    }
};


struct transform_box
{
    template <typename Box1, typename Box2, typename Strategy>
    static inline bool apply(Box1 const& b1, Box2& b2,
                Strategy const& strategy)
    {
        using point_type1 = point_type_t<Box1>;
        using point_type2 = point_type_t<Box2>;

        point_type1 lower_left, upper_right;
        geometry::detail::assign::assign_box_2d_corner<min_corner, min_corner>(
                    b1, lower_left);
        geometry::detail::assign::assign_box_2d_corner<max_corner, max_corner>(
                    b1, upper_right);

        point_type2 p1, p2;
        if (strategy.apply(lower_left, p1) && strategy.apply(upper_right, p2))
        {
            // Create a valid box and therefore swap if necessary
            using coordinate_type = coordinate_type_t<point_type2>;
            coordinate_type x1 = geometry::get<0>(p1);
            coordinate_type y1 = geometry::get<1>(p1);
            coordinate_type x2 = geometry::get<0>(p2);
            coordinate_type y2 = geometry::get<1>(p2);

            if (x1 > x2) { std::swap(x1, x2); }
            if (y1 > y2) { std::swap(y1, y2); }

            geometry::set<min_corner, 0>(b2, x1);
            geometry::set<min_corner, 1>(b2, y1);
            geometry::set<max_corner, 0>(b2, x2);
            geometry::set<max_corner, 1>(b2, y2);

            return true;
        }
        return false;
    }
};

struct transform_box_or_segment
{
    template <typename Geometry1, typename Geometry2, typename Strategy>
    static inline bool apply(Geometry1 const& source, Geometry2& target,
                Strategy const& strategy)
    {
        point_type_t<Geometry1> source_point[2];
        geometry::detail::assign_point_from_index<0>(source, source_point[0]);
        geometry::detail::assign_point_from_index<1>(source, source_point[1]);

        point_type_t<Geometry2> target_point[2];
        if (strategy.apply(source_point[0], target_point[0])
            && strategy.apply(source_point[1], target_point[1]))
        {
            geometry::detail::assign_point_to_index<0>(target_point[0], target);
            geometry::detail::assign_point_to_index<1>(target_point[1], target);
            return true;
        }
        return false;
    }
};


template
<
    typename PointOut,
    typename OutputIterator,
    typename Range,
    typename Strategy
>
inline bool transform_range_out(Range const& range,
    OutputIterator out, Strategy const& strategy)
{
    PointOut point_out;
    for (auto it = boost::begin(range); it != boost::end(range); ++it)
    {
        if (! transform_point::apply(*it, point_out, strategy))
        {
            return false;
        }
        *out++ = point_out;
    }
    return true;
}


struct transform_polygon
{
    template <typename Polygon1, typename Polygon2, typename Strategy>
    static inline bool apply(Polygon1 const& poly1, Polygon2& poly2,
                Strategy const& strategy)
    {
        using point2_type = point_type_t<Polygon2>;

        geometry::clear(poly2);

        if (!transform_range_out<point2_type>(geometry::exterior_ring(poly1),
                    range::back_inserter(geometry::exterior_ring(poly2)), strategy))
        {
            return false;
        }

        // Note: here a resizeable container is assumed.
        traits::resize
            <
                typename std::remove_reference
                <
                    typename traits::interior_mutable_type<Polygon2>::type
                >::type
            >::apply(geometry::interior_rings(poly2),
                     geometry::num_interior_rings(poly1));

        auto const& rings1 = geometry::interior_rings(poly1);
        auto&& rings2 = geometry::interior_rings(poly2);

        auto it1 = boost::begin(rings1);
        auto it2 = boost::begin(rings2);
        for ( ; it1 != boost::end(rings1); ++it1, ++it2)
        {
            if ( ! transform_range_out<point2_type>(*it1,
                                                    range::back_inserter(*it2),
                                                    strategy) )
            {
                return false;
            }
        }

        return true;
    }
};


template <typename Geometry1, typename Geometry2>
struct select_strategy
{
    using type = typename strategy::transform::services::default_strategy
        <
            cs_tag_t<Geometry1>,
            cs_tag_t<Geometry2>,
            coordinate_system_t<Geometry1>,
            coordinate_system_t<Geometry2>,
            dimension<Geometry1>::value,
            dimension<Geometry2>::value,
            point_type_t<Geometry1>,
            point_type_t<Geometry2>
        >::type;
};

struct transform_range
{
    template <typename Range1, typename Range2, typename Strategy>
    static inline bool apply(Range1 const& range1,
            Range2& range2, Strategy const& strategy)
    {
        // "clear" should NOT be done here!
        // geometry::clear(range2);
        return transform_range_out<point_type_t<Range2>>(range1,
                range::back_inserter(range2), strategy);
    }
};


/*!
    \brief Is able to transform any multi-geometry, calling the single-version as policy
*/
template <typename Policy>
struct transform_multi
{
    template <typename Multi1, typename Multi2, typename S>
    static inline bool apply(Multi1 const& multi1, Multi2& multi2, S const& strategy)
    {
        traits::resize<Multi2>::apply(multi2, boost::size(multi1));

        auto it1 = boost::begin(multi1);
        auto it2 = boost::begin(multi2);

        for (; it1 != boost::end(multi1); ++it1, ++it2)
        {
            if (! Policy::apply(*it1, *it2, strategy))
            {
                return false;
            }
        }

        return true;
    }
};


}} // namespace detail::transform
#endif // DOXYGEN_NO_DETAIL


#ifndef DOXYGEN_NO_DISPATCH
namespace dispatch
{

template <concepts::ConstGeometry Geometry1,
          concepts::MutableGeometry Geometry2,
          typename Strategy>
    requires (concepts::ConstPoint<Geometry1> && concepts::Point<Geometry2>)
          || (concepts::ConstLinestring<Geometry1> && concepts::Linestring<Geometry2>)
          || (concepts::ConstRing<Geometry1> && concepts::Ring<Geometry2>)
          || (concepts::ConstPolygon<Geometry1> && concepts::Polygon<Geometry2>)
          || (concepts::ConstBox<Geometry1> && concepts::Box<Geometry2>)
          || (concepts::ConstSegment<Geometry1> && concepts::Segment<Geometry2>)
          || (concepts::ConstMultiPoint<Geometry1> && concepts::MultiPoint<Geometry2>)
          || (concepts::ConstMultiLinestring<Geometry1> && concepts::MultiLinestring<Geometry2>)
          || (concepts::ConstMultiPolygon<Geometry1> && concepts::MultiPolygon<Geometry2>)
inline bool transform(Geometry1 const& geometry1, Geometry2& geometry2,
                      Strategy const& strategy)
{
    if constexpr (concepts::ConstPoint<Geometry1>)
    {
        return detail::transform::transform_point::apply(
            geometry1, geometry2, strategy);
    }
    else if constexpr (concepts::ConstLinestring<Geometry1>
                       || concepts::ConstRing<Geometry1>)
    {
        return detail::transform::transform_range::apply(
            geometry1, geometry2, strategy);
    }
    else if constexpr (concepts::ConstPolygon<Geometry1>)
    {
        return detail::transform::transform_polygon::apply(
            geometry1, geometry2, strategy);
    }
    else if constexpr (concepts::ConstBox<Geometry1>)
    {
        return detail::transform::transform_box::apply(
            geometry1, geometry2, strategy);
    }
    else if constexpr (concepts::ConstSegment<Geometry1>)
    {
        return detail::transform::transform_box_or_segment::apply(
            geometry1, geometry2, strategy);
    }
    else
    {
        range::resize(geometry2, boost::size(geometry1));
        auto out = boost::begin(geometry2);
        for (auto it = boost::begin(geometry1); it != boost::end(geometry1); ++it)
        {
            if (! dispatch::transform(*it, *out++, strategy))
            {
                return false;
            }
        }
        return true;
    }
}


} // namespace dispatch
#endif // DOXYGEN_NO_DISPATCH


namespace resolve_strategy {

template <concepts::ConstGeometry Geometry1,
          concepts::MutableGeometry Geometry2,
          typename Strategy>
inline bool transform(Geometry1 const& geometry1, Geometry2& geometry2,
                      Strategy const& strategy)
{
    if constexpr (std::same_as<Strategy, default_strategy>)
    {
        using strategy_type = typename detail::transform
            ::select_strategy<Geometry1, Geometry2>::type;
        return dispatch::transform(geometry1, geometry2, strategy_type());
    }
    else
    {
        return dispatch::transform(geometry1, geometry2, strategy);
    }
}

} // namespace resolve_strategy


namespace resolve_dynamic {

template <concepts::ConstGeometry Geometry1,
          concepts::MutableGeometry Geometry2,
          typename Strategy>
inline bool transform(Geometry1 const& geometry1, Geometry2& geometry2,
                      Strategy const& strategy)
{
    if constexpr (concepts::ConstDynamicGeometry<Geometry1>)
    {
        bool result = false;
        traits::visit<Geometry1>::apply([&](auto const& source)
        {
            result = resolve_strategy::transform(source, geometry2, strategy);
        }, geometry1);
        return result;
    }
    else
    {
        return resolve_strategy::transform(geometry1, geometry2, strategy);
    }
}

} // namespace resolve_dynamic


/*!
\brief Transforms from one geometry to another geometry  \brief_strategy
\ingroup transform
\tparam Geometry1 \tparam_geometry
\tparam Geometry2 \tparam_geometry
\tparam Strategy strategy
\param geometry1 \param_geometry
\param geometry2 \param_geometry
\param strategy The strategy to be used for transformation
\return True if the transformation could be done

\qbk{distinguish,with strategy}

\qbk{[include reference/algorithms/transform_with_strategy.qbk]}
 */
template <concepts::ConstGeometry Geometry1,
          concepts::MutableGeometry Geometry2,
          typename Strategy>
inline bool transform(Geometry1 const& geometry1, Geometry2& geometry2,
            Strategy const& strategy)
{
    return resolve_dynamic::transform(geometry1, geometry2, strategy);
}


/*!
\brief Transforms from one geometry to another geometry using a strategy
\ingroup transform
\tparam Geometry1 \tparam_geometry
\tparam Geometry2 \tparam_geometry
\param geometry1 \param_geometry
\param geometry2 \param_geometry
\return True if the transformation could be done

\qbk{[include reference/algorithms/transform.qbk]}
 */
template <concepts::ConstGeometry Geometry1,
          concepts::MutableGeometry Geometry2>
inline bool transform(Geometry1 const& geometry1, Geometry2& geometry2)
{
    return geometry::transform(geometry1, geometry2, default_strategy());
}


}} // namespace boost::geometry


#endif // BOOST_GEOMETRY_ALGORITHMS_TRANSFORM_HPP
