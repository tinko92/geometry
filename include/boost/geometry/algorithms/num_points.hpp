// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2007-2014 Barend Gehrels, Amsterdam, the Netherlands.
// Copyright (c) 2008-2014 Bruno Lalande, Paris, France.
// Copyright (c) 2009-2014 Mateusz Loskot, London, UK.
// Copyright (c) 2014 Adam Wulkiewicz, Lodz, Poland.

// This file was modified by Oracle on 2014-2021.
// Modifications copyright (c) 2014-2021, Oracle and/or its affiliates.

// Contributed and/or modified by Menelaos Karavelas, on behalf of Oracle
// Contributed and/or modified by Adam Wulkiewicz, on behalf of Oracle

// Parts of Boost.Geometry are redesigned from Geodan's Geographic Library
// (geolib/GGL), copyright (c) 1995-2010 Geodan, Amsterdam, the Netherlands.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_ALGORITHMS_NUM_POINTS_HPP
#define BOOST_GEOMETRY_ALGORITHMS_NUM_POINTS_HPP

#include <cstddef>

#include <boost/range/size.hpp>
#include <boost/range/value_type.hpp>

#include <boost/geometry/algorithms/detail/counting.hpp>
#include <boost/geometry/algorithms/detail/visit.hpp>
#include <boost/geometry/algorithms/not_implemented.hpp>

#include <boost/geometry/core/closure.hpp>
#include <boost/geometry/core/coordinate_dimension.hpp>
#include <boost/geometry/core/tag_cast.hpp>
#include <boost/geometry/core/tags.hpp>
#include <boost/geometry/core/visit.hpp>

#include <boost/geometry/geometries/adapted/boost_variant.hpp> // For backward compatibility
#include <boost/geometry/geometries/concepts/check.hpp>

#include <boost/geometry/util/type_traits_std.hpp>

namespace boost { namespace geometry
{

// Silence warning C4127: conditional expression is constant
#if defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable : 4127)
#endif


#ifndef DOXYGEN_NO_DETAIL
namespace detail { namespace num_points
{


template <bool AddForOpen>
struct range_count
{
    template <typename Range>
    static inline std::size_t apply(Range const& range)
    {
        std::size_t n = boost::size(range);
        if (AddForOpen
            && n > 0
            && geometry::closure<Range>::value == open
            )
        {
            return n + 1;
        }
        return n;
    }
};

}} // namespace detail::num_points
#endif // DOXYGEN_NO_DETAIL


namespace resolve_dynamic
{

template <concepts::ConstGeometry Geometry>
inline std::size_t num_points(Geometry const& geometry, bool add_for_open)
{
    if constexpr (concepts::ConstDynamicGeometry<Geometry>)
    {
        std::size_t result = 0;
        traits::visit<Geometry>::apply([&](auto const& g)
        {
            result = resolve_dynamic::num_points(g, add_for_open);
        }, geometry);
        return result;
    }
    else if constexpr (concepts::ConstGeometryCollection<Geometry>)
    {
        std::size_t result = 0;
        detail::visit_breadth_first([&](auto const& g)
        {
            result += resolve_dynamic::num_points(g, add_for_open);
            return true;
        }, geometry);
        return result;
    }
    else if constexpr (concepts::ConstPoint<Geometry>)
    {
        return 1;
    }
    else if constexpr (concepts::ConstBox<Geometry>)
    {
        return 1 << geometry::dimension<Geometry>::value;
    }
    else if constexpr (concepts::ConstSegment<Geometry>)
    {
        return 2;
    }
    else if constexpr (concepts::ConstLinestring<Geometry>
                    || concepts::ConstRing<Geometry>)
    {
        return add_for_open
             ? detail::num_points::range_count<true>::apply(geometry)
             : detail::num_points::range_count<false>::apply(geometry);
    }
    else if constexpr (concepts::ConstPolygon<Geometry>)
    {
        std::size_t result = resolve_dynamic::num_points(
            exterior_ring(geometry), add_for_open);
        auto const& rings = interior_rings(geometry);
        for (auto it = boost::begin(rings); it != boost::end(rings); ++it)
        {
            result += resolve_dynamic::num_points(*it, add_for_open);
        }
        return result;
    }
    else
    {
        std::size_t result = 0;
        for (auto it = boost::begin(geometry); it != boost::end(geometry); ++it)
        {
            result += resolve_dynamic::num_points(*it, add_for_open);
        }
        return result;
    }
}

} // namespace resolve_dynamic


/*!
\brief \brief_calc{number of points}
\ingroup num_points
\details \details_calc{num_points, number of points}.
\tparam Geometry \tparam_geometry
\param geometry \param_geometry
\param add_for_open add one for open geometries (i.e. polygon types which are not closed)
\return \return_calc{number of points}

\qbk{[include reference/algorithms/num_points.qbk]}
*/
template <concepts::ConstGeometry Geometry>
inline std::size_t num_points(Geometry const& geometry, bool add_for_open = false)
{
    return resolve_dynamic::num_points(geometry, add_for_open);
}

#if defined(_MSC_VER)
#pragma warning(pop)
#endif

}} // namespace boost::geometry


#endif // BOOST_GEOMETRY_ALGORITHMS_NUM_POINTS_HPP
