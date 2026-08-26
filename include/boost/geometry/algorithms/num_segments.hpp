// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2014-2023, Oracle and/or its affiliates.

// Contributed and/or modified by Vissarion Fysikopoulos, on behalf of Oracle
// Contributed and/or modified by Menelaos Karavelas, on behalf of Oracle
// Contributed and/or modified by Adam Wulkiewicz, on behalf of Oracle

// Licensed under the Boost Software License version 1.0.
// http://www.boost.org/users/license.html

#ifndef BOOST_GEOMETRY_ALGORITHMS_NUM_SEGMENTS_HPP
#define BOOST_GEOMETRY_ALGORITHMS_NUM_SEGMENTS_HPP

#include <cstddef>

#include <boost/range/size.hpp>
#include <boost/range/value_type.hpp>

#include <boost/geometry/algorithms/detail/counting.hpp>
#include <boost/geometry/algorithms/detail/visit.hpp>
#include <boost/geometry/algorithms/not_implemented.hpp>

#include <boost/geometry/core/closure.hpp>
#include <boost/geometry/core/tag.hpp>
#include <boost/geometry/core/tags.hpp>
#include <boost/geometry/core/visit.hpp>

#include <boost/geometry/geometries/adapted/boost_variant.hpp> // For backward compatibility
#include <boost/geometry/geometries/concepts/check.hpp>

namespace boost { namespace geometry
{

#ifndef DOXYGEN_NO_DETAIL
namespace detail { namespace num_segments
{


struct range_count
{
    template <typename Range>
    static inline std::size_t apply(Range const& range)
    {
        std::size_t n = boost::size(range);
        if ( n <= 1 )
        {
            return 0;
        }

        return
            geometry::closure<Range>::value == open
            ?
            n
            :
            static_cast<std::size_t>(n - 1);
    }
};

}} // namespace detail::num_segments
#endif // DOXYGEN_NO_DETAIL



namespace resolve_dynamic
{

template <concepts::ConstGeometry Geometry>
inline std::size_t num_segments(Geometry const& geometry)
{
    if constexpr (concepts::ConstDynamicGeometry<Geometry>)
    {
        std::size_t result = 0;
        traits::visit<Geometry>::apply([&](auto const& g)
        {
            result = resolve_dynamic::num_segments(g);
        }, geometry);
        return result;
    }
    else if constexpr (concepts::ConstGeometryCollection<Geometry>)
    {
        std::size_t result = 0;
        detail::visit_breadth_first([&](auto const& g)
        {
            result += resolve_dynamic::num_segments(g);
            return true;
        }, geometry);
        return result;
    }
    else if constexpr (concepts::ConstPoint<Geometry>
                    || concepts::ConstMultiPoint<Geometry>)
    {
        return 0;
    }
    else if constexpr (concepts::ConstBox<Geometry>)
    {
        constexpr auto dimensions = geometry::dimension<Geometry>::value;
        return dimensions * (1 << (dimensions - 1));
    }
    else if constexpr (concepts::ConstSegment<Geometry>)
    {
        return 1;
    }
    else if constexpr (concepts::ConstLinestring<Geometry>
                    || concepts::ConstRing<Geometry>)
    {
        return detail::num_segments::range_count::apply(geometry);
    }
    else if constexpr (concepts::ConstPolygon<Geometry>)
    {
        std::size_t result = resolve_dynamic::num_segments(
            exterior_ring(geometry));
        auto const& rings = interior_rings(geometry);
        for (auto it = boost::begin(rings); it != boost::end(rings); ++it)
        {
            result += resolve_dynamic::num_segments(*it);
        }
        return result;
    }
    else
    {
        std::size_t result = 0;
        for (auto it = boost::begin(geometry); it != boost::end(geometry); ++it)
        {
            result += resolve_dynamic::num_segments(*it);
        }
        return result;
    }
}

} // namespace resolve_dynamic



/*!
\brief \brief_calc{number of segments}
\ingroup num_segments
\details \details_calc{num_segments, number of segments}.
\tparam Geometry \tparam_geometry
\param geometry \param_geometry
\return \return_calc{number of segments}

\qbk{[include reference/algorithms/num_segments.qbk]}
*/
template <concepts::ConstGeometry Geometry>
inline std::size_t num_segments(Geometry const& geometry)
{
    return resolve_dynamic::num_segments(geometry);
}



}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_ALGORITHMS_NUM_SEGMENTS_HPP
