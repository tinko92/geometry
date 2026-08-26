// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2014-2021, Oracle and/or its affiliates.

// Contributed and/or modified by Menelaos Karavelas, on behalf of Oracle
// Contributed and/or modified by Adam Wulkiewicz, on behalf of Oracle

// Licensed under the Boost Software License version 1.0.
// http://www.boost.org/users/license.html

#ifndef BOOST_GEOMETRY_ALGORITHMS_DETAIL_IS_SIMPLE_AREAL_HPP
#define BOOST_GEOMETRY_ALGORITHMS_DETAIL_IS_SIMPLE_AREAL_HPP

#include <boost/range/begin.hpp>
#include <boost/range/empty.hpp>
#include <boost/range/end.hpp>
#include <boost/range/value_type.hpp>

#include <boost/geometry/core/closure.hpp>
#include <boost/geometry/core/exterior_ring.hpp>
#include <boost/geometry/core/interior_rings.hpp>
#include <boost/geometry/core/ring_type.hpp>
#include <boost/geometry/core/tags.hpp>

#include <boost/geometry/algorithms/detail/is_simple/failure_policy.hpp>
#include <boost/geometry/algorithms/detail/is_valid/has_duplicates.hpp>

#include <boost/geometry/geometries/concepts/check.hpp>


namespace boost { namespace geometry
{


#ifndef DOXYGEN_NO_DETAIL
namespace detail { namespace is_simple
{


template <typename Ring, typename Strategy>
inline bool is_simple_ring(Ring const& ring, Strategy const& strategy)
{
    simplicity_failure_policy policy;
    return ! boost::empty(ring)
        && ! detail::is_valid::has_duplicates<Ring>::apply(ring, policy, strategy);
}

template <typename InteriorRings, typename Strategy>
inline bool are_simple_interior_rings(InteriorRings const& interior_rings,
                                      Strategy const& strategy)
{
    return std::all_of(boost::begin(interior_rings),
                       boost::end(interior_rings),
                       [&](auto const& r)
                       {
                           return is_simple_ring(r, strategy);
                       }); // non-simple ring not found
    // allow empty ring
}

template <typename Polygon, typename Strategy>
inline bool is_simple_polygon(Polygon const& polygon, Strategy const& strategy)
{
    return is_simple_ring(geometry::exterior_ring(polygon), strategy)
        && are_simple_interior_rings(geometry::interior_rings(polygon), strategy);
}

template <concepts::ConstRing Ring, typename Strategy>
inline bool apply(Ring const& ring, Strategy const& strategy)
{
    return is_simple_ring(ring, strategy);
}

template <concepts::ConstPolygon Polygon, typename Strategy>
inline bool apply(Polygon const& polygon, Strategy const& strategy)
{
    return is_simple_polygon(polygon, strategy);
}

template <concepts::ConstMultiPolygon MultiPolygon, typename Strategy>
inline bool apply(MultiPolygon const& multipolygon, Strategy const& strategy)
{
    return std::none_of(boost::begin(multipolygon), boost::end(multipolygon),
        [&](auto const& polygon)
        {
            return ! is_simple_polygon(polygon, strategy);
        });
}


}} // namespace detail::is_simple
#endif // DOXYGEN_NO_DETAIL

}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_ALGORITHMS_DETAIL_IS_SIMPLE_AREAL_HPP
