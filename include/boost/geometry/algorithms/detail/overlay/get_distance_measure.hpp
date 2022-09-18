// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2019-2021 Barend Gehrels, Amsterdam, the Netherlands.

// This file was modified by Oracle on 2022.
// Modifications copyright (c) 2022 Oracle and/or its affiliates.
// Contributed and/or modified by Adam Wulkiewicz, on behalf of Oracle

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_ALGORITHMS_DETAIL_OVERLAY_GET_DISTANCE_MEASURE_HPP
#define BOOST_GEOMETRY_ALGORITHMS_DETAIL_OVERLAY_GET_DISTANCE_MEASURE_HPP

#include <boost/geometry/core/coordinate_type.hpp>
#include <boost/geometry/util/select_coordinate_type.hpp>

namespace boost { namespace geometry
{

#ifndef DOXYGEN_NO_DETAIL
namespace detail
{

template <typename T>
struct distance_measure
{
    T measure;

    distance_measure()
        : measure(T())
    {}

    // Returns true if the distance measure is absolutely zero
    bool is_zero() const
    {
      return ! is_positive() && ! is_negative();
    }

    // Returns true if the distance measure is positive. Distance measure
    // algorithm returns positive value if it is located on the left side.
    bool is_positive() const { return measure > T(0); }

    // Returns true if the distance measure is negative. Distance measure
    // algorithm returns negative value if it is located on the right side.
    bool is_negative() const { return measure < T(0); }
};

} // detail

namespace detail
{

// Returns a (often very tiny) value to indicate its side, and distance,
// 0 (absolutely 0, not even an epsilon) means collinear. Like side,
// a negative means that p is to the right of p1-p2. And a positive value
// means that p is to the left of p1-p2.

template <typename SegmentPoint, typename Point, typename Strategies>
inline auto get_distance_measure(SegmentPoint const& p1, SegmentPoint const& p2, Point const& p,
                                 Strategies const& strategy)
{
    using calculation_type = typename select_coordinate_type<SegmentPoint, Point>::type;
    detail::distance_measure<calculation_type> result;
    result.measure = strategy.side_value_line().apply(p1, p2, p);
    return result;
}

} // namespace detail
#endif // DOXYGEN_NO_DETAIL

}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_ALGORITHMS_DETAIL_OVERLAY_GET_DISTANCE_MEASURE_HPP
