// Boost.Geometry

// Copyright (c) 2021, Oracle and/or its affiliates.

// Contributed and/or modified by Adam Wulkiewicz, on behalf of Oracle

// Licensed under the Boost Software License version 1.0.
// http://www.boost.org/users/license.html

#ifndef BOOST_GEOMETRY_STRATEGIES_DISTANCE_DETAIL_HPP
#define BOOST_GEOMETRY_STRATEGIES_DISTANCE_DETAIL_HPP


#include <boost/geometry/util/type_traits.hpp>


namespace boost { namespace geometry
{

namespace strategies { namespace distance
{

#ifndef DOXYGEN_NO_DETAIL
namespace detail
{

template <typename Geometry1, typename Geometry2>
concept point_point = util::pointlike<Geometry1> && util::pointlike<Geometry2>;

template <typename Geometry1, typename Geometry2>
concept point_segment = (util::pointlike<Geometry1> && util::segmental<Geometry2>)
                     || (util::segmental<Geometry1> && util::pointlike<Geometry2>)
                     || (util::segmental<Geometry1> && util::segmental<Geometry2>);

template <typename Geometry1, typename Geometry2>
concept point_box = util::pointlike<Geometry1> && util::box<Geometry2>;

template <typename Geometry1, typename Geometry2>
concept segment_box = util::segmental<Geometry1> && util::box<Geometry2>;

template <typename Geometry1, typename Geometry2>
concept box_box = util::box<Geometry1> && util::box<Geometry2>;

template <typename Geometry1, typename Geometry2>
concept geometry_pair = point_point<Geometry1, Geometry2>
                     || point_segment<Geometry1, Geometry2>
                     || point_box<Geometry1, Geometry2>
                     || segment_box<Geometry1, Geometry2>
                     || box_box<Geometry1, Geometry2>;

template <typename Geometry1, typename Geometry2>
inline constexpr bool is_pp_v = point_point<Geometry1, Geometry2>;
template <typename Geometry1, typename Geometry2>
inline constexpr bool is_ps_v = point_segment<Geometry1, Geometry2>;
template <typename Geometry1, typename Geometry2>
inline constexpr bool is_pb_v = point_box<Geometry1, Geometry2>;
template <typename Geometry1, typename Geometry2>
inline constexpr bool is_sb_v = segment_box<Geometry1, Geometry2>;
template <typename Geometry1, typename Geometry2>
inline constexpr bool is_bb_v = box_box<Geometry1, Geometry2>;

} // namespace detail
#endif // DOXYGEN_NO_DETAIL

}} // namespace strategies::distance

}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_STRATEGIES_DISTANCE_DETAIL_HPP
