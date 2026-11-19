// Boost.Geometry (aka GGL, Generic Geometry Library)
// Copyright (c) 2007-2015 Barend Gehrels, Amsterdam, the Netherlands.
// Modifications copyright (c) 2015-2024 Oracle and/or its affiliates.
// Contributed and/or modified by Vissarion Fysikopoulos, on behalf of Oracle
// Contributed and/or modified by Adam Wulkiewicz, on behalf of Oracle
// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_STRATEGIES_DETAIL_BOX_INTERSECTION_HPP
#define BOOST_GEOMETRY_STRATEGIES_DETAIL_BOX_INTERSECTION_HPP

#include <cstddef>
#include <type_traits>

#include <boost/geometry/core/access.hpp>
#include <boost/geometry/core/coordinate_dimension.hpp>

namespace boost { namespace geometry { namespace strategies { namespace detail
{

namespace box_intersection_detail
{

template <std::size_t Dimension, std::size_t DimensionCount>
struct intersection_loop
{
    template <typename Box1, typename Box2, typename BoxOut>
    static inline bool apply(Box1 const& box1, Box2 const& box2, BoxOut& box_out)
    {
        auto max1 = get<max_corner, Dimension>(box1);
        auto min2 = get<min_corner, Dimension>(box2);
        if (max1 < min2)
        {
            return false;
        }
        auto max2 = get<max_corner, Dimension>(box2);
        auto min1 = get<min_corner, Dimension>(box1);
        if (max2 < min1)
        {
            return false;
        }

        // Preserve earlier dimensions in the output if a later one is disjoint.
        set<min_corner, Dimension>(box_out, min1 < min2 ? min2 : min1);
        set<max_corner, Dimension>(box_out, max1 > max2 ? max2 : max1);
        return intersection_loop<Dimension + 1, DimensionCount>::apply(box1, box2, box_out);
    }
};

template <std::size_t DimensionCount>
struct intersection_loop<DimensionCount, DimensionCount>
{
    template <typename Box1, typename Box2, typename BoxOut>
    static inline bool apply(Box1 const&, Box2 const&, BoxOut&) { return true; }
};

} // namespace box_intersection_detail

// No arithmetic constructions are required: output bounds are selected from
// input coordinates. Exact-point traits consequently retain their exact values.
struct legacy_box_intersection
{
    template <typename Box1, typename Box2, typename BoxOut>
    static inline bool apply(Box1 const& box1, Box2 const& box2, BoxOut& box_out)
    {
        return box_intersection_detail::intersection_loop<0, geometry::dimension<Box1>::value>
            ::apply(box1, box2, box_out);
    }
};

namespace box_intersection_detail
{

template <typename Strategy>
inline auto get_box_intersection(Strategy const& strategy, int)
    -> decltype(strategy.box_intersection())
{
    return strategy.box_intersection();
}

template <typename Strategy>
inline legacy_box_intersection get_box_intersection(Strategy const&, long) { return {}; }

} // namespace box_intersection_detail

// An override's apply(box1, box2, out) must preserve dimension-ordered partial
// writes on failure. It may delegate to legacy_box_intersection without recursion.
template <typename Strategy>
inline decltype(auto) get_box_intersection(Strategy const& strategy)
{
    return box_intersection_detail::get_box_intersection(strategy, 0);
}

template <std::size_t Dimension, std::size_t DimensionCount>
struct intersection_box_box
{
private:
    template <typename Box1, typename Box2, typename BoxOut, typename Strategy>
    static inline bool apply_impl(Box1 const& box1, Box2 const& box2, BoxOut& box_out,
                                  Strategy const& strategy, std::true_type)
    {
        return get_box_intersection(strategy).apply(box1, box2, box_out);
    }

    template <typename Box1, typename Box2, typename BoxOut, typename Strategy>
    static inline bool apply_impl(Box1 const& box1, Box2 const& box2, BoxOut& box_out,
                                  Strategy const&, std::false_type)
    {
        return box_intersection_detail::intersection_loop<Dimension, DimensionCount>
            ::apply(box1, box2, box_out);
    }

public:
    template <typename Box1, typename Box2, typename BoxOut, typename Strategy>
    static inline bool apply(Box1 const& box1, Box2 const& box2, BoxOut& box_out,
                             Strategy const& strategy)
    {
        // Prefix/suffix compatibility calls must not invoke a full-box override.
        using full_box = std::integral_constant<bool,
            Dimension == 0 && DimensionCount == geometry::dimension<Box1>::value>;
        return apply_impl(box1, box2, box_out, strategy, full_box());
    }
};

template <std::size_t DimensionCount>
struct intersection_box_box<DimensionCount, DimensionCount>
{
    template <typename Box1, typename Box2, typename BoxOut, typename Strategy>
    static inline bool apply(Box1 const&, Box2 const&, BoxOut&, Strategy const&) { return true; }
};

}}}} // namespace boost::geometry::strategies::detail

#endif // BOOST_GEOMETRY_STRATEGIES_DETAIL_BOX_INTERSECTION_HPP
