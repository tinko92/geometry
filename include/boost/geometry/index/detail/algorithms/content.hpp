// Boost.Geometry Index
//
// n-dimensional content (hypervolume) - 2d area, 3d volume, ...
//
// Copyright (c) 2011-2014 Adam Wulkiewicz, Lodz, Poland.
//
// This file was modified by Oracle on 2020-2023.
// Modifications copyright (c) 2020-2023 Oracle and/or its affiliates.
// Contributed and/or modified by Vissarion Fysikopoulos, on behalf of Oracle
// Contributed and/or modified by Adam Wulkiewicz, on behalf of Oracle
//
// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_INDEX_DETAIL_ALGORITHMS_CONTENT_HPP
#define BOOST_GEOMETRY_INDEX_DETAIL_ALGORITHMS_CONTENT_HPP

#include <limits>
#include <type_traits>

#include <boost/geometry/core/access.hpp>
#include <boost/geometry/core/coordinate_dimension.hpp>
#include <boost/geometry/core/coordinate_type.hpp>
#include <boost/geometry/core/static_assert.hpp>
#include <boost/geometry/core/tag.hpp>
#include <boost/geometry/core/tags.hpp>
#include <boost/geometry/util/select_most_precise.hpp>

namespace boost { namespace geometry { namespace index { namespace detail {

template <typename Indexable>
struct default_content_result
{
    using type = typename select_most_precise
        <
            typename coordinate_type<Indexable>::type,
            double
        >::type;
};

namespace dispatch {

template <typename Box,
          std::size_t CurrentDimension = dimension<Box>::value>
struct content_box
{
    BOOST_STATIC_ASSERT(0 < CurrentDimension);

    static inline typename detail::default_content_result<Box>::type apply(Box const& b)
    {
        return content_box<Box, CurrentDimension - 1>::apply(b) *
            ( get<max_corner, CurrentDimension - 1>(b) - get<min_corner, CurrentDimension - 1>(b) );
    }
};

template <typename Box>
struct content_box<Box, 1>
{
    static inline typename detail::default_content_result<Box>::type apply(Box const& b)
    {
        return get<max_corner, 0>(b) - get<min_corner, 0>(b);
    }
};

template <typename Indexable, typename Tag>
struct content
{
    BOOST_GEOMETRY_STATIC_ASSERT_FALSE(
        "Not implemented for this Indexable and Tag.",
        Indexable, Tag);
};

template <typename Indexable>
struct content<Indexable, point_tag>
{
    static typename detail::default_content_result<Indexable>::type apply(Indexable const&)
    {
        return 0;
    }
};

template <typename Indexable>
struct content<Indexable, box_tag>
{
    static typename default_content_result<Indexable>::type apply(Indexable const& b)
    {
        return dispatch::content_box<Indexable>::apply(b);
    }
};

template
<
    typename Box1,
    typename Box2,
    bool is_float = std::is_floating_point<typename default_content_result<Box1>::type>::value
>
struct content_diff
{
    static auto apply(Box1 const& b1, Box2 const& b2)
    {
        return content_box<Box1>::apply(b1) - content_box<Box2>::apply(b2);
    }

    static auto apply(Box1 const& b1, Box2 const& b2, Box2 const& b3)
    {
        return content_box<Box1>::apply(b1) - content_box<Box2>::apply(b2)
            - content_box<Box2>::apply(b3);
    }
};

template <typename Box1, typename Box2>
struct content_diff<Box1, Box2, true>
{
    using return_type = typename default_content_result<Box1>::type;

    static auto apply(Box1 const& b1, Box2 const& b2)
    {
        auto const content1 = content_box<Box1>::apply(b1);
        auto const content2 = content_box<Box2>::apply(b2);
        auto const diff = content1 - content2;
        // Analogous to page 348 in https://doi.org/10.1007/PL00009321
        auto constexpr box_dim = dimension<Box1>::value;
        auto constexpr static_coeff =
              (box_dim * 2 - 1) * (std::numeric_limits<return_type>::epsilon() / 2.);
        return diff < static_coeff * (content1 + content2) ? return_type{0.} : diff;
    }

    static auto apply(Box1 const& b1, Box2 const& b2, Box2 const& b3)
    {
        auto const content1 = content_box<Box1>::apply(b1);
        auto const content2 = content_box<Box2>::apply(b2);
        auto const content3 = content_box<Box2>::apply(b2);
        auto const diff = content1 - content2 - content3;
        // Analogous to page 348 in https://doi.org/10.1007/PL00009321
        auto constexpr box_dim = dimension<Box1>::value;
        auto constexpr static_coeff =
              (box_dim * 2) * (std::numeric_limits<return_type>::epsilon() / 2.);
        return diff < static_coeff * (content1 + content2 + content3) ? return_type{0.} : diff;
    }
};

} // namespace dispatch

template <typename Indexable>
typename default_content_result<Indexable>::type content(Indexable const& b)
{
    return dispatch::content
            <
                Indexable,
                tag_t<Indexable>
            >::apply(b);
}

template <typename Box1, typename Box2>
auto content_diff(Box1 const& b1, Box2 const& b2)
{
    return dispatch::content_diff<Box1, Box2>::apply(b1, b2);
}

template <typename Box1, typename Box2>
auto content_diff(Box1 const& b1, Box2 const& b2, Box2 const& b3)
{
    return dispatch::content_diff<Box1, Box2>::apply(b1, b2, b3);
}

}}}} // namespace boost::geometry::index::detail

#endif // BOOST_GEOMETRY_INDEX_DETAIL_ALGORITHMS_CONTENT_HPP
