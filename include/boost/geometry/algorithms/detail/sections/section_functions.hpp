// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2015 Barend Gehrels, Amsterdam, the Netherlands.

// This file was modified by Oracle on 2015-2021.
// Modifications copyright (c) 2015-2021, Oracle and/or its affiliates.
// Contributed and/or modified by Adam Wulkiewicz, on behalf of Oracle

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_ALGORITHMS_DETAIL_SECTIONS_FUNCTIONS_HPP
#define BOOST_GEOMETRY_ALGORITHMS_DETAIL_SECTIONS_FUNCTIONS_HPP


namespace boost { namespace geometry
{

#ifndef DOXYGEN_NO_DETAIL
namespace detail { namespace section
{

template
<
    std::size_t Dimension,
    typename Point,
    typename Box,
    typename Strategy
>
inline bool preceding(int dir,
                      Point const& point,
                      Box const& point_box,
                      Box const& other_box,
                      Strategy const& strategy)
{
    return strategy.preceding().template apply<Dimension>(dir, point,
                                                          point_box,
                                                          other_box);
}

template
<
    std::size_t Dimension,
    typename Point,
    typename Box,
    typename Strategy
>
inline bool exceeding(int dir,
                      Point const& point,
                      Box const& point_box,
                      Box const& other_box,
                      Strategy const& strategy)
{
    return preceding<Dimension>(-dir, point, point_box, other_box, strategy);
}


}} // namespace detail::section
#endif


}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_ALGORITHMS_DETAIL_SECTIONS_FUNCTIONS_HPP
