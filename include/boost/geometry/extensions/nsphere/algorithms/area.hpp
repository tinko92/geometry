// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2007-2012 Barend Gehrels, Amsterdam, the Netherlands.
// Copyright (c) 2008-2012 Bruno Lalande, Paris, France.
// Copyright (c) 2009-2012 Mateusz Loskot, London, UK.

// This file was modified by Oracle on 2020.
// Modifications copyright (c) 2020, Oracle and/or its affiliates.
// Contributed and/or modified by Adam Wulkiewicz, on behalf of Oracle

// Parts of Boost.Geometry are redesigned from Geodan's Geographic Library
// (geolib/GGL), copyright (c) 1995-2010 Geodan, Amsterdam, the Netherlands.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_EXTENSIONS_NSPHERE_ALGORITHMS_AREA_HPP
#define BOOST_GEOMETRY_EXTENSIONS_NSPHERE_ALGORITHMS_AREA_HPP


#include <type_traits>

#include <boost/math/constants/constants.hpp>

#include <boost/geometry/algorithms/area.hpp>
#include <boost/geometry/extensions/nsphere/core/radius.hpp>
#include <boost/geometry/extensions/nsphere/core/tags.hpp>
#include <boost/geometry/extensions/nsphere/geometries/concepts/nsphere_concept.hpp>



namespace boost { namespace geometry
{

#ifndef DOXYGEN_NO_DETAIL
namespace detail
{

template <concepts::ConstNsphere NSphere>
inline auto area_nsphere(NSphere const& nsphere)
{
    using coordinate_type = coordinate_type_t<NSphere>;

    using return_type = std::conditional_t
        <
            std::is_integral_v<coordinate_type>,
            double,
            coordinate_type
        >;

    assert_dimension<NSphere, 2>();

    return_type radius = get_radius<0>(nsphere);
    radius *= radius * boost::math::constants::pi<return_type>();
    return radius;
}

} // namespace detail

#endif // DOXYGEN_NO_DETAIL

template <concepts::ConstGeometry NSphere>
    requires concepts::ConstNsphere<NSphere>
inline auto area(NSphere const& nsphere)
{
    return detail::area_nsphere(nsphere);
}

template <concepts::ConstGeometry NSphere, typename Strategy>
    requires concepts::ConstNsphere<NSphere>
inline auto area(NSphere const& nsphere, Strategy const&)
{
    return detail::area_nsphere(nsphere);
}


}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_EXTENSIONS_NSPHERE_ALGORITHMS_AREA_HPP
