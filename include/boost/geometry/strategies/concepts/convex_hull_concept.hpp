// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2007-2012 Barend Gehrels, Amsterdam, the Netherlands.
// Copyright (c) 2008-2012 Bruno Lalande, Paris, France.
// Copyright (c) 2009-2012 Mateusz Loskot, London, UK.

// This file was modified by Oracle on 2014.
// Modifications copyright (c) 2014 Oracle and/or its affiliates.

// Contributed and/or modified by Adam Wulkiewicz, on behalf of Oracle

// Parts of Boost.Geometry are redesigned from Geodan's Geographic Library
// (geolib/GGL), copyright (c) 1995-2010 Geodan, Amsterdam, the Netherlands.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_STRATEGIES_CONCEPTS_CONVEX_HULL_CONCEPT_HPP
#define BOOST_GEOMETRY_STRATEGIES_CONCEPTS_CONVEX_HULL_CONCEPT_HPP


#include <vector>

#include <concepts>


namespace boost { namespace geometry { namespace concepts
{


/*!
    \brief Checks strategy for convex_hull
    \ingroup convex_hull
*/
template <typename Strategy>
concept ConvexHullStrategy =
    requires(Strategy const& strategy,
             typename Strategy::state_type& state,
             typename Strategy::geometry_type const& geometry,
             std::vector<typename Strategy::point_type>& points)
    {
        typename Strategy::state_type;
        typename Strategy::point_type;
        typename Strategy::geometry_type;
        strategy.apply(geometry, state);
        strategy.result(state, std::back_inserter(points), true, true);
    };


}}} // namespace boost::geometry::concepts


#endif // BOOST_GEOMETRY_STRATEGIES_CONCEPTS_CONVEX_HULL_CONCEPT_HPP
