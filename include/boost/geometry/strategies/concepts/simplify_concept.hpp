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

#ifndef BOOST_GEOMETRY_STRATEGIES_CONCEPTS_SIMPLIFY_CONCEPT_HPP
#define BOOST_GEOMETRY_STRATEGIES_CONCEPTS_SIMPLIFY_CONCEPT_HPP

#include <concepts>
#include <iterator>
#include <vector>

#include <boost/geometry/strategies/concepts/distance_concept.hpp>


namespace boost { namespace geometry { namespace concepts
{


/*!
    \brief Checks strategy for simplify
    \ingroup simplify
*/
template <typename Strategy, typename Point>
concept SimplifyStrategy =
    requires
    {
        typename Strategy::distance_strategy_type;
        requires PointSegmentDistanceStrategy
            <typename Strategy::distance_strategy_type, Point, Point>;
    }
    && requires(Strategy& strategy,
                std::vector<Point> const& input,
                std::vector<Point>& output)
    {
        strategy.apply(input, std::back_inserter(output), 1.0);
    };



}}} // namespace boost::geometry::concepts

#endif // BOOST_GEOMETRY_STRATEGIES_CONCEPTS_SIMPLIFY_CONCEPT_HPP
