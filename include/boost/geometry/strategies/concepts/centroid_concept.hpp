// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2007-2012 Barend Gehrels, Amsterdam, the Netherlands.
// Copyright (c) 2008-2012 Bruno Lalande, Paris, France.
// Copyright (c) 2009-2012 Mateusz Loskot, London, UK.

// Parts of Boost.Geometry are redesigned from Geodan's Geographic Library
// (geolib/GGL), copyright (c) 1995-2010 Geodan, Amsterdam, the Netherlands.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_STRATEGIES_CONCEPTS_CENTROID_CONCEPT_HPP
#define BOOST_GEOMETRY_STRATEGIES_CONCEPTS_CENTROID_CONCEPT_HPP



#include <concepts>

namespace boost { namespace geometry { namespace concepts
{


/*!
    \brief Checks strategy for centroid
    \ingroup centroid
*/
template <typename Strategy>
concept CentroidStrategy =
    requires(Strategy& strategy,
             typename Strategy::state_type& state,
             typename Strategy::point_type& point,
             typename Strategy::segment_point_type const& segment_point)
    {
        typename Strategy::state_type;
        typename Strategy::point_type;
        typename Strategy::segment_point_type;
        strategy.apply(segment_point, segment_point, state);
        { strategy.result(state, point) } -> std::convertible_to<bool>;
    };


}}} // namespace boost::geometry::concepts


#endif // BOOST_GEOMETRY_STRATEGIES_CONCEPTS_CENTROID_CONCEPT_HPP
