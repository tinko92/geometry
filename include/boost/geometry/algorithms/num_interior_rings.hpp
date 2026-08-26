// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2007-2014 Barend Gehrels, Amsterdam, the Netherlands.
// Copyright (c) 2008-2014 Bruno Lalande, Paris, France.
// Copyright (c) 2009-2014 Mateusz Loskot, London, UK.

// This file was modified by Oracle on 2014-2020.
// Modifications copyright (c) 2014-2020, Oracle and/or its affiliates.
// Contributed and/or modified by Adam Wulkiewicz, on behalf of Oracle
// Contributed and/or modified by Menelaos Karavelas, on behalf of Oracle

// Parts of Boost.Geometry are redesigned from Geodan's Geographic Library
// (geolib/GGL), copyright (c) 1995-2010 Geodan, Amsterdam, the Netherlands.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_ALGORITHMS_NUM_INTERIOR_RINGS_HPP
#define BOOST_GEOMETRY_ALGORITHMS_NUM_INTERIOR_RINGS_HPP

#include <cstddef>

#include <boost/range/size.hpp>
#include <boost/range/value_type.hpp>

#include <boost/variant/apply_visitor.hpp>
#include <boost/variant/static_visitor.hpp>
#include <boost/variant/variant_fwd.hpp>

#include <boost/geometry/core/tag.hpp>
#include <boost/geometry/core/tags.hpp>
#include <boost/geometry/core/visit.hpp>

#include <boost/geometry/core/interior_rings.hpp>

#include <boost/geometry/algorithms/detail/counting.hpp>
#include <boost/geometry/algorithms/detail/visit.hpp>

#include <boost/geometry/geometries/concepts/check.hpp>


namespace boost { namespace geometry
{

/*!
\brief \brief_calc{number of interior rings}
\ingroup num_interior_rings
\details \details_calc{num_interior_rings, number of interior rings}.
\tparam Geometry \tparam_geometry
\param geometry \param_geometry
\return \return_calc{number of interior rings}

\qbk{[include reference/algorithms/num_interior_rings.qbk]}

\note Defined by OGC as "numInteriorRing". To be consistent with "numPoints"
    letter "s" is appended
*/
template <concepts::ConstGeometry Geometry>
inline std::size_t num_interior_rings(Geometry const& geometry)
{
    if constexpr (concepts::ConstDynamicGeometry<Geometry>)
    {
        std::size_t result = 0;
        traits::visit<Geometry>::apply([&](auto const& g)
        {
            result = geometry::num_interior_rings(g);
        }, geometry);
        return result;
    }
    else if constexpr (concepts::ConstGeometryCollection<Geometry>)
    {
        std::size_t result = 0;
        detail::visit_breadth_first([&](auto const& g)
        {
            result += geometry::num_interior_rings(g);
            return true;
        }, geometry);
        return result;
    }
    else if constexpr (concepts::ConstPolygon<Geometry>)
    {
        return boost::size(geometry::interior_rings(geometry));
    }
    else if constexpr (concepts::ConstMultiPolygon<Geometry>)
    {
        std::size_t result = 0;
        for (auto it = boost::begin(geometry); it != boost::end(geometry); ++it)
        {
            result += geometry::num_interior_rings(*it);
        }
        return result;
    }
    else
    {
        return 0;
    }
}


}} // namespace boost::geometry


#endif // BOOST_GEOMETRY_ALGORITHMS_NUM_INTERIOR_RINGS_HPP
