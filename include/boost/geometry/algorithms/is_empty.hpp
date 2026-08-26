// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2015-2023, Oracle and/or its affiliates.

// Contributed and/or modified by Vissarion Fysikopoulos, on behalf of Oracle
// Contributed and/or modified by Menelaos Karavelas, on behalf of Oracle
// Contributed and/or modified by Adam Wulkiewicz, on behalf of Oracle

// Licensed under the Boost Software License version 1.0.
// http://www.boost.org/users/license.html

#ifndef BOOST_GEOMETRY_ALGORITHMS_IS_EMPTY_HPP
#define BOOST_GEOMETRY_ALGORITHMS_IS_EMPTY_HPP

#include <boost/range/begin.hpp>
#include <boost/range/empty.hpp>
#include <boost/range/end.hpp>

#include <boost/geometry/algorithms/not_implemented.hpp>
#include <boost/geometry/algorithms/detail/visit.hpp>

#include <boost/geometry/core/exterior_ring.hpp>
#include <boost/geometry/core/interior_rings.hpp>
#include <boost/geometry/core/tag.hpp>
#include <boost/geometry/core/tags.hpp>
#include <boost/geometry/core/visit.hpp>

#include <boost/geometry/geometries/adapted/std_variant.hpp>
#include <boost/geometry/geometries/concepts/check.hpp>

#include <boost/geometry/util/type_traits_std.hpp>

namespace boost { namespace geometry
{


#ifndef DOXYGEN_NO_DETAIL
namespace detail { namespace is_empty
{

struct always_not_empty
{
    template <typename Geometry>
    static inline bool apply(Geometry const&)
    {
        return false;
    }
};

struct range_is_empty
{
    template <typename Range>
    static inline bool apply(Range const& range)
    {
        return boost::empty(range);
    }
};

class polygon_is_empty
{
    template <typename InteriorRings>
    static inline bool check_interior_rings(InteriorRings const& interior_rings)
    {
        return std::all_of(boost::begin(interior_rings), boost::end(interior_rings),
                           []( auto const& range ){ return boost::empty(range); });
    }

public:
    template <typename Polygon>
    static inline bool apply(Polygon const& polygon)
    {
        return boost::empty(exterior_ring(polygon))
            && check_interior_rings(interior_rings(polygon));
    }
};

template <typename Policy = range_is_empty>
struct multi_is_empty
{
    template <typename MultiGeometry>
    static inline bool apply(MultiGeometry const& multigeometry)
    {
        return std::all_of(boost::begin(multigeometry),
                           boost::end(multigeometry),
                           []( auto const& range ){ return Policy::apply(range); });
    }
};

}} // namespace detail::is_empty
#endif // DOXYGEN_NO_DETAIL


namespace resolve_dynamic
{

template <concepts::ConstGeometry Geometry>
inline bool is_empty(Geometry const& geometry)
{
    if constexpr (concepts::ConstDynamicGeometry<Geometry>)
    {
        bool result = true;
        traits::visit<Geometry>::apply([&](auto const& g)
        {
            result = resolve_dynamic::is_empty(g);
        }, geometry);
        return result;
    }
    else if constexpr (concepts::ConstGeometryCollection<Geometry>)
    {
        bool result = true;
        detail::visit_breadth_first([&](auto const& g)
        {
            result = resolve_dynamic::is_empty(g);
            return result;
        }, geometry);
        return result;
    }
    else if constexpr (concepts::ConstPoint<Geometry>
                    || concepts::ConstBox<Geometry>
                    || concepts::ConstSegment<Geometry>)
    {
        return false;
    }
    else if constexpr (concepts::ConstLinestring<Geometry>
                    || concepts::ConstRing<Geometry>
                    || concepts::ConstMultiPoint<Geometry>)
    {
        return boost::empty(geometry);
    }
    else if constexpr (concepts::ConstPolygon<Geometry>)
    {
        return detail::is_empty::polygon_is_empty::apply(geometry);
    }
    else if constexpr (concepts::ConstMultiPolygon<Geometry>)
    {
        return detail::is_empty::multi_is_empty
            <detail::is_empty::polygon_is_empty>::apply(geometry);
    }
    else
    {
        return detail::is_empty::multi_is_empty<>::apply(geometry);
    }
}

} // namespace resolve_dynamic


/*!
\brief \brief_check{is the empty set}
\ingroup is_empty
\tparam Geometry \tparam_geometry
\param geometry \param_geometry
\return \return_check{is the empty set}

\qbk{[include reference/algorithms/is_empty.qbk]}
*/
template <concepts::ConstGeometry Geometry>
inline bool is_empty(Geometry const& geometry)
{
    return resolve_dynamic::is_empty(geometry);
}


}} // namespace boost::geometry


#endif // BOOST_GEOMETRY_ALGORITHMS_IS_EMPTY_HPP
