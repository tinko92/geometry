// Boost.Geometry

// Copyright (c) 2023 Adam Wulkiewicz, Lodz, Poland.

// Copyright (c) 2017-2023, Oracle and/or its affiliates.
// Contributed and/or modified by Vissarion Fysikopoulos, on behalf of Oracle
// Contributed and/or modified by Adam Wulkiewicz, on behalf of Oracle

// Licensed under the Boost Software License version 1.0.
// http://www.boost.org/users/license.html

#ifndef BOOST_GEOMETRY_ALGORITHMS_DENSIFY_HPP
#define BOOST_GEOMETRY_ALGORITHMS_DENSIFY_HPP


#include <boost/range/size.hpp>
#include <boost/range/value_type.hpp>
#include <boost/throw_exception.hpp>

#include <boost/geometry/algorithms/clear.hpp>
#include <boost/geometry/algorithms/convert.hpp>
#include <boost/geometry/algorithms/detail/convert_point_to_point.hpp>
#include <boost/geometry/algorithms/detail/visit.hpp>
#include <boost/geometry/algorithms/not_implemented.hpp>
#include <boost/geometry/core/closure.hpp>
#include <boost/geometry/core/exception.hpp>
#include <boost/geometry/core/tag.hpp>
#include <boost/geometry/core/tags.hpp>
#include <boost/geometry/core/visit.hpp>
#include <boost/geometry/geometries/adapted/boost_variant.hpp> // For backward compatibility
#include <boost/geometry/strategies/default_strategy.hpp>
#include <boost/geometry/strategies/densify/cartesian.hpp>
#include <boost/geometry/strategies/densify/geographic.hpp>
#include <boost/geometry/strategies/densify/spherical.hpp>
#include <boost/geometry/strategies/detail.hpp>
#include <boost/geometry/util/range.hpp>


namespace boost { namespace geometry
{


#ifndef DOXYGEN_NO_DETAIL
namespace detail { namespace densify
{

template <typename Range>
struct push_back_policy
{
    typedef typename boost::range_value<Range>::type point_type;

    inline explicit push_back_policy(Range & rng)
        : m_rng(rng)
    {}

    inline void apply(point_type const& p)
    {
        range::push_back(m_rng, p);
    }

private:
    Range & m_rng;
};

template <typename Range, typename Point>
inline void convert_and_push_back(Range & range, Point const& p)
{
    typename boost::range_value<Range>::type p2;
    geometry::detail::conversion::convert_point_to_point(p, p2);
    range::push_back(range, p2);
}

template <bool AppendLastPoint = true>
struct densify_range
{
    template <typename FwdRng, typename MutRng, typename T, typename Strategies>
    static inline void apply(FwdRng const& rng, MutRng & rng_out,
                             T const& len, Strategies const& strategies)
    {
        typedef typename boost::range_value<FwdRng>::type point_t;

        auto it = boost::begin(rng);
        auto const end = boost::end(rng);

        if (it == end) // empty(rng)
        {
            return;
        }

        auto strategy = strategies.densify(rng);
        push_back_policy<MutRng> policy(rng_out);

        auto prev = it;
        for ( ++it ; it != end ; prev = it++)
        {
            point_t const& p0 = *prev;
            point_t const& p1 = *it;

            convert_and_push_back(rng_out, p0);

            strategy.apply(p0, p1, policy, len);
        }

        if constexpr (AppendLastPoint)
        {
            convert_and_push_back(rng_out, *prev); // back(rng)
        }
    }
};

template <bool IsClosed1, bool IsClosed2> // false, X
struct densify_ring
{
    template <typename Geometry, typename GeometryOut, typename T, typename Strategies>
    static inline void apply(Geometry const& ring, GeometryOut & ring_out,
                             T const& len, Strategies const& strategies)
    {
        geometry::detail::densify::densify_range<true>
            ::apply(ring, ring_out, len, strategies);

        if (boost::size(ring) <= 1)
            return;

        auto const& p0 = range::back(ring);
        auto const& p1 = range::front(ring);

        auto strategy = strategies.densify(ring);
        push_back_policy<GeometryOut> policy(ring_out);

        strategy.apply(p0, p1, policy, len);

        if constexpr (IsClosed2)
        {
            convert_and_push_back(ring_out, p1);
        }
    }
};

template <>
struct densify_ring<true, true>
    : densify_range<true>
{};

template <>
struct densify_ring<true, false>
    : densify_range<false>
{};

struct densify_convert
{
    template <typename GeometryIn, typename GeometryOut, typename T, typename Strategy>
    static void apply(GeometryIn const& in, GeometryOut &out,
                      T const& , Strategy const& )
    {
        geometry::convert(in, out);
    }
};

}} // namespace detail::densify
#endif // DOXYGEN_NO_DETAIL


#ifndef DOXYGEN_NO_DISPATCH
namespace dispatch
{


template <concepts::ConstGeometry Geometry,
          concepts::MutableGeometry GeometryOut,
          typename T,
          typename Strategies>
    requires (concepts::ConstPoint<Geometry> && concepts::Point<GeometryOut>)
          || (concepts::ConstSegment<Geometry> && concepts::Segment<GeometryOut>)
          || (concepts::ConstBox<Geometry> && concepts::Box<GeometryOut>)
          || (concepts::ConstMultiPoint<Geometry> && concepts::MultiPoint<GeometryOut>)
          || (concepts::ConstLinestring<Geometry> && concepts::Linestring<GeometryOut>)
          || (concepts::ConstMultiLinestring<Geometry> && concepts::MultiLinestring<GeometryOut>)
          || (concepts::ConstRing<Geometry> && concepts::Ring<GeometryOut>)
          || (concepts::ConstPolygon<Geometry> && concepts::Polygon<GeometryOut>)
          || (concepts::ConstMultiPolygon<Geometry> && concepts::MultiPolygon<GeometryOut>)
inline void densify(Geometry const& geometry, GeometryOut& out,
                    T const& len, Strategies const& strategies)
{
    if constexpr (concepts::ConstPoint<Geometry>
                  || concepts::ConstSegment<Geometry>
                  || concepts::ConstBox<Geometry>
                  || concepts::ConstMultiPoint<Geometry>)
    {
        geometry::detail::densify::densify_convert::apply(
            geometry, out, len, strategies);
    }
    else if constexpr (concepts::ConstLinestring<Geometry>)
    {
        geometry::detail::densify::densify_range<>::apply(
            geometry, out, len, strategies);
    }
    else if constexpr (concepts::ConstMultiLinestring<Geometry>)
    {
        std::size_t const count = boost::size(geometry);
        range::resize(out, count);
        for (std::size_t i = 0; i < count; ++i)
        {
            geometry::detail::densify::densify_range<>::apply(
                range::at(geometry, i), range::at(out, i), len, strategies);
        }
    }
    else if constexpr (concepts::ConstRing<Geometry>)
    {
        geometry::detail::densify::densify_ring
            <
                geometry::closure<Geometry>::value != geometry::open,
                geometry::closure<GeometryOut>::value != geometry::open
            >::apply(geometry, out, len, strategies);
    }
    else if constexpr (concepts::ConstPolygon<Geometry>)
    {
        dispatch::densify(exterior_ring(geometry), exterior_ring(out),
                          len, strategies);

        std::size_t const count = boost::size(interior_rings(geometry));
        range::resize(interior_rings(out), count);
        for (std::size_t i = 0; i < count; ++i)
        {
            dispatch::densify(range::at(interior_rings(geometry), i),
                              range::at(interior_rings(out), i),
                              len, strategies);
        }
    }
    else
    {
        std::size_t const count = boost::size(geometry);
        range::resize(out, count);
        for (std::size_t i = 0; i < count; ++i)
        {
            dispatch::densify(range::at(geometry, i), range::at(out, i),
                              len, strategies);
        }
    }
}


} // namespace dispatch
#endif // DOXYGEN_NO_DISPATCH


namespace resolve_strategy
{

template <concepts::ConstGeometry Geometry,
          concepts::MutableGeometry GeometryOut,
          typename Distance,
          typename Strategy>
inline void densify(Geometry const& geometry, GeometryOut& out,
                    Distance const& max_distance, Strategy const& strategy)
{
    if constexpr (std::same_as<Strategy, default_strategy>)
    {
        using strategies_type = typename strategies::densify::services
            ::default_strategy<Geometry>::type;
        dispatch::densify(geometry, out, max_distance, strategies_type());
    }
    else if constexpr (strategies::detail::is_umbrella_strategy<Strategy>::value)
    {
        dispatch::densify(geometry, out, max_distance, strategy);
    }
    else
    {
        using strategies::densify::services::strategy_converter;
        dispatch::densify(geometry, out, max_distance,
                          strategy_converter<Strategy>::get(strategy));
    }
}

} // namespace resolve_strategy


namespace resolve_dynamic {

template <concepts::MutableGeometry Geometry,
          typename Distance,
          typename Strategy>
inline void densify(Geometry const& geometry, Geometry& out,
                    Distance const& max_distance, Strategy const& strategy)
{
    if constexpr (concepts::ConstDynamicGeometry<Geometry>)
    {
        traits::visit<Geometry>::apply([&](auto const& g)
        {
            using geometry_type = util::remove_cref_t<decltype(g)>;
            geometry_type result;
            resolve_dynamic::densify(
                g, result, max_distance, strategy);
            out = std::move(result);
        }, geometry);
    }
    else if constexpr (concepts::ConstGeometryCollection<Geometry>)
    {
        detail::visit_breadth_first([&](auto const& g)
        {
            using geometry_type = util::remove_cref_t<decltype(g)>;
            geometry_type result;
            resolve_dynamic::densify(
                g, result, max_distance, strategy);
            traits::emplace_back<Geometry>::apply(out, std::move(result));
            return true;
        }, geometry);
    }
    else
    {
        resolve_strategy::densify(
            geometry, out, max_distance, strategy);
    }
}

} // namespace resolve_dynamic


/*!
\brief Densify a geometry using a specified strategy
\ingroup densify
\tparam Geometry \tparam_geometry
\tparam Distance A numerical distance measure
\tparam Strategy A type fulfilling a DensifyStrategy concept
\param geometry Input geometry, to be densified
\param out Output geometry, densified version of the input geometry
\param max_distance Distance threshold (in units depending on strategy)
\param strategy Densify strategy to be used for densification

\qbk{distinguish,with strategy}
\qbk{[include reference/algorithms/densify.qbk]}

\qbk{
[heading Available Strategies]
\* [link geometry.reference.strategies.strategy_densify_cartesian Cartesian]
\* [link geometry.reference.strategies.strategy_densify_spherical Spherical]
\* [link geometry.reference.strategies.strategy_densify_geographic Geographic]

[heading Example]
[densify_strategy]
[densify_strategy_output]

[heading See also]
\* [link geometry.reference.algorithms.line_interpolate line_interpolate]
}
*/
template <concepts::MutableGeometry Geometry,
          typename Distance,
          typename Strategy>
inline void densify(Geometry const& geometry,
                    Geometry& out,
                    Distance const& max_distance,
                    Strategy const& strategy)
{
    if (max_distance <= Distance(0))
    {
        BOOST_THROW_EXCEPTION(geometry::invalid_input_exception());
    }

    geometry::clear(out);

    resolve_dynamic::densify(geometry, out, max_distance, strategy);
}


/*!
\brief Densify a geometry
\ingroup densify
\tparam Geometry \tparam_geometry
\tparam Distance A numerical distance measure
\param geometry Input geometry, to be densified
\param out Output geometry, densified version of the input geometry
\param max_distance Distance threshold (in units depending on coordinate system)

\qbk{[include reference/algorithms/densify.qbk]}

\qbk{
[heading Example]
[densify]
[densify_output]

[heading See also]
\* [link geometry.reference.algorithms.line_interpolate line_interpolate]
}
*/
template <concepts::MutableGeometry Geometry, typename Distance>
inline void densify(Geometry const& geometry,
                    Geometry& out,
                    Distance const& max_distance)
{
    densify(geometry, out, max_distance, default_strategy());
}


}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_ALGORITHMS_DENSIFY_HPP
