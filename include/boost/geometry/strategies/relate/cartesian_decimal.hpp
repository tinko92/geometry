// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.

#ifndef BOOST_GEOMETRY_STRATEGIES_RELATE_CARTESIAN_DECIMAL_HPP
#define BOOST_GEOMETRY_STRATEGIES_RELATE_CARTESIAN_DECIMAL_HPP

#include <algorithm>
#include <cstddef>
#include <type_traits>

#include <boost/geometry/strategies/cartesian/exact/decimal_point.hpp>
#include <boost/geometry/strategies/relate/cartesian.hpp>

namespace boost { namespace geometry
{
namespace detail { namespace filtered_decimal
{

namespace ec = strategies::detail::exact_cartesian;
namespace fe = ec::forward_error;

inline approximation const* enclosure(point const& p) { return p.certified_approximation(); }
template <typename P> inline approximation const* enclosure(P const&) { return nullptr; }

struct predicates
{
    using cs_tag = cartesian_tag;

    template <typename A, typename B, typename P>
    static int apply(A const& a, B const& b, P const& p)
    {
        auto const* x = enclosure(a);
        auto const* y = enclosure(b);
        auto const* z = enclosure(p);
        int result;
        if (x && y && z && fe::try_orientation(*x, *y, *z, result))
        {
#ifdef BOOST_GEOMETRY_EXACT_FILTER_STATISTICS
            ++ec::filtering_statistics().decimal_side_hits;
#endif
            return result;
        }
#ifdef BOOST_GEOMETRY_EXACT_FILTER_STATISTICS
        ++ec::filtering_statistics().decimal_side_fallbacks;
#endif
        using fallback = strategy::side::services::default_strategy<cartesian_tag, rational>::type;
        return fallback::apply(a, b, p);
    }

    template <typename A, typename B>
    static bool try_equal(A const& a, B const& b, bool& result)
    {
        auto const* x = enclosure(a);
        auto const* y = enclosure(b);
        if (!x || !y) return false;
        // Distinct nearest-even views imply unequal exact coordinates.
        // Equal views do not certify equality unless both were exact.
        if (x->value != y->value) result = false;
        else if (fe::detail::exact_coordinates(*x) && fe::detail::exact_coordinates(*y)) result = true;
        else return false;
#ifdef BOOST_GEOMETRY_EXACT_FILTER_STATISTICS
        ++ec::filtering_statistics().decimal_equality_hits;
#endif
        return true;
    }

    template <std::size_t D, typename P, typename Q>
    static bool try_disjoint_by_range(P const& p0, P const& p1, Q const& q0, Q const& q1, bool& result)
    {
        auto const* a = enclosure(p0);
        auto const* b = enclosure(p1);
        auto const* c = enclosure(q0);
        auto const* d = enclosure(q1);
        if (!a || !b || !c || !d) return false;
        double const pmin = (std::min)(a->value[D], b->value[D]);
        double const pmax = (std::max)(a->value[D], b->value[D]);
        double const qmin = (std::min)(c->value[D], d->value[D]);
        double const qmax = (std::max)(c->value[D], d->value[D]);
        // Nearest-even is monotone and commutes with min/max. Strict
        // comparisons survive unrounding, ties require rational order.
        if (pmax < qmin || qmax < pmin) result = true;
        else if (pmax > qmin && qmax > pmin) result = false;
        else return false;
#ifdef BOOST_GEOMETRY_EXACT_FILTER_STATISTICS
        ++ec::filtering_statistics().decimal_range_hits;
#endif
        return true;
    }
};

struct segment_intersection : strategy::intersection::cartesian_segments<rational>
{
    template <typename A, typename B, typename Policy>
    static auto apply(A const& a, B const& b, Policy const& policy)
    {
        return strategy::intersection::cartesian_segments<rational>::apply(a, b, policy, predicates());
    }
};

}} // namespace detail::filtered_decimal

namespace strategies { namespace relate
{

// Explicit strategy for prepare_filtered_decimal_polygon; rational
// constructions and output. Same strict FP contract as Cartesian exact
// predicates (round-to-nearest and no fast-math or implicit contraction).
class cartesian_decimal : public cartesian<geometry::detail::filtered_decimal::rational>
{
    using rational = geometry::detail::filtered_decimal::rational;
    using predicates = geometry::detail::filtered_decimal::predicates;
    using base = cartesian<rational>;

public:
    using base::relate;
    static auto side() { return predicates(); }
    static auto relate() { return geometry::detail::filtered_decimal::segment_intersection(); }

    template <typename A, typename B>
    static auto relate(A const&, B const&,
        std::enable_if_t<util::is_pointlike<A>::value
            && (util::is_linear<B>::value || util::is_polygonal<B>::value)>* = nullptr)
    {
        return strategy::within::detail::cartesian_winding_base<predicates, rational>();
    }
};

}} // namespace strategies::relate
}} // namespace boost::geometry

#endif
