// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.

#ifndef BOOST_GEOMETRY_STRATEGIES_RELATE_SIMPLE_CARTESIAN_HPP
#define BOOST_GEOMETRY_STRATEGIES_RELATE_SIMPLE_CARTESIAN_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <type_traits>
#include <boost/geometry/algorithms/area.hpp>
#include <boost/geometry/core/access.hpp>
#include <boost/geometry/core/coordinate_dimension.hpp>
#include <boost/geometry/core/coordinate_system.hpp>
#include <boost/geometry/core/coordinate_type.hpp>
#include <boost/geometry/geometries/segment.hpp>
#include <boost/geometry/policies/relate/intersection_policy.hpp>
#include <boost/geometry/strategies/cartesian/clip.hpp>
#include <boost/geometry/strategies/detail/box_intersection.hpp>
#include <boost/geometry/strategies/detail/box_topology.hpp>
#include <boost/geometry/strategies/detail/distance_measure.hpp>
#include <boost/geometry/strategies/detail/ring_properties.hpp>
#include <boost/geometry/strategies/detail/section.hpp>
#include <boost/geometry/strategies/relate/cartesian.hpp>
#include <boost/geometry/strategies/relate/intersection_types.hpp>

namespace boost { namespace geometry { namespace strategies { namespace detail
{ namespace simple_cartesian
{

template <typename Geometry>
inline void check_geometry()
{
    static_assert(std::is_same<geometry::coordinate_type_t<Geometry>, double>::value,
                  "simple_cartesian requires native double coordinates");
    static_assert(geometry::dimension<Geometry>::value == 2
        && std::is_same<geometry::coordinate_system_t<Geometry>, cs::cartesian>::value,
        "simple_cartesian requires two-dimensional Cartesian geometries");
}

template <std::size_t D, typename Point>
inline double coordinate(Point const& point)
{
    check_geometry<Point>();
    return geometry::get<D>(point);
}

inline int sign(double value) { return value > 0 ? 1 : value < 0 ? -1 : 0; }

template <std::size_t D, typename A, typename B>
inline int compare_coordinate(A const& a, B const& b)
{
    double const x = coordinate<D>(a), y = coordinate<D>(b);
    return x < y ? -1 : x > y ? 1 : 0;
}

template <typename A, typename B, typename P>
inline double determinant(A const& a, B const& b, P const& p)
{
    double const ax = coordinate<0>(a), ay = coordinate<1>(a);
    return (coordinate<0>(b) - ax) * (coordinate<1>(p) - ay)
         - (coordinate<1>(b) - ay) * (coordinate<0>(p) - ax);
}

struct point_equal
{
    using cs_tag = cartesian_tag;
    template <typename A, typename B>
    static bool apply(A const& a, B const& b)
    { return coordinate<0>(a) == coordinate<0>(b) && coordinate<1>(a) == coordinate<1>(b); }
};

struct side
{
    using cs_tag = cartesian_tag;
    template <typename A, typename B, typename P>
    static int apply(A const& a, B const& b, P const& p) { return sign(determinant(a, b, p)); }
    template <typename A, typename B, typename P>
    static int direction_code(A const& a, B const& b, P const& p)
    {
        return sign((coordinate<0>(b) - coordinate<0>(a)) * (coordinate<0>(p) - coordinate<0>(b))
                  + (coordinate<1>(b) - coordinate<1>(a)) * (coordinate<1>(p) - coordinate<1>(b)));
    }
};

template <typename ComparePolicy>
struct point_compare
{
    template <typename A, typename B>
    static bool apply(A const& a, B const& b)
    {
        int const x = compare_coordinate<0>(a, b);
        int const order = x == 0 ? compare_coordinate<1>(a, b) : x;
        return std::is_same<ComparePolicy, strategy::compare::equal_to>::value
            ? order == 0 : ComparePolicy::apply(order, 0);
    }
};

struct winding
{
    using cs_tag = cartesian_tag;
    struct state_type { int count = 0; bool boundary = false; };
    template <typename P, typename A, typename B = A>
    static bool apply(P const& point, A const& a, B const& b, state_type& state)
    {
        int const ay = compare_coordinate<1>(a, point), by = compare_coordinate<1>(b, point);
        if (ay * by > 0) return true;
        int const turn = side::apply(a, b, point);
        if (turn == 0 && compare_coordinate<0>(a, point) * compare_coordinate<0>(b, point) <= 0)
        { state.boundary = true; return false; }
        if (ay <= 0 && by > 0 && turn > 0) ++state.count;
        else if (by <= 0 && ay > 0 && turn < 0) --state.count;
        return true;
    }
    static int result(state_type const& state) { return state.boundary ? 0 : state.count == 0 ? -1 : 1; }
};

struct squared_distance
{
    template <typename A, typename B>
    static double apply(A const& a, B const& b)
    {
        double const dx = coordinate<0>(a) - coordinate<0>(b), dy = coordinate<1>(a) - coordinate<1>(b);
        return dx * dx + dy * dy;
    }
};

struct topology
{
    template <typename A, typename B>
    static bool cluster_equal(A const& a, B const& b) { return point_equal::apply(a, b); }
    template <typename A, typename B>
    static int compare_y(A const& a, B const& b) { return compare_coordinate<1>(a, b); }
    template <typename A, typename B>
    static bool cluster_window_exceeded(A const& upper, B const& lower) { return compare_y(upper, lower) > 0; }
    template <typename Indices, typename Turns>
    static void colocate_cluster(Indices const& indices, Turns& turns)
    {
        auto const point = turns[*indices.begin()].point;
        for (auto const index : indices) turns[index].point = point;
    }
    template <typename A, typename B>
    static bool same_direction_point(A const& a, B const& b) { return point_equal::apply(a, b); }
    template <typename A, typename B>
    static double squared_distance(A const& a, B const& b) { return simple_cartesian::squared_distance::apply(a, b); }
    static bool distance_is_zero(double value) { return value == 0; }
    template <typename A, typename B, typename P>
    static strategies::detail::distance_measure<double> distance_measure(A const& a, B const& b, P const& p)
    {
        strategies::detail::distance_measure<double> result;
        result.measure = determinant(a, b, p);
        return result;
    }
};

struct segment_position
{
    template <typename A, typename B, typename P, typename Q>
    static int apply(A const& a, B const& b, P const& p, Q const& q)
    {
        int const x = compare_coordinate<0>(a, b);
        int const order = x == 0 ? compare_coordinate<1>(p, q) : compare_coordinate<0>(p, q);
        int const direction = x == 0 ? compare_coordinate<1>(a, b) : x;
        return direction > 0 ? -order : order;
    }
    template <typename A, typename B, typename P>
    static int classify_collinear(A const& a, B const& b, P const& p)
    {
        int const from = apply(a, b, p, a);
        if (from == 0) return 1;
        if (from < 0) return 0;
        int const to = apply(a, b, p, b);
        return to < 0 ? 2 : to == 0 ? 3 : 4;
    }
    template <typename A, typename B, typename P>
    static int classify(A const& a, B const& b, P const& p)
    {
        // No-fraction relate helpers also distinguish degenerate/off-line points.
        if (point_equal::apply(a, b)) return point_equal::apply(a, p) ? 5 : 6;
        if (point_equal::apply(a, p)) return 1;
        if (point_equal::apply(b, p)) return 3;
        return side::apply(a, b, p) == 0 ? classify_collinear(a, b, p) : 7;
    }
};

struct position_tag
{
    using no_fraction = void;
    template <typename A, typename B, typename P, typename Q>
    static int compare(A const& a, B const& b, P const& p, Q const& q)
    { return segment_position::apply(a, b, p, q); }
    template <typename A, typename B>
    static bool equal(A const& a, B const& b) { return point_equal::apply(a, b); }
};

template <typename Point>
struct intersection_result
{
    using point_type = Point;
    std::size_t count = 0;
    Point intersections[2];
};

template <typename Source, typename Destination>
inline void copy_point(Source const& source, Destination& destination)
{
    check_geometry<Destination>();
    geometry::set<0>(destination, coordinate<0>(source));
    geometry::set<1>(destination, coordinate<1>(source));
}

template <typename Policy>
struct endpoint_policy : Policy {};

template <typename Result>
struct endpoint_policy<policies::relate::segments_intersection_points<Result>>
    : policies::relate::segments_intersection_points<Result>
{
    using return_type = Result;
    template <typename Info, typename A, typename B>
    static Result segments_crosses(side_info const&, Info const& info, A const& a, B const& b)
    {
        Result result;
        result.count = 1;
        info.calculate(result.intersections[0], a, b);
        return result;
    }
    template <typename Info, typename Point>
    static Result segments_share_common_point(side_info const&, Info const&, Point const& point)
    {
        Result result;
        result.count = 1;
        copy_point(point, result.intersections[0]);
        return result;
    }
    template <typename A, typename B, typename Token>
    static Result segments_collinear(A const& a, B const& b, bool,
        int a1, int a2, int b1, int b2, Token const&, Token const&, Token const&, Token const&)
    {
        Result result;
        std::array<unsigned, 2> selected;
        auto append = [&](unsigned endpoint) { selected[result.count++] = endpoint; };
        if (a1 >= 1 && a1 <= 3 && result.count < 2) append(0);
        if (b1 == 2 && result.count < 2) append(2);
        if (a2 >= 1 && a2 <= 3 && result.count < 2) append(1);
        if (b2 == 2 && result.count < 2) append(3);
        auto visit = [&](unsigned endpoint, auto apply)
        {
            switch (endpoint)
            {
            case 0: return apply(a.first);
            case 1: return apply(a.second);
            case 2: return apply(b.first);
            default: return apply(b.second);
            }
        };
        if (result.count == 2 && visit(selected[1], [&](auto const& second)
            { return visit(selected[0], [&](auto const& first)
                { return segment_position::apply(a.first, a.second, second, first); }); }) < 0)
            std::swap(selected[0], selected[1]);
        for (std::size_t i = 0; i != result.count; ++i)
            visit(selected[i], [&](auto const& point) { copy_point(point, result.intersections[i]); });
        return result;
    }
    template <typename Segment>
    static Result degenerate(Segment const& segment, bool)
    {
        Result result;
        result.count = 1;
        copy_point(segment.first, result.intersections[0]);
        return result;
    }
    template <typename Segment, typename Token>
    static Result one_degenerate(Segment const& segment, Token const&, bool a_degenerate)
    { return degenerate(segment, a_degenerate); }
};

template <typename Result>
struct endpoint_policy<policies::relate::segments_intersection_policy<Result>>
    : policies::relate::segments_intersection_policy<Result>
{
    using original = policies::relate::segments_intersection_policy<Result>;
    using return_type = typename original::return_type;
    using points = endpoint_policy<policies::relate::segments_intersection_points<Result>>;
    using direction = policies::relate::segments_direction;
    template <typename Info, typename A, typename B>
    static return_type segments_crosses(side_info const& sides, Info const& info, A const& a, B const& b)
    { return {points::segments_crosses(sides, info, a, b), direction::segments_crosses(sides, info, a, b)}; }
    template <typename Info, typename Point>
    static return_type segments_share_common_point(side_info const& sides, Info const& info, Point const& point)
    { return {points::segments_share_common_point(sides, info, point), direction::segments_share_common_point(sides, info, point)}; }
    template <typename A, typename B, typename Token>
    static return_type segments_collinear(A const& a, B const& b, bool opposite,
        int a1, int a2, int b1, int b2, Token const& t1, Token const& t2, Token const& t3, Token const& t4)
    {
        return {points::segments_collinear(a, b, opposite, a1, a2, b1, b2, t1, t2, t3, t4),
                direction::segments_collinear(a, b, opposite, a1, a2, b1, b2, t1, t2, t3, t4)};
    }
    template <typename Segment>
    static return_type degenerate(Segment const& segment, bool a_degenerate)
    { return {points::degenerate(segment, a_degenerate), direction::degenerate(segment, a_degenerate)}; }
    template <typename Segment, typename Token>
    static return_type one_degenerate(Segment const& segment, Token const& token, bool a_degenerate)
    { return {points::one_degenerate(segment, token, a_degenerate), direction::one_degenerate(segment, token, a_degenerate)}; }
};

struct segments
{
    using cs_tag = cartesian_tag;
    struct construction_info
    {
        int endpoint;
        template <typename Output, typename A, typename B>
        void calculate(Output& output, A const& a, B const& b) const
        {
            switch (endpoint)
            {
            case 0: copy_point(a.first, output); return;
            case 1: copy_point(a.second, output); return;
            case 2: copy_point(b.first, output); return;
            case 3: copy_point(b.second, output); return;
            }
            double const ax = coordinate<0>(a.first), ay = coordinate<1>(a.first);
            double const dx = coordinate<0>(a.second) - ax, dy = coordinate<1>(a.second) - ay;
            double const ex = coordinate<0>(b.second) - coordinate<0>(b.first);
            double const ey = coordinate<1>(b.second) - coordinate<1>(b.first);
            double const t = ((coordinate<0>(b.first) - ax) * ey - (coordinate<1>(b.first) - ay) * ex)
                           / (dx * ey - dy * ex);
            check_geometry<Output>();
            geometry::set<0>(output, ax + t * dx);
            geometry::set<1>(output, ay + t * dy);
        }
    };

    template <typename Range1, typename Range2, typename Policy>
    static typename Policy::return_type apply(Range1 const& p, Range2 const& q, Policy const&)
    {
        using selected = endpoint_policy<Policy>;
        auto const& p0 = p.at(0); auto const& p1 = p.at(1);
        auto const& q0 = q.at(0); auto const& q1 = q.at(1);
        using P = typename std::decay<decltype(p0)>::type;
        using Q = typename std::decay<decltype(q0)>::type;
        model::referring_segment<P const> const a(p0, p1);
        model::referring_segment<Q const> const b(q0, q1);
        int const sa0 = side::apply(q0, q1, p0), sa1 = side::apply(q0, q1, p1);
        if (sa0 * sa1 > 0) return selected::disjoint();
        int const sb0 = side::apply(p0, p1, q0), sb1 = side::apply(p0, p1, q1);
        if (sb0 * sb1 > 0) return selected::disjoint();
        if (sa0 == 0 && sa1 == 0 && sb0 == 0 && sb1 == 0)
        {
            bool const ap = point_equal::apply(p0, p1), bp = point_equal::apply(q0, q1);
            position_tag const token;
            if (ap && bp) return point_equal::apply(p0, q0) ? selected::degenerate(a, true) : selected::disjoint();
            if (ap || bp)
            {
                int const location = ap ? segment_position::classify_collinear(q0, q1, p0)
                                        : segment_position::classify_collinear(p0, p1, q0);
                if (location < 1 || location > 3) return selected::disjoint();
                return ap ? selected::one_degenerate(a, token, true) : selected::one_degenerate(b, token, false);
            }
            int const a1 = segment_position::classify_collinear(q0, q1, p0);
            int const a2 = segment_position::classify_collinear(q0, q1, p1);
            if ((a1 < 1 && a2 < 1) || (a1 > 3 && a2 > 3)) return selected::disjoint();
            int const b1 = segment_position::classify_collinear(p0, p1, q0);
            int const b2 = segment_position::classify_collinear(p0, p1, q1);
            return selected::segments_collinear(a, b, segment_position::apply(p0, p1, q0, q1) > 0,
                                                a1, a2, b1, b2, token, token, token, token);
        }
        int const endpoint = sa0 == 0 ? 0 : sa1 == 0 ? 1 : sb0 == 0 ? 2 : sb1 == 0 ? 3 : -1;
        if (endpoint < 0)
        {
            double const dx = coordinate<0>(p1) - coordinate<0>(p0), dy = coordinate<1>(p1) - coordinate<1>(p0);
            double const ex = coordinate<0>(q1) - coordinate<0>(q0), ey = coordinate<1>(q1) - coordinate<1>(q0);
            if (dx * ey - dy * ex == 0) return selected::disjoint();
        }
        construction_info const info{endpoint};
        side_info sides;
        sides.template set<0>(sa0, sa1); sides.template set<1>(sb0, sb1);
        if (sa0 == 0 && (point_equal::apply(p0, q0) || point_equal::apply(p0, q1)))
            return selected::segments_share_common_point(sides, info, p0);
        if (sa1 == 0 && (point_equal::apply(p1, q0) || point_equal::apply(p1, q1)))
            return selected::segments_share_common_point(sides, info, p1);
        return selected::segments_crosses(sides, info, a, b);
    }
};

struct segment_box
{
    using cs_tag = cartesian_tag;
    static strategy::covered_by::cartesian_point_box get_disjoint_point_box_strategy() { return {}; }
    template <typename Segment, typename Box>
    static bool apply(Segment const& segment, Box const& box)
    {
        check_geometry<Segment>(); check_geometry<Box>();
        double const ax = geometry::get<0, 0>(segment), ay = geometry::get<0, 1>(segment);
        double const bx = geometry::get<1, 0>(segment), by = geometry::get<1, 1>(segment);
        double const lx = geometry::get<min_corner, 0>(box), ly = geometry::get<min_corner, 1>(box);
        double const hx = geometry::get<max_corner, 0>(box), hy = geometry::get<max_corner, 1>(box);
        if (lx > hx || ly > hy || (std::max)(ax, bx) < lx || (std::min)(ax, bx) > hx
            || (std::max)(ay, by) < ly || (std::min)(ay, by) > hy) return true;
        auto turn = [&](double x, double y) { return sign((bx - ax) * (y - ay) - (by - ay) * (x - ax)); };
        int const first = turn(lx, ly);
        return first != 0 && turn(lx, hy) == first && turn(hx, ly) == first && turn(hx, hy) == first;
    }
};

template <typename Box, typename Point>
struct clipping : strategy::intersection::liang_barsky<Box, Point>
{
    using equals_point_point_strategy_type = point_equal;
    static point_equal get_equals_point_point_strategy() { return {}; }
};

struct ring_topology
{
    template <typename Geometry>
    static double signed_measure(Geometry const& geometry)
    {
        check_geometry<Geometry>();
        return geometry::area(geometry, strategies::relate::cartesian<double>());
    }
    static double absolute_measure(double value) { return std::abs(value); }
    static bool measure_less(double a, double b) { return a < b; }
    static int measure_sign(double value) { return sign(value); }
    static bool measure_is_positive(double value) { return value > 0; }
    static bool measure_is_negative(double value) { return value < 0; }
    static bool measure_is_zero(double value) { return value == 0; }
    template <typename Box> static void expand_envelope(Box&) {}
};

}}}}}

namespace boost { namespace geometry
{
template <typename Point>
struct segment_intersection_points<Point, strategies::detail::simple_cartesian::position_tag>
    : strategies::detail::simple_cartesian::intersection_result<Point> {};
}}

namespace boost { namespace geometry { namespace strategies { namespace relate
{

// Plain double arithmetic, with no tolerance or certified-predicate fallback.
// Requires finite 2D double input and finite arithmetic intermediates. Ring
// parenting uses stock signed area, not the exact strategy's boundary topology.
// Scope: relate and set operations; public equals() legacy tolerance shortcuts
// are not covered. Primitive point equality and coordinate ordering are strict.
class simple_cartesian : public cartesian<double>
{
    using base = cartesian<double>;
    using primitives = strategies::detail::simple_cartesian::point_equal;
public:
    static constexpr bool exact_predicates = false;
    using segment_position_type = strategies::detail::simple_cartesian::position_tag;
    template <typename Point> using intersection_result_type = strategies::detail::simple_cartesian::intersection_result<Point>;
    template <typename Point> using comparable_distance_type = double;
    template <typename Geometry> using ring_measure_type = double;
    template <typename ComparePolicy, typename EqualsPolicy>
    using compare_type = strategies::detail::simple_cartesian::point_compare<ComparePolicy>;
    using exact_cluster_compare = strategies::detail::simple_cartesian::point_compare<strategy::compare::less>;

    template <typename Geometry>
    static auto area(Geometry const& geometry)
    { strategies::detail::simple_cartesian::check_geometry<Geometry>(); return base::area(geometry); }
    template <typename Geometry, typename Box>
    static auto envelope(Geometry const& geometry, Box const& box)
    {
        strategies::detail::simple_cartesian::check_geometry<Geometry>();
        strategies::detail::simple_cartesian::check_geometry<Box>();
        return base::envelope(geometry, box);
    }
    template <typename Box, typename Geometry>
    static auto expand(Box const& box, Geometry const& geometry)
    {
        strategies::detail::simple_cartesian::check_geometry<Box>();
        strategies::detail::simple_cartesian::check_geometry<Geometry>();
        return base::expand(box, geometry);
    }
    template <typename A, typename B>
    static auto within(A const& a, B const& b,
        std::enable_if_t<(util::is_pointlike<A>::value || util::is_box<A>::value) && util::is_box<B>::value>* = nullptr)
    {
        strategies::detail::simple_cartesian::check_geometry<A>();
        strategies::detail::simple_cartesian::check_geometry<B>();
        return base::within(a, b);
    }
    template <typename A, typename B>
    static auto covered_by(A const& a, B const& b,
        std::enable_if_t<(util::is_pointlike<A>::value || util::is_box<A>::value) && util::is_box<B>::value>* = nullptr)
    {
        strategies::detail::simple_cartesian::check_geometry<A>();
        strategies::detail::simple_cartesian::check_geometry<B>();
        return base::covered_by(a, b);
    }
    template <typename A, typename B>
    static auto disjoint(A const& a, B const& b,
        std::enable_if_t<util::is_box<A>::value && util::is_box<B>::value>* = nullptr)
    {
        strategies::detail::simple_cartesian::check_geometry<A>();
        strategies::detail::simple_cartesian::check_geometry<B>();
        return base::disjoint(a, b);
    }
    template <typename A, typename B>
    static auto disjoint(A const&, B const&,
        std::enable_if_t<util::is_segment<A>::value && util::is_box<B>::value>* = nullptr)
    { return strategies::detail::simple_cartesian::segment_box(); }
    static auto side() { return strategies::detail::simple_cartesian::side(); }
    static auto topology() { return strategies::detail::simple_cartesian::topology(); }
    static auto ring_topology() { return strategies::detail::simple_cartesian::ring_topology(); }
    static auto section() { return strategies::detail::strict_cartesian_section(); }
    static auto segment_position() { return strategies::detail::simple_cartesian::segment_position(); }
    static auto box_intersection() { return strategies::detail::legacy_box_intersection(); }
    static auto box_topology() { return strategies::detail::legacy_box_topology(); }
    template <typename A, typename B>
    static auto comparable_distance(A const&, B const&) { return strategies::detail::simple_cartesian::squared_distance(); }
    template <typename Point, typename Box>
    static auto clipping(Box const&) { return strategies::detail::simple_cartesian::clipping<Box, Point>(); }
    static auto relate() { return strategies::detail::simple_cartesian::segments(); }
    template <typename A, typename B>
    static auto relate(A const&, B const&,
        std::enable_if_t<util::is_pointlike<A>::value && util::is_pointlike<B>::value>* = nullptr)
    { return primitives(); }
    template <typename A, typename B>
    static auto relate(A const&, B const&,
        std::enable_if_t<util::is_pointlike<A>::value && (util::is_linear<B>::value || util::is_polygonal<B>::value)>* = nullptr)
    { return strategies::detail::simple_cartesian::winding(); }
};

}}}}

#endif
