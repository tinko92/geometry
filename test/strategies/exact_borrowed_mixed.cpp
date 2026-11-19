// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.

#define BOOST_GEOMETRY_EXACT_FILTER_STATISTICS
#include <boost/core/lightweight_test.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/strategies/relate/cartesian_exact.hpp>
#include <array>
#include <cstdint>
#include <limits>

struct native_point
{
    double x, y;
    native_point(double x, double y) : x(x), y(y) {}
    native_point(native_point const&) = delete;
    native_point& operator=(native_point const&) = delete;
};

struct borrowed_box
{
    native_point const& lower;
    native_point const& upper;
};

namespace boost { namespace geometry { namespace traits
{
template <> struct tag<native_point> { using type = point_tag; };
template <> struct coordinate_type<native_point> { using type = double; };
template <> struct coordinate_system<native_point> { using type = cs::cartesian; };
template <> struct dimension<native_point> : std::integral_constant<std::size_t, 2> {};
template <std::size_t D> struct access<native_point, D>
{
    static double get(native_point const& point) { return D == 0 ? point.x : point.y; }
};
template <> struct tag<borrowed_box> { using type = box_tag; };
template <> struct point_type<borrowed_box> { using type = native_point; };
template <std::size_t I, std::size_t D> struct indexed_access<borrowed_box, I, D>
{
    static double get(borrowed_box const& box)
    { return geometry::get<D>(I == min_corner ? box.lower : box.upper); }
};
}}}

namespace bg = boost::geometry;
namespace ec = bg::strategies::detail::exact_cartesian;
using rational = ec::exact_number;

struct borrowed_range
{
    native_point const& first;
    native_point const& second;
    native_point const& at(std::size_t i) const { return i == 0 ? first : second; }
};

void segment_classification()
{
    native_point const a(0, 0), b(4, 0), c(2, -1), d(2, 1);
    native_point const e(1, 0), f(3, 0), g(4, 1);
    auto check = [](borrowed_range const& a, borrowed_range const& b)
    {
        ec::filtering_statistics() = {};
        auto const result = ec::segment_intersection<double>::apply(a, b,
            bg::policies::relate::segments_direction());
        BOOST_TEST(result.how != 'E');
        BOOST_TEST_EQ(ec::filtering_statistics().native_point_constructions, 0u);
    };
    check({a, b}, {c, d});
    check({a, b}, {e, f});
    check({a, b}, {f, e});
    check({a, b}, {b, g});
    check({a, b}, {e, e});
    check({a, a}, {a, a});
}

void box_access()
{
    native_point const low(0, 0), high(4, 4), inner_low(1, 1), inner_high(3, 3);
    borrowed_box const outer{low, high}, inner{inner_low, inner_high};
    using predicates = ec::static_box_predicates<double>;
    using output_point = bg::model::d2::point_xy<float>;
    bg::model::box<output_point> output;
    std::array<native_point, 3> const points{{{0, 2}, {4, 1}, {1, 3}}};
    ec::filtering_statistics() = {};
    BOOST_TEST(predicates::point_covered_by(inner_low, outer));
    BOOST_TEST(predicates::within(inner, outer));
    BOOST_TEST(!predicates::disjoint(outer, inner));
    BOOST_TEST(!predicates::touches(outer, inner));
    ec::static_envelope_box<double>::apply(outer, output);
    BOOST_TEST_EQ((bg::get<bg::max_corner, 0>(output)), 4);
    ec::static_envelope_range<double>::apply(points, output);
    BOOST_TEST_EQ((bg::get<bg::min_corner, 0>(output)), 0);
    BOOST_TEST_EQ((bg::get<bg::min_corner, 1>(output)), 1);
    BOOST_TEST_EQ((bg::get<bg::max_corner, 0>(output)), 4);
    BOOST_TEST_EQ((bg::get<bg::max_corner, 1>(output)), 3);
    ec::static_expand_box<double>::apply(output, outer);
    BOOST_TEST_EQ((bg::get<bg::max_corner, 1>(output)), 4);
    BOOST_TEST_EQ(ec::filtering_statistics().native_point_constructions, 0u);
}

template <typename A, typename B, typename P>
void orientation(A const& a, B const& b, P const& p)
{
    auto coordinate = [](auto const& point)
    { return std::array<rational, 2>{{rational(bg::get<0>(point)), rational(bg::get<1>(point))}}; };
    auto const x = coordinate(a), y = coordinate(b), z = coordinate(p);
    rational const determinant = (y[0] - x[0]) * (z[1] - x[1]) - (y[1] - x[1]) * (z[0] - x[0]);
    ec::filtering_statistics() = {};
    BOOST_TEST_EQ(ec::static_orientation<double>(a, b, p), determinant.sign());
    BOOST_TEST_EQ(ec::static_orientation<double>(b, a, p), -determinant.sign());
    BOOST_TEST_EQ(ec::static_orientation<double>(b, p, a), determinant.sign());
    BOOST_TEST_EQ(ec::filtering_statistics().native_point_constructions, 0u);
}

void native_support_orientation()
{
    using point = ec::exact_point<double>;
    point const a(0, 0), b(1, 1), c(0, 1), d(1, -1);
    auto const crossing = ec::construct_intersection(a, b, c, d);
    native_point const first(0, 0), second(1, 1), beyond(2, 2), off(2, 3);
    auto const cache_size = ec::point_cache_detail::native_cache<double>::size();
    ec::filtering_statistics() = {};
    BOOST_TEST(crossing.native_points_on_producing_line(first, beyond));
    BOOST_TEST(!crossing.native_points_on_producing_line(first, off));
    BOOST_TEST_EQ(ec::static_orientation<double>(first, second, crossing), 0);
    BOOST_TEST_EQ(ec::static_orientation<double>(second, first, crossing), 0);
    BOOST_TEST_EQ(ec::static_orientation<double>(first, crossing, second), 0);
    BOOST_TEST_EQ(ec::static_orientation<double>(second, crossing, first), 0);
    BOOST_TEST_EQ(ec::static_orientation<double>(crossing, first, second), 0);
    BOOST_TEST_EQ(ec::static_orientation<double>(crossing, second, first), 0);
    BOOST_TEST_EQ(ec::point_cache_detail::native_cache<double>::size(), cache_size);
    BOOST_TEST_EQ(ec::filtering_statistics().native_point_constructions, 0u);
    orientation(first, off, crossing);

    // An exact concurrency on a third, unrelated support still needs a predicate.
    point const h0(0, 1), h1(2, 1), v0(1, 0), v1(1, 2);
    auto const concurrent = ec::construct_intersection(h0, h1, v0, v1);
    BOOST_TEST(!concurrent.native_points_on_producing_line(first, second));
    orientation(first, beyond, concurrent);

    using integer = std::int64_t;
    using int_point = ec::exact_point<integer>;
    using native_int = bg::model::d2::point_xy<integer>;
    native_int const mixed0(0, 0), mixed1(9007199254740993LL, 9007199254740992LL);
    BOOST_TEST(!crossing.native_points_on_producing_line(mixed0, mixed1));
    BOOST_TEST_EQ(ec::static_orientation<>(mixed0, mixed1, crossing), 1);
    ec::point_cache_scope<integer> integer_operation;
    integer const low = (std::numeric_limits<integer>::min)();
    integer const high = (std::numeric_limits<integer>::max)();
    int_point const i0(low, low), i1(high, high), i2(low, high), i3(high, low);
    auto const center = ec::construct_intersection(i0, i1, i2, i3);
    native_int const n0(low, low), n1(high, high), n2(high, high - 1);
    auto const integer_cache_size = ec::point_cache_detail::native_cache<integer>::size();
    BOOST_TEST_EQ(ec::static_orientation<integer>(n0, n1, center), 0);
    BOOST_TEST_EQ(ec::point_cache_detail::native_cache<integer>::size(), integer_cache_size);
    BOOST_TEST_EQ(ec::static_orientation<integer>(n0, n2, center), 1);
}

int main()
{
    ec::point_cache_scope<double> operation;
    segment_classification();
    box_access();
    native_support_orientation();
    using point = ec::exact_point<double>;
    point const a(0, 0), b(1, 1), c(0, 1), d(1, -1), e(0, 2);
    auto const p = ec::construct_intersection(a, b, c, d);
    auto const q = ec::construct_intersection(a, b, e, d);
    native_point const first(0, 0), second(1, 1), query(0.5, 0.25);
    orientation(first, second, query);
    orientation(first, second, p);
    orientation(first, p, query);
    orientation(p, first, second);
    orientation(first, p, q);
    orientation(p, q, query);
    orientation(p, q, a);

    ec::filtering_statistics() = {};
    BOOST_TEST_EQ(ec::static_direction<double>(first, second, query), -1);
    BOOST_TEST_EQ(ec::static_direction<double>(first, second, p), -1);
    BOOST_TEST(!ec::static_equal_points<double>(first, p));
    BOOST_TEST(ec::static_equal_points<double>(first, a));
    ec::point_compare<double, bg::strategy::compare::less> order;
    BOOST_TEST(order.apply(first, p));
    BOOST_TEST(order.apply(p, second));
    ec::winding<double>::state_type state;
    BOOST_TEST(ec::winding<double>::apply(query, first, p, state));
    BOOST_TEST_EQ(ec::filtering_statistics().native_point_constructions, 0u);
    return boost::report_errors();
}
