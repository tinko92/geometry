// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.

#include <boost/geometry.hpp>
#include <boost/geometry/strategies/relate/cartesian_exact.hpp>
#include <boost/core/lightweight_test.hpp>
#include <cstdint>

namespace bg = boost::geometry;

template <typename Bound, typename Working>
void check_exact_output(Bound const& bound, Working const& working, char const* expected_wkt)
{
    using rational = bg::strategies::detail::exact_cartesian::exact_number;
    using point = bg::model::d2::point_xy<rational>;
    using multi = bg::model::multi_polygon<bg::model::polygon<point>>;
    multi actual, expected;
    bound.convert_geometry(working, actual);
    bg::read_wkt(expected_wkt, expected);
    BOOST_TEST_EQ(bg::relation(actual, expected).str(), "2FFF1FFF2");
}

char const* const union_wkt = "MULTIPOLYGON(((0 0,0 6,1 5,7 5,7 -1,1 -1,1 0,0 0)))";
char const* const intersection_wkt = "MULTIPOLYGON(((1 0,1 5,6 0,1 0)))";
char const* const difference_wkt = "MULTIPOLYGON(((0 0,0 6,1 5,1 0,0 0)))";
char const* const symmetric_wkt = "MULTIPOLYGON(((0 0,0 6,1 5,1 0,0 0)),((1 -1,1 0,6 0,1 5,7 5,7 -1,1 -1)))";

template <typename C, typename Policy>
void check()
{
    using point = bg::model::d2::point_xy<C>;
    using polygon = bg::model::polygon<point>;
    using multi = bg::model::multi_polygon<polygon>;
    polygon a, b;
    bg::read_wkt("POLYGON((0 0,0 6,6 0,0 0))", a);
    bg::read_wkt("POLYGON((1 -1,1 5,7 5,7 -1,1 -1))", b);
    bg::strategies::relate::cartesian_exact<C, Policy> strategy;
    BOOST_TEST_EQ(bg::relation(a, b, strategy).str(), "212101212");
    multi u, i, d, s;
    bg::union_(a, b, u, strategy);
    bg::intersection(a, b, i, strategy);
    bg::difference(a, b, d, strategy);
    bg::sym_difference(a, b, s, strategy);
    BOOST_TEST_EQ(u.size(), 1u);
    BOOST_TEST_EQ(i.size(), 1u);
    BOOST_TEST_EQ(d.size(), 1u);
    BOOST_TEST_EQ(s.size(), 2u);
    if constexpr (!std::is_same<Policy, bg::strategies::relate::approximate_output>::value)
    {
        BOOST_TEST_EQ(bg::area(u), 41.5);
        BOOST_TEST_EQ(bg::area(i), 12.5);
        BOOST_TEST_EQ(bg::area(d), 5.5);
        BOOST_TEST_EQ(bg::area(s), 29.0);
    }

    strategy.with_sources(a, b, [&](auto const& bound)
    {
        using bound_type = std::decay_t<decltype(bound)>;
        using working = typename bound_type::template overlay_output_type<multi>;
        working first, second, combined, joined, common, symmetric;
        bg::union_(a, b, joined, bound);
        bg::intersection(a, b, common, bound);
        bg::sym_difference(a, b, symmetric, bound);
        bg::difference(a, b, first, bound);
        // Exact internal geometry, not rounded native coordinates, defines the
        // set result under every policy. Ordinary emission has no area promise.
        check_exact_output(bound, joined, union_wkt);
        check_exact_output(bound, common, intersection_wkt);
        check_exact_output(bound, first, difference_wkt);
        check_exact_output(bound, symmetric, symmetric_wkt);
        bg::difference(b, a, second, bound);
        bg::union_(first, second, combined, bound);
        check_exact_output(bound, combined, symmetric_wkt);
        multi exported;
        bound.append_output(combined, exported);
        if constexpr (!std::is_same<Policy, bg::strategies::relate::approximate_output>::value)
            BOOST_TEST_EQ(bg::area(exported), 29.0);
    });
}

template <typename C, typename Policy>
void check_linear()
{
    using point = bg::model::d2::point_xy<C>;
    using line = bg::model::linestring<point>;
    using lines = bg::model::multi_linestring<line>;
    using points = bg::model::multi_point<point>;
    using box = bg::model::box<point>;
    line const a{{0, 0}, {3, 3}}, b{{0, 3}, {3, 0}};
    bg::strategies::relate::cartesian_exact<C, Policy> strategy;
    points crossing;
    bg::intersection(a, b, crossing, strategy);
    BOOST_TEST_EQ(crossing.size(), 1u);
    BOOST_TEST_EQ(bg::get<0>(crossing.front()), C(1.5));
    BOOST_TEST_EQ(bg::relation(a, b, strategy).str(), "0F1FF0102");
    lines u, d, s;
    bg::union_(a, b, u, strategy);
    bg::difference(a, b, d, strategy);
    bg::sym_difference(a, b, s, strategy);
    BOOST_TEST(!u.empty());
    BOOST_TEST_EQ(d.size(), 1u);
    BOOST_TEST(!s.empty());
    box const bounds(point(1, 0), point(2, 3));
    lines clipped;
    bg::intersection(a, bounds, clipped, strategy);
    BOOST_TEST_EQ(clipped.size(), 1u);
    BOOST_TEST_EQ(bg::get<0>(clipped.front().front()), C(1));
    BOOST_TEST_EQ(bg::get<1>(clipped.front().back()), C(2));
    bg::model::segment<point> const s1(point(0, 0), point(3, 3));
    bg::model::segment<point> const s2(point(0, 3), point(3, 0));
    crossing.clear();
    bg::intersection(s1, s2, crossing, strategy);
    BOOST_TEST_EQ(crossing.size(), 1u);
    BOOST_TEST_EQ(bg::get<0>(crossing.front()), C(1.5));
}

template <bool Clockwise, bool Closed>
void check_ring_traits()
{
    using point = bg::model::d2::point_xy<double>;
    using polygon = bg::model::polygon<point, Clockwise, Closed>;
    polygon a, b;
    bg::read_wkt("POLYGON((0 0,0 6,6 0,0 0))", a);
    bg::read_wkt("POLYGON((1 -1,1 5,7 5,7 -1,1 -1))", b);
    bg::correct(a);
    bg::correct(b);
    bg::model::multi_polygon<polygon> result;
    bg::strategies::relate::cartesian_exact<double> strategy;
    bg::sym_difference(a, b, result, strategy);
    BOOST_TEST_EQ(result.size(), 2u);
    strategy.with_sources(a, b, [&](auto const& bound)
    {
        using native_multi = decltype(result);
        using working = typename std::decay_t<decltype(bound)>::template overlay_output_type<native_multi>;
        working exact;
        bg::sym_difference(a, b, exact, bound);
        check_exact_output(bound, exact, symmetric_wkt);
    });
}

int main()
{
    check<double, bg::strategies::relate::approximate_output>();
    check<double, bg::strategies::relate::eager_nearest_output>();
    check<double, bg::strategies::relate::deferred_nearest_output>();
    check<std::int64_t, bg::strategies::relate::approximate_output>();
    check_linear<double, bg::strategies::relate::approximate_output>();
    check_linear<double, bg::strategies::relate::eager_nearest_output>();
    check_linear<double, bg::strategies::relate::deferred_nearest_output>();
    check_linear<std::int64_t, bg::strategies::relate::approximate_output>();
    check_ring_traits<true, false>();
    check_ring_traits<false, true>();
    check_ring_traits<false, false>();
    return boost::report_errors();
}
