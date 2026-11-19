// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Licensed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt

#define BOOST_GEOMETRY_EXACT_FILTER_STATISTICS
#include <boost/geometry.hpp>
#include <boost/geometry/algorithms/detail/sections/sectionalize.hpp>
#include <boost/geometry/strategies/relate/cartesian_exact.hpp>
#include "exact_test_context.hpp"
#include <boost/core/lightweight_test.hpp>
#include <cstdlib>
#include <cstdint>
#include <cmath>
#include <limits>
#include <new>
#include <type_traits>

namespace bg = boost::geometry;

namespace
{
bool count_adapter_allocations = false;
std::size_t adapter_allocations = 0;
}

void* operator new(std::size_t size)
{
    if (count_adapter_allocations) ++adapter_allocations;
    if (void* result = std::malloc(size ? size : 1)) return result;
    throw std::bad_alloc();
}
void operator delete(void* value) noexcept { std::free(value); }
void* operator new[](std::size_t size) { return ::operator new(size); }
void operator delete[](void* value) noexcept { ::operator delete(value); }
#if defined(__cpp_sized_deallocation)
void operator delete(void* value, std::size_t) noexcept { ::operator delete(value); }
void operator delete[](void* value, std::size_t) noexcept { ::operator delete(value); }
#endif

template <typename Function>
bool without_adapter_materialization(Function function)
{
    namespace ec = bg::strategies::detail::exact_cartesian;
    ec::point_cache_scope operation;
    ec::filtering_statistics() = {};
    function(); // Permit indexed precision-cache growth before counting.
    prepare_exact_cache_repeat();
    adapter_allocations = 0;
    count_adapter_allocations = true;
    bool const result = function();
    count_adapter_allocations = false;
    BOOST_TEST_EQ(adapter_allocations,0u);
    BOOST_TEST_EQ(ec::filtering_statistics().native_rationals,0u);
    BOOST_TEST_EQ(ec::filtering_statistics().point_coordinate_values,0u);
    BOOST_TEST_EQ(ec::filtering_statistics().homogeneous_values,0u);
    return result;
}

template <typename Geometry>
auto validation_area(Geometry const& geometry)
{
    using rational = bg::strategies::detail::exact_cartesian::exact_number;
    return bg::area(geometry, bg::strategies::relate::cartesian<rational>());
}

template <typename Geometry1, typename Geometry2, typename Strategy>
void test_native_predicate_case(Geometry1 const& a, Geometry2 const& b, Strategy const& strategy,
    bool covered_by, bool crosses, bool disjoint, bool intersects,
    bool overlaps, bool touches, bool relate, bool within)
{
    BOOST_TEST_EQ(bg::covered_by(a, b, strategy), covered_by);
    BOOST_TEST_EQ(bg::crosses(a, b, strategy), crosses);
    BOOST_TEST_EQ(bg::disjoint(a, b, strategy), disjoint);
    BOOST_TEST_EQ(bg::intersects(a, b, strategy), intersects);
    BOOST_TEST_EQ(bg::overlaps(a, b, strategy), overlaps);
    BOOST_TEST_EQ(bg::touches(a, b, strategy), touches);
    BOOST_TEST_EQ(bg::relate(a, b, bg::de9im::mask("F0F******"), strategy), relate);
    BOOST_TEST_EQ(bg::within(a, b, strategy), within);
}

template <typename Coordinate>
void test_native_predicates()
{
    using point = bg::model::d2::point_xy<Coordinate>;
    using ring = bg::model::ring<point>;
    using line = bg::model::linestring<point>;
    bg::strategies::relate::cartesian_exact<Coordinate> strategy;
    ring const rectangle{{2, 2}, {2, 4}, {4, 4}, {4, 2}, {2, 2}};
    ring const triangle{{1, 2}, {4, 5}, {5, 4}, {1, 2}};
    ring const touching{{4, 2}, {4, 4}, {6, 4}, {6, 2}, {4, 2}};
    ring const separated{{5, 2}, {5, 4}, {7, 4}, {7, 2}, {5, 2}};
    line const diagonal{{0, 0}, {6, 6}};
    line const horizontal{{0, 3}, {6, 3}};
    test_native_predicate_case(rectangle, triangle, strategy,
        false, false, false, true, true, false, false, false);
    test_native_predicate_case(rectangle.front(), triangle, strategy,
        false, false, true, false, false, false, false, false);
    test_native_predicate_case(rectangle.front(), diagonal, strategy,
        true, false, false, true, false, false, false, true);
    test_native_predicate_case(horizontal, diagonal, strategy,
        false, true, false, true, false, false, false, false);
    test_native_predicate_case(rectangle, touching, strategy,
        false, false, false, true, false, true, false, false);
    test_native_predicate_case(rectangle, separated, strategy,
        false, false, true, false, false, false, false, false);
    test_native_predicate_case(point(2, 3), rectangle, strategy,
        true, false, false, true, false, true, true, false);
    test_native_predicate_case(point(3, 3), rectangle, strategy,
        true, false, false, true, false, false, false, true);
    test_native_predicate_case(diagonal, rectangle, strategy,
        false, true, false, true, false, false, false, false);

    bg::model::polygon<point> polygon1, polygon2;
    bg::exterior_ring(polygon1) = rectangle;
    bg::exterior_ring(polygon2) = triangle;
    test_native_predicate_case(polygon1, polygon2, strategy,
        false, false, false, true, true, false, false, false);

    using segment = bg::model::segment<point>;
    segment const first(point(0, 0), point(6, 6));
    segment const crossing(point(0, 6), point(6, 0));
    segment const remote(point(0, 7), point(1, 7));
    BOOST_TEST(!bg::disjoint(first, crossing, strategy));
    BOOST_TEST(bg::intersects(first, crossing, strategy));
    BOOST_TEST(bg::disjoint(first, remote, strategy));
    BOOST_TEST(!bg::intersects(first, remote, strategy));
}

template <typename Coordinate>
void test_geometry()
{
    using point = bg::model::d2::point_xy<Coordinate>;
    using polygon = bg::model::polygon<point>;
    bg::strategies::relate::cartesian_exact<Coordinate> strategy;
    polygon a, b;
    bg::read_wkt("POLYGON((0 0,0 4,4 4,4 0,0 0))", a);
    bg::read_wkt("POLYGON((2 2,2 6,6 6,6 2,2 2))", b);
    bg::model::multi_polygon<polygon> out;
    bg::union_(a, b, out, strategy);
    BOOST_TEST_EQ(bg::area(out), 28);
    out.clear();
    bg::intersection(a, b, out, strategy);
    BOOST_TEST_EQ(bg::area(out), 4);
    out.clear();
    bg::difference(a, b, out, strategy);
    BOOST_TEST_EQ(bg::area(out), 12);
    out.clear();
    bg::sym_difference(a, b, out, strategy);
    BOOST_TEST_EQ(bg::area(out), 24);
    BOOST_TEST_EQ(bg::relation(a, b, strategy).str(), "212101212");

    using line = bg::model::linestring<point>;
    line x, y;
    bg::read_wkt("LINESTRING(0 0,4 4)", x);
    bg::read_wkt("LINESTRING(0 4,4 0)", y);
    BOOST_TEST_EQ(bg::relation(x, y, strategy).str(), "0F1FF0102");
    BOOST_TEST_EQ(bg::relation(x, a, strategy).str(), "1FFF0F212");
    BOOST_TEST_EQ(bg::relation(point(2, 2), a, strategy).str(), "0FFFFF212");

    using multi_line = bg::model::multi_linestring<line>;
    multi_line lines;
    bg::intersection(x, a, lines, strategy);
    BOOST_TEST_EQ(lines.size(), 1u);
    lines.clear();
    bg::difference(x, b, lines, strategy);
    BOOST_TEST_EQ(lines.size(), 1u);

    bg::model::box<point> const box(point(1, 1), point(3, 3));
    lines.clear();
    bg::intersection(x, box, lines, strategy);
    BOOST_TEST_EQ(lines.size(), 1u);
    bg::model::box<point> clipped;
    BOOST_TEST(bg::intersection(bg::model::box<point>(point(0, 0), point(2, 2)),
        box, clipped, strategy));
    BOOST_TEST_EQ((bg::get<bg::min_corner, 0>(clipped)), 1);
    BOOST_TEST_EQ((bg::get<bg::max_corner, 1>(clipped)), 2);

    using multi_point = bg::model::multi_point<point>;
    using multi_polygon = bg::model::multi_polygon<polygon>;
    std::tuple<multi_point, multi_line, multi_polygon> tuple;
    bg::intersection(x, y, tuple, strategy);
    BOOST_TEST_EQ(std::get<0>(tuple).size(), 1u);
    BOOST_TEST(std::get<1>(tuple).empty());
    BOOST_TEST(std::get<2>(tuple).empty());
    std::get<0>(tuple).clear();
    bg::intersection(x, a, tuple, strategy);
    BOOST_TEST_EQ(std::get<1>(tuple).size(), 1u);
}

void test_small_features()
{
    bg::strategies::detail::exact_cartesian::point_cache_scope<double> operation;
    using point = bg::model::d2::point_xy<double>;
    using polygon = bg::model::polygon<point>;
    bg::strategies::relate::cartesian_exact<double> strategy;
    double const height = std::ldexp(1.0, -60);
    polygon triangle;
    bg::exterior_ring(triangle) = {{0, 0}, {0, height}, {1, 0}, {0, 0}};
    bg::model::multi_polygon<polygon> out;
    bg::intersection(triangle, triangle, out, strategy);
    BOOST_TEST_EQ(out.size(), 1u);
    BOOST_TEST_EQ(bg::relation(triangle, triangle, strategy).str(), "2FFF1FFF2");
    BOOST_TEST_EQ(bg::relation(point(0, 0), point(height, 0), strategy).str(), "FF0FFF0F2");

    using line = bg::model::linestring<point>;
    line const a{{0, 0}, {1, 1}};
    line const b{{0, 1}, {1, -1}};
    bg::model::multi_point<point> points;
    bg::intersection(a, b, points, strategy);
    BOOST_TEST_EQ(points.size(), 1u);
    if (!points.empty())
    {
        BOOST_TEST_EQ(bg::get<0>(points.front()), 1.0 / 3.0);
        BOOST_TEST_EQ(bg::get<1>(points.front()), 1.0 / 3.0);
    }

    using rational_coordinate = bg::strategies::detail::exact_cartesian::exact_number;
    bg::model::multi_point<bg::model::d2::point_xy<rational_coordinate>> rational_points;
    bg::intersection(a, b, rational_points, strategy);
    BOOST_TEST_EQ(rational_points.size(), 1u);
    if (!rational_points.empty())
    {
        BOOST_TEST_EQ(bg::get<0>(rational_points.front()), rational_coordinate(1) / 3);
        BOOST_TEST_EQ(bg::get<1>(rational_points.front()), rational_coordinate(1) / 3);
    }

    polygon notched, bridge;
    bg::exterior_ring(notched) = {{0, 0}, {0, 3}, {3, 3}, {3, 1.25},
        {1, 1}, {3, 0.75}, {3, 0}, {0, 0}};
    double const left = std::nextafter(1.0, 2.0);
    bg::exterior_ring(bridge) = {{left, 0}, {left, 3}, {4, 3}, {4, 0}, {left, 0}};
    using implicit = bg::strategies::detail::exact_cartesian::exact_point<double>;
    using rational = bg::strategies::detail::exact_cartesian::exact_number;
    bg::model::multi_polygon<bg::model::polygon<implicit>> exact;
    bg::union_(notched, bridge, exact, strategy);
    BOOST_TEST_EQ(exact.size(), 1u);
    BOOST_TEST_EQ(bg::num_interior_rings(exact), 1u);
    rational const hole_area = rational(1) / rational(boost::multiprecision::cpp_int(1) << 107);
    BOOST_TEST_EQ(validation_area(exact), rational(12) - hole_area);

    // The true hole survives topology processing even though its two constructed
    // vertices round to the same output point, leaving a collapsed output ring.
    out.clear();
    bg::union_(notched, bridge, out, strategy);
    BOOST_TEST_EQ(bg::num_interior_rings(out), 1u);
}

void test_large_integers()
{
    using coordinate = std::int64_t;
    namespace ec = bg::strategies::detail::exact_cartesian;
    ec::point_cache_scope<coordinate> operation;
    using point = bg::model::d2::point_xy<coordinate>;
    using polygon = bg::model::polygon<point>;
    bg::strategies::relate::cartesian_exact<coordinate> strategy;
    coordinate const origin = coordinate(1) << 62;
    polygon a, b;
    bg::exterior_ring(a) = {{origin, origin}, {origin, origin + 8},
        {origin + 8, origin + 8}, {origin + 8, origin}, {origin, origin}};
    bg::exterior_ring(b) = {{origin + 4, origin + 4}, {origin + 4, origin + 12},
        {origin + 12, origin + 12}, {origin + 12, origin + 4}, {origin + 4, origin + 4}};
    bg::model::multi_polygon<polygon> out;
    bg::union_(a, b, out, strategy);
    BOOST_TEST_EQ(out.size(), 1u);
    bg::model::multi_polygon<bg::model::polygon<ec::exact_point<coordinate>>> exact;
    bg::union_(a, b, exact, strategy);
    BOOST_TEST_EQ(validation_area(exact), 112);
    // Binary64 approximations need not preserve unit-scale increments here.
    // Export retains the selected structure, not an integer area guarantee.
    BOOST_TEST_EQ(bg::num_points(out), bg::num_points(exact));
    BOOST_TEST_EQ(bg::relation(a, b, strategy).str(), "212101212");

    static_assert(!ec::lossless_coordinate_conversion<coordinate, double>::value,
                  "Lossy mixed integer/floating inputs are outside the exact strategy contract");
    using narrow_point = bg::model::d2::point_xy<std::int32_t>;
    using double_point = bg::model::d2::point_xy<double>;
    bg::strategies::relate::cartesian_exact<double> mixed;
    BOOST_TEST_EQ(bg::relation(narrow_point(123456789, 0),
        double_point(123456789.0, 0), mixed).str(), "0FFFFFFF2");
}

template <typename Coordinate>
void test_box_failure_paths()
{
    using point = bg::model::d2::point_xy<Coordinate>;
    using box = bg::model::box<point>;
    bg::strategies::relate::cartesian_exact<Coordinate> strategy;
    using exact_box = bg::model::box<bg::strategies::detail::exact_cartesian::exact_point<Coordinate>>;
    using search_box = bg::strategies::detail::exact_cartesian::bounding_box;
    static_assert(std::is_same<typename decltype(strategy)::template overlay_output_type<box>,
        box>::value, "Public box output must retain partial writes");
    static_assert(std::is_same<typename decltype(strategy)::template overlay_geometry_type<box>,
        exact_box>::value, "Intermediate box geometry must remain exact");
    static_assert(std::is_same<typename decltype(strategy)::template box_type<point>,
        search_box>::value, "Candidate searches use conservative native enclosures");
    box const first(point(0, 0), point(2, 2));
    box const others[] = {
        box(point(3, 0), point(4, 2)),
        box(point(-2, 0), point(-1, 2)),
        box(point(1, 3), point(3, 4)),
        box(point(1, -2), point(3, -1)),
        box(point(1, 1), point(3, 3)),
        box(point(2, 1), point(3, 3))
    };
    for (auto const& second : others)
    {
        box expected(point(-9, -8), point(9, 8));
        box actual = expected;
        bool const expected_result = bg::intersection(first, second, expected);
        bool const actual_result = bg::intersection(first, second, actual, strategy);
        BOOST_TEST_EQ(actual_result, expected_result);
        BOOST_TEST_EQ((bg::get<bg::min_corner, 0>(actual)),
            (bg::get<bg::min_corner, 0>(expected)));
        BOOST_TEST_EQ((bg::get<bg::max_corner, 0>(actual)),
            (bg::get<bg::max_corner, 0>(expected)));
        BOOST_TEST_EQ((bg::get<bg::min_corner, 1>(actual)),
            (bg::get<bg::min_corner, 1>(expected)));
        BOOST_TEST_EQ((bg::get<bg::max_corner, 1>(actual)),
            (bg::get<bg::max_corner, 1>(expected)));
    }

    using polygon = bg::model::polygon<point>;
    polygon a, b;
    bg::read_wkt("POLYGON((0 0,0 2,2 2,2 0,0 0))", a);
    bg::read_wkt("POLYGON((1 1,1 3,3 3,3 1,1 1))", b);
    bg::model::multi_polygon<polygon> output;
    output.push_back(a);
    bg::intersection(a, b, output, strategy);
    BOOST_TEST_EQ(output.size(), 2u);
    BOOST_TEST_EQ(validation_area(output.front()), 4);
    BOOST_TEST_EQ(validation_area(output.back()), 1);
}

void test_box_output_state_and_rounding()
{
    using point = bg::model::d2::point_xy<double>;
    using box = bg::model::box<point>;
    bg::strategies::relate::cartesian_exact<double> strategy;
    box const first(point(0, 0), point(2, 2));
    double const nan = std::numeric_limits<double>::quiet_NaN();
    box output(point(nan, nan), point(nan, nan));
    BOOST_TEST(!bg::intersection(first, box(point(3, 0), point(4, 2)), output, strategy));
    BOOST_TEST(std::isnan((bg::get<bg::min_corner, 0>(output))));
    BOOST_TEST(std::isnan((bg::get<bg::max_corner, 0>(output))));
    BOOST_TEST(std::isnan((bg::get<bg::min_corner, 1>(output))));
    BOOST_TEST(std::isnan((bg::get<bg::max_corner, 1>(output))));
    BOOST_TEST(!bg::intersection(first, box(point(1, 3), point(3, 4)), output, strategy));
    BOOST_TEST_EQ((bg::get<bg::min_corner, 0>(output)), 1);
    BOOST_TEST_EQ((bg::get<bg::max_corner, 0>(output)), 2);
    BOOST_TEST(std::isnan((bg::get<bg::min_corner, 1>(output))));
    BOOST_TEST(std::isnan((bg::get<bg::max_corner, 1>(output))));
    BOOST_TEST(bg::intersection(first, box(point(1, 1), point(3, 3)), output, strategy));
    BOOST_TEST_EQ((bg::get<bg::min_corner, 0>(output)), 1);
    BOOST_TEST_EQ((bg::get<bg::max_corner, 0>(output)), 2);
    BOOST_TEST_EQ((bg::get<bg::min_corner, 1>(output)), 1);
    BOOST_TEST_EQ((bg::get<bg::max_corner, 1>(output)), 2);

    using integer_point = bg::model::d2::point_xy<std::int64_t>;
    using integer_box = bg::model::box<integer_point>;
    integer_box rounded;
    BOOST_TEST(bg::intersection(first, box(point(0.5, 0.5), point(1.5, 1.5)), rounded, strategy));
    BOOST_TEST_EQ((bg::get<bg::min_corner, 0>(rounded)), 0);
    BOOST_TEST_EQ((bg::get<bg::max_corner, 0>(rounded)), 1);
    BOOST_TEST_EQ((bg::get<bg::min_corner, 1>(rounded)), 0);
    BOOST_TEST_EQ((bg::get<bg::max_corner, 1>(rounded)), 1);

    using narrow_point = bg::model::d2::point_xy<std::int32_t>;
    using narrow_box = bg::model::box<narrow_point>;
    std::int32_t const n = 123456789;
    narrow_box const integer(narrow_point(n + 1, 0), narrow_point(n + 3, 2));
    box const disjoint(point(double(n), 0), point(double(n), 2));
    rounded = integer_box(integer_point(-9, -8), integer_point(9, 8));
    BOOST_TEST(!bg::intersection(integer, disjoint, rounded, strategy));
    BOOST_TEST_EQ((bg::get<bg::min_corner, 0>(rounded)), -9);
    BOOST_TEST_EQ((bg::get<bg::max_corner, 0>(rounded)), 9);
    BOOST_TEST_EQ((bg::get<bg::min_corner, 1>(rounded)), -8);
    BOOST_TEST_EQ((bg::get<bg::max_corner, 1>(rounded)), 8);

    using rational = bg::strategies::detail::exact_cartesian::exact_number;
    bg::model::box<bg::model::d2::point_xy<rational>> exact;
    box const overlap(point(double(n), 1), point(double(n + 2), 3));
    BOOST_TEST(bg::intersection(integer, overlap, exact, strategy));
    BOOST_TEST_EQ((bg::get<bg::min_corner, 0>(exact)), rational(n + 1));
    BOOST_TEST_EQ((bg::get<bg::max_corner, 0>(exact)), rational(n + 2));
    BOOST_TEST_EQ((bg::get<bg::min_corner, 1>(exact)), 1);
    BOOST_TEST_EQ((bg::get<bg::max_corner, 1>(exact)), 2);
}

template <typename Coordinate>
void test_native_box_predicates()
{
    using point = bg::model::d2::point_xy<Coordinate>;
    using box = bg::model::box<point>;
    bg::strategies::relate::cartesian_exact<Coordinate> strategy;
    box const outer(point(0, 0), point(4, 4));
    point const points[] = {point(2, 2), point(0, 2), point(4, 4), point(5, 2)};
    for (auto const& p : points)
    {
        BOOST_TEST_EQ(bg::within(p, outer, strategy), bg::within(p, outer));
        BOOST_TEST_EQ(bg::covered_by(p, outer, strategy), bg::covered_by(p, outer));
        BOOST_TEST_EQ(bg::disjoint(p, outer, strategy), bg::disjoint(p, outer));
    }
    box const boxes[] = {outer, box(point(1, 1), point(3, 3)),
        box(point(1, 1), point(1, 3)), box(point(4, 1), point(5, 3)),
        box(point(5, 1), point(6, 3))};
    for (auto const& first : boxes)
    {
        for (auto const& second : boxes)
        {
            BOOST_TEST_EQ(bg::within(first, second, strategy), bg::within(first, second));
            BOOST_TEST_EQ(bg::covered_by(first, second, strategy), bg::covered_by(first, second));
            BOOST_TEST_EQ(bg::disjoint(first, second, strategy), bg::disjoint(first, second));
            BOOST_TEST_EQ(bg::touches(first, second, strategy), bg::touches(first, second));
            BOOST_TEST_EQ(bg::overlaps(first, second, strategy), bg::overlaps(first, second));
        }
    }
}

void test_mixed_box_predicates()
{
    using integer_point = bg::model::d2::point_xy<std::int32_t>;
    using integer_box = bg::model::box<integer_point>;
    using floating_point = bg::model::d2::point_xy<double>;
    using floating_box = bg::model::box<floating_point>;
    bg::strategies::relate::cartesian_exact<double> strategy;
    std::int32_t const n = 123456789;
    integer_point const p(n + 1, 1);
    floating_box const outer(floating_point(double(n), 0), floating_point(double(n + 2), 2));
    floating_box const left(floating_point(double(n - 2), 0), floating_point(double(n), 2));
    BOOST_TEST(bg::within(p, outer, strategy));
    BOOST_TEST(bg::covered_by(p, outer, strategy));
    BOOST_TEST(!bg::disjoint(p, outer, strategy));
    BOOST_TEST(!bg::within(p, left, strategy));
    BOOST_TEST(!bg::covered_by(p, left, strategy));
    BOOST_TEST(bg::disjoint(p, left, strategy));
    BOOST_TEST(bg::disjoint(left, p, strategy));

    bg::model::multi_point<integer_point> points;
    points.push_back(p);
    BOOST_TEST(bg::covered_by(points, outer, strategy));
    BOOST_TEST(!bg::covered_by(points, left, strategy));

    integer_box const right(integer_point(n + 1, 0), integer_point(n + 2, 2));
    BOOST_TEST(bg::within(right, outer, strategy));
    BOOST_TEST(bg::covered_by(right, outer, strategy));
    BOOST_TEST(bg::disjoint(right, left, strategy));
    BOOST_TEST(bg::disjoint(left, right, strategy));
    BOOST_TEST(!bg::touches(right, left, strategy));
    BOOST_TEST(!bg::touches(left, right, strategy));
    floating_box const touching(floating_point(double(n + 2), 0),
        floating_point(double(n + 4), 2));
    BOOST_TEST(bg::touches(right, touching, strategy));
    BOOST_TEST(bg::touches(touching, right, strategy));
    integer_box const overlap_integer(integer_point(n + 3, 0), integer_point(n + 5, 2));
    floating_box const overlap_floating(floating_point(double(n + 2), 1),
        floating_point(double(n + 4), 3));
    BOOST_TEST(bg::overlaps(overlap_integer, overlap_floating, strategy));
    BOOST_TEST(bg::overlaps(overlap_floating, overlap_integer, strategy));
    BOOST_TEST(!bg::overlaps(right, outer, strategy));
    BOOST_TEST(!bg::overlaps(outer, right, strategy));

    integer_box const crossing(integer_point(n - 1, 0), integer_point(n + 1, 2));
    BOOST_TEST(!bg::within(crossing, left, strategy));
    BOOST_TEST(!bg::covered_by(crossing, left, strategy));
    floating_box const crossing_double(floating_point(double(n), 0.5),
        floating_point(double(n + 2), 1.5));
    BOOST_TEST(!bg::within(crossing_double, right, strategy));
    BOOST_TEST(!bg::covered_by(crossing_double, right, strategy));

    bg::model::segment<integer_point> const segment(p, integer_point(n + 2, 1));
    auto const point_box = strategy.disjoint(segment, left).get_disjoint_point_box_strategy();
    BOOST_TEST(!point_box.apply(p, left));
    BOOST_TEST(point_box.apply(p, outer));
    BOOST_TEST(bg::disjoint(segment, left, strategy));
    BOOST_TEST(!bg::disjoint(segment, outer, strategy));

    using rational = bg::strategies::detail::exact_cartesian::exact_number;
    bg::model::d2::point_xy<rational> const rational_point(rational(n + 1), 1);
    BOOST_TEST(bg::within(rational_point, outer, strategy));
    BOOST_TEST(!bg::covered_by(rational_point, left, strategy));

    using unsigned_point = bg::model::d2::point_xy<std::uint32_t>;
    bg::model::box<unsigned_point> const positive(unsigned_point(0, 0), unsigned_point(4, 4));
    integer_box const negative(integer_point(-2, 0), integer_point(-1, 4));
    BOOST_TEST(bg::disjoint(negative, positive, strategy));
    BOOST_TEST(!bg::covered_by(integer_point(-1, 1), positive, strategy));
}

void test_wide_integer_box_union()
{
    bg::strategies::detail::exact_cartesian::point_cache_scope<int> operation;
    using rational = bg::strategies::detail::exact_cartesian::exact_number;
    using point = bg::model::d2::point_xy<int>;
    using box = bg::model::box<point>;
    bg::strategies::relate::cartesian_exact<int> strategy;
    box const wide(point((std::numeric_limits<int>::min)(), 0),
        point((std::numeric_limits<int>::max)(), 1));
    rational const expected(4294967295LL);

    bg::model::polygon<point> disjoint;
    bg::read_wkt("POLYGON((0 2,0 3,1 3,1 2,0 2))", disjoint);
    using implicit = bg::strategies::detail::exact_cartesian::exact_point<int>;
    bg::model::multi_polygon<bg::model::polygon<implicit>> output;
    bg::union_(wide, disjoint, output, strategy);
    BOOST_TEST_EQ(output.size(), 2u);
    BOOST_TEST_EQ(validation_area(output), expected + 1);
    output.clear();
    bg::union_(disjoint, wide, output, strategy);
    BOOST_TEST_EQ(output.size(), 2u);
    BOOST_TEST_EQ(validation_area(output), expected + 1);
}

void test_retained_section_points()
{
    namespace ec = bg::strategies::detail::exact_cartesian;
    ec::point_cache_scope operation;
    using point = ec::exact_point<double>;
    using rational = ec::exact_number;
    using strategy_type = bg::strategies::relate::cartesian_exact<double>;
    using search_box = bg::strategies::relate::services::box_type<strategy_type,point>::type;
    strategy_type const strategy;
    double const delta = std::ldexp(1.0,-54);
    // Intersect y=x-1 with y=delta: the exact x=1+delta rounds to 1.
    point const shifted = ec::construct_intersection(point(0,-1),point(2,1),
                                                     point(0,delta),point(2,delta));
    point const lower = point::splice_coordinates(shifted,0,point(0,0),1);
    point const upper = point::splice_coordinates(shifted,0,point(0,1),1);
    BOOST_TEST(bg::get<0>(lower) == rational(1)+rational(delta));
    bg::model::d2::point_xy<double> rounded;
    ec::assign_output(lower,rounded);
    BOOST_TEST_EQ(bg::get<0>(rounded),1.0);

    bg::model::polygon<point> polygon;
    auto& ring = bg::exterior_ring(polygon);
    ring = {{1,0},{1,1},upper,lower,{1,0}};
    bg::sections<search_box,2> sections;
    bg::sectionalize<false,std::integer_sequence<std::size_t,0,1>>(polygon,sections,strategy);
    BOOST_TEST_EQ(sections.size(),4u);
    int const directions[4][2] = {{0,1},{1,0},{0,-1},{-1,0}};
    std::size_t duplicate_count=0,segment_count=0;
    for(std::size_t i=0;i<sections.size();++i)
    {
        auto const& section = sections[i];
        duplicate_count+=section.duplicate ? 1u : 0u;
        segment_count+=section.count;
        BOOST_TEST(!section.duplicate);
        BOOST_TEST_EQ(section.count,1u);
        if(i<4)
        {
            BOOST_TEST_EQ(section.begin_index,static_cast<bg::signed_size_type>(i));
            BOOST_TEST_EQ(section.end_index,static_cast<bg::signed_size_type>(i+1));
            BOOST_TEST_EQ(section.directions[0],directions[i][0]);
            BOOST_TEST_EQ(section.directions[1],directions[i][1]);
        }
        for(auto index=section.begin_index;index<=section.end_index;++index)
        {
            auto const& vertex = ring[static_cast<std::size_t>(index)];
            auto const x=bg::get<0>(vertex),y=bg::get<1>(vertex);
            auto const& box=section.bounding_box;
            BOOST_TEST((rational(bg::get<bg::min_corner,0>(box))<=x));
            BOOST_TEST((rational(bg::get<bg::max_corner,0>(box))>=x));
            BOOST_TEST((rational(bg::get<bg::min_corner,1>(box))<=y));
            BOOST_TEST((rational(bg::get<bg::max_corner,1>(box))>=y));
        }
    }
    BOOST_TEST_EQ(duplicate_count,0u);
    BOOST_TEST_EQ(segment_count,4u);
}

void test_implicit_segment_adapters()
{
    namespace ec = bg::strategies::detail::exact_cartesian;
    using point = ec::exact_point<double>;
    using rational = ec::exact_number;
    using qpoint = bg::model::d2::point_xy<rational>;
    bg::strategies::relate::cartesian_exact<double> const strategy;
    bg::strategies::relate::cartesian<rational> const validation;
    auto third=[](double x,double y)
    {
        return ec::construct_intersection(point(x,y),point(x+1,y+1),
                                          point(x,y+1),point(x+1,y-1));
    };
    auto lift=[](point const& p) { return qpoint(bg::get<0>(p),bg::get<1>(p)); };
    point const a=third(0,0),b=third(2,2),middle=third(1,1),off=third(1,2);
    bg::model::segment<point> const segment(a,b),same(a,b),reversed(b,a),different(a,off),degenerate(a,a);
    bg::model::referring_segment<point const> const reference(a,b),reverse_reference(b,a);
    bg::model::segment<qpoint> const oracle(lift(a),lift(b)),qsame(lift(a),lift(b)),
        qreversed(lift(b),lift(a)),qdifferent(lift(a),lift(off)),qdegenerate(lift(a),lift(a));
    auto const qa=lift(a),qmiddle=lift(middle),qoff=lift(off);
    BOOST_TEST_EQ(without_adapter_materialization([&]{return bg::within(middle,segment,strategy);}),
                  bg::within(qmiddle,oracle,validation));
    BOOST_TEST_EQ(without_adapter_materialization([&]{return bg::within(middle,reference,strategy);}),
                  bg::within(qmiddle,oracle,validation));
    BOOST_TEST_EQ(without_adapter_materialization([&]{return bg::within(a,segment,strategy);}),
                  bg::within(qa,oracle,validation));
    BOOST_TEST_EQ(without_adapter_materialization([&]{return bg::covered_by(a,reverse_reference,strategy);}),
                  bg::covered_by(qa,qreversed,validation));
    BOOST_TEST_EQ(without_adapter_materialization([&]{return bg::disjoint(off,segment,strategy);}),
                  bg::disjoint(qoff,oracle,validation));
    BOOST_TEST_EQ(without_adapter_materialization([&]{return bg::within(a,degenerate,strategy);}),
                  bg::within(qa,qdegenerate,validation));
    BOOST_TEST_EQ(without_adapter_materialization([&]{return bg::covered_by(a,degenerate,strategy);}),
                  bg::covered_by(qa,qdegenerate,validation));
    BOOST_TEST_EQ(without_adapter_materialization([&]{return bg::equals(segment,same,strategy);}),
                  bg::equals(oracle,qsame,validation));
    BOOST_TEST_EQ(without_adapter_materialization([&]{return bg::equals(segment,reversed,strategy);}),
                  bg::equals(oracle,qreversed,validation));
    BOOST_TEST_EQ(without_adapter_materialization([&]{return bg::equals(reference,reverse_reference,strategy);}),
                  bg::equals(oracle,qreversed,validation));
    BOOST_TEST_EQ(without_adapter_materialization([&]{return bg::equals(segment,different,strategy);}),
                  bg::equals(oracle,qdifferent,validation));
    BOOST_TEST_EQ(without_adapter_materialization([&]{return bg::equals(degenerate,degenerate,strategy);}),
                  bg::equals(qdegenerate,qdegenerate,validation));

    bg::model::polygon<point> inside,outside;
    bg::exterior_ring(inside)={{-1,-1},{-1,3},{3,3},{3,-1},{-1,-1}};
    bg::exterior_ring(outside)={{4,4},{4,5},{5,5},{5,4},{4,4}};
    bg::model::polygon<qpoint> qinside,qoutside;
    bg::convert(inside,qinside);
    bg::convert(outside,qoutside);
    for(bool remote:{false,true})
    {
        auto const& polygon=remote ? outside : inside;
        auto const& qpolygon=remote ? qoutside : qinside;
        bool const expected=bg::disjoint(oracle,qpolygon,validation);
        BOOST_TEST_EQ(without_adapter_materialization([&]{return bg::disjoint(segment,polygon,strategy);}),expected);
        BOOST_TEST_EQ(without_adapter_materialization([&]{return bg::disjoint(reverse_reference,polygon,strategy);}),expected);
        BOOST_TEST_EQ(without_adapter_materialization([&]{return bg::disjoint(segment,bg::exterior_ring(polygon),strategy);}),expected);
    }
}

int main()
{
    test_native_predicates<double>();
    test_native_predicates<float>();
    test_native_predicates<std::int64_t>();
    test_geometry<double>();
    test_geometry<float>();
    test_geometry<std::int64_t>();
    test_small_features();
    test_large_integers();
    test_box_failure_paths<double>();
    test_box_failure_paths<float>();
    test_box_failure_paths<std::int64_t>();
    test_box_output_state_and_rounding();
    test_native_box_predicates<double>();
    test_native_box_predicates<float>();
    test_native_box_predicates<std::int64_t>();
    test_mixed_box_predicates();
    test_wide_integer_box_union();
    test_retained_section_points();
    test_implicit_segment_adapters();
    return boost::report_errors();
}
