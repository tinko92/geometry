// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.

#include <boost/geometry/algorithms/prepare_rational_polygon.hpp>
#include <boost/geometry/algorithms/prepare_filtered_decimal_polygon.hpp>

#include <boost/core/lightweight_test.hpp>
#include <boost/geometry/algorithms/area.hpp>
#include <boost/geometry/algorithms/difference.hpp>
#include <boost/geometry/algorithms/intersection.hpp>
#include <boost/geometry/algorithms/relate.hpp>
#include <boost/geometry/algorithms/relation.hpp>
#include <boost/geometry/algorithms/sym_difference.hpp>
#include <boost/geometry/algorithms/union.hpp>
#include <boost/geometry/io/wkt/read.hpp>
#include <boost/geometry/strategies/relate/cartesian_decimal.hpp>
#include <boost/geometry/strategies/relate/cartesian_exact.hpp>

#include <chrono>
#include <cstddef>
#include <iostream>
#include <string>
#include <utility>

namespace bg = boost::geometry;
namespace fd = bg::detail::filtered_decimal;
namespace ec = bg::strategies::detail::exact_cartesian;
using rational = boost::multiprecision::cpp_rational;
using boost::multiprecision::cpp_int;
using decimal = bg::exact_coordinate<std::string>;

template <typename Coordinate>
using point = bg::model::point<Coordinate, 2, bg::cs::cartesian>;
template <typename Coordinate>
using polygon = bg::model::polygon<point<Coordinate>, false>;
using rational_polygons = bg::model::multi_polygon<polygon<rational>>;

bg::strategies::relate::cartesian<rational> const rational_strategy;
bg::strategies::relate::cartesian_decimal const filtered_strategy;
bg::strategies::relate::cartesian_exact<double> const binary_strategy;

template <typename Geometry>
Geometry read(char const* wkt)
{
    Geometry result;
    bg::read_wkt(wkt, result);
    return result;
}

auto prepare(char const* wkt)
{
    return bg::prepare_rational_polygon(read<polygon<decimal>>(wkt));
}

rational area(rational_polygons const& polygons)
{
    return bg::area(polygons, rational_strategy);
}

// The original CCW triangle pair from the local #663 fixtures. Their exact
// decimal intersection is tiny but is not empty, nor is it the binary set.
char const issue_663_a[] =
    "POLYGON((140.10975222826562 135.99989776999999,140.10975222826562 126.97244753050222,"
    "140.16255693232961 135.99989776999999,140.10975222826562 135.99989776999999))";
char const issue_663_b[] =
    "POLYGON((151.99995129999999 123.99989777000000,139.99994018999999 126.99990054000000,"
    "129.99995129999999 124.99989777000000,151.99995129999999 123.99989777000000))";

void issue_663()
{
    auto const a = prepare(issue_663_a);
    auto const b = prepare(issue_663_b);
    rational_polygons out;
    BOOST_TEST(bg::intersection(a, b, out, rational_strategy));
    BOOST_TEST_EQ(out.size(), 1u);

    // Independent triangle-area oracle: its vertical edge has length h,
    // and the two sloping edges meet h / (slope_a - slope_b) to the right.
    auto const expected_area = [](auto const& first, auto const& second, rational const& h)
    {
        rational const x = bg::get<0>(first.outer()[1]), y = bg::get<1>(first.outer()[1]);
        rational const sa = (rational(bg::get<1>(first.outer()[2])) - y)
                          / (rational(bg::get<0>(first.outer()[2])) - x);
        rational const sb = (rational(bg::get<1>(second.outer()[0])) - rational(bg::get<1>(second.outer()[1])))
                          / (rational(bg::get<0>(second.outer()[0])) - rational(bg::get<0>(second.outer()[1])));
        return rational(h * h / (2 * (sa - sb)));
    };
    rational const expected = expected_area(a, b, rational(595245407, cpp_int("120000111100000000000000")));
    BOOST_TEST_EQ(area(out), expected);
    BOOST_TEST_EQ(bg::relation(a, b, rational_strategy).str(), "212101212");

    auto const fa = bg::prepare_filtered_decimal_polygon(read<polygon<decimal>>(issue_663_a));
    auto const fb = bg::prepare_filtered_decimal_polygon(read<polygon<decimal>>(issue_663_b));
    rational_polygons filtered;
    BOOST_TEST(bg::intersection(fa, fb, filtered, filtered_strategy));
    BOOST_TEST_EQ(filtered.size(), 1u);
    BOOST_TEST_EQ(area(filtered), expected);
    BOOST_TEST_EQ(bg::relation(fa, fb, filtered_strategy).str(), "212101212");

    auto const ba = read<polygon<double>>(issue_663_a);
    auto const bb = read<polygon<double>>(issue_663_b);
    rational_polygons binary;
    BOOST_TEST(bg::intersection(ba, bb, binary, binary_strategy));
    BOOST_TEST_EQ(binary.size(), 1u);
    rational const binary_h(cpp_int("11804726901675"), cpp_int("9903529483292268689467899904"));
    BOOST_TEST_EQ(area(binary), expected_area(ba, bb, binary_h));
    BOOST_TEST(area(binary) < expected);
    BOOST_TEST_EQ(bg::relation(ba, bb, binary_strategy).str(), "212101212");
}

char const square[] = "POLYGON((0 0,1 0,1 1,0 1,0 0))";
char const near_square[] =
    "POLYGON((1.00000000000000000001 0,2 0,2 1,1.00000000000000000001 1,1.00000000000000000001 0))";

void decimal_gap()
{
    auto const input = read<polygon<decimal>>(near_square);
    BOOST_TEST_EQ(bg::get<0>(input.outer().front()).scalar(), "1.00000000000000000001");
    BOOST_TEST_EQ(bg::get<0>(input.outer().front()).binary64(), 1);
    auto const a = prepare(square);
    auto const b = bg::prepare_rational_polygon(input);
    BOOST_TEST_EQ(bg::relation(a, b, rational_strategy).str(), "FF2FF1212");
    BOOST_TEST(bg::relate(a, b, bg::de9im::mask("FF*FF****"), rational_strategy));

    auto const ba = read<polygon<double>>(square);
    auto const bb = read<polygon<double>>(near_square);
    BOOST_TEST_EQ(bg::relation(ba, bb, binary_strategy).str(), "FF2F11212");

    rational_polygons decimal_union, binary_union;
    bg::union_(a, b, decimal_union, rational_strategy);
    bg::union_(ba, bb, binary_union, binary_strategy);
    BOOST_TEST_EQ(decimal_union.size(), 2u);
    BOOST_TEST_EQ(binary_union.size(), 1u);

    auto const fa = bg::prepare_filtered_decimal_polygon(read<polygon<decimal>>(square));
    auto const fb = bg::prepare_filtered_decimal_polygon(input);
    BOOST_TEST_EQ(bg::relation(fa, fb, filtered_strategy).str(), "FF2FF1212");
    BOOST_TEST(bg::relate(fa, fb, bg::de9im::mask("FF*FF****"), filtered_strategy));
    rational_polygons filtered_union;
    bg::union_(fa, fb, filtered_union, filtered_strategy);
    BOOST_TEST_EQ(filtered_union.size(), 2u);
}

void hole_and_reuse()
{
    char const wkt[] = "MULTIPOLYGON(((0 0,4 0,4 4,0 4,0 0),(1 1,1 2,2 2,2 1,1 1)))";
    auto const source_donut = read<bg::model::multi_polygon<polygon<decimal>>>(wkt);
    auto const source_cutter = read<polygon<decimal>>("POLYGON((1.5 0,3 0,3 3,1.5 3,1.5 0))");
    auto const donut = bg::prepare_rational_polygon(source_donut);
    auto const cutter = bg::prepare_rational_polygon(source_cutter);
    BOOST_TEST_EQ(donut.size(), 1u);
    BOOST_TEST_EQ(donut.front().inners().size(), 1u);
    BOOST_TEST_EQ(bg::get<0>(donut.front().inners().front().front()), 1);

    rational_polygons common, rest, either;
    BOOST_TEST(bg::intersection(donut, cutter, common, rational_strategy));
    BOOST_TEST_EQ(area(common), 4);
    bg::difference(donut, cutter, rest, rational_strategy);
    BOOST_TEST_EQ(area(rest), 11);
    bg::sym_difference(donut, cutter, either, rational_strategy);
    BOOST_TEST_EQ(area(either), rational(23, 2));
    BOOST_TEST(bg::relate(donut, cutter, bg::de9im::mask("T********"), rational_strategy));
    auto const fd_donut = bg::prepare_filtered_decimal_polygon(source_donut);
    auto const fd_cutter = bg::prepare_filtered_decimal_polygon(source_cutter);
    auto const check_public_filters = [](char const* name, auto const& run)
    {
#ifdef BOOST_GEOMETRY_EXACT_FILTER_STATISTICS
        ec::filtering_statistics() = {};
#endif
        run();
#ifdef BOOST_GEOMETRY_EXACT_FILTER_STATISTICS
        auto const& hits = ec::filtering_statistics();
        BOOST_TEST(hits.decimal_side_hits > 0);
        BOOST_TEST(hits.decimal_range_hits > 0);
        BOOST_TEST(hits.decimal_equality_hits > 0);
        std::cout << name << ": " << hits.decimal_side_hits << " side, "
                  << hits.decimal_range_hits << " range, " << hits.decimal_equality_hits
                  << " equality hits; " << hits.decimal_side_fallbacks << " side fallbacks\n";
#else
        (void)name;
#endif
    };
    rational_polygons filtered_common, filtered_rest, filtered_either, filtered_union, plain_union;
    check_public_filters("intersection", [&]
    { BOOST_TEST(bg::intersection(fd_donut, fd_cutter, filtered_common, filtered_strategy)); });
    check_public_filters("difference", [&]
    { bg::difference(fd_donut, fd_cutter, filtered_rest, filtered_strategy); });
    check_public_filters("sym_difference", [&]
    { bg::sym_difference(fd_donut, fd_cutter, filtered_either, filtered_strategy); });
    check_public_filters("union", [&]
    { bg::union_(fd_donut, fd_cutter, filtered_union, filtered_strategy); });
    bg::union_(donut, cutter, plain_union, rational_strategy);
    auto const same_set = [](auto const& a, auto const& b)
    {
        BOOST_TEST_EQ(bg::relation(a, b, rational_strategy).str(), "2FFF1FFF2");
    };
    same_set(filtered_common, common);
    same_set(filtered_rest, rest);
    same_set(filtered_either, either);
    same_set(filtered_union, plain_union);
    check_public_filters("relation", [&]
    {
        BOOST_TEST_EQ(bg::relation(fd_donut, fd_cutter, filtered_strategy).str(),
                      bg::relation(donut, cutter, rational_strategy).str());
    });

    bg::model::polygon<point<decimal>, true, false> open;
    open.outer().emplace_back(decimal("0"), decimal("0"));
    open.outer().emplace_back(decimal("0"), decimal("1"));
    open.outer().emplace_back(decimal("1"), decimal("0"));
    auto const prepared = bg::prepare_rational_polygon(open);
    static_assert(bg::closure<decltype(prepared)>::value == bg::open, "Preserve open rings");
    static_assert(bg::point_order<decltype(prepared)>::value == bg::clockwise, "Preserve ring order");
    BOOST_TEST_EQ(prepared.outer().size(), 3u);
    auto const filtered_open = bg::prepare_filtered_decimal_polygon(open);
    static_assert(bg::closure<decltype(filtered_open)>::value == bg::open, "Preserve open rings");
    static_assert(bg::point_order<decltype(filtered_open)>::value == bg::clockwise, "Preserve ring order");
    BOOST_TEST_EQ(filtered_open.outer().size(), 3u);

    using wrapped = polygon<decimal>;
    BOOST_TEST_THROWS(read<wrapped>("POLYGON((0 0,1e+ 0,1 1,0 0))"), bg::read_wkt_exception);
}

void unbounded_rational()
{
    // Both determinants in a ratio can exceed binary64's range, or
    // underflow it. The ratios and ordered turns must still be exact.
    cpp_int const huge = boost::multiprecision::pow(cpp_int(10), 800u);
    for (rational const& unit : {rational(huge), rational(1, huge)})
    {
        bg::segment_ratio<rational> const one_third(unit, 3 * unit);
        bg::segment_ratio<rational> const same(2 * unit, 6 * unit);
        bg::segment_ratio<rational> const two_thirds(2 * unit, 3 * unit);
        BOOST_TEST(one_third == same);
        BOOST_TEST(one_third < two_thirds);
        BOOST_TEST(!(two_thirds < one_third));
    }

    for (char const* exponent : {"0", "400", "-400"})
    {
        std::string const e = std::string("e") + exponent;
        auto const source_a = read<polygon<decimal>>(("POLYGON((0 0,3" + e + " 0,0 1,0 0))").c_str());
        auto const source_b = read<polygon<decimal>>(
            ("POLYGON((1" + e + " 0,4" + e + " 0,4" + e + " 1,1" + e + " 1,1" + e + " 0))").c_str());
        auto const a = bg::prepare_rational_polygon(source_a);
        auto const b = bg::prepare_rational_polygon(source_b);
        rational_polygons common;
        BOOST_TEST(bg::intersection(a, b, common, rational_strategy));
        BOOST_TEST_EQ(common.size(), 1u);
        rational const scale = decimal("1" + e).rational();
        BOOST_TEST_EQ(area(common), 2 * scale / 3);
        BOOST_TEST_EQ(bg::relation(a, b, rational_strategy).str(), "212111212");

        bool found = false;
        for (auto const& part : common)
        for (auto const& vertex : part.outer())
        {
            if (bg::get<0>(vertex) == scale && bg::get<1>(vertex) == rational(2, 3))
            {
                found = true;
            }
        }
        BOOST_TEST(found);
        auto const fa = bg::prepare_filtered_decimal_polygon(source_a);
        auto const fb = bg::prepare_filtered_decimal_polygon(source_b);
        BOOST_TEST(bool(fb.outer().front().certified_approximation()) == (e == "e0"));
        rational_polygons filtered;
        BOOST_TEST(bg::intersection(fa, fb, filtered, filtered_strategy));
        BOOST_TEST_EQ(bg::relation(fa, fb, filtered_strategy).str(), "212111212");
        BOOST_TEST_EQ(bg::relation(filtered, common, rational_strategy).str(), "2FFF1FFF2");
        BOOST_TEST_EQ(area(filtered), 2 * scale / 3);
        found = false;
        for (auto const& part : filtered)
        for (auto const& vertex : part.outer())
            if (bg::get<0>(vertex) == scale && bg::get<1>(vertex) == rational(2, 3))
                found = true;
        BOOST_TEST(found);
    }
}

void certificate_lifetime()
{
    fd::point const a(decimal("0.1"), decimal("0.2")), b(decimal("2.1"), decimal("0.2"));
    fd::point p(decimal("0.6"), decimal("1.2"));
    auto const* approx = p.certified_approximation();
    BOOST_TEST(approx != nullptr);
    if (approx)
    {
        BOOST_TEST_EQ(approx->value[0], decimal("0.6").binary64());
        BOOST_TEST(approx->error[0] > 0);
        BOOST_TEST(abs(bg::get<0>(p) - rational(approx->value[0])) <= rational(approx->error[0]));
    }
    BOOST_TEST_EQ(fd::predicates::apply(a, b, p), 1);
    fd::point moved(std::move(p));
    BOOST_TEST(p.certified_approximation() == nullptr);
    BOOST_TEST(moved.certified_approximation() != nullptr);
    p = moved;
    BOOST_TEST(p.certified_approximation() != nullptr);
#ifdef BOOST_GEOMETRY_EXACT_FILTER_STATISTICS
    ec::filtering_statistics() = {};
#endif
    bg::set<1>(p, rational(0));
    BOOST_TEST(p.certified_approximation() == nullptr);
    BOOST_TEST_EQ(fd::predicates::apply(a, b, p), -1);
    BOOST_TEST_EQ(fd::predicates::apply(a, b, moved), 1);
#ifdef BOOST_GEOMETRY_EXACT_FILTER_STATISTICS
    BOOST_TEST_EQ(ec::filtering_statistics().decimal_side_fallbacks, 1u);
    BOOST_TEST_EQ(ec::filtering_statistics().decimal_side_hits, 1u);
#endif

    fd::point const one(decimal("1"), decimal("0"));
    fd::point const close(decimal("1.00000000000000000001"), decimal("0"));
    bool answer = false;
    BOOST_TEST(!fd::predicates::try_equal(one, close, answer));
    BOOST_TEST(!fd::predicates::try_disjoint_by_range<0>(one, one, close, close, answer));
    BOOST_TEST(fd::predicates::try_equal(one, one, answer) && answer);
}

template <typename F>
void measure(char const* fixture, char const* stage, F&& run, int iterations = 1000)
{
    using clock = std::chrono::steady_clock;
    std::size_t checksum = run();
    auto const start = clock::now();
    for (int i = 0; i < iterations; ++i)
    {
        checksum += run();
    }
    auto const elapsed = std::chrono::duration<double, std::micro>(clock::now() - start).count();
    std::cout << fixture << '\t' << stage << '\t' << iterations << '\t'
              << elapsed / iterations << '\t' << checksum << '\n';
}

void benchmark(char const* fixture, char const* first, char const* second, int iterations = 1000)
{
    auto const a = read<polygon<decimal>>(first);
    auto const b = read<polygon<decimal>>(second);
    auto const ra = bg::prepare_rational_polygon(a);
    auto const rb = bg::prepare_rational_polygon(b);
    auto const fa = bg::prepare_filtered_decimal_polygon(a);
    auto const fb = bg::prepare_filtered_decimal_polygon(b);

    measure(fixture, "parse_decimal_wrappers", [&]
    {
        return read<polygon<decimal>>(first).outer().size() + read<polygon<decimal>>(second).outer().size();
    }, iterations);
    measure(fixture, "prepare_rational_pair", [&]
    {
        return bg::prepare_rational_polygon(a).outer().size() + bg::prepare_rational_polygon(b).outer().size();
    }, iterations);
    measure(fixture, "prepare_filtered_pair", [&]
    {
        return bg::prepare_filtered_decimal_polygon(a).outer().size()
            + bg::prepare_filtered_decimal_polygon(b).outer().size();
    }, iterations);
    auto const intersect = [](auto const& x, auto const& y, auto const& strategy)
    {
        rational_polygons out;
        bg::intersection(x, y, out, strategy);
        return out.size();
    };
    measure(fixture, "rational_intersection", [&] { return intersect(ra, rb, rational_strategy); }, iterations);
    measure(fixture, "annotated_unfiltered_intersection", [&] { return intersect(fa, fb, rational_strategy); }, iterations);
    measure(fixture, "filtered_intersection", [&] { return intersect(fa, fb, filtered_strategy); }, iterations);
    measure(fixture, "rational_relation", [&] { return bg::relation(ra, rb, rational_strategy).str()[0]; }, iterations);
    measure(fixture, "annotated_unfiltered_relation", [&] { return bg::relation(fa, fb, rational_strategy).str()[0]; }, iterations);
    measure(fixture, "filtered_relation", [&] { return bg::relation(fa, fb, filtered_strategy).str()[0]; }, iterations);
}

std::string corridor(int x_offset, int y_offset, int phase)
{
    auto decimal_text = [](int tenths)
    { return std::to_string(tenths / 10) + '.' + std::to_string(tenths % 10); };
    auto vertex = [&](int x, int y)
    { return decimal_text(x + x_offset) + ' ' + decimal_text(y + y_offset); };
    std::string result = "POLYGON((";
    auto bottom = [&](int i) { return vertex(10 * i, 10 + (i + phase) % 5 * 7); };
    for (int i = 0; i <= 48; ++i)
        result += bottom(i) + ',';
    return result + vertex(480, 92) + ',' + vertex(0, 92) + ',' + bottom(0) + "))";
}

int main(int argc, char** argv)
{
    if (argc == 2 && std::string(argv[1]) == "--benchmark")
    {
        std::cout << "fixture\tstage\titerations\tus_per_pair\tchecksum\n";
        benchmark("crossing",
            "POLYGON((0.1 0.1,4.1 0.1,0.1 4.1,0.1 0.1))",
            "POLYGON((1.2 1.3,4.2 1.3,4.2 3.3,1.2 3.3,1.2 1.3))");
        benchmark("issue_663", issue_663_a, issue_663_b);
        benchmark("decimal_gap", square, near_square);
        benchmark("corridors", corridor(1, 1, 0).c_str(), corridor(23, 22, 2).c_str(), 60);
        return 0;
    }
    issue_663();
    decimal_gap();
    hole_and_reuse();
    unbounded_rational();
    certificate_lifetime();
    return boost::report_errors();
}
