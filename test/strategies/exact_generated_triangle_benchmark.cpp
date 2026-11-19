// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.

#ifdef BG_TRIANGLE_PROBE
namespace triangle_probe
{
unsigned candidates = 0, certified = 0, boundary_queries = 0, raw_kept = 0, direct = 0, split = 0;
}
#endif

#include <boost/core/lightweight_test.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/strategies/relate/cartesian_exact.hpp>

#include <chrono>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <string>

namespace bg = boost::geometry;
using rational = bg::strategies::detail::exact_cartesian::exact_number;
template <typename C> using polygon = bg::model::polygon<bg::model::d2::point_xy<C>>;
template <typename C> using multi = bg::model::multi_polygon<polygon<C>>;
using strategy = bg::strategies::relate::cartesian_exact<double>;

struct fixture
{
    char const* name;
    polygon<double> a;
    multi<double> b;
    multi<rational> expected;
    int area;
};

fixture fixture_named(std::string const& name)
{
    fixture f;
    if (name == "touching-triangles")
    {
        f.name = "touching-triangles";
        bg::read_wkt("POLYGON((-2 0,-2 4,2 4,2 0,-2 0))", f.a);
        bg::read_wkt("MULTIPOLYGON(((0 2,4 1,4 0,0 2)),((0 2,4 4,4 3,0 2)))", f.b);
        f.expected.resize(2);
        f.expected[0].outer() = {{0, 2}, {2, rational(3, 2)}, {2, 1}, {0, 2}};
        f.expected[1].outer() = {{0, 2}, {2, 3}, {2, rational(5, 2)}, {0, 2}};
        f.area = 1;
    }
    else
    {
        f.name = "convex-clip";
        bg::read_wkt("POLYGON((0 0,0 4,5 4,5 0,0 0))", f.a);
        bg::read_wkt("MULTIPOLYGON(((0 -1,1 5,4 5,3 -1,0 -1)))", f.b);
        f.expected.resize(1);
        f.expected[0].outer() = {{rational(1, 6), 0}, {rational(5, 6), 4},
                                 {rational(23, 6), 4}, {rational(19, 6), 0}, {rational(1, 6), 0}};
        f.area = 12;
    }
    return f;
}

void verify(fixture const& f)
{
    polygon<rational> exact_a;
    multi<rational> exact_b, result;
    bg::convert(f.a, exact_a);
    bg::convert(f.b, exact_b);
    bg::strategies::relate::cartesian<rational> const validation;
    BOOST_TEST(bg::is_valid(exact_a, validation));
    BOOST_TEST(bg::is_valid(exact_b, validation));
    BOOST_TEST(bg::is_valid(f.expected, validation));
#ifdef BG_TRIANGLE_PROBE
    namespace p = triangle_probe;
    p::candidates = p::certified = p::boundary_queries = p::raw_kept = p::direct = p::split = 0;
#endif
    BOOST_TEST(bg::intersection(f.a, f.b, result, strategy{}));
#ifdef BG_TRIANGLE_PROBE
    std::cout << f.name << " candidates=" << p::candidates << " certified=" << p::certified
              << " boundary_queries=" << p::boundary_queries << " raw_kept=" << p::raw_kept
              << " direct=" << p::direct << " split=" << p::split << '\n';
#endif
    BOOST_TEST_EQ(result.size(), f.expected.size());
    BOOST_TEST(bg::is_valid(result, validation));
    BOOST_TEST_EQ(bg::area(result, validation), f.area);
    BOOST_TEST_EQ(bg::relation(result, f.expected, validation).str(), "2FFF1FFF2");
    std::cout << f.name << ' ' << bg::wkt(result) << '\n';
}

BOOST_NOINLINE std::size_t run(fixture const& f, multi<double>& result)
{
    result.clear();
    bg::intersection(f.a, f.b, result, strategy{});
    return result.size();
}

void measure(fixture const& f, double seconds)
{
    using clock = std::chrono::steady_clock;
    multi<double> result;
    auto const expected = run(f, result);
    BOOST_TEST_EQ(expected, f.expected.size());
    auto const cpu_begin = std::clock();
    auto const begin = clock::now();
    std::size_t calls = 0, checksum = 0;
    do
    {
        for (int i = 0; i < 64; ++i) checksum += run(f, result);
        calls += 64;
    } while (std::chrono::duration<double>(clock::now() - begin).count() < seconds);
    double const wall = std::chrono::duration<double, std::micro>(clock::now() - begin).count() / calls;
    double const cpu = 1e6 * double(std::clock() - cpu_begin) / CLOCKS_PER_SEC / calls;
    BOOST_TEST_EQ(checksum, expected * calls);
    std::cout << f.name << '\t' << std::setprecision(9) << cpu << '\t' << wall << '\t'
              << calls << '\t' << checksum << '\n';
}

int main(int argc, char** argv)
{
    if (argc == 1)
    {
        verify(fixture_named("touching-triangles"));
        verify(fixture_named("convex-clip"));
    }
    else if (argc == 3 && (std::string(argv[1]) == "touching-triangles" || std::string(argv[1]) == "convex-clip"))
        measure(fixture_named(argv[1]), std::stod(argv[2]));
    else return 2;
    return boost::report_errors();
}
