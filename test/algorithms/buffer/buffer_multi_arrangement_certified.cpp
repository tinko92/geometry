// Boost.Geometry
// Unit Test
// Copyright (c) 2026
// Distributed under the Boost Software License, Version 1.0.

#include <boost/geometry/algorithms/buffer_multi_arrangement_certified.hpp>
#include <boost/geometry/geometries/multi_linestring.hpp>
#include <boost/geometry/geometries/multi_point.hpp>
#include <boost/geometry/geometries/multi_polygon.hpp>
#include <boost/geometry/io/wkt/read.hpp>
#include <boost/core/lightweight_test.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>

namespace bg = boost::geometry;
namespace ec = bg::strategies::detail::exact_cartesian;
using number = ec::exact_number;
using status = bg::certified_buffer_status;
template <typename T> using point = bg::model::point<T, 2, bg::cs::cartesian>;
using p = point<double>;
using mp = bg::model::multi_point<p>;
using line = bg::model::linestring<p>;
using ml = bg::model::multi_linestring<line>;
using polygon = bg::model::polygon<p>;
using multi = bg::model::multi_polygon<polygon>;
using input_polygon = bg::model::polygon<p, false>;
using input_multi = bg::model::multi_polygon<input_polygon>;

struct xy { number x, y; };
number q(int n, int d) { return number(n) / d; }
number const displacement = number(1) / (boost::multiprecision::cpp_int(1) << 42);

xy read(p const& p) { return {ec::exact_value(bg::get<0>(p)), ec::exact_value(bg::get<1>(p))}; }

template <typename Ring, typename Visit>
void edges(Ring const& ring, Visit const& visit)
{
    for (std::size_t i = 1; i < ring.size(); ++i) { visit(read(ring[i - 1]), read(ring[i])); }
}

int locate(multi const& m, xy const& p)
{
    // Independent exact oracle on the stored binary coordinates.
    bool inside = false, boundary = false;
    auto const crossing = [&](xy const& a, xy const& b)
    {
        if ((p.x - a.x) * (p.x - b.x) <= 0 && (p.y - a.y) * (p.y - b.y) <= 0
            && (p.x - a.x) * (b.y - a.y) == (p.y - a.y) * (b.x - a.x)) { boundary = true; }
        if ((a.y > p.y) != (b.y > p.y) && a.x + (b.x - a.x) * (p.y - a.y) / (b.y - a.y) > p.x)
        { inside = !inside; }
    };
    for (auto const& poly : m)
    {
        edges(poly.outer(), crossing);
        for (auto const& h : poly.inners()) { edges(h, crossing); }
    }
    return boundary ? 0 : inside ? 1 : -1;
}

number area2(multi const& m)
{
    number sum = 0;
    auto const add = [&](xy const& a, xy const& b) { sum += a.y * b.x - a.x * b.y; };
    for (auto const& poly : m)
    {
        edges(poly.outer(), add);
        for (auto const& h : poly.inners()) { edges(h, add); }
    }
    return sum;
}

std::size_t holes(multi const& m)
{
    std::size_t n = 0;
    for (auto const& poly : m) { n += poly.inners().size(); }
    return n;
}

input_polygon rectangle(double x, double y, double u, double v)
{
    input_polygon poly;
    poly.outer() = {p(x, y), p(u, y), p(u, v), p(x, v), p(x, y)};
    return poly;
}

void test_points()
{
    // Quarter chords are unit diamonds. The two overlapping diamonds
    // have an intersection of area 1/8, plus one separate unit diamond.
    mp const members{p(0, 0), p(1.5, 0), p(40, 0), p(0, 0)};
    multi out;
    BOOST_TEST(bg::buffer_multi_arrangement_certified(members, out, 1, q(1, 2), 0) == status::success);
    BOOST_TEST_EQ(out.size(), 2u);
    BOOST_TEST_EQ(area2(out), q(47, 4));
    BOOST_TEST_EQ(locate(out, {q(3, 4), 0}), 1);
    BOOST_TEST_EQ(locate(out, {20, 0}), -1);
    BOOST_TEST_EQ(locate(out, {40, 0}), 1);

    // At exactly two units the closed diamonds touch at a single point.
    mp const touch{p(0, 0), p(2, 0)};
    BOOST_TEST(bg::buffer_multi_arrangement_certified(touch, out, 1, q(1, 2), 0) == status::cannot_certify);
    BOOST_TEST_EQ(area2(out), q(47, 4));
    using rp = point<number>;
    // The ideal diamond gap is positive, but both extreme x coordinates
    // round to 1. The displacement itself fits; native contact must fail.
    number const tiny = number(1) / (boost::multiprecision::cpp_int(1) << 60);
    bg::model::multi_point<rp> const near{rp(0, 0), rp(2 + tiny, 0)};
    BOOST_TEST(bg::buffer_multi_arrangement_certified(near, out, 1, q(1, 2), displacement)
               == status::cannot_certify);
    BOOST_TEST_EQ(area2(out), q(47, 4));
    mp too_many;
    for (int i = 0; i < 17; ++i) { too_many.emplace_back(4 * i, 0); }
    BOOST_TEST(bg::buffer_multi_arrangement_certified(too_many, out, 1, q(1, 2), 0) == status::cannot_certify);
    BOOST_TEST_EQ(area2(out), q(47, 4));
    BOOST_TEST(bg::buffer_multi_arrangement_certified(touch, out, -1, q(1, 2), 0) == status::invalid_input);
    double const infinity = (std::numeric_limits<double>::infinity)();
    BOOST_TEST(bg::buffer_multi_arrangement_certified(mp{p(0, infinity)}, out, 1, q(1, 2), 0)
               == status::invalid_input);
    BOOST_TEST_EQ(area2(out), q(47, 4));

    // Rounding by (1/3,1/3) meets a 2/5 per-axis tolerance but fails the
    // shared Euclidean check. This must not publish a proposed polygon.
    bg::model::multi_point<rp> const fractional{rp(q(13, 3), q(13, 3))};
    bg::model::multi_polygon<bg::model::polygon<point<int>>> integral;
    BOOST_TEST(bg::buffer_multi_arrangement_certified(fractional, integral, 1, q(1, 2), q(2, 5))
               == status::cannot_certify);
    BOOST_TEST(integral.empty());
}

void test_lines_and_holes()
{
    // Members meet at (0,12); their other ends nearly close the irrational
    // loop. A repeated-point member contributes a disjoint disk.
    ml const input{line{p(0, 0), p(8, 8), p(0, 12)}, line{p(0, 12), p(-8, 4), p(-0.5, 0.25)},
                   line{p(20, 20), p(20, 20)}, line{}};
    multi out;
    BOOST_TEST(bg::buffer_multi_arrangement_certified(input, out, 1, q(1, 8), displacement) == status::success);
    BOOST_TEST_EQ(out.size(), 2u);
    BOOST_TEST_EQ(holes(out), 1u);
    BOOST_TEST_EQ(locate(out, {0, 6}), -1);
    BOOST_TEST_EQ(locate(out, {4, 4}), 1);
    BOOST_TEST_EQ(locate(out, {20, 20}), 1);
    BOOST_TEST_EQ(locate(out, {q(-1, 4), q(1, 8)}), 1);
}

void test_1233_original_multi()
{
    multi out;
    double constexpr length = 0x1p52;
    ml const reduced{line{p(0, 100), p(-length, 0)}, line{p(-length, 0), p(0, -100)},
                     line{p(1000, 0), p(1002, 0)}};
    BOOST_TEST(bg::buffer_multi_arrangement_certified(reduced, out, 1, q(1, 8), q(1, 2)) == status::success);
    BOOST_TEST_EQ(out.size(), 2u);
    BOOST_TEST_EQ(locate(out, {-length, 0}), 1);
    BOOST_TEST_EQ(locate(out, {-length / 2, 0}), -1);
    BOOST_TEST_EQ(locate(out, {1001, 0}), 1);
    auto const before = area2(out);

    // Full published WKT: https://github.com/boostorg/geometry/issues/1233
    ml original;
    bg::read_wkt("MULTILINESTRING((11393 -19530,-5728 30983,-9223372036854775807 15322,27608 14504),"
                 "(-32468 -20811,-27782 -14559,-30822 11140,-12550 27458),"
                 "(24176 -9103,4178 -868,7332 11648,21439 30360,-121 11529,8803 -10835,5660 20213))",
                 original);
    BOOST_TEST_EQ(original.size(), 3u);
    BOOST_TEST_EQ(original.front().size(), 4u);
    double constexpr min_x = -0x1p63;
    BOOST_TEST_EQ(bg::get<0>(original.front()[2]), min_x);
    double const predecessor = std::nextafter(min_x, -(std::numeric_limits<double>::infinity)());
    // Any valid stored polygon would need a minimum x in [x-3/2,x-3/8],
    // from the outer/inner neighborhood bounds. There is no binary64 value
    // in this interval. The next smaller x is x-2048.
    BOOST_TEST_EQ(ec::exact_value(min_x) - ec::exact_value(predecessor), 2048);
    BOOST_TEST(bg::buffer_multi_arrangement_certified(original, out, 1, q(1, 8), q(1, 2)) == status::cannot_certify);
    BOOST_TEST_EQ(out.size(), 2u);
    BOOST_TEST_EQ(area2(out), before);
}

void test_polygons()
{
    auto annulus = rectangle(0, 0, 20, 20);
    annulus.inners().resize(1);
    annulus.inners().front() = rectangle(5, 5, 15, 15).outer();
    std::reverse(annulus.inners().front().begin(), annulus.inners().front().end());
    // This member is nested in the hole, not in the other's material.
    input_multi island{annulus, rectangle(8, 8, 12, 12), input_polygon{}};
    multi out;
    BOOST_TEST(bg::buffer_multi_arrangement_certified(island, out, 1, q(1, 2), 0) == status::success);
    BOOST_TEST_EQ(out.size(), 2u);
    BOOST_TEST_EQ(holes(out), 1u);
    BOOST_TEST_EQ(locate(out, {10, 10}), 1);
    BOOST_TEST_EQ(locate(out, {q(13, 2), 10}), -1);
    BOOST_TEST_EQ(locate(out, {q(11, 2), 10}), 1);
    BOOST_TEST(bg::buffer_multi_arrangement_certified(island, out, -1, q(1, 2), 0) == status::success);
    BOOST_TEST_EQ(out.size(), 2u);
    BOOST_TEST_EQ(holes(out), 1u);
    BOOST_TEST_EQ(locate(out, {10, 10}), 1);
    BOOST_TEST_EQ(locate(out, {7, 10}), -1);
    BOOST_TEST_EQ(locate(out, {2, 10}), 1);

    input_multi disjoint{rectangle(0, 0, 4, 4), rectangle(5.5, 0, 9.5, 4)};
    BOOST_TEST(bg::buffer_multi_arrangement_certified(disjoint, out, 1, q(1, 2), 0) == status::success);
    BOOST_TEST_EQ(out.size(), 1u);
    BOOST_TEST_EQ(locate(out, {q(19, 4), 2}), 1);
    auto const before = area2(out);
    input_multi const overlap{annulus, rectangle(2, 2, 3, 3)};
    BOOST_TEST(bg::buffer_multi_arrangement_certified(overlap, out, 1, q(1, 2), 0) == status::invalid_input);
    input_multi const touch{rectangle(0, 0, 4, 4), rectangle(4, 4, 6, 6)};
    BOOST_TEST(bg::buffer_multi_arrangement_certified(touch, out, 1, q(1, 2), 0) == status::invalid_input);
    BOOST_TEST_EQ(area2(out), before);
}

void test_empties()
{
    multi out;
    out.resize(1);
    BOOST_TEST(bg::buffer_multi_arrangement_certified(mp{}, out, 1, q(1, 2), 0) == status::success);
    BOOST_TEST(out.empty());
    out.resize(1);
    BOOST_TEST(bg::buffer_multi_arrangement_certified(ml{line{}}, out, 1, q(1, 2), 0) == status::success);
    BOOST_TEST(out.empty());
    out.resize(1);
    BOOST_TEST(bg::buffer_multi_arrangement_certified(input_multi{input_polygon{}}, out, -1, q(1, 2), 0)
               == status::success);
    BOOST_TEST(out.empty());
    out.resize(1);
    BOOST_TEST(bg::buffer_multi_arrangement_certified(mp{}, out, 1, 0, 0) == status::cannot_certify);
    BOOST_TEST_EQ(out.size(), 1u);
}

void test_irrational_polygon_members()
{
    auto triangle = [](double x)
    {
        input_polygon p;
        p.outer().emplace_back(x, 0);
        p.outer().emplace_back(x + 3, 0);
        p.outer().emplace_back(x, 2);
        p.outer().push_back(p.outer().front());
        return p;
    };
    input_multi const members{triangle(0), triangle(5.5)};
    multi out;
    BOOST_TEST(bg::buffer_multi_arrangement_certified(members, out, 1, q(1, 5), displacement) == status::success);
    BOOST_TEST_EQ(out.size(), 2u);
    BOOST_TEST_EQ(locate(out, {q(17, 4), q(1, 4)}), -1);
    BOOST_TEST(bg::buffer_multi_arrangement_certified(members, out, q(-1, 2), q(1, 10), displacement)
               == status::success);
    BOOST_TEST_EQ(out.size(), 2u);
    BOOST_TEST_EQ(locate(out, {q(3, 5), q(3, 5)}), 1);
    BOOST_TEST_EQ(locate(out, {q(61, 10), q(3, 5)}), 1);
}

int main()
{
    test_points();
    test_lines_and_holes();
    test_1233_original_multi();
    test_polygons();
    test_empties();
    test_irrational_polygon_members();
    return boost::report_errors();
}
