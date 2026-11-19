// Boost.Geometry
// Copyright (c) 2026
// Distributed under the Boost Software License, Version 1.0.

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

#include <boost/core/lightweight_test.hpp>
#include <boost/geometry/algorithms/erode_convex_polygon_certified.hpp>
#include <boost/geometry.hpp>

#include "exact_issue_1281_fixture.hpp"

namespace bg = boost::geometry;
namespace erosion = bg::detail::certified_convex_erosion;
using q = erosion::exact_number;
using status = bg::certified_buffer_status;
using point = bg::model::d2::point_xy<double>;
using polygon = bg::model::polygon<point, false>;
using multi = bg::model::multi_polygon<polygon>;
using rational_polygon = bg::model::polygon<bg::model::d2::point_xy<q>, false>;
using rational_strategy = bg::strategies::relate::cartesian<q>;

template <typename Polygon = polygon>
Polygon source(char const* wkt)
{
    Polygon result;
    bg::read_wkt(wkt, result);
    return result;
}

bool exact_valid(polygon const& p)
{
    rational_polygon lifted;
    bg::convert(p, lifted);
    return bg::is_valid(lifted, rational_strategy());
}

void test_irrational()
{
    auto const a = source("POLYGON((0 0,4 0,0 4,0 0))");
    multi out;
    BOOST_TEST(bg::erode_convex_polygon_certified(a, out, 1, 1e-14) == status::success);
    BOOST_TEST_EQ(out.size(), 1u);
    if (out.size() != 1) { return; }
    BOOST_TEST(exact_valid(out.front()));
    BOOST_TEST_EQ(out.front().outer().size(), 4u);
    double const width = 2 - std::sqrt(2.0);
    BOOST_TEST(std::abs(bg::area(out) - width*width/2) < 1e-14);
    bool corner = false;
    for (auto const& p : out.front().outer())
        corner |= bg::get<0>(p) == 1 && bg::get<1>(p) == 1;
    BOOST_TEST(corner);
    double const area = bg::area(out);
    BOOST_TEST(bg::erode_convex_polygon_certified(a, out, 1, 0) == status::cannot_certify);
    BOOST_TEST_EQ(bg::area(out), area);
}

void test_three_root_sign()
{
    erosion::side const a({0,0},{1,1},1), b({0,0},{2,1},1), c({0,0},{3,2},1);
    // 2, 5, 13 are independent square classes. At the rational endpoints
    // of this enclosure, only the exact field can decide the strict sign.
    std::array<erosion::term, 3> const terms{{{1,&a}, {1,&b}, {1,&c}}};
    erosion::enclosure bound{0,0};
    for (auto const& t : terms) { bound.add_scaled(t.support->root, t.coefficient); }
    BOOST_TEST_EQ(erosion::sign(-bound.low, terms), 1);
    BOOST_TEST_EQ(erosion::sign(-bound.high, terms), -1);
    erosion::side const d({0,0},{2,2},1);
    BOOST_TEST_EQ(erosion::sign(0, {{{2,&a}, {-1,&d}, {0,&c}}}), 0);

    // Finite inputs whose normals/products cannot be enclosed in binary64
    // still have a rational result. A lost filter is not an empty erosion.
    auto box = source("POLYGON((0 0,1 0,1 1,0 1,0 0))");
    double const w = 0x1p600, r = 0x1p598;
    for (auto& p : box.outer())
    {
        bg::set<0>(p, bg::get<0>(p)*w);
        bg::set<1>(p, bg::get<1>(p)*w);
    }
    multi out;
    BOOST_TEST(bg::erode_convex_polygon_certified(box, out, r, 0) == status::success);
    BOOST_TEST_EQ(out.size(), 1u);
    if (out.size() != 1) { return; }
    BOOST_TEST_EQ(out.front().outer().size(), 5u);
    for (auto const& p : out.front().outer())
    {
        BOOST_TEST(bg::get<0>(p) == r || bg::get<0>(p) == w-r);
        BOOST_TEST(bg::get<1>(p) == r || bg::get<1>(p) == w-r);
    }
}

void test_redundancy()
{
    auto const a = source("POLYGON((0 0,2 0,2 0,5 0,6 1,6 6,0 6,0 0))");
    multi out;
    BOOST_TEST(bg::erode_convex_polygon_certified(a, out, 2, 0) == status::success);
    BOOST_TEST_EQ(out.size(), 1u);
    if (!out.empty()) { BOOST_TEST_EQ(out.front().outer().size(), 5u); }
    BOOST_TEST_EQ(bg::area(out), 4.0);

    // The shifted 3-4-5 bevel passes through (7,3), an inset corner.
    // Rotate input storage to exercise both ways of entering an on-line vertex.
    auto b = source("POLYGON((0 0,8.5 0,10 2,10 9,0 9,0 0))");
    for (int i = 0; i < 5; ++i)
    {
        BOOST_TEST(bg::erode_convex_polygon_certified(b, out, 3, 0) == status::success);
        BOOST_TEST_EQ(bg::area(out), 12.0);
        BOOST_TEST_EQ(bg::num_points(out), 5u);
        auto& ring = b.outer();
        std::rotate(ring.begin(), ring.begin()+1, ring.end()-1);
        ring.back() = ring.front();
    }

    using cw_open = bg::model::polygon<point, true, false>;
    auto const c = source<cw_open>("POLYGON((0 0,0 4,6 4,6 0))");
    bg::model::multi_polygon<cw_open> open;
    BOOST_TEST(bg::erode_convex_polygon_certified(c, open, 1, 0) == status::success);
    BOOST_TEST_EQ(bg::area(open), 8.0);
    BOOST_TEST_EQ(bg::num_points(open), 4u);
    BOOST_TEST(bg::is_valid(open));
}

void test_collapse()
{
    multi out;
    out.push_back(source("POLYGON((0 0,10 0,0 10,0 0))"));
    auto const square = source("POLYGON((0 0,4 0,4 4,0 4,0 0))");
    auto const rect = source("POLYGON((0 0,6 0,6 4,0 4,0 0))");
    for (auto const& a : {square, rect})
    {
        // The square collapses to a point; the rectangle to a segment.
        BOOST_TEST(bg::erode_convex_polygon_certified(a, out, 2, 0) == status::success);
        BOOST_TEST(out.empty());
        BOOST_TEST(bg::erode_convex_polygon_certified(a, out, std::nextafter(2.0, 3.0), 0)
                   == status::success);
        BOOST_TEST(out.empty());
    }

    q const below = q(1171572875253809LL) / 1000000000000000LL;
    q const above = q(1171572875253810LL) / 1000000000000000LL;
    auto const local = source("POLYGON((0 0,4 0,0 4,0 0))");
    // The local near-collapse is also too small for a strictly convex
    // binary32 result; the ideal polygon is nonetheless nonempty.
    using float_polygon = bg::model::polygon<bg::model::d2::point_xy<float>, false>;
    bg::model::multi_polygon<float_polygon> single;
    BOOST_TEST(bg::erode_convex_polygon_certified(local, single, below, 1) == status::cannot_certify);
    BOOST_TEST(bg::erode_convex_polygon_certified(local, out, below, 1e-15) == status::success);
    BOOST_TEST_EQ(out.size(), 1u);

    auto shifted = local;
    for (auto& p : shifted.outer())
    {
        bg::set<0>(p, bg::get<0>(p) + 0x1p40);
        bg::set<1>(p, bg::get<1>(p) + 0x1p40);
    }
    auto const count = bg::num_points(out);
    BOOST_TEST(bg::erode_convex_polygon_certified(shifted, out, below, 1) == status::cannot_certify);
    BOOST_TEST_EQ(bg::num_points(out), count);
    BOOST_TEST(bg::erode_convex_polygon_certified(shifted, out, above, 0) == status::success);
    BOOST_TEST(out.empty());
}

void test_invalid()
{
    multi out;
    out.push_back(source("POLYGON((0 0,10 0,0 10,0 0))"));
    auto const reject = [&](polygon const& p)
    {
        BOOST_TEST(bg::erode_convex_polygon_certified(p, out, 1, 1) == status::invalid_input);
        BOOST_TEST_EQ(bg::area(out), 50.0);
    };
    reject(source("POLYGON((0 0,4 0,2 1,4 4,0 4,0 0))"));
    reject(source("POLYGON((0 0,4 0,2 0,4 4,0 4,0 0))"));
    auto star = source("POLYGON((0 3,2 -3,-3 1,3 1,-2 -3,0 3))");
    bg::correct(star);
    // All local turns are left, but they make two sweeps.
    auto const& ring = star.outer();
    for (std::size_t i = 0; i < 5; ++i)
    {
        auto const p = erosion::common::read_point(ring[i]);
        auto const a = erosion::subtract(erosion::common::read_point(ring[(i+1)%5]),p);
        auto const b = erosion::subtract(erosion::common::read_point(ring[(i+2)%5]),p);
        BOOST_TEST(erosion::cross(a,b) > 0);
    }
    reject(star);
    reject(source("POLYGON((0 0,4 0,4 4,0 4,0 0),(1 1,1 3,3 3,3 1,1 1))"));
    reject(source("POLYGON((0 0,0 4,4 0,0 0))"));
    auto closed = source("POLYGON((0 0,4 0,0 4,0 0))");
    closed.outer().pop_back();
    reject(closed);
    bg::set<0>(closed.outer().front(),std::numeric_limits<double>::infinity());
    reject(closed);
    closed.outer().resize(514, point(0,0));
    bg::set<0>(closed.outer().front(),0.0);
    BOOST_TEST(bg::erode_convex_polygon_certified(closed, out, 1, 1) == status::cannot_certify);
    q const huge(boost::multiprecision::cpp_int(1) << 4096);
    BOOST_TEST(bg::erode_convex_polygon_certified(closed, out, huge, 1) == status::cannot_certify);
    BOOST_TEST_EQ(bg::area(out), 50.0);
}

void replay_1281(char const* wkt)
{
    auto const a = source(wkt);
    BOOST_TEST(exact_valid(a));
    auto strict_source = [](polygon const& p)
    {
        std::vector<erosion::rational_xy> corners;
        BOOST_TEST(erosion::read_source(p, corners) == status::success);
        BOOST_TEST_EQ(corners.size(), p.outer().size()-1);
    };
    strict_source(a);
    BOOST_TEST_EQ(a.outer().size(), 104u);
    multi positive;
    bg::buffer(a, positive, bg::strategy::buffer::distance_symmetric<double>(9),
               bg::strategy::buffer::side_straight(), bg::strategy::buffer::join_round(36),
               bg::strategy::buffer::end_round(36), bg::strategy::buffer::point_circle(36));
    BOOST_TEST_EQ(positive.size(), 1u);
    if (positive.size() != 1) { return; }
    auto const& b = positive.front();
    BOOST_TEST(exact_valid(b));
    BOOST_TEST(b.inners().empty());
    strict_source(b);
    BOOST_TEST_EQ(b.outer().size(), 185u);
    multi out;
    double constexpr budget = 1e-9;
    auto const start = std::chrono::steady_clock::now();
    BOOST_TEST(bg::erode_convex_polygon_certified(b, out, 9, budget) == status::success);
    double const ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start).count();
    BOOST_TEST_EQ(out.size(), 1u);
    if (out.size() != 1) { return; }
    BOOST_TEST(exact_valid(out.front()));
    BOOST_TEST(out.front().inners().empty());
    BOOST_TEST_EQ(out.front().outer().size(), 162u);
    double const area = bg::area(out);
    BOOST_TEST(std::abs(area - 764.98706637614) < 1e-5);

    // Independent audit of the actual stored contour: each exported vertex
    // clears every intermediate side by at least 9 - the displacement budget.
    auto const& ring = b.outer();
    q const r = 9 - erosion::common::rational_value(budget);
    bool clears = true;
    for (std::size_t i = 1; i < ring.size(); ++i)
    {
        auto const p = erosion::common::read_point(ring[i-1]);
        auto const d = erosion::subtract(erosion::common::read_point(ring[i]),p);
        q const minimum2 = r*r*erosion::dot(d,d);
        for (auto const& v : out.front().outer())
        {
            auto const offset = erosion::subtract(erosion::common::read_point(v),p);
            q const left = erosion::cross(d,offset);
            clears &= left >= 0 && left*left >= minimum2;
        }
    }
    BOOST_TEST(clears);
    std::cout << std::setprecision(14) << "1281: " << a.outer().size()-1 << " original sides, "
              << b.outer().size()-1 << " intermediate sides, " << out.front().outer().size()-1
              << " output sides, area=" << area << ", displacement <= " << budget << ", " << ms << " ms\n";
}

int main(int argc, char** argv)
{
    test_irrational();
    test_three_root_sign();
    test_redundancy();
    test_collapse();
    test_invalid();
    if (argc == 2)
    {
        std::ifstream input(argv[1]);
        std::string wkt;
        if (!std::getline(input,wkt)) { BOOST_ERROR("Cannot read the supplied 1281 WKT"); }
        else { replay_1281(wkt.c_str()); }
    }
    else { replay_1281(issue_1281_wkt); }
    return boost::report_errors();
}
