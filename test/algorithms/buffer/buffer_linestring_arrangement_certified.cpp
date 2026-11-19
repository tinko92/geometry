// Boost.Geometry
// Unit Test
// Copyright (c) 2026
// Distributed under the Boost Software License, Version 1.0.

#include <boost/geometry/algorithms/buffer_linestring_arrangement_certified.hpp>
#include <boost/geometry/geometries/linestring.hpp>
#include <boost/core/lightweight_test.hpp>

#include <cstddef>
#include <vector>

namespace bg = boost::geometry;
namespace ec = bg::strategies::detail::exact_cartesian;
using number = ec::exact_number;
using status = bg::certified_buffer_status;
template <typename T> using point = bg::model::point<T, 2, bg::cs::cartesian>;
using p = point<double>;
using line = bg::model::linestring<p>;
using polygon = bg::model::polygon<p>;

struct xy { number x, y; };
struct edge { xy a, b; };

number q(int a, int b) { return number(a) / b; }
xy operator-(xy const& a, xy const& b) { return {a.x - b.x, a.y - b.y}; }
number det(xy const& a, xy const& b) { return a.x * b.y - a.y * b.x; }

template <typename T>
xy read(point<T> const& p) { return {ec::exact_value(bg::get<0>(p)), ec::exact_value(bg::get<1>(p))}; }

template <typename Range>
std::vector<edge> edges(Range const& range)
{
    std::vector<edge> result;
    for (std::size_t i = 1; i < range.size(); ++i) { result.push_back({read(range[i - 1]), read(range[i])}); }
    return result;
}

template <typename Range>
number area2(Range const& ring)
{
    number area = 0;
    for (auto const& e : edges(ring)) { area += det(e.a, e.b); }
    return area;
}

bool on(xy const& p, edge const& e)
{
    return det(p - e.a, e.b - e.a) == 0 && (p.x - e.a.x) * (p.x - e.b.x) <= 0
        && (p.y - e.a.y) * (p.y - e.b.y) <= 0;
}

bool meet(edge const& a, edge const& b)
{
    xy const u = a.b - a.a, v = b.b - b.a, w = b.a - a.a;
    number const d = det(u, v);
    if (d == 0) { return on(a.a, b) || on(a.b, b) || on(b.a, a) || on(b.b, a); }
    number const t = det(w, v) / d, s = det(w, u) / d;
    return 0 <= t && t <= 1 && 0 <= s && s <= 1;
}

template <typename Polygon>
std::vector<edge> boundary(Polygon const& polygon)
{
    auto result = edges(polygon.outer());
    for (auto const& h : polygon.inners())
    {
        auto const ring = edges(h);
        result.insert(result.end(), ring.begin(), ring.end());
    }
    return result;
}

int locate(std::vector<edge> const& boundary, xy const& p)
{
    bool inside = false;
    for (auto const& e : boundary)
    {
        if (on(p, e)) { return 0; }
        auto const& a = e.a;
        auto const& b = e.b;
        if ((a.y > p.y) != (b.y > p.y) && a.x + (b.x - a.x) * (p.y - a.y) / (b.y - a.y) > p.x)
        {
            inside = !inside;
        }
    }
    return inside ? 1 : -1;
}

void check_closed_hole(polygon const& poly)
{
    BOOST_TEST_EQ(poly.inners().size(), 1u);
    BOOST_TEST(area2(poly.outer()) < 0);
    if (poly.inners().size() == 1) { BOOST_TEST(area2(poly.inners().front()) > 0); }
    auto const check = [](auto const& ring)
    {
        BOOST_TEST(ring.size() >= 4);
        if (ring.empty()) { return; }
        BOOST_TEST_EQ(bg::get<0>(ring.front()), bg::get<0>(ring.back()));
        BOOST_TEST_EQ(bg::get<1>(ring.front()), bg::get<1>(ring.back()));
    };
    check(poly.outer());
    for (auto const& h : poly.inners()) { check(h); }
}

void test_open_irrational_loop()
{
    line const trace{p(0, 0), p(8, 8), p(0, 12), p(-8, 4), p(-0.5, 0.25)};
    auto const support = edges(trace);
    BOOST_TEST_EQ(support.size(), 4u);
    // Opposite parallel supports are sqrt(72) apart, so even their full
    // radius-1 capsules cannot share a witness. The open endpoint gap is
    // nevertheless covered, even at the smaller certified inner radius.
    xy const d = support[0].b - support[0].a, shift = support[2].a - support[0].a;
    number const cross = det(d, shift);
    BOOST_TEST(cross * cross > 4 * (d.x * d.x + d.y * d.y));
    BOOST_TEST(q(5, 16) < q(25, 64)); // gap^2 < (3/4 - 1/8)^2

    ec::round_buffer::bounded_radical_field f;
    ec::round_buffer::bounded_radical_field::value root;
    BOOST_TEST(f.try_sqrt(128, root));
    BOOST_TEST(f.try_sqrt(80, root));
    BOOST_TEST_EQ(f.rank(), 2u);

    number const e = q(1, 8);
    number const displacement = number(1) / (boost::multiprecision::cpp_int(1) << 42);
    polygon inner, outer;
    BOOST_TEST(bg::buffer_linestring_arrangement_certified(trace, outer, 1, e, displacement) == status::success);
    BOOST_TEST(bg::buffer_linestring_arrangement_certified(trace, inner, q(3, 4), e, displacement) == status::success);
    check_closed_hole(outer);
    check_closed_hole(inner);
    auto const a = boundary(inner), b = boundary(outer);
    BOOST_TEST_EQ(locate(b, {0, 6}), -1);  // hole
    BOOST_TEST_EQ(locate(b, {4, 4}), 1);   // on the source trace, buffered interior
    BOOST_TEST_EQ(locate(b, {15, 0}), -1);
    BOOST_TEST_EQ(locate(a, {q(-1, 4), q(1, 8)}), 1); // actual open-end gap is bridged

    // Sufficient inclusion with separate native displacement budgets;
    // also test directly against the stored rings, including their holes.
    BOOST_TEST(q(3, 4) + displacement < 1 - e - displacement);
    for (auto const& i : a)
    {
        BOOST_TEST_EQ(locate(b, i.a), 1);
        for (auto const& j : b) { BOOST_TEST_NOT(meet(i, j)); }
    }
    for (auto const& j : b) { BOOST_TEST_EQ(locate(a, j.a), -1); }

    auto const size = outer.outer().size();
    BOOST_TEST(bg::buffer_linestring_arrangement_certified(trace, outer, 1, e, 0) == status::cannot_certify);
    BOOST_TEST_EQ(outer.outer().size(), size);
    BOOST_TEST_EQ(outer.inners().size(), 1u);
    bg::model::polygon<point<float>> coarse;
    coarse.outer().emplace_back(1, 2);
    BOOST_TEST(bg::buffer_linestring_arrangement_certified(trace, coarse, 1, e, displacement)
               == status::cannot_certify);
    BOOST_TEST_EQ(coarse.outer().size(), 1u);
}

void test_exact_native_and_work_limit()
{
    // Quarter charts yield integer vertices; matching and opposing
    // collinear atoms near each joint must be noded without perturbation.
    line const square{p(0, 0), p(20, 0), p(20, 20), p(0, 20), p(0, 0)};
    polygon out;
    BOOST_TEST(bg::buffer_linestring_arrangement_certified(square, out, 5, 4, 0) == status::success);
    check_closed_hole(out);
    BOOST_TEST_EQ(area2(out.outer()), -1700);
    if (out.inners().size() == 1) { BOOST_TEST_EQ(area2(out.inners().front()), 200); }

    line many;
    for (int i = 0; i < 18; ++i) { many.emplace_back(i % 2, i % 2); }
    BOOST_TEST(bg::buffer_linestring_arrangement_certified(many, out, 1, q(1, 4), q(1, 1000))
               == status::cannot_certify);
    BOOST_TEST_EQ(area2(out.outer()), -1700);
}

void test_two_holes()
{
    line const trace{p(0, 0), p(20, 0), p(20, 20), p(0, 20), p(0, 0),
                     p(0, -20), p(-20, -20), p(-20, 0), p(0, 0)};
    polygon out;
    BOOST_TEST(bg::buffer_linestring_arrangement_certified(trace, out, 5, 4, 0) == status::success);
    BOOST_TEST_EQ(out.inners().size(), 2u);
    for (auto const& h : out.inners()) { BOOST_TEST_EQ(area2(h), 200); }
    auto const ring = boundary(out);
    BOOST_TEST_EQ(locate(ring, {10, 10}), -1);
    BOOST_TEST_EQ(locate(ring, {-10, -10}), -1);
    BOOST_TEST_EQ(locate(ring, {0, 0}), 1);
}

void test_touching_hole_refusal()
{
    line trace{p(0, 0), p(-10, 0), p(-10, 10), p(12, 10), p(12, 0), p(1.5, 0)};
    polygon out;
    BOOST_TEST(bg::buffer_linestring_arrangement_certified(trace, out, 1, q(1, 2), 0) == status::success);
    check_closed_hole(out);
    auto const area = area2(out.outer());
    // At gap two, the quarter-chart caps meet only at (1,0), pinching
    // the exterior and hole boundaries together in the exact arrangement.
    trace.back() = p(2, 0);
    BOOST_TEST(bg::buffer_linestring_arrangement_certified(trace, out, 1, q(1, 2), 0)
               == status::cannot_certify);
    BOOST_TEST_EQ(out.inners().size(), 1u);
    BOOST_TEST_EQ(area2(out.outer()), area);
}

int main()
{
    test_open_irrational_loop();
    test_exact_native_and_work_limit();
    test_two_holes();
    test_touching_hole_refusal();
    return boost::report_errors();
}
