// Boost.Geometry
// Unit Test
// Copyright (c) 2026
// Distributed under the Boost Software License, Version 1.0.

#include <boost/geometry/algorithms/buffer_polygon_certified.hpp>

#include <boost/core/lightweight_test.hpp>
#include <boost/geometry/geometries/multi_polygon.hpp>

#include <algorithm>
#include <cstddef>
#include <initializer_list>

namespace bg = boost::geometry;
using status = bg::certified_buffer_status;
using mode = bg::certified_polygon_buffer_mode;
using number = boost::multiprecision::cpp_rational;
template <typename T> using point = bg::model::point<T, 2, bg::cs::cartesian>;
using p = point<number>;
using source_polygon = bg::model::polygon<p, false, true>;
using mp = bg::model::multi_polygon<bg::model::polygon<p>>;
using open_mp = bg::model::multi_polygon<bg::model::polygon<p, false, false>>;

number q(int n, int d) { return number(n) / d; }

struct xy { number x, y; };

template <typename T>
xy read(point<T> const& a) { return {number(bg::get<0>(a)), number(bg::get<1>(a))}; }

template <typename Ring>
void set_closed(Ring& ring, std::initializer_list<xy> vertices)
{
    for (auto const& v : vertices) { ring.emplace_back(v.x, v.y); }
    ring.push_back(ring.front());
}

source_polygon polygon(std::initializer_list<xy> vertices)
{
    source_polygon out;
    set_closed(out.outer(), vertices);
    return out;
}

template <typename Ring, typename Visitor>
void visit_edges(Ring const& ring, Visitor const& visit)
{
    bool constexpr closed = bg::closure<Ring>::value == bg::closed;
    if (ring.size() < (closed ? 4u : 3u)) { return; }
    std::size_t const size = ring.size() - (closed ? 1 : 0);
    for (std::size_t i = 0; i < size; ++i)
    {
        visit(read(ring[i]), read(ring[(i + 1) % size]));
    }
}

template <typename Ring>
int locate(Ring const& ring, xy const& r)
{
    // Independent test oracle: compute the ray intersection coordinate,
    // rather than interpreting the production determinant sign.
    bool inside = false, boundary = false;
    visit_edges(ring, [&](xy const& a, xy const& b)
    {
        if ((r.x - a.x) * (r.x - b.x) <= 0 && (r.y - a.y) * (r.y - b.y) <= 0
            && (r.x - a.x) * (b.y - a.y) == (r.y - a.y) * (b.x - a.x))
        {
            boundary = true;
        }
        if ((a.y > r.y) != (b.y > r.y)
            && a.x + (b.x - a.x) * (r.y - a.y) / (b.y - a.y) > r.x)
        {
            inside = !inside;
        }
    });
    return boundary ? 0 : inside ? 1 : -1;
}

template <typename Multi>
int locate_multi(Multi const& result, xy const& r)
{
    for (auto const& poly : result)
    {
        auto const in = locate(poly.outer(), r);
        if (in == 0) { return 0; }
        if (in < 0) { continue; }
        bool inside = true;
        for (auto const& h : poly.inners())
        {
            auto const hole = locate(h, r);
            if (hole == 0) { return 0; }
            if (hole > 0) { inside = false; }
        }
        if (inside) { return 1; }
    }
    return -1;
}

template <typename Ring>
number area2(Ring const& ring)
{
    number area = 0;
    visit_edges(ring, [&](xy const& a, xy const& b) { area += a.x * b.y - a.y * b.x; });
    return area;
}

template <typename Multi>
void check_orientation(Multi const& result)
{
    using poly = typename Multi::value_type;
    bool constexpr closed = bg::closure<poly>::value == bg::closed;
    int constexpr sign = bg::point_order<poly>::value == bg::clockwise ? -1 : 1;
    auto const check = [&](auto const& ring, int order)
    {
        BOOST_TEST(ring.size() >= (closed ? 4u : 3u));
        BOOST_TEST(area2(ring) * order > 0);
        if (ring.empty()) { return; }
        bool const endpoints_equal = bg::get<0>(ring.front()) == bg::get<0>(ring.back())
                                  && bg::get<1>(ring.front()) == bg::get<1>(ring.back());
        BOOST_TEST_EQ(endpoints_equal, closed);
    };
    for (auto const& poly : result)
    {
        check(poly.outer(), sign);
        for (auto const& h : poly.inners()) { check(h, -sign); }
    }
}

source_polygon square_with_hole()
{
    auto out = polygon({{0, 0}, {20, 0}, {20, 20}, {0, 20}});
    out.inners().resize(1);
    set_closed(out.inners().back(), {{6, 6}, {6, 14}, {14, 14}, {14, 6}});
    return out;
}

void test_shell_and_hole_dilation()
{
    auto const input = square_with_hole();
    mp out;
    BOOST_TEST(bg::buffer_polygon_certified(input, out, 2, q(1, 2)) == status::success);
    check_orientation(out);
    BOOST_TEST_EQ(out.size(), 1u);
    if (out.size() != 1) { return; }
    BOOST_TEST_EQ(out.front().inners().size(), 1u);
    if (out.front().inners().size() == 1)
    {
        // Straight exact supports give the shrunken [8,12]^2 hole.
        BOOST_TEST_EQ(area2(out.front().inners().front()), 32);
    }
    BOOST_TEST_EQ(locate_multi(out, {-q(3, 2), 10}), 1);
    BOOST_TEST_EQ(locate_multi(out, {-q(9, 4), 10}), -1);
    BOOST_TEST_EQ(locate_multi(out, {8, 10}), 0);
    BOOST_TEST_EQ(locate_multi(out, {10, 10}), -1);
    BOOST_TEST_EQ(locate_multi(out, {7, 10}), 1);

    // At exactly the inradius the overlapping band leaves no zero-area hole.
    BOOST_TEST(bg::buffer_polygon_certified(input, out, 4, 1) == status::success);
    BOOST_TEST_EQ(out.size(), 1u);
    if (out.size() == 1) { BOOST_TEST(out.front().inners().empty()); }
    BOOST_TEST_EQ(locate_multi(out, {10, 10}), 1);
}

void test_polygon_boundary_band()
{
    auto const input = square_with_hole();
    open_mp out;
    BOOST_TEST(bg::buffer_polygon_certified(input, out, 2, q(1, 2), mode::boundary_band) == status::success);
    check_orientation(out);
    // The outer-shell band and the hole band are disjoint nested components,
    // each with a hole. The inner hole belongs to the innermost shell.
    BOOST_TEST_EQ(out.size(), 2u);
    for (auto const& poly : out) { BOOST_TEST_EQ(poly.inners().size(), 1u); }
    BOOST_TEST_EQ(locate_multi(out, {0, 10}), 1);
    BOOST_TEST_EQ(locate_multi(out, {6, 10}), 1);
    BOOST_TEST_EQ(locate_multi(out, {3, 10}), -1);
    BOOST_TEST_EQ(locate_multi(out, {10, 10}), -1);
}

void test_erosion_with_hole_and_split()
{
    auto const input = square_with_hole();
    mp out;
    BOOST_TEST(bg::buffer_polygon_certified(input, out, 2, q(1, 2), mode::erosion) == status::success);
    BOOST_TEST_EQ(out.size(), 1u);
    if (out.size() == 1)
    {
        BOOST_TEST_EQ(out.front().inners().size(), 1u);
        BOOST_TEST_EQ(area2(out.front().outer()), -512); // [2,18]^2, exact straight sides
    }
    check_orientation(out);
    BOOST_TEST_EQ(locate_multi(out, {2, 10}), 0);
    BOOST_TEST_EQ(locate_multi(out, {3, 10}), 1);
    BOOST_TEST_EQ(locate_multi(out, {5, 10}), -1);

    auto const dumbbell = polygon({{0, 0}, {6, 0}, {6, 2}, {10, 2}, {10, 0}, {16, 0},
                                  {16, 6}, {10, 6}, {10, 4}, {6, 4}, {6, 6}, {0, 6}});
    // The neck has half width 1, less than r-e = 5/4; both chamber
    // centers have clearance 3 > r. The certificate must separate them.
    BOOST_TEST(bg::buffer_polygon_certified(dumbbell, out, q(3, 2), q(1, 4), mode::erosion)
              == status::success);
    BOOST_TEST_EQ(out.size(), 2u);
    check_orientation(out);
    for (auto const& poly : out)
    {
        BOOST_TEST(poly.inners().empty());
        if (poly.outer().empty()) { continue; }
        bool const left = bg::get<0>(poly.outer().front()) < 8;
        for (auto const& v : poly.outer()) { BOOST_TEST(left ? bg::get<0>(v) < 8 : bg::get<0>(v) > 8); }
    }
    BOOST_TEST_EQ(locate_multi(out, {3, 3}), 1);
    BOOST_TEST_EQ(locate_multi(out, {13, 3}), 1);
    BOOST_TEST_EQ(locate_multi(out, {8, 3}), -1);

    auto const flat = polygon({{0, 0}, {4, 0}, {4, 2}, {0, 2}});
    BOOST_TEST(bg::buffer_polygon_certified(flat, out, 1, q(1, 4), mode::erosion) == status::success);
    BOOST_TEST(out.empty()); // true radius-1 closed erosion is only a segment
}

void test_oblique_support_and_refusals()
{
    auto const triangle = polygon({{0, 0}, {4, 0}, {0, 3}});
    mp out;
    BOOST_TEST(bg::buffer_polygon_certified(triangle, out, q(1, 2), q(1, 8)) == status::success);
    BOOST_TEST_EQ(out.size(), 1u);
    // At the hypotenuse tangent the exact outward normal is (3/5,4/5).
    BOOST_TEST_EQ(locate_multi(out, {q(43, 10), q(2, 5)}), 0);
    BOOST_TEST_EQ(locate_multi(out, {1, 1}), 1);
    BOOST_TEST_EQ(locate_multi(out, {5, 0}), -1);

    bg::model::multi_polygon<bg::model::polygon<point<int>>> integral;
    integral.resize(1);
    integral.front().outer().emplace_back(9, 7);
    BOOST_TEST(bg::buffer_polygon_certified(triangle, integral, q(1, 2), q(1, 8)) == status::cannot_certify);
    BOOST_TEST_EQ(integral.size(), 1u);
    BOOST_TEST_EQ(integral.front().outer().size(), 1u);

    auto const irrational = polygon({{0, 0}, {3, 0}, {0, 2}});
    BOOST_TEST(bg::buffer_polygon_certified(irrational, out, 1, q(1, 4)) == status::cannot_certify);
    auto touching = polygon({{0, 0}, {20, 0}, {20, 20}, {0, 20}});
    touching.inners().resize(1);
    set_closed(touching.inners().back(), {{0, 4}, {0, 8}, {4, 8}, {4, 4}});
    BOOST_TEST(bg::buffer_polygon_certified(touching, out, 1, q(1, 4)) == status::invalid_input);
    auto const crossing = polygon({{0, 0}, {5, 0}, {1, 5}, {1, -1}, {5, 4}, {0, 5}});
    BOOST_TEST(bg::buffer_polygon_certified(crossing, out, 1, q(1, 4)) == status::invalid_input);
    auto not_closed = triangle;
    not_closed.outer().pop_back();
    BOOST_TEST(bg::buffer_polygon_certified(not_closed, out, 1, q(1, 4)) == status::invalid_input);
    BOOST_TEST_EQ(locate_multi(out, {q(43, 10), q(2, 5)}), 0);
}

void test_exact_convex_erosion()
{
    auto const triangle = polygon({{0, 0}, {4, 0}, {0, 3}});
    mp out;
    BOOST_TEST(bg::buffer_polygon_certified(triangle, out, q(1, 2), q(1, 8), mode::erosion)
              == status::success);
    BOOST_TEST_EQ(out.size(), 1u);
    if (out.size() == 1)
    {
        BOOST_TEST(out.front().inners().empty());
        BOOST_TEST_EQ(area2(out.front().outer()), -3);
        // Independent inward halfplanes: x>=1/2, y>=1/2, 3x+4y<=19/2.
        for (auto const& v : out.front().outer())
        {
            auto const p = read(v);
            BOOST_TEST(p.x >= q(1, 2) && p.y >= q(1, 2) && 3 * p.x + 4 * p.y <= q(19, 2));
        }
    }
    BOOST_TEST_EQ(locate_multi(out, {q(5, 2), q(1, 2)}), 0);
    BOOST_TEST_EQ(locate_multi(out, {q(1, 2), 2}), 0);
    BOOST_TEST_EQ(locate_multi(out, {1, 1}), 1);
    BOOST_TEST(bg::buffer_polygon_certified(triangle, out, 1, q(1, 4), mode::erosion)
              == status::success);
    BOOST_TEST(out.empty()); // at the inradius the regularized result is empty
}

void test_rational_nodes()
{
    namespace arr = bg::detail::rational_chord_arrangement;
    arr::edge const diagonal({0, 0}, {3, 3});
    arr::edge const horizontal({0, q(4, 3)}, {3, q(4, 3)});
    arr::edge const overlap({2, 2}, {1, 1});
    arr::intersection hit;
    BOOST_TEST(arr::intersect(diagonal, horizontal, hit));
    BOOST_TEST_EQ(hit.t0, q(4, 9));
    BOOST_TEST_EQ(hit.u1, q(4, 9));
    BOOST_TEST(arr::intersect(diagonal, overlap, hit));
    BOOST_TEST_EQ(hit.t0, q(1, 3));
    BOOST_TEST_EQ(hit.t1, q(2, 3));
    BOOST_TEST_EQ(hit.u0, 1);
    BOOST_TEST_EQ(hit.u1, 0);
}

int main()
{
    test_shell_and_hole_dilation();
    test_polygon_boundary_band();
    test_erosion_with_hole_and_split();
    test_oblique_support_and_refusals();
    test_exact_convex_erosion();
    test_rational_nodes();
    return boost::report_errors();
}
