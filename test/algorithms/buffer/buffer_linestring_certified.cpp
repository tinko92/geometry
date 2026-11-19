// Boost.Geometry
// Unit Test
// Copyright (c) 2026
// Distributed under the Boost Software License, Version 1.0.

#include <boost/geometry/algorithms/buffer_linestring_certified.hpp>

#include <boost/core/lightweight_test.hpp>
#include <boost/geometry/geometries/linestring.hpp>
#include <boost/geometry/geometries/segment.hpp>

#include <cstddef>

namespace bg = boost::geometry;
using status = bg::certified_buffer_status;
using number = boost::multiprecision::cpp_rational;
template <typename T> using point = bg::model::point<T, 2, bg::cs::cartesian>;
using p = point<number>;
using line = bg::model::linestring<p>;

number q(int n, int d) { return number(n) / d; }

struct xy { number x, y; };

template <typename T>
xy read(point<T> const& pt) { return {number(bg::get<0>(pt)), number(bg::get<1>(pt))}; }

xy operator-(xy const& a, xy const& b) { return {a.x - b.x, a.y - b.y}; }

number cross(xy const& a, xy const& b) { return a.x * b.y - a.y * b.x; }

bool on_segment(xy const& r, xy const& a, xy const& b)
{
    return cross(r - a, b - a) == 0
        && (r.x - a.x) * (r.x - b.x) <= 0
        && (r.y - a.y) * (r.y - b.y) <= 0;
}

bool segments_meet(xy const& a, xy const& b, xy const& c, xy const& d)
{
    xy const u = b - a, v = d - c, w = c - a;
    number const determinant = cross(u, v);
    if (determinant == 0)
    {
        return on_segment(a, c, d) || on_segment(b, c, d)
            || on_segment(c, a, b) || on_segment(d, a, b);
    }
    number const t = cross(w, v) / determinant, s = cross(w, u) / determinant;
    return t >= 0 && t <= 1 && s >= 0 && s <= 1;
}

template <typename Ring, typename Visitor>
void for_edges(Ring const& ring, Visitor const& visit)
{
    bool constexpr closed = bg::closure<Ring>::value == bg::closed;
    if (ring.size() < (closed ? 4u : 3u)) { return; }
    std::size_t const n = ring.size() - (closed ? 1 : 0);
    for (std::size_t i = 0; i < n; ++i)
    {
        visit(read(ring[i]), read(ring[(i + 1) % n]));
    }
}

template <typename Ring>
bool strictly_inside(Ring const& ring, xy const& r)
{
    bool inside = false, boundary = false;
    for_edges(ring, [&](xy const& a, xy const& b)
    {
        if (on_segment(r, a, b)) { boundary = true; }
        if ((a.y > r.y) != (b.y > r.y)
            && a.x + (b.x - a.x) * (r.y - a.y) / (b.y - a.y) > r.x)
        {
            inside = !inside;
        }
    });
    return inside && !boundary;
}

void test_578_nested_traces()
{
    // A direct reproduction of the inclusion property, with noncollinear
    // rational-length input edges (including the shortened Pythagorean legs).
    // This is a purpose-built fixture, not the issue's historical input.
    line const outer{p(-6, -8), p(0, 0), p(9, 0), p(12, 4), p(20, 4)};
    line const inner{p(-3, -4), p(0, 0), p(9, 0), p(q(21, 2), 2)};
    bg::model::polygon<p> larger;
    bg::model::polygon<p, false, false> smaller;
    number const outer_error = q(1, 16);
    number const inner_radius = 4 - outer_error;
    BOOST_TEST(bg::buffer_linestring_certified(outer, larger, 4, outer_error) == status::success);
    BOOST_TEST(bg::buffer_linestring_certified(inner, smaller, inner_radius, q(1, 8)) == status::success);
    BOOST_TEST(larger.outer().size() > 12 && smaller.outer().size() > 12);
    BOOST_TEST(larger.inners().empty() && smaller.inners().empty());
    if (larger.outer().empty() || smaller.outer().empty()) { return; }

    // Exact independent ray/crossing oracle for strict containment, not a
    // comparison of areas, floating vertex probes or per-capsule results.
    for (auto const& vertex : smaller.outer())
    {
        BOOST_TEST(strictly_inside(larger.outer(), read(vertex)));
    }
    for_edges(smaller.outer(), [&](xy const& a, xy const& b)
    {
        for_edges(larger.outer(), [&](xy const& c, xy const& d)
        {
            BOOST_TEST_NOT(segments_meet(a, b, c, d));
        });
    });
    BOOST_TEST(strictly_inside(larger.outer(), {17, 4}));
    BOOST_TEST_NOT(strictly_inside(smaller.outer(), {17, 4}));
}

void test_lossless_and_refusals()
{
    line const bent{p(0, 0), p(20, 0), p(20, 20)};
    bg::model::polygon<point<int>> integral;
    // An integral chord union is representable losslessly on the native grid.
    BOOST_TEST(bg::buffer_linestring_certified(bent, integral, 5, 4) == status::success);
    BOOST_TEST(strictly_inside(integral.outer(), {20, 10}));
    BOOST_TEST_NOT(strictly_inside(integral.outer(), {10, 10}));
    auto const size = integral.outer().size();
    BOOST_TEST(bg::buffer_linestring_certified(bent, integral, q(4, 3), q(1, 2)) == status::cannot_certify);
    BOOST_TEST_EQ(integral.outer().size(), size);

    bg::model::polygon<p> rational;
    rational.outer().emplace_back(123, 321);
    line const irrational{p(0, 0), p(2, 0), p(3, 1)};
    BOOST_TEST(bg::buffer_linestring_certified(irrational, rational, 1, q(1, 8)) == status::cannot_certify);
    line const loop{p(0, 0), p(20, 0), p(20, 20), p(0, 20), p(0, 0)};
    BOOST_TEST(bg::buffer_linestring_certified(loop, rational, 1, q(1, 8)) == status::cannot_certify);
    BOOST_TEST_EQ(rational.outer().size(), 1u);
    BOOST_TEST_EQ(bg::get<0>(rational.outer().front()), 123);

    line const backtracking{p(0, 0), p(0, 0), p(5, 0), p(0, 0), p(5, 0), p(5, 4)};
    BOOST_TEST(bg::buffer_linestring_certified(backtracking, rational, 1, q(1, 8)) == status::success);
    BOOST_TEST(strictly_inside(rational.outer(), {5, 2}));

    line const singleton{p(2, 3), p(2, 3)};
    BOOST_TEST(bg::buffer_linestring_certified(singleton, rational, 1, q(1, 4)) == status::success);
    BOOST_TEST(strictly_inside(rational.outer(), {2, 3}));
    BOOST_TEST_NOT(strictly_inside(rational.outer(), {4, 3}));
    auto const point_size = rational.outer().size();
    BOOST_TEST(bg::buffer_linestring_certified(line{}, rational, 1, q(1, 4)) == status::invalid_input);
    line many;
    for (int i = 0; i < 34; ++i) { many.emplace_back(2 * (i % 2), 0); }
    BOOST_TEST(bg::buffer_linestring_certified(many, rational, 5, 4) == status::cannot_certify);
    BOOST_TEST_EQ(rational.outer().size(), point_size);
}

void test_rounding_can_change_contact()
{
    // The exact capsules are separated by 1/12. Both individual integer
    // exports meet their geometric error budgets, yet share (4,1).
    bg::model::segment<p> const a(p(0, 0), p(8, 0)), b(p(0, q(11, 4)), p(8, q(11, 4)));
    number const radius = q(4, 3), error = q(2, 3);
    BOOST_TEST_EQ(q(11, 4) - 2 * radius, q(1, 12));
    bg::model::polygon<point<int>> ra, rb;
    BOOST_TEST(bg::buffer_segment_certified(a, ra, radius, error) == status::success);
    BOOST_TEST(bg::buffer_segment_certified(b, rb, radius, error) == status::success);
    for (auto const* poly : {&ra, &rb})
    {
        bool touches = false;
        for_edges(poly->outer(), [&](xy const& u, xy const& v)
        {
            touches = touches || on_segment({4, 1}, u, v);
        });
        BOOST_TEST(touches);
    }
    auto const size = ra.outer().size();
    BOOST_TEST(bg::buffer_linestring_certified(line{p(0, 0), p(8, 0)}, ra, radius, error)
              == status::cannot_certify);
    BOOST_TEST_EQ(ra.outer().size(), size);
}

void test_exact_atom_cases()
{
    namespace impl = bg::detail::certified_linestring_buffer;
    impl::piece const a{{0, 0}, {2, 0}, {2, 2}, {0, 2}};
    impl::piece const abutting{{2, 0}, {3, 0}, {3, 2}, {2, 2}};
    impl::piece const overlapping{{1, 0}, {4, 0}, {4, 2}, {1, 2}};
    impl::piece const point_contact{{2, 2}, {3, 2}, {3, 3}, {2, 3}};
    impl::piece const crossing{{1, -1}, {3, 1}, {1, 2}};
    auto area2 = [](impl::piece const& ring)
    {
        number result = 0;
        for (std::size_t i = 0; i < ring.size(); ++i)
        {
            auto const& a = ring[i];
            auto const& b = ring[(i + 1) % ring.size()];
            result += a.x * b.y - a.y * b.x;
        }
        return result;
    };
    impl::piece ring;
    BOOST_TEST(impl::union_one_ring({a, abutting}, ring));
    BOOST_TEST_EQ(area2(ring), 12);
    ring.clear();
    BOOST_TEST(impl::union_one_ring({a, overlapping, a}, ring));
    BOOST_TEST_EQ(area2(ring), 16);
    ring.clear();
    BOOST_TEST_NOT(impl::union_one_ring({a, point_contact}, ring));
    ring.clear();
    BOOST_TEST(impl::union_one_ring({a, crossing}, ring));
    bool rational_node = false;
    for (auto const& v : ring) { rational_node = rational_node || (v.x == 2 && v.y == q(3, 2)); }
    BOOST_TEST(rational_node);
}

int main()
{
    test_578_nested_traces();
    test_lossless_and_refusals();
    test_rounding_can_change_contact();
    test_exact_atom_cases();
    return boost::report_errors();
}
