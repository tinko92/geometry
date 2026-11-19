// Boost.Geometry
// Unit Test
// Copyright (c) 2026
// Distributed under the Boost Software License, Version 1.0.

#include <boost/geometry/algorithms/buffer_segment_certified.hpp>

#include <boost/core/lightweight_test.hpp>
#include <boost/geometry/geometries/linestring.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/segment.hpp>

#include <algorithm>
#include <cmath>
#include <limits>

namespace bg = boost::geometry;
using status = bg::certified_buffer_status;
using number = boost::multiprecision::cpp_rational;

template <typename T>
using point = bg::model::point<T, 2, bg::cs::cartesian>;

struct xy
{
    number x, y;
};

number q(int n, int d)
{
    return number(n) / d;
}

template <typename T>
xy read(point<T> const& p)
{
    return {number(bg::get<0>(p)), number(bg::get<1>(p))};
}

xy operator-(xy const& a, xy const& b)
{
    return {a.x - b.x, a.y - b.y};
}

number length2(xy const& a)
{
    return a.x * a.x + a.y * a.y;
}

number determinant(xy const& a, xy const& b)
{
    return a.x * b.y - a.y * b.x;
}

number oracle_point_segment(xy const& p, xy const& a, xy const& b)
{
    xy const v = b - a, w = p - a;
    number const d = length2(v);
    if (d == 0) { return length2(w); }
    number const t = (w.x * v.x + w.y * v.y) / d;
    number const s = (std::max)(number(0), (std::min)(number(1), t));
    return length2({w.x - s * v.x, w.y - s * v.y});
}

number oracle_segment_segment(xy const& a, xy const& b, xy const& p, xy const& q)
{
    // Solve for ray parameters; this does not use the production side tests.
    xy const u = b - a, v = q - p, w = p - a;
    number const d = determinant(u, v);
    if (d != 0)
    {
        number const t = determinant(w, v) / d, s = determinant(w, u) / d;
        if (0 <= t && t <= 1 && 0 <= s && s <= 1) { return 0; }
    }
    return (std::min)({oracle_point_segment(a, p, q), oracle_point_segment(b, p, q),
                       oracle_point_segment(p, a, b), oracle_point_segment(q, a, b)});
}

template <typename Polygon>
void check_export(Polygon const& out, xy const& a, xy const& b, number const& r, number const& error)
{
    auto const& ring = out.outer();
    bool constexpr closed = bg::closure<Polygon>::value == bg::closed;
    BOOST_TEST(out.inners().empty());
    BOOST_TEST(ring.size() >= (closed ? 4u : 3u));
    if (ring.size() < (closed ? 4u : 3u)) { return; }
    if constexpr (closed)
    {
        BOOST_TEST_EQ(length2(read(ring.front()) - read(ring.back())), 0);
    }
    std::size_t const size = ring.size() - (closed ? 1 : 0);
    number const inner = r > error ? number(r - error) : number(0);
    number const outer = r + error;
    number area2 = 0;
    bool inside_a = false, inside_b = false;
    auto const crosses_ray = [](xy const& p, xy const& u, xy const& v)
    {
        return (u.y > p.y) != (v.y > p.y)
            && u.x + (v.x - u.x) * (p.y - u.y) / (v.y - u.y) > p.x;
    };
    for (std::size_t i = 0; i < size; ++i)
    {
        xy const u = read(ring[i]), v = read(ring[(i + 1) % size]);
        BOOST_TEST(length2(v - u) > 0);
        BOOST_TEST(oracle_point_segment(u, a, b) < outer * outer);
        BOOST_TEST(oracle_segment_segment(a, b, u, v) > inner * inner);
        area2 += determinant(u, v);
        inside_a = inside_a != crosses_ray(a, u, v);
        inside_b = inside_b != crosses_ray(b, u, v);
    }
    BOOST_TEST(bg::point_order<Polygon>::value == bg::clockwise ? area2 < 0 : area2 > 0);
    BOOST_TEST(inside_a && inside_b);
}

void check_nearest_sqrt_half(double stored, int base, int sign)
{
    double const inf = (std::numeric_limits<double>::infinity)();
    number const lo = (number(std::nextafter(stored, -inf)) + number(stored)) / 2;
    number const hi = (number(std::nextafter(stored, inf)) + number(stored)) / 2;
    auto const cmp = [&](number const& limit)
    {
        number const diff = sign * (limit - base);
        return sign * (diff < 0 || diff * diff < q(1, 2) ? -1 : 1);
    };
    BOOST_TEST(cmp(lo) < 0 && cmp(hi) > 0);
}

void test_irrational_sides()
{
    bg::model::segment<point<int>> const segment(point<int>(0, 0), point<int>(7, 7));
    bg::model::polygon<point<double>, false, false> out;
    number const tolerance = q(1, 200);
    BOOST_TEST(bg::buffer_segment_certified(segment, out, 1, tolerance) == status::success);
    check_export(out, {0, 0}, {7, 7}, 1, tolerance);

    // The two straight side chords meet the semicircles at exact tangencies.
    // Verify all four emitted vertices against independent binary64 rounding
    // cells for base +/- sqrt(1/2), not against trig or sampled normals.
    auto const& ring = out.outer();
    if (ring.size() < 6) { return; }
    auto const check = [&](std::size_t i, int base, int sign)
    {
        check_nearest_sqrt_half(bg::get<0>(ring[i]), base, sign);
        check_nearest_sqrt_half(bg::get<1>(ring[i]), base, -sign);
    };
    check(0, 7, 1);                     // b + right normal
    check(ring.size() / 2 - 1, 7, -1);  // b + left normal
    check(ring.size() / 2, 0, -1);      // a + left normal
    check(ring.size() - 1, 0, 1);       // a + right normal

    bg::model::polygon<point<number>> rational;
    BOOST_TEST(bg::buffer_segment_certified(segment, rational, 1, tolerance) == status::success);
    check_export(rational, {0, 0}, {7, 7}, 1, tolerance);
}

void test_pythagorean_sides()
{
    point<number> const a(q(2, 11), q(-1, 13));
    point<number> const b(bg::get<0>(a) + q(3, 7), bg::get<1>(a) + q(4, 7));
    bg::model::segment<point<number>> const segment(a, b);
    bg::model::polygon<point<number>> clockwise;
    number const r = q(5, 6), error = q(1, 100);
    BOOST_TEST(bg::buffer_segment_certified(segment, clockwise, r, error) == status::success);
    check_export(clockwise, read(a), read(b), r, error);
    auto const& ring = clockwise.outer();
    if (ring.size() < 7) { return; }
    auto const check = [&](std::size_t i, point<number> const& base, int sign)
    {
        xy const delta = read(ring[i]) - read(base);
        BOOST_TEST_EQ(delta.x, sign * q(2, 3));
        BOOST_TEST_EQ(delta.y, -sign * q(1, 2));
    };
    std::size_t const n = ring.size() - 1;
    check(0, a, 1);
    check(n / 2 - 1, a, -1);
    check(n / 2, b, -1);
    check(n - 1, b, 1);
}

void test_reversed_oblique_segment()
{
    using p = point<number>;
    p const a(q(1, 7), q(-1, 3)), b(q(-11, 4), q(1, 8));
    bg::model::segment<p> const segment(a, b), reverse(b, a);
    bg::model::polygon<point<double>, false, false> forward, backward;
    number const r = q(2, 9), error = q(1, 300);
    BOOST_TEST(bg::buffer_segment_certified(segment, forward, r, error) == status::success);
    BOOST_TEST(bg::buffer_segment_certified(reverse, backward, r, error) == status::success);
    check_export(forward, read(a), read(b), r, error);
    check_export(backward, read(a), read(b), r, error);
    std::size_t const n = forward.outer().size();
    BOOST_TEST_EQ(n, backward.outer().size());
    if (n != backward.outer().size()) { return; }
    for (std::size_t i = 0; i < n; ++i)
    {
        BOOST_TEST_EQ(length2(read(forward.outer()[i]) - read(backward.outer()[(i + n / 2) % n])), 0);
    }
}

void test_enclosure_refinement()
{
    using boost::multiprecision::cpp_int;
    namespace impl = bg::detail::certified_segment_buffer;
    cpp_int p = 1, q = 1;
    for (unsigned i = 0; i < 90; ++i)
    {
        cpp_int const next_q = 2 * p + q;
        p += q;
        q = next_q;
    }
    for (unsigned i = 0; i < 2; ++i)
    {
        // Alternating Pell approximants to 1/sqrt(2) put the true value on
        // either side of a rounding midpoint by much less than 2^-96.
        impl::inverse_length_enclosure scale(2);
        int rounded = -1;
        BOOST_TEST(scale.round(number(1) / 2 - number(p) / q, 1, rounded));
        BOOST_TEST_EQ(rounded, 2 * p * p < q * q ? 1 : 0);
        cpp_int const next_q = 2 * p + q;
        p += q;
        q = next_q;
    }
}

void test_export_failures()
{
    bg::model::segment<point<int>> const horizontal(point<int>(0, 0), point<int>(4, 0));
    bg::model::polygon<point<int>, false, false> integral;
    number const r = q(4, 3), exact_error = q(1, 3), above = q(1, 3) + q(1, 1000);
    BOOST_TEST(bg::buffer_segment_certified(horizontal, integral, r, above) == status::success);
    check_export(integral, {0, 0}, {4, 0}, r, above);
    BOOST_TEST_EQ(integral.outer().size(), 10u);
    BOOST_TEST(bg::buffer_segment_certified(horizontal, integral, r, exact_error) == status::cannot_certify);
    BOOST_TEST_EQ(integral.outer().size(), 10u);

    double const x = 16777216;
    bg::model::segment<point<double>> const translated(point<double>(x, x), point<double>(x + 4, x + 4));
    bg::model::polygon<point<float>> coarse_grid;
    coarse_grid.outer().emplace_back(5, 6);
    BOOST_TEST(bg::buffer_segment_certified(translated, coarse_grid, q(1, 4), q(1, 16)) == status::cannot_certify);
    BOOST_TEST_EQ(coarse_grid.outer().size(), 1u);
    BOOST_TEST_EQ(bg::get<0>(coarse_grid.outer().front()), 5);

    BOOST_TEST(bg::buffer_segment_certified(horizontal, integral, r, 0) == status::cannot_certify);
    BOOST_TEST(bg::buffer_segment_certified(horizontal, integral, -1, above) == status::invalid_input);
    double const inf = (std::numeric_limits<double>::infinity)();
    bg::model::segment<point<double>> const nonfinite(point<double>(inf, 0), point<double>(2, 1));
    BOOST_TEST(bg::buffer_segment_certified(nonfinite, integral, r, above) == status::invalid_input);
    BOOST_TEST_EQ(integral.outer().size(), 10u);
}

void test_proper_crossings()
{
    namespace impl = bg::detail::certified_segment_buffer;
    bg::model::polygon<point<int>, false, false> ring;
    ring.outer().emplace_back(-1, -1);
    ring.outer().emplace_back(1, -1);
    ring.outer().emplace_back(1, 1);
    ring.outer().emplace_back(-1, 1);
    BOOST_TEST(impl::certify_export(ring, {-q(1, 2), 0}, {q(1, 2), 0}, 1, 2, 4));
    // No vertex or source endpoint touches the other segment. Only the
    // interior crossings reveal that the source exits the polygon.
    BOOST_TEST_NOT(impl::certify_export(ring, {-3, 0}, {3, 0}, 1, 2, 4));
}

void test_collinear_trace()
{
    using p = point<number>;
    bg::model::linestring<p> const backtrack{p(1, 1), p(1, 1), p(2, 2), p(-3, -3), p(5, 5), p(0, 0)};
    bg::model::segment<p> const span(p(-3, -3), p(5, 5));
    bg::model::polygon<p> from_line, from_span;
    number const r = q(3, 5), e = q(1, 25);
    BOOST_TEST(bg::buffer_collinear_linestring_certified(backtrack, from_line, r, e) == status::success);
    BOOST_TEST(bg::buffer_segment_certified(span, from_span, r, e) == status::success);
    check_export(from_line, {-3, -3}, {5, 5}, r, e);
    BOOST_TEST_EQ(from_span.outer().size(), from_line.outer().size());
    for (std::size_t i = 0; i < (std::min)(from_span.outer().size(), from_line.outer().size()); ++i)
    {
        BOOST_TEST_EQ(length2(read(from_span.outer()[i]) - read(from_line.outer()[i])), 0);
    }
    auto const size = from_line.outer().size();
    bg::model::linestring<p> const bent{p(0, 0), p(2, 2), p(2, 3)};
    number const tiny = number(1) / (boost::multiprecision::cpp_int(1) << 100);
    bg::model::linestring<p> const almost{p(0, 0), p(1, 1), p(2, number(2) + tiny)};
    BOOST_TEST(bg::buffer_collinear_linestring_certified(bent, from_line, r, e) == status::cannot_certify);
    BOOST_TEST(bg::buffer_collinear_linestring_certified(almost, from_line, r, e) == status::cannot_certify);
    BOOST_TEST(bg::buffer_collinear_linestring_certified(bg::model::linestring<p>{}, from_line, r, e)
              == status::invalid_input);
    BOOST_TEST_EQ(from_line.outer().size(), size);

    bg::model::linestring<p> const singleton{p(2, 3), p(2, 3)};
    bg::model::segment<p> const zero(p(2, 3), p(2, 3));
    BOOST_TEST(bg::buffer_collinear_linestring_certified(singleton, from_line, r, e) == status::success);
    BOOST_TEST(bg::buffer_segment_certified(zero, from_span, r, e) == status::success);
    check_export(from_line, {2, 3}, {2, 3}, r, e);
    BOOST_TEST_EQ(from_line.outer().size(), from_span.outer().size());
}

int main()
{
    test_irrational_sides();
    test_pythagorean_sides();
    test_reversed_oblique_segment();
    test_enclosure_refinement();
    test_export_failures();
    test_proper_crossings();
    test_collinear_trace();
    return boost::report_errors();
}
