// Boost.Geometry
// Unit Test
// Copyright (c) 2026
// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include <boost/geometry/algorithms/buffer_point_certified.hpp>

#include <boost/core/lightweight_test.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/polygon.hpp>

#include <limits>

namespace bg = boost::geometry;
using status = bg::certified_point_buffer_status;
using number = boost::multiprecision::cpp_rational;

template <typename T>
using point = bg::model::point<T, 2, bg::cs::cartesian>;

number q(int numerator, int denominator)
{
    return number(numerator) / denominator;
}

template <typename Polygon, typename T>
bool check_export(Polygon const& out, point<T> const& center, number const& radius, number const& error)
{
    // Independent oracle on stored coordinates: minimize the squared norm
    // using the clamped segment parameter (not the production cross formula).
    auto const& ring = out.outer();
    bool constexpr closed = bg::closure<Polygon>::value == bg::closed;
    BOOST_TEST(ring.size() >= (closed ? 4u : 3u));
    BOOST_TEST(out.inners().empty());
    if (ring.size() < (closed ? 4u : 3u)) { return false; }
    if constexpr (closed)
    {
        BOOST_TEST_EQ(bg::get<0>(ring.front()), bg::get<0>(ring.back()));
        BOOST_TEST_EQ(bg::get<1>(ring.front()), bg::get<1>(ring.back()));
    }
    number const cx(bg::get<0>(center)), cy(bg::get<1>(center));
    number const inner = radius - error, outer = radius + error;
    BOOST_TEST(inner > 0);
    std::size_t const size = ring.size() - (closed ? 1 : 0);
    bool rounded = false;
    number twice_area = 0;
    for (std::size_t i = 0; i < size; ++i)
    {
        auto const& a = ring[i];
        auto const& b = ring[(i + 1) % size];
        number const x = number(bg::get<0>(a)) - cx, y = number(bg::get<1>(a)) - cy;
        number const u = number(bg::get<0>(b)) - cx, v = number(bg::get<1>(b)) - cy;
        number const dx = u - x, dy = v - y, len2 = dx * dx + dy * dy;
        BOOST_TEST(len2 > 0);
        if (len2 == 0) { continue; }
        number t = -(x * dx + y * dy) / len2;
        if (t < 0) { t = 0; }
        if (t > 1) { t = 1; }
        number const px = x + t * dx, py = y + t * dy;
        BOOST_TEST(px * px + py * py > inner * inner);
        number const vertex2 = x * x + y * y;
        BOOST_TEST(vertex2 < outer * outer);
        rounded = rounded || vertex2 != radius * radius;
        twice_area += x * v - y * u;
    }
    BOOST_TEST(bg::point_order<Polygon>::value == bg::clockwise ? twice_area < 0 : twice_area > 0);
    return rounded;
}

void test_exact_and_subdivision()
{
    point<number> const center(q(2, 7), q(-3, 11));
    bg::model::polygon<point<number>, false, false> open;
    bg::model::polygon<point<number>> coarse, fine;
    coarse.inners().resize(1);
    number const r = q(5, 3), loose = q(1, 2), tight = q(1, 200);

    BOOST_TEST(bg::buffer_point_certified(center, coarse, r, loose) == status::success);
    BOOST_TEST(bg::buffer_point_certified(center, fine, r, tight) == status::success);
    BOOST_TEST(bg::buffer_point_certified(center, open, r, tight) == status::success);
    BOOST_TEST(fine.outer().size() > coarse.outer().size());
    BOOST_TEST_EQ(fine.outer().size(), open.outer().size() + 1);
    BOOST_TEST_NOT(check_export(coarse, center, r, loose));
    BOOST_TEST_NOT(check_export(fine, center, r, tight));
    BOOST_TEST_NOT(check_export(open, center, r, tight));
}

void test_rounded_success()
{
    point<double> const center(1073741824.125, -1073741824.25);
    bg::model::polygon<point<double>> out;
    number const r = q(7, 3), error = q(1, 40);
    BOOST_TEST(bg::buffer_point_certified(center, out, r, error) == status::success);
    BOOST_TEST(check_export(out, center, r, error));

    point<int> const origin(0, 0);
    bg::model::polygon<point<int>, false, false> integral;
    // The rounded polygon is a square with extra axis vertices. Its exact
    // maximum inward error is 1/3. The strict certificate must reject equality
    // and smaller requests, and accept a slightly larger one.
    number const radius = q(4, 3), equality = q(1, 3);
    number const above = equality + q(1, 1000), below = equality - q(1, 1000);
    BOOST_TEST(bg::buffer_point_certified(origin, integral, radius, above) == status::success);
    BOOST_TEST_EQ(integral.outer().size(), 8u);
    BOOST_TEST(check_export(integral, origin, radius, above));
    BOOST_TEST(bg::buffer_point_certified(origin, integral, radius, equality) == status::cannot_certify);
    BOOST_TEST(bg::buffer_point_certified(origin, integral, radius, below) == status::cannot_certify);
    BOOST_TEST_EQ(integral.outer().size(), 8u);
}

void test_edge_interior_failure()
{
    point<double> const center(0.5, 0);
    bg::model::polygon<point<int>, false, false> out;
    number const radius = q(9, 4), loose = q(7, 10), tight = q(1, 2);
    // Both tolerances select the same eight exact chords. Every rounded
    // vertex passes the tighter annulus, but the edge (0,-2) -> (2,-1) has
    // minimum squared radius 49/20, below (r - tight)^2 = 49/16.
    BOOST_TEST(bg::buffer_point_certified(center, out, radius, loose) == status::success);
    BOOST_TEST_EQ(out.outer().size(), 8u);
    BOOST_TEST(check_export(out, center, radius, loose));
    for (auto const& p : out.outer())
    {
        number const x = number(bg::get<0>(p)) - q(1, 2), y = bg::get<1>(p);
        BOOST_TEST(x * x + y * y > q(49, 16));
        BOOST_TEST(x * x + y * y < q(121, 16));
    }
    BOOST_TEST(bg::buffer_point_certified(center, out, radius, tight) == status::cannot_certify);
    BOOST_TEST_EQ(out.outer().size(), 8u);
}

void test_export_winding_guard()
{
    namespace impl = bg::detail::certified_point_buffer;
    bg::model::polygon<point<int>, false, false> out;
    out.outer().emplace_back(1, 0);
    out.outer().emplace_back(0, 1);
    out.outer().emplace_back(-1, 0);
    out.outer().emplace_back(0, -1);
    BOOST_TEST(impl::certify_export(out, {0, 0}, 1, q(1, 2), 4));

    // All consecutive gaps and radial bounds still pass, but a repeated
    // traversal of the diamond is not a simple boundary.
    auto const once = out.outer();
    out.outer().insert(out.outer().end(), once.begin(), once.end());
    BOOST_TEST_NOT(impl::certify_export(out, {0, 0}, 1, q(1, 2), 8));
}

template <typename T>
bg::model::polygon<point<T>> sentinel()
{
    bg::model::polygon<point<T>> out;
    out.outer().emplace_back(8, -3);
    out.inners().resize(1);
    out.inners().back().emplace_back(-2, 7);
    return out;
}

template <typename Polygon>
void check_unchanged(Polygon const& out)
{
    BOOST_TEST_EQ(out.outer().size(), 1u);
    BOOST_TEST_EQ(out.inners().size(), 1u);
    if (out.outer().empty() || out.inners().empty()) { return; }
    BOOST_TEST_EQ(out.inners().front().size(), 1u);
    BOOST_TEST_EQ(bg::get<0>(out.outer().front()), 8);
    BOOST_TEST_EQ(bg::get<1>(out.outer().front()), -3);
    if (out.inners().front().empty()) { return; }
    BOOST_TEST_EQ(bg::get<0>(out.inners().front().front()), -2);
    BOOST_TEST_EQ(bg::get<1>(out.inners().front().front()), 7);
}

void test_grid_failures()
{
    auto integer = sentinel<int>();
    point<int> const origin(0, 0);
    BOOST_TEST(bg::buffer_point_certified(origin, integer, q(1, 2), q(1, 8)) == status::cannot_certify);
    check_unchanged(integer);
    point<int> const near_limit((std::numeric_limits<int>::max)(), 0);
    BOOST_TEST(bg::buffer_point_certified(near_limit, integer, 2, q(1, 2)) == status::cannot_certify);
    check_unchanged(integer);

    auto floating = sentinel<float>();
    point<double> const translated(16777216.0, 16777216.0);
    BOOST_TEST(bg::buffer_point_certified(translated, floating, q(1, 4), q(1, 16)) == status::cannot_certify);
    check_unchanged(floating);

    auto doubles = sentinel<double>();
    double const largest = (std::numeric_limits<double>::max)();
    BOOST_TEST(bg::buffer_point_certified(point<double>(largest, 0), doubles, largest, largest)
              == status::cannot_certify);
    check_unchanged(doubles);
}

void test_invalid_and_limits()
{
    point<double> const center(0, 0);
    auto out = sentinel<double>();
    double const inf = (std::numeric_limits<double>::infinity)();
    number const tiny = number(1) / (boost::multiprecision::cpp_int(1) << 100);
    BOOST_TEST(bg::buffer_point_certified(center, out, 0, 1) == status::invalid_input);
    BOOST_TEST(bg::buffer_point_certified(center, out, 1, -1) == status::invalid_input);
    BOOST_TEST(bg::buffer_point_certified(center, out, inf, 1) == status::invalid_input);
    BOOST_TEST(bg::buffer_point_certified(center, out, 1, inf) == status::invalid_input);
    BOOST_TEST(bg::buffer_point_certified(point<double>(0, inf), out, 1, 1) == status::invalid_input);
    BOOST_TEST(bg::buffer_point_certified(center, out, 1, 0) == status::cannot_certify);
    BOOST_TEST(bg::buffer_point_certified(center, out, 1, tiny) == status::cannot_certify);
    BOOST_TEST(bg::buffer_point_certified(center, out, 1, q(1, 10000000)) == status::cannot_certify);
    check_unchanged(out);
}

int main()
{
    test_exact_and_subdivision();
    test_rounded_success();
    test_edge_interior_failure();
    test_export_winding_guard();
    test_grid_failures();
    test_invalid_and_limits();
    return boost::report_errors();
}
