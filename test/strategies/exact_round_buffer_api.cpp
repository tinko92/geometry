// Boost.Geometry
// Unit Test
// Copyright (c) 2026
// Distributed under the Boost Software License, Version 1.0.

#include "exact_issue_1281_fixture.hpp"

#include <boost/geometry/algorithms/buffer_round_certified.hpp>
#include <boost/geometry/geometries/multi_linestring.hpp>
#include <boost/geometry/geometries/multi_point.hpp>
#include <boost/geometry/geometries/segment.hpp>
#include <boost/geometry/io/wkt/read.hpp>
#include <boost/core/lightweight_test.hpp>

#include <algorithm>
#include <limits>
#include <vector>

namespace bg = boost::geometry;
namespace ec = bg::strategies::detail::exact_cartesian;
using number = ec::exact_number;
using status = bg::certified_buffer_status;
template <typename T> using point = bg::model::point<T, 2, bg::cs::cartesian>;
using p = point<double>;
using mp = bg::model::multi_point<p>;
using line = bg::model::linestring<p>;
using ml = bg::model::multi_linestring<line>;
using segment = bg::model::segment<p>;
using polygon = bg::model::polygon<p>;
using multi = bg::model::multi_polygon<polygon>;
using rational_polygon = bg::model::polygon<point<number>, true, false>;

number lift(number const& value) { return value; }
number lift(double value) { return ec::exact_value(value); }

template <typename Polygon>
bool same_storage(Polygon const& a, Polygon const& b)
{
    auto const same_ring = [](auto const& a, auto const& b)
    {
        return std::equal(a.begin(), a.end(), b.begin(), b.end(), [](auto const& a, auto const& b)
        {
            return bg::get<0>(a) == bg::get<0>(b) && bg::get<1>(a) == bg::get<1>(b);
        });
    };
    return same_ring(a.outer(), b.outer())
        && std::equal(a.inners().begin(), a.inners().end(), b.inners().begin(), b.inners().end(), same_ring);
}

template <typename Multi>
bool same_multi(Multi const& a, Multi const& b)
{
    return std::equal(a.begin(), a.end(), b.begin(), b.end(), [](auto const& u, auto const& v)
    { return same_storage(u, v); });
}

template <typename Ring>
number area2(Ring const& ring)
{
    number result = 0;
    for (std::size_t i = 0; i < ring.size(); ++i)
    {
        auto const& a = ring[i];
        auto const& b = ring[(i + 1) % ring.size()];
        result += lift(bg::get<0>(a)) * lift(bg::get<1>(b)) - lift(bg::get<1>(a)) * lift(bg::get<0>(b));
    }
    return result;
}

void test_point_and_segment()
{
    number const sagitta = number(1) / 8, displacement = number(1) / (1LL << 48);
    rational_polygon ideal;
    multi native;
    BOOST_TEST(bg::buffer_round_certified(p(0, 0), ideal, 1, sagitta, 0) == status::success);
    BOOST_TEST(bg::buffer_round_certified(p(0, 0), native, 1, sagitta, displacement) == status::success);
    BOOST_TEST_EQ(native.size(), 1u);
    if (native.size() != 1) { return; }
    auto const& ring = native.front().outer();
    BOOST_TEST(area2(ring) < 0);
    BOOST_TEST_EQ(ring.size(), ideal.outer().size() + 1);
    if (ring.size() != ideal.outer().size() + 1) { return; }
    BOOST_TEST_EQ(bg::get<0>(ring.front()), bg::get<0>(ring.back()));
    BOOST_TEST_EQ(bg::get<1>(ring.front()), bg::get<1>(ring.back()));
    for (std::size_t i = 0; i < ideal.outer().size(); ++i)
    {
        auto const& a = ideal.outer()[i];
        auto const& b = ideal.outer()[(i + 1) % ideal.outer().size()];
        number const x = bg::get<0>(a), y = bg::get<1>(a);
        number const sx = x + bg::get<0>(b), sy = y + bg::get<1>(b);
        number const dx = ec::exact_value(bg::get<0>(ring[i])) - x;
        number const dy = ec::exact_value(bg::get<1>(ring[i])) - y;
        BOOST_TEST(x*x + y*y == 1);
        BOOST_TEST(sx*sx + sy*sy >= 4 * (1 - sagitta) * (1 - sagitta));
        BOOST_TEST(dx*dx + dy*dy <= displacement * displacement);
    }
    auto const before = native.front();
    BOOST_TEST(bg::buffer_round_certified(p(0, 0), native, 1, sagitta, 0) == status::cannot_certify);
    BOOST_TEST_EQ(native.size(), 1u);
    if (native.size() == 1) { BOOST_TEST(same_storage(native.front(), before)); }

    // Each axis rounds by 1/3, but their Euclidean displacement exceeds 2/5.
    point<number> const fractional(number(1)/3, number(1)/3);
    bg::model::polygon<point<int>> integer;
    integer.outer().emplace_back(7, 8);
    auto const integer_before = integer;
    BOOST_TEST(bg::buffer_round_certified(fractional, integer, 1, 0.75, 0.4) == status::cannot_certify);
    BOOST_TEST(same_storage(integer, integer_before));
    BOOST_TEST(bg::buffer_round_certified(fractional, integer, 1, 0.75, 0.5) == status::success);
    BOOST_TEST_EQ(integer.outer().size(), 5u);

    BOOST_TEST(bg::buffer_round_certified(segment(p(0, 0), p(3, 4)), ideal, 1, sagitta, 0)
               == status::success);
    auto const capsule_area = area2(ideal.outer());
    BOOST_TEST(capsule_area < -20); // The 5-by-2 straight strip plus round ends.
    line const backtracking{p(0, 0), p(3, 4), p(1.5, 2), p(1.5, 2)};
    BOOST_TEST(bg::buffer_round_certified(backtracking, ideal, 1, sagitta, 0) == status::success);
    BOOST_TEST_EQ(area2(ideal.outer()), capsule_area);

    polygon out;
    BOOST_TEST(bg::buffer_round_certified(segment(p(0, 0), p(1, 1)), out, 1, sagitta, displacement)
               == status::success);
    BOOST_TEST(area2(out.outer()) < -5);
    BOOST_TEST(bg::buffer_round_certified(segment(p(0, 0), p(1, 1)), ideal, 1, sagitta, displacement)
               == status::cannot_certify);
    BOOST_TEST(bg::buffer_round_certified(segment(p(0, 0), p(0, 0)), ideal, 1, sagitta, 0)
               == status::success);
    BOOST_TEST(area2(ideal.outer()) > -8);
}

void test_star_and_nonstar()
{
    polygon out;
    line const bent{p(0, 0), p(1, 0), p(1, 1)};
    number const e = number(1) / 8, d = number(1) / (1LL << 42);
    BOOST_TEST(bg::buffer_round_certified(bent, out, 2, e, d) == status::success);
    polygon star;
    BOOST_TEST(bg::buffer_linestring_star_certified(bent, p(0.5, 0.5), star, 2, e, d) == status::success);
    BOOST_TEST_EQ(area2(out.outer()), area2(star.outer()));
    rational_polygon exact;
    BOOST_TEST(bg::buffer_round_certified(bent, exact, 2, e, 0) == status::success);

    // Opposite supports are sqrt(72) apart, so no common star witness exists.
    line const loop{p(0, 0), p(8, 8), p(0, 12), p(-8, 4), p(-0.5, 0.25)};
    BOOST_TEST(bg::buffer_round_certified(loop, out, 1, e, d) == status::success);
    BOOST_TEST_EQ(out.inners().size(), 1u);
    BOOST_TEST(area2(out.outer()) < 0);
    if (out.inners().size() == 1) { BOOST_TEST(area2(out.inners().front()) > 0); }
    auto const before = out;
    line const pinch{p(0, 0), p(-10, 0), p(-10, 10), p(12, 10), p(12, 0), p(2, 0)};
    BOOST_TEST(bg::buffer_round_certified(pinch, out, 1, 0.5, 0) == status::cannot_certify);
    BOOST_TEST(same_storage(out, before));

    // Integrated Euclidean export certificate and loose public budget.
    line const hairpin{p(0, 100), p(-0x1p52, 0), p(0, -100)};
    BOOST_TEST(bg::buffer_round_certified(hairpin, out, 1, e, 2) == status::success);
    BOOST_TEST(out.inners().empty());
    BOOST_TEST(area2(out.outer()) < 0);
}

void test_polygon()
{
    polygon box;
    box.outer() = {p(0, 0), p(0, 4), p(4, 4), p(4, 0), p(0, 0)};
    rational_polygon exact;
    BOOST_TEST(bg::buffer_round_certified(box, exact, 1, 0.125, 0) == status::success);
    BOOST_TEST(area2(exact.outer()) < -64);
    polygon out;
    BOOST_TEST(bg::buffer_round_certified(box, out, -1, 0, 0) == status::success);
    BOOST_TEST_EQ(area2(out.outer()), -8);
    auto const before = out;

    polygon concave;
    concave.outer() = {p(0, 0), p(0, 3), p(1, 3), p(1, 1), p(3, 1), p(3, 0), p(0, 0)};
    BOOST_TEST(bg::buffer_round_certified(concave, out, -10, 0, 0) == status::cannot_certify);

    polygon irrational;
    irrational.outer() = {p(0, 0), p(0, 4), p(4, 0), p(0, 0)};
    BOOST_TEST(bg::buffer_round_certified(irrational, out, -0.5, 0, 0) == status::cannot_certify);
    BOOST_TEST(same_storage(out, before));
    BOOST_TEST(bg::buffer_round_certified(box, out, 1, 0.125, 1) == status::success);
    BOOST_TEST(area2(out.outer()) < -64);
    BOOST_TEST(bg::buffer_round_certified(irrational, out, -0.5, 0, 1e-12) == status::success);
    // The exact eroded triangle has leg length 3 - sqrt(2)/2.
    auto const triangle_area = area2(out.outer());
    BOOST_TEST(triangle_area > number(-53)/10 && triangle_area < number(-52)/10);

    polygon translated;
    double const x = 0x1p52;
    translated.outer() = {p(x, 0), p(x, 4), p(x + 4, 4), p(x + 4, 0), p(x, 0)};
    auto const prior_to_collapse = out;
    // The exact erosion is 1/2 wide, but both sides round to x+2.
    BOOST_TEST(bg::buffer_round_certified(translated, out, -1.75, 0.125, 1) == status::cannot_certify);
    BOOST_TEST(same_storage(out, prior_to_collapse));

    BOOST_TEST(bg::buffer_round_certified(box, out, -2, 0, 0) == status::success);
    BOOST_TEST(out.outer().empty());
    box.outer()[2] = p(6, 4);
    box.outer()[3] = p(6, 0);
    multi collapsed{box};
    BOOST_TEST(bg::buffer_round_certified(box, collapsed, -2, 0, 0) == status::success);
    BOOST_TEST(collapsed.empty());
    BOOST_TEST(bg::buffer_round_certified(box, collapsed, -3, 0, 0) == status::success);
    BOOST_TEST(collapsed.empty());

    polygon invalid;
    invalid.outer() = {p(0, 0), p(4, 4), p(0, 4), p(4, 0), p(0, 0)};
    BOOST_TEST(bg::buffer_round_certified(invalid, out, -1, 0, 0) == status::invalid_input);
    BOOST_TEST(out.outer().empty());
}

struct xy
{
    number x, y;
};

xy sub(xy const& a, xy const& b) { return {a.x-b.x, a.y-b.y}; }
number cross(xy const& a, xy const& b) { return a.x*b.y-a.y*b.x; }
number dot(xy const& a, xy const& b) { return a.x*b.x+a.y*b.y; }

template <typename Polygon>
std::vector<xy> ccw_vertices(Polygon const& polygon)
{
    std::vector<xy> points;
    for (auto const& v : polygon.outer()) { points.push_back({lift(bg::get<0>(v)), lift(bg::get<1>(v))}); }
    if (bg::closure<Polygon>::value == bg::closed) { points.pop_back(); }
    if (bg::point_order<Polygon>::value == bg::clockwise) { std::reverse(points.begin(), points.end()); }
    return points;
}

bool clear_supports(std::vector<xy> const& shell, std::vector<xy> const& points, number const& radius)
{
    for (std::size_t j = 0; j < shell.size(); ++j)
    {
        auto const edge = sub(shell[(j+1) % shell.size()], shell[j]);
        number const length2 = dot(edge, edge), min2 = radius*radius*length2;
        if (length2 == 0) { return false; }
        for (auto const& v : points)
        {
            number const side = cross(edge, sub(v, shell[j]));
            if (side < 0 || side*side < min2) { return false; }
        }
    }
    return true;
}

bool offset_endpoint_matches(xy const& q, xy const& p, xy const& edge,
                             number const& r, number const& d)
{
    // Independent oracle in the edge's tangent/normal frame. The outward
    // projection, in unnormalized units, should be r*sqrt(length2).
    auto const delta = sub(q, p);
    number const length2 = dot(edge, edge);
    number const normal = -cross(edge, delta);
    number const tangent = dot(edge, delta);
    number const lhs = normal*normal + tangent*tangent + (r*r-d*d)*length2;
    if (lhs >= 0 && normal < 0) { return false; }
    if (lhs <= 0 && normal >= 0) { return true; }
    number const rhs2 = 4*r*r*normal*normal*length2;
    return lhs > 0 ? lhs*lhs <= rhs2 : lhs*lhs >= rhs2;
}

bool audit_convex_buffer(std::vector<xy> const& input, std::vector<xy> const& stored,
                         number const& r, number const& e, number const& d)
{
    // Identify every straight offset by its two adjacent stored endpoints.
    // A support merely parallel to the input, shifted inward to a disk
    // chord, cannot pass this test at small d.
    std::vector<std::size_t> start;
    auto const n = stored.size();
    for (std::size_t i = 0; i < input.size(); ++i)
    {
        auto const& a = input[i];
        auto const& b = input[(i+1) % input.size()];
        auto const edge = sub(b,a);
        std::size_t j = 0;
        for (; j < n; ++j)
        {
            if (offset_endpoint_matches(stored[j], a, edge, r, d)
                && offset_endpoint_matches(stored[(j+1) % n], b, edge, r, d)) { break; }
        }
        if (j == n) { return false; }
        start.push_back(j);
    }

    // Between identified offsets, every vertex belongs to that source
    // corner's circle up to export error. Account for every stored edge.
    std::size_t accounted = 0;
    number const lower = r-d, upper = r+d;
    for (std::size_t i = 0; i < input.size(); ++i)
    {
        std::size_t const from = (start[(i+input.size()-1) % input.size()]+1) % n;
        std::size_t const to = start[i];
        accounted += (to+n-from) % n + 1;
        for (std::size_t j = from;; j = (j+1) % n)
        {
            auto const v = sub(stored[j], input[i]);
            number const distance2 = dot(v,v);
            if (distance2 < lower*lower || distance2 > upper*upper) { return false; }
            if (j == to) { break; }
        }
    }
    return accounted == n && clear_supports(stored, stored, 0)
        && clear_supports(stored, input, r-e-d);
}

void test_support_exact_convex_dilation()
{
    using ccw_open = bg::model::polygon<p, false, false>;
    ccw_open source;
    // 184 distinct parabola secants with irrational lengths, and one
    // closing edge. No single small radical field represents this contour.
    for (int x = -92; x <= 92; ++x) { source.outer().emplace_back(x, x*x); }
    number const e = number(1)/16, d = number(1)/(1LL<<35);
    polygon out;
    if (!BOOST_TEST(bg::buffer_round_certified(source, out, 2, e, d) == status::success)) { return; }
    BOOST_TEST(out.inners().empty());
    BOOST_TEST(area2(out.outer()) < 0);
    BOOST_TEST_EQ(bg::get<0>(out.outer().front()), bg::get<0>(out.outer().back()));
    BOOST_TEST_EQ(bg::get<1>(out.outer().front()), bg::get<1>(out.outer().back()));
    BOOST_TEST(audit_convex_buffer(ccw_vertices(source), ccw_vertices(out), 2, e, d));
    auto const before = out;
    BOOST_TEST(bg::buffer_round_certified(source, out, 2, e, 0) == status::cannot_certify);
    BOOST_TEST(same_storage(before, out));
    BOOST_TEST(bg::buffer_round_certified(source, out, 2, 0, d) == status::cannot_certify);
    BOOST_TEST(same_storage(before, out));

    // Displacement <= 1 is loose enough for rounding on this grid, but
    // the result still cannot publish a collapsed or reversed edge.
    auto translated = source;
    for (auto& v : translated.outer())
    {
        bg::set<0>(v, bg::get<0>(v) + 0x1p52);
        bg::set<1>(v, bg::get<1>(v) + 0x1p52);
    }
    BOOST_TEST(bg::buffer_round_certified(translated, out, 2, e, 1) == status::cannot_certify);
    BOOST_TEST(same_storage(before, out));

    ccw_open subdivided;
    subdivided.outer() = {p(0,0), p(2,0), p(2,0), p(4,0), p(4,4), p(0,4)};
    ccw_open open;
    BOOST_TEST(bg::buffer_round_certified(subdivided, open, 1, 0.5, 0) == status::success);
    BOOST_TEST_EQ(open.outer().size(), 9u); // 8 arc ends + the retained subdivision.
    BOOST_TEST(std::count_if(open.outer().begin(), open.outer().end(), [](auto const& v)
    { return bg::get<0>(v) == 2 && bg::get<1>(v) == -1; }) == 1);
    std::reverse(subdivided.outer().begin(), subdivided.outer().end());
    BOOST_TEST(bg::buffer_round_certified(subdivided, out, 1, 0.5, 0) == status::invalid_input);
    BOOST_TEST(same_storage(before, out));
    bg::set<0>(subdivided.outer().front(), std::numeric_limits<double>::infinity());
    BOOST_TEST(bg::buffer_round_certified(subdivided, out, 1, 0.5, 0) == status::invalid_input);
    BOOST_TEST(same_storage(before, out));

    ccw_open star;
    star.outer() = {p(0,0), p(4,2), p(0,3), p(3,0), p(2,4)};
    auto const cycle = ccw_vertices(star);
    for (std::size_t i = 0; i < cycle.size(); ++i)
    {
        BOOST_TEST(cross(sub(cycle[(i+1)%5], cycle[i]),
                         sub(cycle[(i+2)%5], cycle[(i+1)%5])) > 0);
    }
    BOOST_TEST(bg::buffer_round_certified(star, out, 1, e, d) == status::invalid_input);
    BOOST_TEST(same_storage(before, out));
}

void test_support_exact_1281()
{
    using ccw_closed = bg::model::polygon<p, false>;
    ccw_closed source, intermediate, eroded;
    bg::read_wkt(issue_1281_wkt, source);
    BOOST_TEST_EQ(source.outer().size(), 104u);
    number const e = number(1)/128, d = number(1)/1000000000;
    if (!BOOST_TEST(bg::buffer_round_certified(source, intermediate, 9, e, d) == status::success)) { return; }
    if (!BOOST_TEST(bg::buffer_round_certified(intermediate, eroded, -9, 0, d) == status::success)) { return; }
    BOOST_TEST(intermediate.inners().empty());
    BOOST_TEST(eroded.inners().empty());
    BOOST_TEST(intermediate.outer().size() > 184u);
    BOOST_TEST(area2(eroded.outer()) > 0);
    auto const stored = ccw_vertices(intermediate);
    BOOST_TEST(audit_convex_buffer(ccw_vertices(source), stored, 9, e, d));
    BOOST_TEST(clear_supports(stored, ccw_vertices(eroded), 9-d));
}

void test_irrational_polygon_arrangement()
{
    using ccw_polygon = bg::model::polygon<p, false>;
    ccw_polygon holed;
    holed.outer() = {p(0, 0), p(24, 0), p(24, 24), p(15, 24), p(14, 20), p(10, 20),
                     p(9, 24), p(0, 24), p(0, 0)};
    holed.inners() = {{p(4, 6), p(5, 12), p(10, 7), p(4, 6)},
                     {p(16, 6), p(17, 12), p(21, 7), p(16, 6)}};
    number const e = number(1) / 5, d = number(1) / (1LL << 40);
    multi out;
    for (int distance : {1, -1})
    {
        BOOST_TEST(bg::buffer_round_certified(holed, out, distance, e, d) == status::success);
        BOOST_TEST_EQ(out.size(), 1u);
        if (out.size() != 1) { return; }
        auto const before = out.front();
        BOOST_TEST_EQ(before.inners().size(), 2u);
        auto const source_area = area2(holed.outer());
        number const shell_area = -area2(before.outer());
        BOOST_TEST(distance > 0 ? shell_area > source_area : shell_area < source_area);
        for (auto const& hole : before.inners()) { BOOST_TEST(area2(hole) > 0); }

        BOOST_TEST(bg::buffer_round_certified(holed, out, distance, e, 0) == status::cannot_certify);
        BOOST_TEST_EQ(out.size(), 1u);
        if (out.size() == 1) { BOOST_TEST(same_storage(out.front(), before)); }
    }

    rational_polygon exact;
    BOOST_TEST(bg::buffer_round_certified(holed, exact, 0, 0, 0) == status::success);
    BOOST_TEST_EQ(exact.outer().size() + 1, holed.outer().size());
    BOOST_TEST_EQ(area2(exact.outer()), -area2(holed.outer()));
    BOOST_TEST_EQ(exact.inners().size(), 2u);
    if (exact.inners().size() == 2)
    {
        BOOST_TEST_EQ(exact.inners()[0].size() + 1, holed.inners()[0].size());
        BOOST_TEST_EQ(area2(exact.inners()[0]), -area2(holed.inners()[0]));
    }
    ccw_polygon roundtrip;
    BOOST_TEST(bg::buffer_round_certified(exact, roundtrip, 0, 0, 0) == status::success);
    BOOST_TEST(same_storage(roundtrip, holed));
    auto const exact_before = exact;
    for (int distance : {1, -1})
    {
        BOOST_TEST(bg::buffer_round_certified(holed, exact, distance, e, d) == status::cannot_certify);
        BOOST_TEST(same_storage(exact, exact_before));
    }

    auto fractional = exact;
    bg::set<0>(fractional.outer().front(), number(-1) / 3);
    polygon native = out.front();
    auto const native_before = native;
    BOOST_TEST(bg::buffer_round_certified(fractional, native, 0, 0, 1) == status::cannot_certify);
    BOOST_TEST(same_storage(native, native_before));
    auto unclosed = holed;
    unclosed.outer().pop_back();
    BOOST_TEST(bg::buffer_round_certified(unclosed, native, 0, 0, 0) == status::invalid_input);
    BOOST_TEST(same_storage(native, native_before));

    ccw_polygon dumbbell;
    dumbbell.outer() = {p(0, 0), p(6, 0), p(6, 2), p(10, 2), p(10, 0), p(16, 0),
                        p(16, 6), p(10, 6), p(10, 4), p(6, 4), p(6, 6), p(0, 6), p(0, 0)};
    for (auto& p : dumbbell.outer()) { bg::set<1>(p, bg::get<1>(p) + bg::get<0>(p) / 4); }
    BOOST_TEST(bg::buffer_round_certified(dumbbell, out, -1.5, 0.25, d) == status::success);
    BOOST_TEST_EQ(out.size(), 2u);
    unsigned left = 0, right = 0;
    for (auto const& component : out)
    {
        BOOST_TEST(area2(component.outer()) < 0);
        BOOST_TEST(component.inners().empty());
        left += std::all_of(component.outer().begin(), component.outer().end(),
                            [](p const& p) { return bg::get<0>(p) < 8; });
        right += std::all_of(component.outer().begin(), component.outer().end(),
                             [](p const& p) { return bg::get<0>(p) > 8; });
    }
    BOOST_TEST_EQ(left, 1u);
    BOOST_TEST_EQ(right, 1u);
    BOOST_TEST(bg::buffer_round_certified(dumbbell, native, -1.5, 0.25, d) == status::cannot_certify);
    BOOST_TEST(same_storage(native, native_before));
}

void test_chord_erosion_empty_implication()
{
    polygon elbow;
    elbow.outer() = {p(0, 0), p(0, 3), p(1, 3), p(1, 1), p(3, 1), p(3, 0), p(0, 0)};
    // At r=5/8, neither unit-width arm can contain a surviving interior
    // point. In the unit joint, x,y>r gives distance to (1,1) < sqrt(2)*(1-r) < r.
    // The chord cap instead cuts off x+y>=2-r, leaving a triangle of leg 1/8.
    multi native;
    rational_polygon exact;
    BOOST_TEST(bg::buffer_round_certified(elbow, native, -0.625, 0.5, 0) == status::success);
    BOOST_TEST_EQ(native.size(), 1u);
    if (native.size() == 1) { BOOST_TEST_EQ(area2(native.front().outer()), number(-1) / 64); }
    BOOST_TEST(bg::buffer_round_certified(elbow, exact, -0.625, 0.5, 0) == status::success);
    BOOST_TEST_EQ(area2(exact.outer()), number(-1) / 64);

    // At r=3/4 there is no positive-area ideal chord face, so the
    // outer-bound certificate can report true regularized empty.
    BOOST_TEST(bg::buffer_round_certified(elbow, native, -0.75, 0.5, 0) == status::success);
    BOOST_TEST(native.empty());
    BOOST_TEST(bg::buffer_round_certified(elbow, exact, -0.75, 0.5, 0) == status::success);
    BOOST_TEST(exact.outer().empty());
}

void test_empty_vs_failure()
{
    polygon sentinel;
    sentinel.outer().emplace_back(7, 8);
    polygon out = sentinel;
    BOOST_TEST(bg::buffer_round_certified(p(0, 0), out, -1, 0, 0) == status::success);
    BOOST_TEST(out.outer().empty());
    out = sentinel;
    BOOST_TEST(bg::buffer_round_certified(segment(p(0, 0), p(1, 0)), out, -1, 0, 0) == status::success);
    BOOST_TEST(out.outer().empty());
    line const trace{p(0, 0), p(1, 0), p(1, 1)};
    multi all{sentinel, sentinel};
    BOOST_TEST(bg::buffer_round_certified(trace, all, -1, 0, 0) == status::success);
    BOOST_TEST(all.empty());
    all.push_back(sentinel);
    BOOST_TEST(bg::buffer_round_certified(trace, all, 0, 0, 0) == status::success);
    BOOST_TEST(all.empty());

    auto const nan = std::numeric_limits<double>::quiet_NaN();
    auto const inf = std::numeric_limits<double>::infinity();
    out = sentinel;
    BOOST_TEST(bg::buffer_round_certified(line{}, out, -1, 0, 0) == status::invalid_input);
    BOOST_TEST(bg::buffer_round_certified(p(nan, 0), out, -1, 0, 0) == status::invalid_input);
    BOOST_TEST(bg::buffer_round_certified(trace, out, inf, 0, 0) == status::invalid_input);
    BOOST_TEST(bg::buffer_round_certified(trace, out, 1, -1, 0) == status::invalid_input);
    BOOST_TEST(bg::buffer_round_certified(trace, out, 1, 1, -1) == status::invalid_input);
    BOOST_TEST(bg::buffer_round_certified(trace, out, 1, 0, 0) == status::cannot_certify);
    BOOST_TEST_EQ(out.outer().size(), 1u);
    if (!out.outer().empty()) { BOOST_TEST_EQ(bg::get<0>(out.outer().front()), 7); }
}

void test_public_multipoint()
{
    mp const overlap{p(0, 0), p(1.5, 0)};
    polygon single;
    // A loose public sagitta is tightened for the arrangement. Quarter
    // chords give two unit diamonds whose overlap has area 1/8.
    BOOST_TEST(bg::buffer_round_certified(overlap, single, 1, 1, 0) == status::success);
    BOOST_TEST_EQ(area2(single.outer()), number(-31) / 4);
    auto const before = single;
    mp const disjoint{p(0, 0), p(20, 0)};
    multi out;
    BOOST_TEST(bg::buffer_round_certified(disjoint, out, 1, 0.5, 0) == status::success);
    BOOST_TEST_EQ(out.size(), 2u);
    BOOST_TEST(bg::buffer_round_certified(disjoint, single, 1, 0.5, 0) == status::cannot_certify);
    BOOST_TEST(same_storage(single, before));
    auto const prior = out;
    BOOST_TEST(bg::buffer_round_certified(mp{p(0, 0), p(2, 0)}, out, 1, 0.5, 0) == status::cannot_certify);
    BOOST_TEST(same_multi(out, prior));
    BOOST_TEST(bg::buffer_round_certified(mp{}, out, 1, 0, 0) == status::success);
    BOOST_TEST(out.empty());
    BOOST_TEST(bg::buffer_round_certified(overlap, single, 0, 0, 0) == status::success);
    BOOST_TEST(single.outer().empty());
}

void test_public_multiline()
{
    number const e = number(1) / 8, d = number(1) / (1LL << 42);
    ml const loop{line{p(0, 0), p(8, 8), p(0, 12)}, line{p(0, 12), p(-8, 4), p(-0.5, 0.25)}, line{}};
    polygon single;
    BOOST_TEST(bg::buffer_round_certified(loop, single, 1, e, d) == status::success);
    BOOST_TEST_EQ(single.inners().size(), 1u);
    auto const before = single;
    ml const disjoint{line{p(0, 0), p(2, 0)}, line{p(20, 0), p(22, 0)}};
    multi out;
    BOOST_TEST(bg::buffer_round_certified(disjoint, out, 1, 0.5, 0) == status::success);
    BOOST_TEST_EQ(out.size(), 2u);
    BOOST_TEST(bg::buffer_round_certified(disjoint, single, 1, 0.5, 0) == status::cannot_certify);
    BOOST_TEST(same_storage(single, before));

    // This fraction's exposed endpoint cannot export losslessly on its
    // own. The closed member's exact chord band hides the entire capsule.
    // A union of member exports could not succeed at zero displacement.
    using rp = point<number>;
    using rl = bg::model::linestring<rp>;
    rl const hidden{rp(number(13) / 3, 4), rp(number(14) / 3, 4)};
    BOOST_TEST(bg::buffer_round_certified(hidden, single, 5, 4, 0) == status::cannot_certify);
    BOOST_TEST(same_storage(single, before));
    bg::model::multi_linestring<rl> const combination{
        rl{rp(0, 0), rp(8, 0), rp(8, 8), rp(0, 8), rp(0, 0)}, hidden};
    BOOST_TEST(bg::buffer_round_certified(combination, single, 5, 4, 0) == status::success);
    BOOST_TEST(single.inners().empty());
    BOOST_TEST_EQ(area2(single.outer()), -548);

    auto const original = out;
    double const inf = (std::numeric_limits<double>::infinity)();
    ml const nonfinite{line{p(0, 0), p(inf, 0)}};
    BOOST_TEST(bg::buffer_round_certified(nonfinite, out, -1, 0, 0) == status::invalid_input);
    BOOST_TEST(same_multi(out, original));
    BOOST_TEST(bg::buffer_round_certified(disjoint, out, -1, 0, 0) == status::success);
    BOOST_TEST(out.empty());
    BOOST_TEST(bg::buffer_round_certified(ml{line{}}, single, 1, 0, 0) == status::success);
    BOOST_TEST(single.outer().empty());
}

polygon box(double x, double y, double u, double v)
{
    polygon result;
    result.outer() = {p(x, y), p(x, v), p(u, v), p(u, y), p(x, y)};
    return result;
}

number area2(multi const& m)
{
    number area = 0;
    for (auto const& p : m)
    {
        area += area2(p.outer());
        for (auto const& h : p.inners()) { area += area2(h); }
    }
    return area;
}

void test_public_multipolygon()
{
    auto annulus = box(0, 0, 20, 20);
    annulus.inners().push_back(box(5, 5, 15, 15).outer());
    std::reverse(annulus.inners().front().begin(), annulus.inners().front().end());
    // The island is inside the hole, hence not in the other member's fill.
    multi const input{annulus, box(8, 8, 12, 12)};
    multi out;
    BOOST_TEST(bg::buffer_round_certified(input, out, 1, 0.5, 0) == status::success);
    BOOST_TEST_EQ(out.size(), 2u);
    BOOST_TEST(bg::buffer_round_certified(input, out, -1, 0.5, 0) == status::success);
    BOOST_TEST_EQ(out.size(), 2u);
    // Square erosion [1,19]^2, minus the [5,15]^2 hole grown by the
    // unit diamond (area 142), plus the eroded island [9,11]^2.
    BOOST_TEST_EQ(area2(out), -372);
    std::size_t holes = 0;
    for (auto const& p : out) { holes += p.inners().size(); }
    BOOST_TEST_EQ(holes, 1u);
    auto const before = out;

    // Sharing the seam of two rectangles does not make that seam an
    // erosion boundary of their union. This source is refused, not eroded
    // componentwise. Contained overlap is refused even with no ring crossing.
    multi const shared{box(0, 0, 2, 2), box(2, 0, 4, 2)};
    for (double r : {-0.5, 0.0, 1.0})
    {
        BOOST_TEST(bg::buffer_round_certified(shared, out, r, 0.25, 0) == status::invalid_input);
        BOOST_TEST(same_multi(out, before));
    }
    multi const overlap{annulus, box(2, 2, 3, 3)};
    BOOST_TEST(bg::buffer_round_certified(overlap, out, -1, 0.5, 0) == status::invalid_input);
    BOOST_TEST(same_multi(out, before));

    using exact_multi = bg::model::multi_polygon<rational_polygon>;
    exact_multi exact;
    BOOST_TEST(bg::buffer_round_certified(input, exact, 0, 0, 0) == status::success);
    BOOST_TEST_EQ(exact.size(), 2u);
    if (exact.size() != 2) { return; }
    auto const before_exact = exact;
    BOOST_TEST(bg::buffer_round_certified(input, exact, 1, 0.5, 0) == status::cannot_certify);
    BOOST_TEST(same_multi(exact, before_exact));
    BOOST_TEST_EQ(exact.front().outer().size() + 1, input.front().outer().size());
    BOOST_TEST(bg::buffer_round_certified(exact, out, 0, 0, 0) == status::success);
    BOOST_TEST(same_multi(out, input));
    for (auto& p : exact.back().outer()) { bg::set<0>(p, number(bg::get<0>(p) + number(1) / 3)); }
    BOOST_TEST(bg::buffer_round_certified(exact, out, 0, 0, 1) == status::cannot_certify);
    BOOST_TEST(same_multi(out, input));

    polygon single = input.front();
    auto const old_single = single;
    BOOST_TEST(bg::buffer_round_certified(input, single, 0, 0, 0) == status::cannot_certify);
    BOOST_TEST(same_storage(single, old_single));
    BOOST_TEST(bg::buffer_round_certified(exact.back(), single, 6, 5, 0) == status::cannot_certify);
    BOOST_TEST(same_storage(single, old_single));
    // Dilation fills the annulus hole and hides the island's non-dyadic
    // coordinates. Only the final exact union has to export.
    BOOST_TEST(bg::buffer_round_certified(exact, single, 6, 5, 0) == status::success);
    BOOST_TEST(single.inners().empty());
    BOOST_TEST_EQ(area2(single.outer()), -1904);
    BOOST_TEST(bg::buffer_round_certified(multi{polygon{}}, out, 1, 0, 0) == status::success);
    BOOST_TEST(out.empty());
}

int main()
{
    test_point_and_segment();
    test_star_and_nonstar();
    test_polygon();
    test_support_exact_convex_dilation();
    test_support_exact_1281();
    test_irrational_polygon_arrangement();
    test_chord_erosion_empty_implication();
    test_empty_vs_failure();
    test_public_multipoint();
    test_public_multiline();
    test_public_multipolygon();
    return boost::report_errors();
}
