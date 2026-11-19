// Boost.Geometry
// Unit Test
// Copyright (c) 2026
// Distributed under the Boost Software License, Version 1.0.

#include "issue_578_fixture.hpp"

#include <boost/core/lightweight_test.hpp>
#include <boost/geometry/algorithms/buffer_linestring_certified.hpp>
#include <boost/geometry/algorithms/buffer_linestring_star_certified.hpp>
#include <boost/geometry/geometries/linestring.hpp>
#include <boost/geometry/io/wkt/read.hpp>
#include <boost/geometry/io/wkt/write.hpp>
#include <boost/geometry/strategies/cartesian/exact/radical_buffer_noding.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <type_traits>
#include <vector>

namespace bg = boost::geometry;
namespace ec = bg::strategies::detail::exact_cartesian;
namespace rb = ec::round_buffer;
namespace rn = rb::radical_noding;
namespace seg = bg::detail::certified_segment_buffer;
using number = ec::exact_number;
using xy = rb::rational_xy;
using native_point = bg::model::point<double, 2, bg::cs::cartesian>;
using linestring = bg::model::linestring<native_point>;

std::vector<xy> exact_trace(linestring const& line)
{
    std::vector<xy> out;
    for (auto const& p : line) { out.push_back(seg::read_point(p)); }
    return out;
}

std::vector<seg::rational_segment> edges(std::vector<xy> const& trace)
{
    std::vector<seg::rational_segment> out;
    for (std::size_t i = 1; i < trace.size(); ++i) { out.emplace_back(trace[i - 1], trace[i]); }
    return out;
}

rn::value dot(rn::field const& f, rn::point const& a, rn::point const& b)
{
    return f.add(f.multiply(a.x, b.x), f.multiply(a.y, b.y));
}

bool below_distance(rn::field const& f, rn::point const& p,
                    seg::rational_segment const& s, number const& radius2)
{
    auto const v = rn::subtract(f, p, rn::lift(f, s.from));
    auto const d = rn::lift(f, s.direction);
    auto const projection = dot(f, v, d);
    if (f.sign(projection) <= 0)
    {
        return f.sign(f.subtract(dot(f, v, v), f.rational(radius2))) < 0;
    }
    if (f.sign(f.subtract(projection, f.rational(s.norm2))) >= 0)
    {
        auto const w = rn::subtract(f, p, rn::lift(f, s.to));
        return f.sign(f.subtract(dot(f, w, w), f.rational(radius2))) < 0;
    }
    auto const area = rn::cross(f, d, v);
    return f.sign(f.subtract(f.multiply(area, area), f.rational(radius2 * s.norm2))) < 0;
}

struct cover_certificate
{
    std::vector<seg::rational_segment> const& outer;
    number const inner_radius, outer_radius2;
    unsigned arcs = 0;
    unsigned static constexpr max_arcs = 4096;
    unsigned static constexpr max_depth = 12;

    bool cover_arc(xy const& base, unsigned quadrant, number const& a, number const& b, unsigned depth)
    {
        if (++arcs > max_arcs) { return false; }
        auto const u = rb::rational_arc_sample(quadrant, a), v = rb::rational_arc_sample(quadrant, b);
        number const scale = inner_radius / (1 + seg::dot(u, v));
        // The arc lies in the triangle between its chord and the two
        // tangents. The tangents meet at r*(u+v)/(1+u.v), all rational.
        std::array<xy, 3> const triangle{{{base.x + inner_radius * u.x, base.y + inner_radius * u.y},
            {base.x + scale * (u.x + v.x), base.y + scale * (u.y + v.y)},
            {base.x + inner_radius * v.x, base.y + inner_radius * v.y}}};
        for (auto const& s : outer)
        {
            auto const inside = [&](xy const& p) { return s.squared_distance(p) < outer_radius2; };
            if (std::all_of(triangle.begin(), triangle.end(), inside)) { return true; }
        }
        if (depth == max_depth) { return false; }
        number const mid = (a + b) / 2;
        return cover_arc(base, quadrant, a, mid, depth + 1) && cover_arc(base, quadrant, mid, b, depth + 1);
    }

    bool cover_side(seg::rational_segment const& edge, rb::offset_side side) const
    {
        rn::field f;
        rn::segment precise;
        if (!rn::lift_offset(f, rb::straight_offset(edge.from, edge.to, inner_radius, side), precise))
        {
            return false;
        }
        for (auto const& s : outer)
        {
            if (below_distance(f, precise.from, s, outer_radius2) && below_distance(f, precise.to, s, outer_radius2))
            {
                return true; // The entire true algebraic side is in this convex capsule.
            }
        }
        return false;
    }
};

void check_exposed_irrational_side(std::vector<xy> const& source, std::size_t i, number const& radius)
{
    // Every source point is on or to the right of the producer. Its true
    // left offset reaches the global extremum in that normal direction.
    seg::rational_segment const edge(source[i], source[i + 1]);
    for (auto const& p : source) { BOOST_TEST(edge.side(p) <= 0); }
    rn::field f;
    rn::value root;
    BOOST_TEST(f.try_sqrt(edge.norm2, root));
    BOOST_TEST_EQ(f.rank(), 1u);
    number rational = 0;
    // cross(d,p-a)=r*sqrt(d.d) on the entire exposed side. No rational
    // point satisfies that equality, so the ideal union cannot be exported
    // losslessly as any finite polygon with native or rational coordinates.
    BOOST_TEST_NOT(f.try_rational_quotient(f.scale(root, radius), f.rational(1), rational));
}

void check_original_mixed_chord_node(std::vector<seg::rational_segment> const& edges, number const& error)
{
    // A proper crossing of two permitted radius-4 cap chords from the first
    // two actual outer source edges. These are not historical join samples:
    // the charts and the sagitta budget define the support-exact construction.
    rn::field f;
    auto const cap_chord = [&](seg::rational_segment const& edge, unsigned quadrant,
                               number const& from, number const& to)
    {
        rn::value root;
        BOOST_TEST(f.try_sqrt(edge.norm2, root));
        auto const u = rb::rational_arc_sample(quadrant, from);
        auto const v = rb::rational_arc_sample(quadrant, to);
        BOOST_TEST(rb::minor_arc_radial_error_within(rb::circle_direction(u), rb::circle_direction(v), 4, error));
        auto const& center = quadrant < 2 ? edge.to : edge.from;
        auto const& d = edge.direction;
        auto const at = [&](xy const& sample)
        {
            return rn::point{
                f.add(f.rational(center.x), f.scale(root, 4 * (sample.x * d.y + sample.y * d.x) / edge.norm2)),
                f.add(f.rational(center.y), f.scale(root, 4 * (-sample.x * d.x + sample.y * d.y) / edge.norm2))};
        };
        return rn::segment{at(u), at(v)};
    };
    auto const a = cap_chord(edges[0], 3, number(63) / 64, 1);
    auto const b = cap_chord(edges[1], 0, 0, number(1) / 64);
    BOOST_TEST_EQ(f.rank(), 2u);
    auto const i = rn::intersect(f, a, b);
    BOOST_TEST(i.kind == rn::intersection_kind::contact);
    if (i.kind != rn::intersection_kind::contact) { return; }
    BOOST_TEST(f.sign(i.t0.n) > 0 && f.sign(f.subtract(i.t0.d, i.t0.n)) > 0);
    BOOST_TEST(f.sign(i.u0.n) > 0 && f.sign(f.subtract(i.u0.d, i.u0.n)) > 0);
    auto const p = rn::at(f, a, i.t0), q = rn::at(f, b, i.u0);
    BOOST_TEST_EQ(rn::compare(f, p.x, q.x), 0);
    BOOST_TEST_EQ(rn::compare(f, p.y, q.y), 0);
    number rational;
    BOOST_TEST_NOT(f.try_rational_quotient(p.x.n, p.x.d, rational));
}

using native_polygon = bg::model::polygon<native_point, false, false>;

// Independent exact rational sweep on what was stored. Check both sets of
// angular knots; between them the radial order is the sign of a linear
// functional of the ray, so positive order at both ends excludes a crossing.
template <typename A, typename B>
bool check_stored_radial_order(A const& inner, B const& outer, xy const& witness, bool equal = false)
{
    auto const angle_less = [](xy const& a, xy const& b)
    {
        bool const ha = seg::common::second_half(a), hb = seg::common::second_half(b);
        return ha != hb ? !ha : seg::cross(a, b) > 0;
    };
    auto const directions = [&](auto const& polygon)
    {
        using polygon_type = typename std::decay<decltype(polygon)>::type;
        std::vector<xy> points;
        for (auto const& p : polygon.outer()) { points.push_back(seg::difference(seg::read_point(p), witness)); }
        if constexpr (bg::closure<polygon_type>::value == bg::closed)
        { if (!points.empty()) { points.pop_back(); } }
        if constexpr (bg::point_order<polygon_type>::value == bg::clockwise)
        { std::reverse(points.begin(), points.end()); }
        std::rotate(points.begin(), std::min_element(points.begin(), points.end(), angle_less), points.end());
        return points;
    };
    auto const a = directions(inner), b = directions(outer);
    if (a.size() < 3 || b.size() < 3) { return false; }
    auto const check = [&](std::vector<xy> const& probes, std::vector<xy> const& ring, bool inside)
    {
        std::size_t i = 0;
        for (auto const& p : probes)
        {
            while (i < ring.size() && !angle_less(p, ring[i])) { ++i; }
            auto const& from = i == 0 ? ring.back() : ring[i - 1];
            auto const& to = i == ring.size() ? ring.front() : ring[i];
            auto const side = seg::cross(seg::difference(to, from), seg::difference(p, from));
            if (equal ? side != 0 : inside ? side <= 0 : side >= 0) { return false; }
        }
        return true;
    };
    return check(a, b, true) && check(b, a, false);
}

void test_original(bool emit_wkt)
{
    linestring input_inner, input_outer;
    bg::read_wkt(issue_578_fixture::inner, input_inner);
    bg::read_wkt(issue_578_fixture::outer, input_outer);
    auto const inner = exact_trace(input_inner), outer = exact_trace(input_outer);
    BOOST_TEST_EQ(inner.size(), 16u);
    BOOST_TEST_EQ(outer.size(), 101u);
    if (inner.size() != 16 || outer.size() != 101) { return; }
    auto const outer_edges = edges(outer);
    auto const inner_edges = edges(inner);

    // Re-evaluate the original binary64 formula, including its 1e-6 slack.
    // The exact certificate uses rational directions and conservative radii,
    // not floating evaluation of distance or the historical chord schedule.
    double constexpr pi = 0x1.921fb54442d18p+1;
    double const nominal_sagitta = (1. - std::cos(pi / 180.)) * 4.;
    double const original_inner_radius = 4. - nominal_sagitta - 1e-6;
    number const error = number(60925) / 100000000;
    number const inner_bound = number(19996949) / 5000000;
    number const reduced_outer = 4 - error;
    BOOST_TEST(ec::exact_value(nominal_sagitta) < error);
    BOOST_TEST(ec::exact_value(original_inner_radius) < inner_bound);
    // Also bound the mathematical formula independently of libm, using
    // 3.14159 < pi < 3.14160 and alternating cosine bounds at pi/180.
    number const lo = number(314159) / 18000000, hi = number(314160) / 18000000;
    BOOST_TEST(2 * hi * hi < error);
    BOOST_TEST(4 - 2 * lo * lo + hi * hi * hi * hi / 6 - number(1) / 1000000 < inner_bound);

    // The stated subtrace premise does not hold even approximately at the
    // 1e-6 radial slack. Vertex 5 has clearance >0.39 from every outer edge.
    number const separation2 = number(1521) / 10000;
    for (auto const& s : outer_edges) { BOOST_TEST(s.squared_distance(inner[5]) > separation2); }

    cover_certificate cert{outer_edges, inner_bound, reduced_outer * reduced_outer};
    // A point common to all convex outer capsules makes their union
    // star-shaped. Covering each inner capsule's boundary is then enough:
    // extend the ray from this common point through any interior point
    // to an inner boundary point, then use the covering outer capsule.
    for (auto const& s : outer_edges) { BOOST_TEST(s.squared_distance({85, -62}) < cert.outer_radius2); }
    for (auto const& p : inner)
    {
        for (unsigned q = 0; q < 4; ++q) { BOOST_TEST(cert.cover_arc(p, q, 0, 1, 0)); }
    }
    for (auto const& s : inner_edges)
    {
        BOOST_TEST(cert.cover_side(s, rb::offset_side::left));
        BOOST_TEST(cert.cover_side(s, rb::offset_side::right));
    }

    check_exposed_irrational_side(inner, 0, ec::exact_value(original_inner_radius));
    check_exposed_irrational_side(outer, 29, 4);
    check_original_mixed_chord_node(outer_edges, error);
    bg::detail::certified_linestring_buffer::exact_polygon output;
    output.outer().emplace_back(1, 2);
    using status = bg::certified_buffer_status;
    BOOST_TEST(bg::buffer_linestring_certified(input_inner, output, original_inner_radius, error) == status::cannot_certify);
    BOOST_TEST(bg::buffer_linestring_certified(input_outer, output, 4, error) == status::cannot_certify);
    BOOST_TEST_EQ(output.outer().size(), 1u);
    BOOST_TEST_EQ(bg::get<0>(output.outer().front()), 1);

    // Native export uses a separate contract. Each stored vertex is at
    // Euclidean distance <= displacement from the ideal contour vertex.
    number const displacement = number(1) / (boost::multiprecision::cpp_int(1) << 40);
    BOOST_TEST(ec::exact_value(original_inner_radius) + 2 * displacement < inner_bound);
    native_point const witness(85, -62);
    native_polygon stored_inner, stored_outer;
    BOOST_TEST(bg::buffer_linestring_star_certified(input_inner, witness, stored_inner,
              original_inner_radius, error, displacement) == status::success);
    BOOST_TEST(bg::buffer_linestring_star_certified(input_outer, witness, stored_outer,
              4, error, displacement) == status::success);
    bool const nested = check_stored_radial_order(stored_inner, stored_outer, {85, -62});
    BOOST_TEST(nested);
    if (emit_wkt && nested)
    {
        std::cout << std::setprecision(std::numeric_limits<double>::max_digits10)
                  << bg::wkt(stored_inner) << '\n' << bg::wkt(stored_outer) << '\n';
    }

    auto const count = stored_inner.outer().size();
    BOOST_TEST(bg::buffer_linestring_star_certified(input_inner, witness, stored_inner,
              original_inner_radius, error, 0) == status::cannot_certify);
    BOOST_TEST(bg::buffer_linestring_star_certified(input_inner, native_point(0, 0), stored_inner,
              original_inner_radius, error, displacement) == status::cannot_certify);
    BOOST_TEST_EQ(stored_inner.outer().size(), count);

    using float_point = bg::model::point<float, 2, bg::cs::cartesian>;
    bg::model::polygon<float_point> coarse;
    coarse.outer().emplace_back(9, 4);
    BOOST_TEST(bg::buffer_linestring_star_certified(input_inner, witness, coarse,
              original_inner_radius, error, displacement) == status::cannot_certify);
    BOOST_TEST_EQ(coarse.outer().size(), 1u);
}

void test_against_rational_union()
{
    // The four axial directions suffice at this sagitta. This L-shaped
    // example has coincident and crossing supports but an integral union.
    linestring const trace{native_point(0, 0), native_point(2, 0), native_point(2, 2)};
    bg::model::polygon<native_point> stored; // clockwise and closed
    bg::detail::certified_linestring_buffer::exact_polygon reference;
    using status = bg::certified_buffer_status;
    BOOST_TEST(bg::buffer_linestring_certified(trace, reference, 5, 4) == status::success);
    BOOST_TEST(bg::buffer_linestring_star_certified(trace, native_point(2, 0), stored, 5, 4, 0) == status::success);
    BOOST_TEST(check_stored_radial_order(stored, reference, {2, 0}, true));
}

void test_collapsed_export()
{
    // All x coordinates collapse to x on this binary32 grid. A generous
    // displacement allowance does not bypass the stored-ring topology test.
    double constexpr x = 16777216;
    linestring const trace{native_point(x, 0), native_point(x + 0.125, 0)};
    using float_point = bg::model::point<float, 2, bg::cs::cartesian>;
    bg::model::polygon<float_point> stored;
    stored.outer().emplace_back(9, 4);
    BOOST_TEST(bg::buffer_linestring_star_certified(trace, native_point(x, 0), stored, 0.25, 0.125, 1)
               == bg::certified_buffer_status::cannot_certify);
    BOOST_TEST_EQ(stored.outer().size(), 1u);
    BOOST_TEST_EQ(bg::get<0>(stored.outer().front()), 9);
}

int main(int argc, char** argv)
{
    test_original(argc == 2 && std::string(argv[1]) == "--emit-wkt");
    test_against_rational_union();
    test_collapsed_export();
    return boost::report_errors();
}
