// Boost.Geometry
// Unit Test
// Copyright (c) 2026
// Distributed under the Boost Software License, Version 1.0.

#include <boost/geometry/algorithms/buffer_polygon_arrangement_certified.hpp>
#include <boost/geometry/geometries/multi_polygon.hpp>
#include <boost/core/lightweight_test.hpp>

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <initializer_list>
#include <iostream>
#include <string>
#include <vector>

namespace bg = boost::geometry;
namespace ec = bg::strategies::detail::exact_cartesian;
using number = ec::exact_number;
using status = bg::certified_buffer_status;
template <typename T> using point = bg::model::point<T, 2, bg::cs::cartesian>;
using input_polygon = bg::model::polygon<point<double>, false>;
using output_polygon = bg::model::polygon<point<double>>;
using multi = bg::model::multi_polygon<output_polygon>;

struct xy { number x, y; };
struct edge { xy a, b; };

number q(int a, int b) { return number(a) / b; }
number const displacement = number(1) / (boost::multiprecision::cpp_int(1) << 40);
xy operator-(xy const& a, xy const& b) { return {a.x - b.x, a.y - b.y}; }
number det(xy const& a, xy const& b) { return a.x * b.y - a.y * b.x; }
number squared(xy const& a) { return a.x * a.x + a.y * a.y; }

template <typename T>
xy read(point<T> const& p) { return {ec::exact_value(bg::get<0>(p)), ec::exact_value(bg::get<1>(p))}; }

template <typename Ring>
void close(Ring& ring, std::initializer_list<std::pair<double, double>> points, bool shear = false)
{
    for (auto const& p : points) { ring.emplace_back(p.first, p.second + (shear ? p.first / 4 : 0)); }
    ring.push_back(ring.front());
}

input_polygon mixed_holes()
{
    input_polygon p;
    close(p.outer(), {{0, 0}, {24, 0}, {24, 24}, {15, 24}, {14, 20}, {10, 20}, {9, 24}, {0, 24}});
    p.inners().resize(2);
    close(p.inners()[0], {{4, 6}, {5, 12}, {10, 7}});
    close(p.inners()[1], {{16, 6}, {17, 12}, {21, 7}});
    return p;
}

input_polygon dumbbell()
{
    input_polygon p;
    close(p.outer(), {{0, 0}, {6, 0}, {6, 2}, {10, 2}, {10, 0}, {16, 0},
                      {16, 6}, {10, 6}, {10, 4}, {6, 4}, {6, 6}, {0, 6}}, true);
    return p;
}

template <typename Ring>
void add_edges(std::vector<edge>& out, Ring const& ring)
{
    for (std::size_t i = 1; i < ring.size(); ++i) { out.push_back({read(ring[i - 1]), read(ring[i])}); }
}

template <typename Polygon>
void add_boundary(std::vector<edge>& out, Polygon const& p)
{
    add_edges(out, p.outer());
    for (auto const& h : p.inners()) { add_edges(out, h); }
}

std::vector<edge> boundary(multi const& m)
{
    std::vector<edge> out;
    for (auto const& p : m) { add_boundary(out, p); }
    return out;
}

int locate(std::vector<edge> const& ring, xy const& p)
{
    // Independent rational ray-crossing oracle using the intersection x,
    // not the production algebraic side predicate or a floating strategy.
    bool inside = false;
    for (auto const& e : ring)
    {
        auto const& a = e.a;
        auto const& b = e.b;
        if ((p.x - a.x) * (p.x - b.x) <= 0 && (p.y - a.y) * (p.y - b.y) <= 0
            && det(p - a, b - a) == 0) { return 0; }
        if ((a.y > p.y) != (b.y > p.y) && a.x + (b.x - a.x) * (p.y - a.y) / (b.y - a.y) > p.x)
        { inside = !inside; }
    }
    return inside ? 1 : -1;
}

number distance2(std::vector<edge> const& source, xy const& p)
{
    number best = squared(p - source.front().a);
    for (auto const& e : source)
    {
        xy const v = e.b - e.a, w = p - e.a;
        number const t = (v.x * w.x + v.y * w.y) / squared(v);
        number const s = (std::max)(number(0), (std::min)(t, number(1)));
        number const candidate = squared({w.x - s * v.x, w.y - s * v.y});
        if (candidate < best) { best = candidate; }
    }
    return best;
}

template <typename Ring>
void check_orientation(Ring const& ring, bool ccw)
{
    BOOST_TEST(ring.size() >= 4);
    if (ring.empty()) { return; }
    BOOST_TEST_EQ(bg::get<0>(ring.front()), bg::get<0>(ring.back()));
    BOOST_TEST_EQ(bg::get<1>(ring.front()), bg::get<1>(ring.back()));
    number area = 0;
    for (std::size_t i = 1; i < ring.size(); ++i) { area += det(read(ring[i - 1]), read(ring[i])); }
    BOOST_TEST(ccw ? area > 0 : area < 0);
}

void test_holed_both_signs()
{
    auto const input = mixed_holes();
    // These are source squared lengths. There are five independent square
    // classes globally, so this fixture requires producer-local fields.
    ec::round_buffer::bounded_radical_field field;
    ec::round_buffer::bounded_radical_field::value root;
    for (int d : {17, 37, 50, 41}) { BOOST_TEST(field.try_sqrt(d, root)); }
    BOOST_TEST_EQ(field.rank(), 4u);
    BOOST_TEST_NOT(field.try_sqrt(26, root));

    multi expanded, eroded;
    number const e = q(1, 5);
    BOOST_TEST(bg::buffer_polygon_arrangement_certified(input, expanded, 1, e, displacement) == status::success);
    BOOST_TEST(bg::buffer_polygon_arrangement_certified(input, eroded, -1, e, displacement) == status::success);
    for (auto const* m : {&expanded, &eroded})
    {
        BOOST_TEST_EQ(m->size(), 1u);
        if (m->size() != 1) { return; }
        BOOST_TEST_EQ(m->front().inners().size(), 2u);
        check_orientation(m->front().outer(), false);
        for (auto const& h : m->front().inners()) { check_orientation(h, true); }
    }

    std::vector<edge> source;
    add_boundary(source, input);
    auto const plus = boundary(expanded), minus = boundary(eroded);
    BOOST_TEST_EQ(locate(plus, {12, q(41, 2)}), 1); // dilated concave notch
    BOOST_TEST_EQ(locate(minus, {12, q(41, 2)}), -1);
    for (auto const& hole : {xy{6, 8}, xy{18, 8}})
    {
        BOOST_TEST_EQ(locate(plus, hole), -1);
        BOOST_TEST_EQ(locate(minus, hole), -1);
    }

    number const close_radius = 1 - e - displacement, far_radius = 1 + displacement;
    number const close2 = close_radius * close_radius, far2 = far_radius * far_radius;
    unsigned in_plus = 0, out_plus = 0, in_minus = 0, out_minus = 0;
    for (int y = -2; y <= 26; ++y)
    {
        for (int x = -2; x <= 26; ++x)
        {
            xy const p{q(4 * x + 1, 4), q(2 * y + 1, 2)};
            int const material = locate(source, p);
            number const clearance2 = distance2(source, p);
            if (material >= 0 || clearance2 < close2) { BOOST_TEST_EQ(locate(plus, p), 1); ++in_plus; }
            if (material < 0 && clearance2 > far2) { BOOST_TEST_EQ(locate(plus, p), -1); ++out_plus; }
            if (material > 0 && clearance2 > far2) { BOOST_TEST_EQ(locate(minus, p), 1); ++in_minus; }
            if (material <= 0 || clearance2 < close2) { BOOST_TEST_EQ(locate(minus, p), -1); ++out_minus; }
        }
    }
    BOOST_TEST(in_plus > 500 && out_plus > 100 && in_minus > 300 && out_minus > 200);

    auto const size = expanded.front().outer().size();
    BOOST_TEST(bg::buffer_polygon_arrangement_certified(input, expanded, 1, e, 0) == status::cannot_certify);
    BOOST_TEST_EQ(expanded.front().outer().size(), size);
    bg::model::multi_polygon<bg::model::polygon<point<float>>> coarse;
    coarse.resize(1);
    coarse.front().outer().emplace_back(3, 5);
    BOOST_TEST(bg::buffer_polygon_arrangement_certified(input, coarse, -1, e, displacement) == status::cannot_certify);
    BOOST_TEST_EQ(coarse.size(), 1u);
    BOOST_TEST_EQ(coarse.front().outer().size(), 1u);

    BOOST_TEST(bg::buffer_polygon_arrangement_certified(input, expanded, 3, q(3, 5), displacement) == status::success);
    BOOST_TEST_EQ(expanded.size(), 1u);
    if (expanded.size() == 1) { BOOST_TEST(expanded.front().inners().empty()); }
    BOOST_TEST_EQ(locate(boundary(expanded), {6, 8}), 1);
}

void test_split_and_empty()
{
    auto const input = dumbbell();
    multi out;
    BOOST_TEST(bg::buffer_polygon_arrangement_certified(input, out, q(-3, 2), q(1, 4), displacement) == status::success);
    BOOST_TEST_EQ(out.size(), 2u);
    // The left and right supports are 16 apart, hence no common witness.
    // The cross section x=8, 4<=y<=6 has clearance <=1 (project onto a
    // neck side), below r-e-d. Both chamber samples have clearance >r+d.
    std::vector<edge> src;
    add_boundary(src, input);
    auto const ring = boundary(out);
    number const upper = q(3, 2) + displacement;
    for (auto const& center : {xy{3, q(15, 4)}, xy{13, q(25, 4)}})
    {
        BOOST_TEST(distance2(src, center) > upper * upper);
        BOOST_TEST_EQ(locate(ring, center), 1);
    }
    BOOST_TEST_EQ(locate(ring, {8, 5}), -1);
    for (auto const& poly : out)
    {
        check_orientation(poly.outer(), false);
        BOOST_TEST(poly.inners().empty());
        if (poly.outer().empty()) { continue; }
        bool const left = bg::get<0>(poly.outer().front()) < 8;
        for (auto const& p : poly.outer()) { BOOST_TEST(left ? bg::get<0>(p) < 8 : bg::get<0>(p) > 8); }
    }

    input_polygon flat;
    close(flat.outer(), {{0, 0}, {4, 0}, {4, 2}, {0, 2}}, true);
    BOOST_TEST(bg::buffer_polygon_arrangement_certified(flat, out, -1, q(1, 4), 0) == status::success);
    BOOST_TEST(out.empty());
    BOOST_TEST(bg::buffer_polygon_arrangement_certified(flat, out, 0, q(1, 4), 0) == status::invalid_input);
    flat.outer().pop_back();
    BOOST_TEST(bg::buffer_polygon_arrangement_certified(flat, out, 1, q(1, 4), displacement) == status::invalid_input);
}

void benchmark(char const* label, input_polygon const& p, number const& radius, number const& error)
{
    constexpr unsigned iterations = 20;
    multi out;
    auto const start = std::chrono::steady_clock::now();
    for (unsigned i = 0; i < iterations; ++i)
    {
        BOOST_TEST(bg::buffer_polygon_arrangement_certified(p, out, radius, error, displacement) == status::success);
    }
    auto const elapsed = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start);
    std::size_t vertices = 0;
    for (auto const& component : out)
    {
        vertices += component.outer().size() - 1;
        for (auto const& h : component.inners()) { vertices += h.size() - 1; }
    }
    std::cout << label << ": " << elapsed.count() / iterations << " ms/call, " << vertices << " vertices\n";
}

int main(int argc, char** argv)
{
    if (argc == 2 && std::string(argv[1]) == "--benchmark")
    {
        auto const input = mixed_holes();
        benchmark("mixed holes dilation", input, 1, q(1, 5));
        benchmark("mixed holes erosion", input, -1, q(1, 5));
        benchmark("sheared dumbbell erosion", dumbbell(), q(-3, 2), q(1, 4));
    }
    else
    {
        test_holed_both_signs();
        test_split_and_empty();
    }
    return boost::report_errors();
}
