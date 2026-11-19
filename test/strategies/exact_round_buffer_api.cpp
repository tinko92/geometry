// Boost.Geometry
// Unit Test
// Copyright (c) 2026
// Distributed under the Boost Software License, Version 1.0.

#include <boost/geometry/algorithms/buffer_round_certified.hpp>
#include <boost/geometry/geometries/segment.hpp>
#include <boost/core/lightweight_test.hpp>

#include <algorithm>
#include <limits>

namespace bg = boost::geometry;
namespace ec = bg::strategies::detail::exact_cartesian;
using number = ec::exact_number;
using status = bg::certified_buffer_status;
template <typename T> using point = bg::model::point<T, 2, bg::cs::cartesian>;
using p = point<double>;
using line = bg::model::linestring<p>;
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

int main()
{
    test_point_and_segment();
    test_star_and_nonstar();
    test_polygon();
    test_irrational_polygon_arrangement();
    test_chord_erosion_empty_implication();
    test_empty_vs_failure();
    return boost::report_errors();
}
