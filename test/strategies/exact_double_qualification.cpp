// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.

#include <boost/core/lightweight_test.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/strategies/relate/cartesian_exact.hpp>
#include <array>
#include <cfenv>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <limits>
#include <random>

namespace bg = boost::geometry;
namespace ec = bg::strategies::detail::exact_cartesian;
namespace cr = ec::construction_rounding;
namespace fe = ec::forward_error;
namespace policy = bg::strategies::relate;
using q = boost::multiprecision::cpp_rational;
using point = bg::model::d2::point_xy<double>;
using input = std::array<double, 8>;
std::size_t orientations = 0, intersections = 0, parallels = 0, eft_residuals = 0;

std::uint64_t bits(double value)
{ std::uint64_t result; std::memcpy(&result, &value, sizeof(value)); return result; }
double previous(double value) { return std::nextafter(value, -std::numeric_limits<double>::infinity()); }
double next(double value) { return std::nextafter(value, std::numeric_limits<double>::infinity()); }

// Test the candidate's exact rounding cell, never a floating conversion oracle.
bool nearest_cell(double candidate, q const& exact)
{
    if (!std::isfinite(candidate)) return false;
    q const lower = (q(previous(candidate)) + q(candidate)) / 2;
    q const upper = (q(candidate) + q(next(candidate))) / 2;
    bool const even = (bits(candidate) & 1) == 0;
    return (exact > lower || (even && exact == lower))
        && (exact < upper || (even && exact == upper));
}

void orientation(input const& v)
{
    q const determinant = (q(v[2]) - q(v[0])) * (q(v[5]) - q(v[1]))
                        - (q(v[3]) - q(v[1])) * (q(v[4]) - q(v[0]));
    int const expected = determinant.sign();
    point const a(v[0], v[1]), b(v[2], v[3]), c(v[4], v[5]);
    BOOST_TEST_EQ(ec::native_orientation<double>(a, b, c), expected);
    BOOST_TEST_EQ(ec::native_orientation<double>(b, a, c), -expected);
    BOOST_TEST_EQ(ec::native_orientation<double>(c, a, b), expected);
    fe::point_approximation pa, pb, pc;
    BOOST_TEST(fe::try_native_point(v[0], v[1], pa));
    BOOST_TEST(fe::try_native_point(v[2], v[3], pb));
    BOOST_TEST(fe::try_native_point(v[4], v[5], pc));
    int certified = 17;
    if (fe::try_orientation(pa, pb, pc, certified)) BOOST_TEST_EQ(certified, expected);
    ++orientations;
}

bool intersection_oracle(input const& v, std::array<q, 2>& output)
{
    std::array<q, 8> p;
    for (unsigned i = 0; i != 8; ++i) p[i] = q(v[i]);
    q const a = p[1] - p[3], b = p[2] - p[0], c = p[0] * p[3] - p[1] * p[2];
    q const d = p[5] - p[7], e = p[6] - p[4], f = p[4] * p[7] - p[5] * p[6];
    q const w = a * e - b * d;
    if (w == 0) return false;
    output = {{(b * f - c * e) / w, (c * d - a * f) / w}};
    return true;
}

void enclosure(ec::interval const& bound, q const& exact)
{
    BOOST_TEST(bound.valid());
    BOOST_TEST(std::isfinite(bound.approximate) && std::isfinite(bound.error));
    if (std::isfinite(bound.approximate) && std::isfinite(bound.error))
    {
        BOOST_TEST(q(bound.approximate) - q(bound.error) <= exact);
        BOOST_TEST(q(bound.approximate) + q(bound.error) >= exact);
    }
}

template <typename Policy>
void construction(input const& v, std::array<q, 2> const& exact)
{
    ec::point_cache_scope<double> scope;
    cr::point_state<double, Policy> state;
    state.init(v);
    for (unsigned axis = 0; axis != 2; ++axis)
    {
        double const stored = state.approximate(axis);
        enclosure(state.bounds(axis), exact[axis]);
        if constexpr (cr::policy_traits<Policy>::store_errors)
        {
            BOOST_TEST_EQ(bits(state.template emit<double>(axis, v)), bits(stored));
            BOOST_TEST_EQ(bits(state.approximate(axis)), bits(stored));
        }
        else
        {
            BOOST_TEST(nearest_cell(stored, exact[axis]));
            BOOST_TEST(nearest_cell(state.template emit<double>(axis, v), exact[axis]));
        }
        auto const refined = state.refine(v, axis);
        enclosure(refined, exact[axis]);
        BOOST_TEST(nearest_cell(refined.approximate, exact[axis]));
        BOOST_TEST_EQ(bits(state.template emit<double>(axis, v)), bits(refined.approximate));
    }
}

void check(input const& v)
{
    orientation(v);
    std::array<q, 2> exact;
    fe::point_approximation bound;
    bool const unique = intersection_oracle(v, exact);
    bool const accepted = fe::try_intersection(v, bound);
    BOOST_TEST(!accepted || unique);
    if (!unique) { ++parallels; return; }
    if (accepted) for (unsigned axis = 0; axis != 2; ++axis)
        enclosure({bound.value[axis], bound.error[axis]}, exact[axis]);
    construction<policy::approximate_output>(v, exact);
    construction<policy::nearest_output>(v, exact);
    ++intersections;
}

void eft_and_ties()
{
    double const u = std::ldexp(1.0, -52);
    for (auto const& ab : std::array<std::array<double, 2>, 4>{{
        {{1 + u, 1 - u}}, {{0.1, 0.2}}, {{134217729, 134217727}}, {{1.001, -2.003}}}})
    {
        double const head = ec::arithmetic_detail::rounded_mul(ab[0], ab[1]);
        double const tail = std::fma(ab[0], ab[1], -head);
        BOOST_TEST(q(head) + q(tail) == q(ab[0]) * q(ab[1]));
        if (tail != 0) ++eft_residuals;
        ec::construction_detail::expansion<8> expansion;
        BOOST_TEST(expansion.product(ab[0], ab[1]));
        q total = 0;
        for (std::size_t i = 0; i != expansion.size; ++i) total += q(expansion.components[i]);
        BOOST_TEST(total == q(ab[0]) * q(ab[1]));
        auto const sum = ec::orientation_detail::two_sum(ab[0], ab[1]);
        auto const difference = ec::orientation_detail::two_diff(ab[0], ab[1]);
        BOOST_TEST(q(sum[0]) + q(sum[1]) == q(ab[0]) + q(ab[1]));
        BOOST_TEST(q(difference[0]) + q(difference[1]) == q(ab[0]) - q(ab[1]));
    }
    BOOST_TEST(eft_residuals > 0);
    q const midpoint = (q(1) + q(next(1))) / 2;
    BOOST_TEST(nearest_cell(1, midpoint)); BOOST_TEST(!nearest_cell(next(1), midpoint));
    q const upper_midpoint = (q(next(1)) + q(next(next(1)))) / 2;
    BOOST_TEST(!nearest_cell(next(1), upper_midpoint)); BOOST_TEST(nearest_cell(next(next(1)), upper_midpoint));
    for (double sign : {-1.0, 1.0}) for (unsigned offset : {0u, 1u})
        check({{sign * (1 + offset * u), 0, sign * (1 + (offset + 1) * u), 2, 0, 1, 2, 1}});
}

template <typename Policy>
void public_operations()
{
    using polygon = bg::model::polygon<point>;
    using multi = bg::model::multi_polygon<polygon>;
    polygon a, b;
    bg::read_wkt("POLYGON((0 0,0 2,2 2,2 0,0 0))", a);
    bg::read_wkt("POLYGON((1 1,1 3,3 3,3 1,1 1))", b);
    bg::model::polygon<point, false, false> reversed_open;
    bg::convert(b, reversed_open);
    bg::strategies::relate::cartesian_exact<double, Policy> strategy;
    multi actual[4];
    bg::union_(a, reversed_open, actual[0], strategy);
    bg::intersection(a, reversed_open, actual[1], strategy);
    bg::difference(a, reversed_open, actual[2], strategy);
    bg::sym_difference(a, reversed_open, actual[3], strategy);
    char const* expected[] = {
        "MULTIPOLYGON(((0 0,0 2,1 2,1 3,3 3,3 1,2 1,2 0,0 0)))",
        "MULTIPOLYGON(((1 1,1 2,2 2,2 1,1 1)))",
        "MULTIPOLYGON(((0 0,0 2,1 2,1 1,2 1,2 0,0 0)))",
        "MULTIPOLYGON(((0 0,0 2,1 2,1 1,2 1,2 0,0 0)),((1 2,1 3,3 3,3 1,2 1,2 2,1 2)))"};
    for (unsigned i = 0; i != 4; ++i)
    {
        multi control; bg::read_wkt(expected[i], control);
        BOOST_TEST(bg::is_valid(actual[i]));
        BOOST_TEST_EQ(bg::relation(actual[i], control).str(), "2FFF1FFF2");
    }
    BOOST_TEST_EQ(bg::relation(a, reversed_open, strategy).str(), "212101212");
    bg::model::linestring<point> const first{{-1, 1}, {3, 1}}, second{{1, -1}, {1, 3}};
    BOOST_TEST_EQ(bg::relation(first, second, strategy).str(), "0F1FF0102");
    BOOST_TEST_EQ(bg::relation(first, a, strategy).str(), "101FF0212");
    BOOST_TEST(!bg::within(point(0, 1), a, strategy));
    BOOST_TEST(bg::covered_by(point(0, 1), a, strategy));
}

int main()
{
    BOOST_TEST_EQ(std::fegetround(), FE_TONEAREST);
    eft_and_ties();
    double const n = std::ldexp(1.0, 27);
    check({{0, 0, n, n - 1, n - 1, n - 2, 1, 0}});
    check({{0, 0, 4, 4, 0, 1, 4, 5}});
    check({{0, 0, 4, 4, 1, 1, 3, 3}});
    check({{0, 0, 1, 1, 0, 1, 1, -1}});
    std::mt19937_64 random(0x564f4c4154494c45ULL);
    for (unsigned i = 0; i != 2048; ++i)
    {
        input v;
        for (double& x : v) x = double(std::int64_t(random() % 2000001) - 1000000) / 1000;
        check(v);
        if (i % 4 == 0) { double const e = std::ldexp(1.0, -40 - int(i % 9)); check({{0, 0, 1, 1, 0, e, 1, 1 - e}}); }
        if (i % 8 == 0) check({{0, 0, n + 1, n, 0, 1, n, n}});
        if (i % 16 == 0) check({{0, 0, 4, 4, 2, 2, 6, 6}});
    }
    public_operations<policy::approximate_output>();
    public_operations<policy::nearest_output>();
    std::cout << orientations << " double orientation cases; " << intersections << " unique intersections x2 policies; "
              << parallels << " parallel/collinear cases; " << eft_residuals << " nonzero FMA residuals.\n";
    return boost::report_errors();
}
