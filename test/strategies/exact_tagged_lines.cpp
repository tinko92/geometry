// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.

#include <boost/geometry/strategies/cartesian/exact/tagged_orientation.hpp>

#include <array>
#include <cfenv>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <string>
#include <vector>

#include <boost/core/lightweight_test.hpp>
#include <boost/geometry/strategies/cartesian/exact/arithmetic.hpp>
#include <boost/geometry/strategies/cartesian/exact/fixed_rounding.hpp>

namespace ec = boost::geometry::strategies::detail::exact_cartesian;
namespace et = ec::tagged_orientation;
using rational = ec::exact_number;
using boost::multiprecision::cpp_int;
using triple = std::array<et::producers, 3>;
using tagged = std::array<et::intersection, 3>;
using fixed = std::array<ec::native_fixed_traits<double>::intersection_type, 3>;

static_assert(et::native::capacity == 101 && et::number<1, 2>::capacity == 101);
static_assert(et::number<2, 2>::capacity == 202 && et::number<2, 8>::capacity == 202);
static_assert(et::number<3, 8>::capacity == 303);
static_assert(et::minor_xy::capacity == 505 && et::minor_w::capacity == 605);
static_assert(et::determinant::capacity == 807);
static_assert(2 * et::number<3, 8>::capacity + et::number<2, 8>::capacity == 808);

et::producers at(double x, double y)
{
    return {{x, y, x == 0 ? 1.0 : 0.0, y, x, y, x, y == 0 ? 1.0 : 0.0}};
}

et::producers reverse_line(et::producers p)
{
    std::swap(p[0], p[2]);
    std::swap(p[1], p[3]);
    return p;
}

tagged prepare_tagged(triple const& p)
{
    return {{et::make_intersection(p[0]), et::make_intersection(p[1]), et::make_intersection(p[2])}};
}

fixed prepare_fixed(triple const& p)
{
    return {{ec::direct_native_homogeneous(p[0]), ec::direct_native_homogeneous(p[1]),
             ec::direct_native_homogeneous(p[2])}};
}

int tagged_sign(tagged const& p)
{ return et::orientation(p[0], p[1], p[2]); }

int fixed_sign(fixed const& p)
{
    // Mirror the current predicate's compact-first choice for cached points.
    ec::fixed_homogeneous<128> a, b, c;
    if (ec::try_fixed_narrow(p[0], a) && ec::try_fixed_narrow(p[1], b) && ec::try_fixed_narrow(p[2], c))
        return ec::fixed_orientation(a, b, c);
    return ec::fixed_orientation(p[0], p[1], p[2]);
}

std::array<rational, 2> cartesian(et::producers const& v)
{
    rational const ax(v[0]), ay(v[1]), bx(v[2]), by(v[3]);
    rational const cx(v[4]), cy(v[5]), dx(v[6]), dy(v[7]);
    rational const rx = bx - ax, ry = by - ay, sx = dx - cx, sy = dy - cy;
    rational const t = ((cx - ax) * sy - (cy - ay) * sx) / (rx * sy - ry * sx);
    return {{ax + t * rx, ay + t * ry}};
}

// Read back a radix record independently of carry/product implementation.
template <typename Number>
rational read_number(Number const& n)
{
    if (n.empty()) return 0;
    cpp_int value = 0;
    for (int k = n.last(); k > n.first();)
    {
        value *= et::radix;
        value += static_cast<std::int32_t>(*n.data(--k));
    }
    int const shift = et::digit_bits * n.first();
    if (shift >= 0) return rational(value << shift);
    return rational(value, cpp_int(1) << -shift);
}

int check(triple const& v, bool readback = false)
{
    auto const a = cartesian(v[0]), b = cartesian(v[1]), c = cartesian(v[2]);
    rational const d = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
    int const expected = (d > 0) - (d < 0);
    auto const p = prepare_tagged(v);
    BOOST_TEST_EQ(tagged_sign(p), expected);
    BOOST_TEST_EQ(fixed_sign(prepare_fixed(v)), expected);
    if (readback)
    {
        std::array<std::array<rational, 2>, 3> const rows{{a, b, c}};
        for (std::size_t i = 0; i != p.size(); ++i)
        {
            auto const w = read_number(p[i].w);
            BOOST_TEST_EQ(read_number(p[i].x) / w, rows[i][0]);
            BOOST_TEST_EQ(read_number(p[i].y) / w, rows[i][1]);
        }
    }
    return expected;
}

void boundary_cases()
{
    double const tiny = (std::numeric_limits<double>::denorm_min)();
    double const high = (std::numeric_limits<double>::max)();
    double const below = std::nextafter(high, 0.0);
    double const above_one = std::nextafter(1.0, 2.0);
    et::producers const max_cross{{-high, -high, high, high, -high, high, high, -high}};
    et::producers const half_min_x{{0, -1, tiny, 1, -1, 0, 1, 0}};
    et::producers const half_min_y{{-1, 0, 1, tiny, 0, -1, 0, 1}};
    et::producers const half_ulp_x{{1, 0, above_one, 2, 0, 1, 2, 1}};
    et::producers const outside_double{{0, 0, high, 1, 0, 2, high, 1.5}};

    triple const subnormal{{max_cross, half_min_x, half_min_y}};
    triple const rounding_tie{{at(1, 1), half_ulp_x, at(1, above_one)}};
    triple const cancellation{{at(tiny, 0), at(high, high), at(below, below)}};
    BOOST_TEST_EQ(check(subnormal, true), 1);
    BOOST_TEST_EQ(check(rounding_tie, true), 1);
    BOOST_TEST_EQ(check(cancellation, true), 1);
    BOOST_TEST_EQ(check({{at(tiny, tiny), at(2 * tiny, tiny), at(tiny, 2 * tiny)}}, true), 1);
    BOOST_TEST_EQ(check({{at(high, high), at(-high, high), at(high, -high)}}, true), 1);
    BOOST_TEST_EQ(check({{at(-0.0, 0.0), at(tiny, tiny), at(high, high)}}, true), 0);
    BOOST_TEST_EQ(check({{max_cross, reverse_line(max_cross), at(high, tiny)}}, true), 0);
    BOOST_TEST_EQ(check({{at(0, 0), outside_double, at(0, 1)}}, true), 1);
    BOOST_TEST(cartesian(outside_double)[0] > rational(high));

    // Cartesian nearest-even conversion really collapses the first two x
    // values, but it is never used by either exact predicate.
    auto const f = prepare_fixed(rounding_tie);
    double x1, x2;
    BOOST_TEST(ec::try_fixed_round(f[0], 0, x1));
    BOOST_TEST(ec::try_fixed_round(f[1], 0, x2));
    BOOST_TEST_EQ(x1, x2);

    for (auto v : {subnormal, rounding_tie, cancellation})
    {
        auto r = prepare_tagged(v);
        BOOST_TEST_EQ(et::orientation(r[1], r[2], r[0]), 1);
        BOOST_TEST_EQ(et::orientation(r[0], r[2], r[1]), -1);
        for (auto& p : v)
        {
            p = reverse_line(p);
            std::rotate(p.begin(), p.begin() + 4, p.end());
        }
        BOOST_TEST_EQ(check(v), 1);
    }
    for (unsigned mask = 1; mask != 8; ++mask)
    {
        auto v = rounding_tie;
        for (unsigned i = 0; i != v.size(); ++i)
            if (mask & (1u << i)) v[i] = reverse_line(v[i]);
        BOOST_TEST_EQ(check(v), 1);
    }

    auto const rounding = std::fegetround();
    for (auto mode : {FE_DOWNWARD, FE_UPWARD, FE_TOWARDZERO, FE_TONEAREST})
    {
        BOOST_TEST_EQ(std::fesetround(mode), 0);
        BOOST_TEST_EQ(tagged_sign(prepare_tagged(cancellation)), 1);
        BOOST_TEST_EQ(tagged_sign(prepare_tagged(subnormal)), 1);
    }
    BOOST_TEST_EQ(std::fesetround(rounding), 0);

    et::producers const parallel{{0, 0, 1, 0, 0, 1, 1, 1}};
    et::producers const coincident{{0, 0, 1, 0, 2, 0, 3, 0}};
    et::producers const degenerate{{0, 0, 0, 0, 0, 1, 1, 1}};
    auto invalid = at(0, 0);
    invalid[2] = std::numeric_limits<double>::infinity();
    BOOST_TEST_THROWS(et::make_intersection(parallel), std::domain_error);
    BOOST_TEST_THROWS(et::make_intersection(coincident), std::domain_error);
    BOOST_TEST_THROWS(et::make_intersection(degenerate), std::domain_error);
    BOOST_TEST_THROWS(et::make_intersection(invalid), std::domain_error);
    invalid[2] = std::numeric_limits<double>::quiet_NaN();
    BOOST_TEST_THROWS(et::make_intersection(invalid), std::domain_error);
}

double random_double(std::mt19937_64& rng)
{
    double value;
    do
    {
        std::uint64_t const bits = rng();
        std::memcpy(&value, &bits, sizeof(value));
    }
    while (!std::isfinite(value));
    return value;
}

struct corpus
{
    std::vector<triple> grid, near_tie, sparse_extreme, full_range;
};

corpus make_corpus()
{
    corpus result;
    std::mt19937_64 rng(0x1517);
    double const high = (std::numeric_limits<double>::max)();
    double const tiny = (std::numeric_limits<double>::denorm_min)();
    auto const sample = [](double x, double y)
    {
        return et::producers{{x, y, x + 3, y + 3, x, y + 1, x + 2, y}};
    };
    for (int i = 0; i != 128; ++i)
    {
        double const x = double(rng() % 661), y = double(rng() % 647);
        result.grid.push_back({{sample(x, y), sample(x + (i % 31), y + (i % 37)), sample(x + 9, y + 19)}});
        auto v = triple{{sample(x / 7, y / 3), sample(x / 7 + 2, y / 3 + 2), sample(x / 7 + 4, y / 3 + 4)}};
        v[2][7] = std::nextafter(v[2][7], (i % 2) ? 1024.0 : -1024.0);
        result.near_tie.push_back(v);

        double const a = std::ldexp(high, -i);
        double const b = std::nextafter(a, 0.0);
        result.sparse_extreme.push_back({{at(tiny * (i + 1), 0), at(a, a), at(b, b)}});
        if (i % 2 == 0)
        {
            for (auto& p : v)
                for (auto& coordinate : p) coordinate = random_double(rng);
            result.full_range.push_back(v);
        }
    }
    return result;
}

template <typename F>
void measure(char const* fixture, char const* path, std::size_t count, int repeats, F evaluate)
{
    using clock = std::chrono::steady_clock;
    int checksum = 0;
    auto const start = clock::now();
    for (int pass = 0; pass != repeats; ++pass)
        for (std::size_t i = 0; i != count; ++i) checksum += evaluate(i);
    double const ns = std::chrono::duration<double, std::nano>(clock::now() - start).count();
    std::cout << fixture << '\t' << path << '\t' << count * repeats << '\t'
              << std::fixed << std::setprecision(1) << ns / (count * repeats) << '\t' << checksum << '\n';
}

void benchmark(char const* name, std::vector<triple> const& v)
{
    std::vector<tagged> t;
    std::vector<fixed> f;
    t.reserve(v.size());
    f.reserve(v.size());
    for (auto const& p : v)
    {
        t.push_back(prepare_tagged(p));
        f.push_back(prepare_fixed(p));
    }
    measure(name, "fixed-prepared", v.size(), 32, [&](auto i) { return fixed_sign(f[i]); });
    measure(name, "tagged-prepared", v.size(), 32, [&](auto i) { return tagged_sign(t[i]); });
    measure(name, "fixed-cold", v.size(), 8, [&](auto i) { return fixed_sign(prepare_fixed(v[i])); });
    measure(name, "tagged-cold", v.size(), 8, [&](auto i) { return tagged_sign(prepare_tagged(v[i])); });
}

int main(int argc, char** argv)
{
    boundary_cases();
    auto const cases = make_corpus();
    for (auto const* values : {&cases.grid, &cases.near_tie, &cases.sparse_extreme, &cases.full_range})
        for (std::size_t i = 0; i != values->size(); ++i)
            check((*values)[i], i < 8);
    if (argc == 2 && std::string(argv[1]) == "--benchmark")
    {
        std::cout << "record_bytes\tfixed=" << sizeof(fixed) / 3 << "\ttagged=" << sizeof(et::intersection) << '\n';
        std::cout << "fixture\tpath\tcount\tns_per_sign\tchecksum\n";
        benchmark("integer-grid", cases.grid);
        benchmark("near-collinear", cases.near_tie);
        benchmark("sparse-extreme", cases.sparse_extreme);
        benchmark("full-range", cases.full_range);
    }
    return boost::report_errors();
}
