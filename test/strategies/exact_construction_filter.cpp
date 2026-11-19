// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Licensed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt

#include <boost/core/lightweight_test.hpp>
#include <boost/geometry/strategies/cartesian/exact/construction_filter.hpp>
#include <boost/geometry/strategies/cartesian/exact/point.hpp>

#include <array>
#include <atomic>
#include <cfenv>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <new>
#include <random>
#include <thread>
#include <vector>

namespace ec = boost::geometry::strategies::detail::exact_cartesian;
namespace cd = ec::construction_detail;
using rational = ec::exact_number;
using endpoints = std::array<double, 8>;

namespace
{
thread_local bool counting = false;
thread_local std::size_t allocations = 0;
}

void* operator new(std::size_t size)
{
    if (counting) ++allocations;
    if (void* result = std::malloc(size == 0 ? 1 : size)) return result;
    throw std::bad_alloc();
}
void operator delete(void* value) noexcept { std::free(value); }
void* operator new[](std::size_t size) { return ::operator new(size); }
void operator delete[](void* value) noexcept { ::operator delete(value); }
#if defined(__cpp_sized_deallocation)
void operator delete(void* value, std::size_t) noexcept { std::free(value); }
void operator delete[](void* value, std::size_t) noexcept { std::free(value); }
#endif

bool reference(endpoints const& v, unsigned dimension, rational& result)
{
    std::array<rational, 8> q;
    for (unsigned i = 0; i < 8; ++i) q[i] = ec::exact_value(v[i]);
    rational const ax = q[2] - q[0], ay = q[1] - q[3];
    rational const bx = q[6] - q[4], by = q[5] - q[7];
    rational const a = q[0] * q[3] - q[1] * q[2];
    rational const b = q[4] * q[7] - q[5] * q[6];
    rational const w = ay * bx - ax * by;
    if (w == 0) return false;
    result = (dimension == 0 ? rational(ax * b - a * bx) : rational(a * by - ay * b)) / w;
    return true;
}

template <typename Output>
bool check(endpoints const& v, unsigned dimension)
{
    Output result = Output(17);
    allocations = 0;
    counting = true;
    bool const certain = ec::try_round_intersection(v, dimension, result);
    counting = false;
    BOOST_TEST_EQ(allocations, 0u);
    if (!certain) { BOOST_TEST(result == Output(17)); return false; }
    rational exact;
    BOOST_TEST(reference(v, dimension, exact));
    Output const expected = ec::point_detail::output_rounder<Output>::apply(exact);
    BOOST_TEST(result == expected);
    if (std::is_floating_point<Output>::value && result == Output(0))
        BOOST_TEST_EQ(std::signbit(double(result)), std::signbit(double(expected)));
    return true;
}

void check_interval(endpoints const& v, unsigned dimension)
{
    ec::interval bound{17, 19};
    allocations = 0;
    counting = true;
    bool const certain = ec::try_intersection_interval(v, dimension, bound);
    counting = false;
    BOOST_TEST_EQ(allocations, 0u);
    if (!certain) { BOOST_TEST_EQ(bound.approximate, 17); BOOST_TEST_EQ(bound.error, 19); return; }
    rational q;
    BOOST_TEST(reference(v, dimension, q));
    BOOST_TEST(bound.valid());
    BOOST_TEST(ec::exact_value(bound.lower()) <= q);
    BOOST_TEST(ec::exact_value(bound.upper()) >= q);
    BOOST_TEST_EQ(bound.approximate, ec::point_detail::output_rounder<double>::apply(q));
}

bool check_coordinate_comparison(endpoints const& v, unsigned dimension, double target)
{
    int result = 17;
    allocations = 0;
    counting = true;
    bool const certain = ec::try_compare_intersection_coordinate(v, dimension, target, result);
    counting = false;
    BOOST_TEST_EQ(allocations, 0u);
    if (!certain) { BOOST_TEST_EQ(result, 17); return false; }
    rational q;
    BOOST_TEST(reference(v, dimension, q));
    rational const t = ec::exact_value(target);
    BOOST_TEST_EQ(result, int(q > t) - int(q < t));
    return true;
}

bool check_coordinates_comparison(endpoints const& first, endpoints const& second, unsigned dimension)
{
    int result = 17;
    allocations = 0;
    counting = true;
    bool const certain = ec::try_compare_intersection_coordinates(first, second, dimension, result);
    counting = false;
    BOOST_TEST_EQ(allocations, 0u);
    if (!certain) { BOOST_TEST_EQ(result, 17); return false; }
    rational a, b;
    BOOST_TEST(reference(first, dimension, a));
    BOOST_TEST(reference(second, dimension, b));
    BOOST_TEST_EQ(result, int(a > b) - int(a < b));
    return true;
}

endpoints permute(endpoints v, unsigned mask)
{
    if (mask & 1) { std::swap(v[0], v[2]); std::swap(v[1], v[3]); }
    if (mask & 2) { std::swap(v[4], v[6]); std::swap(v[5], v[7]); }
    if (mask & 4) for (unsigned i = 0; i < 4; ++i) std::swap(v[i], v[i + 4]);
    return v;
}

void ordinary()
{
    for (endpoints const& v : {
        endpoints{{0, 0, 3, 3, 0, 3, 3, 0}},
        endpoints{{-10, -5, 6, 8, 8, -4, -2, 9}},
        endpoints{{0, 0, 1, 1, -1, 0, 1, 0}},
        endpoints{{4710043, -2665034, 4685409, -2666385, 4700000, -2500000, 4700000, -3500000}}})
        for (unsigned mask = 0; mask < 8; ++mask)
            for (unsigned dimension = 0; dimension < 2; ++dimension)
            {
                auto const p = permute(v, mask);
                BOOST_TEST(check<double>(p, dimension));
                BOOST_TEST(check<float>(p, dimension));
                BOOST_TEST(check<std::int64_t>(p, dimension));
                check_interval(p, dimension);
                auto const existing = ec::construct_intersection(
                    ec::exact_point<double>(p[0], p[1]), ec::exact_point<double>(p[2], p[3]),
                    ec::exact_point<double>(p[4], p[5]), ec::exact_point<double>(p[6], p[7]));
                rational q;
                BOOST_TEST(reference(p, dimension, q));
                BOOST_TEST(q == existing.exact_coordinate(dimension));
            }
    double zero = -0.0;
    BOOST_TEST(ec::try_round_intersection(endpoints{{0, 0, 1, 1, -1, 0, 1, 0}}, 0, zero));
    BOOST_TEST(!std::signbit(zero));
}

void ties()
{
    double const one_up = std::nextafter(1.0, 2.0);
    double const two_up = std::nextafter(one_up, 2.0);
    double const one_down = std::nextafter(1.0, 0.0);
    for (double lo : {one_down, 1.0, one_up})
    {
        double const hi = std::nextafter(lo, 2.0);
        for (double height : {0.5, 1.0, 1.5}) for (double sign : {-1.0, 1.0})
        {
            endpoints const v{{sign * lo, 0, sign * hi, 2, 0, height, 3, height}};
            for (unsigned mask = 0; mask < 8; ++mask)
                BOOST_TEST(check<double>(permute(v, mask), 0));
        }
    }
    double d = 17;
    BOOST_TEST(ec::try_round_intersection(endpoints{{1, 0, one_up, 2, 0, 1, 3, 1}}, 0, d));
    BOOST_TEST_EQ(d, 1.0);
    BOOST_TEST(ec::try_round_intersection(endpoints{{one_up, 0, two_up, 2, 0, 1, 3, 1}}, 0, d));
    BOOST_TEST_EQ(d, two_up);

    // Double-rounding adversaries: q is just above/below a float midpoint,
    // although rounding q first to double can hit the midpoint exactly.
    for (float lo : {1.0f, std::nextafter(1.0f, 2.0f)})
        for (double perturb : {-std::ldexp(1.0, -52), 0.0, std::ldexp(1.0, -52)})
            for (double sign : {-1.0, 1.0})
            {
                double const hi = double(std::nextafter(lo, 2.0f)) + perturb;
                endpoints const v{{sign * double(lo), 0, sign * hi, 2, 0, 1, 3, 1}};
                BOOST_TEST(check<float>(v, 0));
            }

    for (double lo : {-129, -128, -3, -2, -1, 0, 1, 2, 126, 127})
    {
        endpoints const v{{lo, 0, lo + 1, 2, 0, 1, 3, 1}};
        BOOST_TEST(check<std::int32_t>(v, 0));
        if (lo >= -128 && lo < 127) BOOST_TEST(check<std::int8_t>(v, 0));
    }
    double const limit = std::ldexp(1.0, 52);
    BOOST_TEST(check<std::int64_t>(endpoints{{limit, 0, limit + 1, 2, 0, 1, 3, 1}}, 0));
    BOOST_TEST(check<std::uint64_t>(endpoints{{limit, 0, limit + 1, 2, 0, 1, 3, 1}}, 0));
    BOOST_TEST(check<std::uint32_t>(endpoints{{0, 0, 1, 2, 0, 1, 3, 1}}, 0));
}

void cancellation_and_guards()
{
    for (int exponent : {-150, -100, -1, 0, 1, 100, 149})
    {
        double const base = std::ldexp(1.0, exponent);
        double const step = std::nextafter(base, 2 * base) - base;
        endpoints const v{{base, base, base + step, base + 2 * step,
                            base, base + 2 * step, base + step, base}};
        for (unsigned mask = 0; mask < 8; ++mask)
            for (unsigned dimension = 0; dimension < 2; ++dimension)
            {
                BOOST_TEST(check<double>(permute(v, mask), dimension));
                check_interval(permute(v, mask), dimension);
            }
    }
    double const minimum = std::ldexp(1.0, -150), maximum = std::ldexp(1.0, 150);
    BOOST_TEST(check<double>(endpoints{{-maximum, -maximum, maximum, maximum,
                                       -maximum, maximum, maximum, -maximum}}, 0));
    BOOST_TEST(check<double>(endpoints{{0, minimum, 1, 1, 0, 1, 1, minimum}}, 0));
    endpoints const regular{{0, 0, 1, 1, 0, 1, 1, 0}};
    for (double unsupported : {std::nextafter(minimum, 0.0), std::nextafter(maximum, 2 * maximum),
          (std::numeric_limits<double>::infinity)(), (std::numeric_limits<double>::quiet_NaN)()})
        for (unsigned i = 0; i < 8; ++i)
        {
            auto v = regular; v[i] = unsupported;
            BOOST_TEST(!check<double>(v, 0));
        }
    // Unsupported quotient range and non-unique lines are fallbacks, not
    // overflow repairs. The caller retains its existing exact path.
    double const large = std::ldexp(1.0, 100);
    BOOST_TEST(!check<double>(endpoints{{0, 0, large, large, 0, large,
        large, std::nextafter(2 * large, 3 * large)}}, 0));
    BOOST_TEST(!check<double>(endpoints{{0, 0, 1, minimum, minimum, 0, 0, 1}}, 1));
    for (endpoints const& v : {endpoints{{0, 0, 1, 1, 0, 1, 1, 2}},
                               endpoints{{0, 0, 1, 1, 2, 2, 3, 3}},
                               endpoints{{0, 0, 0, 0, 1, 0, 1, 1}}})
    {
        BOOST_TEST(!check<double>(v, 0));
        check_interval(v, 0);
    }
    double value = 17;
    BOOST_TEST(!ec::try_round_intersection(regular, 2, value)); BOOST_TEST_EQ(value, 17);
    long double wide = 17;
    BOOST_TEST(!ec::try_round_intersection(regular, 0, wide)); BOOST_TEST(wide == 17);
    rational exact = 17;
    BOOST_TEST(!ec::try_round_intersection(regular, 0, exact)); BOOST_TEST(exact == 17);
    bool boolean = true;
    BOOST_TEST(!ec::try_round_intersection(regular, 0, boolean)); BOOST_TEST(boolean);
    BOOST_TEST(!check<std::int8_t>(endpoints{{128, 0, 128, 2, 0, 1, 3, 1}}, 0));
    BOOST_TEST(!check<std::int64_t>(endpoints{{2 * maximum, 0, 2 * maximum, 2, 0, 1, 3, 1}}, 0));
    // Capacity exhaustion must be visible, with no silently lost low limb.
    cd::expansion<1> short_expansion;
    BOOST_TEST(short_expansion.add(1.0));
    BOOST_TEST(!short_expansion.add(std::ldexp(1.0, -53)));
    BOOST_TEST(!short_expansion.valid);
    BOOST_TEST(!short_expansion.add(2.0));
}

double random_value(std::mt19937_64& generator, bool wide)
{
    if (!wide) return (double(generator() % 2000001) - 1000000) / 1024;
    int const exponent = int(generator() % 300) - 150;
    double const mantissa = 1 + std::ldexp(double(generator() & ((std::uint64_t(1) << 52) - 1)), -52);
    return std::ldexp(generator() & 1 ? mantissa : -mantissa, exponent);
}

void randomized()
{
    std::mt19937_64 generator(429151);
    unsigned ordinary_certified = 0, ordinary_calls = 0, wide_certified = 0;
    for (bool wide : {false, true}) for (unsigned i = 0; i < 1024; ++i)
    {
        endpoints v;
        for (auto& value : v) value = random_value(generator, wide);
        for (unsigned dimension = 0; dimension < 2; ++dimension)
        {
            bool const certified = check<double>(v, dimension);
            if (wide) wide_certified += certified;
            else { ++ordinary_calls; ordinary_certified += certified; }
            check<float>(v, dimension);
            if (!wide) check<std::int64_t>(v, dimension);
            check_interval(v, dimension);
        }
    }
    BOOST_TEST(ordinary_certified > ordinary_calls * 99 / 100);
    BOOST_TEST(wide_certified > 512);
    std::cout << "construction double certified ordinary=" << ordinary_certified << '/' << ordinary_calls
              << " wide=" << wide_certified << "/2048\n";
}

void coordinate_comparisons()
{
    std::array<endpoints, 5> const cases{{
        {{0, 0, 3, 3, 0, 3, 3, 0}},
        {{1.5, -2, 1.5, 5, -1, 1.5, 4, 1.5}},
        {{0, 0, 1, 1, 0, 1, 2, 0}},
        {{-5, 7, 5, 7, 0, -3, 0, 9}},
        {{-10, -5, 6, 8, 8, -4, -2, 9}}
    }};
    // Reversal of either line changes denominator signs. Swapping the lines
    // and using different producer segments must preserve exact equality.
    for (auto const& first : cases) for (auto const& second : cases)
        for (unsigned a = 0; a < 8; ++a) for (unsigned b = 0; b < 8; ++b)
            for (unsigned dimension = 0; dimension < 2; ++dimension)
                BOOST_TEST(check_coordinates_comparison(permute(first, a), permute(second, b), dimension));
    for (auto const& item : cases) for (unsigned mask = 0; mask < 8; ++mask)
        for (unsigned dimension = 0; dimension < 2; ++dimension)
        {
            auto const v = permute(item, mask);
            rational q;
            BOOST_TEST(reference(v, dimension, q));
            double const rounded = ec::point_detail::output_rounder<double>::apply(q);
            for (double target : {0.0, -1.0, 1.0, 1.5, 7.0, rounded,
                    std::nextafter(rounded, -1e100), std::nextafter(rounded, 1e100)})
            {
                // A neighbour of exact zero is outside the conservative guard.
                if (target != 0 && std::abs(target) < std::ldexp(1.0, -150)) continue;
                BOOST_TEST(check_coordinate_comparison(v, dimension, target));
            }
        }

    // These two implicit x coordinates round to the SAME double yet differ
    // by 2^-106. Rounded output/interval equality cannot certify this ordering.
    double const next = std::nextafter(1.0, 2.0);
    endpoints const a{{1, 0, next, 2, 0, 0.5, 3, 0.5}};
    double const height = std::nextafter(0.5, 1.0);
    endpoints const b{{1, 0, next, 2, 0, height, 3, height}};
    rational qa, qb;
    BOOST_TEST(reference(a, 0, qa)); BOOST_TEST(reference(b, 0, qb));
    BOOST_TEST(qa < qb);
    BOOST_TEST_EQ(ec::point_detail::output_rounder<double>::apply(qa),
                  ec::point_detail::output_rounder<double>::apply(qb));
    BOOST_TEST(check_coordinates_comparison(a, b, 0));
    BOOST_TEST(check_coordinates_comparison(b, a, 0));
    BOOST_TEST(check_coordinate_comparison(a, 0, 1));

    for (int exponent : {-50, -25, 0, 25, 49})
    {
        double const base = std::ldexp(1.0, exponent);
        double const step = std::nextafter(base, 2 * base) - base;
        endpoints const first{{base, base, base + step, base + 2 * step,
                                base, base + 2 * step, base + step, base}};
        auto second = first;
        second[7] += step;
        // At exponent -50 the nonzero degree-five difference reaches 2^-510.
        for (unsigned mask = 0; mask < 8; ++mask) for (unsigned dimension = 0; dimension < 2; ++dimension)
        {
            BOOST_TEST(check_coordinates_comparison(first, permute(second, mask), dimension));
            BOOST_TEST(check_coordinates_comparison(second, permute(first, mask), dimension));
            BOOST_TEST(check_coordinate_comparison(permute(first, mask), dimension, base));
        }
    }
    double const maximum = std::ldexp(1.0, 50);
    endpoints const boundary{{-maximum, -maximum, maximum, maximum,
                              -maximum, maximum, maximum, -maximum}};
    BOOST_TEST(check_coordinates_comparison(boundary, boundary, 0));
    BOOST_TEST(check_coordinate_comparison(boundary, 0, 0));
}

void comparison_guards()
{
    endpoints const v{{0, 0, 3, 3, 0, 3, 3, 0}};
    for (endpoints const& invalid : {endpoints{{0, 0, 1, 1, 0, 1, 1, 2}},
                                     endpoints{{0, 0, 1, 1, 2, 2, 3, 3}},
                                     endpoints{{0, 0, 0, 0, 1, 0, 1, 1}}})
    {
        BOOST_TEST(!check_coordinate_comparison(invalid, 0, 0));
        BOOST_TEST(!check_coordinates_comparison(invalid, v, 0));
        BOOST_TEST(!check_coordinates_comparison(v, invalid, 0));
    }
    for (int exponent : {50, 150})
    {
        double const minimum = std::ldexp(1.0, -exponent), maximum = std::ldexp(1.0, exponent);
        for (double outside : {std::nextafter(minimum, 0.0), std::nextafter(maximum, 2 * maximum),
                (std::numeric_limits<double>::infinity)(), (std::numeric_limits<double>::quiet_NaN)()})
        {
            for (unsigned i = 0; i < 8; ++i)
            {
                auto invalid = v; invalid[i] = outside;
                if (exponent == 50)
                {
                    BOOST_TEST(!check_coordinates_comparison(invalid, v, 0));
                    BOOST_TEST(!check_coordinates_comparison(v, invalid, 0));
                }
                else BOOST_TEST(!check_coordinate_comparison(invalid, 0, 0));
            }
            if (exponent == 150) BOOST_TEST(!check_coordinate_comparison(v, 0, outside));
        }
    }
    BOOST_TEST(!check_coordinate_comparison(v, 2, 0));
    BOOST_TEST(!check_coordinates_comparison(v, v, 2));
    cd::expansion<2> limbs, one;
    BOOST_TEST(limbs.add(1)); BOOST_TEST(limbs.add(std::ldexp(1.0, -53))); BOOST_TEST(one.add(1));
    cd::expansion<1> too_short;
    BOOST_TEST(!cd::add_product(too_short, limbs, one));
    BOOST_TEST(!too_short.valid);
}

void randomized_comparisons()
{
    std::mt19937_64 generator(218499);
    unsigned scalar = 0, pair = 0;
    for (bool wide : {false, true}) for (unsigned i = 0; i < 512; ++i)
    {
        endpoints a, b;
        for (auto* item : {&a, &b}) for (auto& value : *item)
        {
            if (!wide) value = random_value(generator, false);
            else
            {
                int const exponent = int(generator() % 100) - 50;
                double const mantissa = 1 + std::ldexp(double(generator() & ((std::uint64_t(1) << 52) - 1)), -52);
                value = std::ldexp(generator() & 1 ? mantissa : -mantissa, exponent);
            }
        }
        for (unsigned dimension = 0; dimension < 2; ++dimension)
        {
            pair += check_coordinates_comparison(a, b, dimension);
            scalar += check_coordinate_comparison(a, dimension, random_value(generator, true));
            BOOST_TEST(check_coordinates_comparison(a, a, dimension));
        }
    }
    BOOST_TEST(pair > 2000); BOOST_TEST(scalar > 2000);
    std::cout << "construction comparisons certified scalar=" << scalar << "/2048 pair=" << pair << "/2048\n";
}

void active_expansion_storage()
{
    using expansion = ec::construction_detail::expansion<512>;
    expansion value;
    value.components.fill(17);
    BOOST_TEST(value.add(1));
    BOOST_TEST_EQ(value.size, 1u);
    BOOST_TEST_EQ(value.components[0], 1);
    BOOST_TEST_EQ(value.components[1], 17);
    BOOST_TEST_EQ(value.components.back(), 17);

    expansion copy;
    copy.components.fill(23);
    copy = value;
    BOOST_TEST_EQ(copy.size, 1u);
    BOOST_TEST_EQ(copy.components[0], 1);
    BOOST_TEST_EQ(copy.components[1], 23);
    BOOST_TEST_EQ(copy.components.back(), 23);
    copy = copy;
    BOOST_TEST_EQ(copy.sign(), 1);
    expansion moved(std::move(copy));
    BOOST_TEST_EQ(moved.sign(), 1);
    BOOST_TEST(moved.add(-1));
    BOOST_TEST_EQ(moved.size, 0u);
    BOOST_TEST_EQ(moved.sign(), 0);
}

void modes_and_threads()
{
    int const original = std::fegetround();
    endpoints const v{{0, 0, 3, 3, 0, 3, 3, 0}};
    for (int mode : {FE_TONEAREST})
    {
        if (std::fesetround(mode) != 0) continue;
        bool const expected = mode == FE_TONEAREST;
        BOOST_TEST_EQ(check<double>(v, 0), expected);
        BOOST_TEST_EQ(check<float>(v, 0), expected);
        BOOST_TEST_EQ(check<std::int64_t>(v, 0), expected);
        BOOST_TEST_EQ(check_coordinate_comparison(v, 0, 1.5), expected);
        BOOST_TEST_EQ(check_coordinates_comparison(v, v, 0), expected);
        check_interval(v, 0);
        BOOST_TEST_EQ(std::fegetround(), mode);
    }
    std::fesetround(original);
    std::atomic<bool> correct{true};
    std::vector<std::thread> workers;
    for (unsigned worker = 0; worker != 4; ++worker)
    {
        int const mode = FE_TONEAREST;
        workers.emplace_back([&, mode]
        {
            if (std::fesetround(mode) != 0) return;
            for (unsigned i = 0; i < 256; ++i)
            {
                double output = 17;
                int scalar = 17, pair = 17;
                counting = true; allocations = 0;
                bool const success = ec::try_round_intersection(v, i % 2, output);
                bool const scalar_success = ec::try_compare_intersection_coordinate(v, i % 2, 1.0, scalar);
                bool const pair_success = ec::try_compare_intersection_coordinates(v, v, i % 2, pair);
                counting = false;
                if (success != (mode == FE_TONEAREST) || allocations != 0
                    || scalar_success != success || pair_success != success
                    || output != (success ? 1.5 : 17)
                    || scalar != (success ? 1 : 17) || pair != (success ? 0 : 17)) correct.store(false);
            }
        });
    }
    for (auto& worker : workers) worker.join();
    BOOST_TEST(correct.load());
    BOOST_TEST_EQ(std::fegetround(), original);
}

int main()
{
    int const original = std::fegetround();
    if (std::fesetround(FE_TONEAREST) != 0) return 1;
    ordinary();
    ties();
    cancellation_and_guards();
    randomized();
    coordinate_comparisons();
    comparison_guards();
    randomized_comparisons();
    active_expansion_storage();
    modes_and_threads();
    std::fesetround(original);
    return boost::report_errors();
}
