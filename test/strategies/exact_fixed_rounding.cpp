// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.

#include <boost/core/lightweight_test.hpp>
#include <boost/geometry/strategies/cartesian/exact/fixed_rounding.hpp>
#include <boost/geometry/strategies/cartesian/exact/point.hpp>

#include <array>
#include <cfenv>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <new>
#include <random>
#include <stdexcept>
#include <type_traits>

namespace ec = boost::geometry::strategies::detail::exact_cartesian;
using rational = ec::exact_number;
using boost::multiprecision::cpp_int;

namespace
{
bool counting = false;
std::size_t allocations = 0;
struct allocation_scope
{
    allocation_scope() { allocations = 0; counting = true; }
    ~allocation_scope() { counting = false; }
};
}

void* operator new(std::size_t size)
{
    if (counting) { ++allocations; }
    if (void* result = std::malloc(size ? size : 1)) { return result; }
    throw std::bad_alloc();
}
void operator delete(void* value) noexcept { std::free(value); }
void* operator new[](std::size_t size) { return ::operator new(size); }
void operator delete[](void* value) noexcept { ::operator delete(value); }
#if defined(__cpp_sized_deallocation)
void operator delete(void* value, std::size_t) noexcept { ::operator delete(value); }
void operator delete[](void* value, std::size_t) noexcept { ::operator delete(value); }
#endif

template <typename Output>
void compare_output(Output actual, Output expected, std::false_type)
{
    BOOST_TEST_EQ(actual, expected);
}

template <typename Output>
void compare_output(Output actual, Output expected, std::true_type)
{
    BOOST_TEST_EQ(actual, expected);
    BOOST_TEST_EQ(std::signbit(actual), std::signbit(expected));
}

template <typename Output, std::size_t N, std::size_t D>
Output check_round(ec::fixed_integer<N> const& numerator, ec::fixed_integer<D> const& denominator)
{
    rational const value = rational(cpp_int(numerator)) / rational(cpp_int(denominator));
    Output expected{};
    bool overflow = false;
    try { expected = ec::point_detail::output_rounder<Output>::apply(value); }
    catch (std::overflow_error const&) { overflow = true; }
    Output actual = Output(17);
    if (overflow)
    {
        BOOST_TEST_THROWS(ec::try_fixed_round<Output>(numerator, denominator, actual), std::overflow_error);
        BOOST_TEST_EQ(actual, Output(17));
        return actual;
    }
    bool supported;
    {
        allocation_scope scope;
        supported = ec::try_fixed_round<Output>(numerator, denominator, actual);
    }
    BOOST_TEST(supported);
    BOOST_TEST_EQ(allocations, 0u);
    compare_output(actual, expected, std::is_floating_point<Output>());
    return actual;
}

template <typename Output>
Output check_round(rational const& value)
{
    // The oracle/input conversion is deliberately outside the allocation guard.
    ec::fixed_integer<4096> const numerator(boost::multiprecision::numerator(value));
    ec::fixed_integer<4096> const denominator(boost::multiprecision::denominator(value));
    return check_round<Output>(numerator, denominator);
}

template <typename Output>
void test_integer_rounding()
{
    for (int numerator = -25; numerator <= 25; ++numerator)
    for (int denominator = 1; denominator <= 7; ++denominator)
    {
        check_round<Output>(rational(numerator) / denominator);
    }
    rational const lo((std::numeric_limits<Output>::min)());
    rational const hi((std::numeric_limits<Output>::max)());
    for (rational const& value : std::array<rational, 12>{{lo, hi, lo - 1, hi + 1,
        lo - rational(1) / 2, lo - rational(1) / 4, lo + rational(1) / 2,
        hi - rational(1) / 2, hi + rational(1) / 4, hi + rational(1) / 2,
        lo - rational(3) / 4, hi + rational(3) / 4}})
    {
        check_round<Output>(value);
    }
}

template <typename Output>
void test_binary_boundaries()
{
    constexpr int digits = std::numeric_limits<Output>::digits;
    rational const ulp = rational(1) / (cpp_int(1) << (digits - 1));
    rational const half = ulp / 2, small = ulp / 1024;
    Output const next = std::nextafter(Output(1), Output(2));
    Output const next_next = std::nextafter(next, Output(2));
    BOOST_TEST_EQ(check_round<Output>(1 + half), Output(1));
    BOOST_TEST_EQ(check_round<Output>(1 + half + small), next);
    BOOST_TEST_EQ(check_round<Output>(1 + half - small), Output(1));
    BOOST_TEST_EQ(check_round<Output>(1 + ulp + half), next_next);
    BOOST_TEST_EQ(check_round<Output>(-1 - half), Output(-1));
    BOOST_TEST_EQ(check_round<Output>(-1 - ulp - half), -next_next);

    Output const tiny = (std::numeric_limits<Output>::denorm_min)();
    rational const exact_tiny = ec::exact_value(tiny);
    BOOST_TEST_EQ(check_round<Output>(exact_tiny / 2), Output(0));
    BOOST_TEST(std::signbit(check_round<Output>(-exact_tiny / 2)));
    BOOST_TEST_EQ(check_round<Output>(exact_tiny * 3 / 2), 2 * tiny);
    BOOST_TEST_EQ(check_round<Output>(exact_tiny * 7 / 4), 2 * tiny);
    BOOST_TEST_EQ(check_round<Output>(exact_tiny / 4), Output(0));
    BOOST_TEST(std::signbit(check_round<Output>(-exact_tiny / 4)));
    rational const normal = ec::exact_value((std::numeric_limits<Output>::min)());
    BOOST_TEST_EQ(check_round<Output>(normal - exact_tiny / 2), (std::numeric_limits<Output>::min)());
    BOOST_TEST_EQ(check_round<Output>(normal - exact_tiny * 3 / 4),
        std::nextafter((std::numeric_limits<Output>::min)(), Output(0)));
    check_round<Output>(normal + exact_tiny / 2);

    Output const maximum = (std::numeric_limits<Output>::max)();
    rational const maximum_ulp(cpp_int(1) << (std::numeric_limits<Output>::max_exponent - digits));
    rational const overflow_tie = ec::exact_value(maximum) + maximum_ulp / 2;
    BOOST_TEST_EQ(check_round<Output>(overflow_tie - maximum_ulp / 16), maximum);
    BOOST_TEST(std::isinf(check_round<Output>(overflow_tie)));
    BOOST_TEST(std::isinf(check_round<Output>(-overflow_tie)));
    check_round<Output>(overflow_tie + maximum_ulp / 16);
    check_round<Output>(rational(cpp_int(1) << (std::numeric_limits<Output>::max_exponent + 100)));
    check_round<Output>(-rational(cpp_int(1) << (std::numeric_limits<Output>::max_exponent + 100)));

    // Exercise both exponent-adjustment branches and rounding carry into 2^e.
    for (int exponent : {-150, -25, -1, 0, 1, 24, 52, 127})
    {
        rational const scale = exponent < 0 ? rational(1) / rational(cpp_int(1) << -exponent)
                                           : rational(cpp_int(1) << exponent);
        check_round<Output>(scale * (1 - half));
        check_round<Output>(scale * (1 + half));
        check_round<Output>(scale * rational(7) / 9);
    }
}

void test_homogeneous_and_contract()
{
    ec::fixed_homogeneous<128, 64> const point{ec::fixed_integer<128>((std::uint64_t(1) << 63) + 1),
                                            ec::fixed_integer<128>(-3), ec::fixed_integer<64>(2)};
    std::int64_t rounded;
    BOOST_TEST(ec::try_fixed_round<std::int64_t>(point, 0, rounded));
    BOOST_TEST_EQ(rounded, std::int64_t(1) << 62);
    BOOST_TEST(ec::try_fixed_round<std::int64_t>(point, 1, rounded));
    BOOST_TEST_EQ(rounded, -2);
    BOOST_TEST_THROWS(ec::try_fixed_round<std::int64_t>(point, 2, rounded), std::out_of_range);
    BOOST_TEST_EQ(rounded, -2);

    ec::fixed_integer<1> const one(1);
    ec::fixed_integer<1> const zero(0), negative(-1);
    double output = 17;
    BOOST_TEST(!ec::try_fixed_round<double>(one, zero, output));
    BOOST_TEST(!ec::try_fixed_round<double>(one, negative, output));
    BOOST_TEST_EQ(output, 17);
    long double unsupported = 23;
    BOOST_TEST(!ec::try_fixed_round<long double>(one, one, unsupported));
    BOOST_TEST_EQ(unsupported, 23);
    rational unchanged = 29;
    BOOST_TEST(!ec::try_fixed_round<rational>(one, one, unchanged));
    BOOST_TEST_EQ(unchanged, 29);
    check_round<double>(zero, one);
    check_round<float>(one, one);
    check_round<std::int8_t>(negative, one);

    // Saturate a word and its denominator, including division's multi-limb path.
    ec::fixed_integer<5000> const wide = (ec::fixed_integer<5000>(1) << 4999) - 1;
    ec::fixed_integer<5000> const other = (ec::fixed_integer<5000>(1) << 4998) + 3;
    check_round<double>(wide, other);
    check_round<float>(-wide, other);
    check_round<std::int64_t>(wide, other);
    check_round<double>(wide, ec::fixed_integer<1>(1));
    check_round<double>(ec::fixed_integer<1>(-1), wide);
    // The actual upper bit, rather than just a power of two minus one.
    ec::fixed_integer<5000> const full = (wide << 1) + 1;
    check_round<double>(full, other);

    std::uint64_t const large = std::uint64_t(1) << 53;
    BOOST_TEST_EQ(check_round<double>(ec::fixed_integer<64>(large + 1), one), double(large));
    BOOST_TEST_EQ(check_round<double>(ec::fixed_integer<64>(large + 3), one), double(large + 4));
    check_round<float>(ec::fixed_integer<64>((std::numeric_limits<std::uint64_t>::max)()), one);
    std::int64_t const base = std::int64_t(large);
    auto const fractional = ec::direct_native_homogeneous<std::int64_t>
        ({{base, 0, base + 4, 0, base + 1, -1, base + 2, 2}});
    BOOST_TEST_EQ(check_round<double>(fractional.X, fractional.W), double(large + 2));
    BOOST_TEST_EQ(check_round<std::int64_t>(fractional.X, fractional.W), base + 1);
    check_round<float>(fractional.X, fractional.W);

    // Also accept the common interchange record's integer format directly.
    ec::fixed_integer<32768> const scalar_n = (ec::fixed_integer<32768>(1) << 32767) + 17;
    ec::fixed_integer<32768> const scalar_d = (ec::fixed_integer<32768>(1) << 32766) - 1;
    check_round<double>(scalar_n, scalar_d);
    check_round<float>(-scalar_n, scalar_d);
}

void test_randomized()
{
    std::mt19937_64 random(0x46524f554e44494e);
    for (unsigned i = 0; i < 128; ++i)
    {
        ec::fixed_integer<3200> numerator(random()), denominator(random() | 1);
        for (unsigned limb = 0; limb < i % 16; ++limb)
        {
            numerator <<= 64;
            numerator += random();
            denominator <<= 64;
            denominator += random();
        }
        numerator <<= random() % 1100;
        denominator <<= random() % 1100;
        if (random() & 1) { numerator = -numerator; }
        check_round<double>(numerator, denominator);
        check_round<float>(numerator, denominator);
        check_round<std::int64_t>(numerator, denominator);
        check_round<std::uint64_t>(numerator, denominator);
    }
    for (unsigned i = 0; i < 128; ++i)
    {
        std::uint64_t bits = random();
        if (((bits >> 52) & 2047) == 2047) { bits ^= std::uint64_t(1) << 52; }
        double input;
        std::memcpy(&input, &bits, sizeof(input));
        auto const value = ec::fixed_input_homogeneous(input, input);
        check_round<double>(value.X, value.W);
        check_round<float>(value.X, value.W);
        check_round<std::int64_t>(value.X, value.W);
    }
    for (unsigned i = 0; i < 16; ++i)
    {
        std::array<std::int64_t, 8> input;
        for (auto& coordinate : input) { coordinate = std::int64_t(random() % 2001) - 1000; }
        try
        {
            auto const intersection = ec::direct_native_homogeneous(input);
            check_round<double>(intersection.X, intersection.W);
            check_round<float>(intersection.Y, intersection.W);
            check_round<std::int64_t>(intersection.X, intersection.W);
        }
        catch (std::domain_error const&) {}
    }
}

int main()
{
    int const original_rounding = std::fegetround();
    for (int mode : {FE_TONEAREST, FE_DOWNWARD, FE_UPWARD, FE_TOWARDZERO})
    {
        if (std::fesetround(mode) != 0) { continue; }
        test_binary_boundaries<float>();
        test_binary_boundaries<double>();
        test_homogeneous_and_contract();
    }
    BOOST_TEST_EQ(std::fesetround(original_rounding), 0);
    test_integer_rounding<bool>();
    test_integer_rounding<std::int8_t>();
    test_integer_rounding<std::uint8_t>();
    test_integer_rounding<std::int16_t>();
    test_integer_rounding<std::uint16_t>();
    test_integer_rounding<std::int32_t>();
    test_integer_rounding<std::uint32_t>();
    test_integer_rounding<std::int64_t>();
    test_integer_rounding<std::uint64_t>();
    test_randomized();
    return boost::report_errors();
}
