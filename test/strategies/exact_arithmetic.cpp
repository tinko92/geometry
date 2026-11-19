// Boost.Geometry
// Unit Test
//
// Copyright (c) 2026 Tinko Bartels, Berlin, Germany.
//
// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include <boost/geometry/strategies/cartesian/exact/arithmetic.hpp>

#include <boost/core/lightweight_test.hpp>

#include <algorithm>
#include <array>
#include <cfenv>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <random>
#include <type_traits>
#include <utility>

namespace ec = boost::geometry::strategies::detail::exact_cartesian;
using boost::multiprecision::cpp_int;

template <typename T, typename = void>
struct has_exact_value : std::false_type {};

template <typename T>
struct has_exact_value<T, decltype(void(ec::exact_value(std::declval<T>())))>
    : std::true_type {};

// Independent binary-value oracle: decode the IEC 60559 double fields rather
// than using the Multiprecision floating conversion under test.
ec::exact_number binary_oracle(double value)
{
    static_assert(sizeof(double) == sizeof(std::uint64_t)
                  && std::numeric_limits<double>::digits == 53, "Binary64 test");
    if (value == 0)
    {
        return ec::exact_number(0);
    }
    std::uint64_t bits;
    std::memcpy(&bits, &value, sizeof(bits));
    unsigned const exponent = static_cast<unsigned>((bits >> 52) & 0x7ff);
    std::uint64_t const fraction = bits & ((std::uint64_t(1) << 52) - 1);
    cpp_int numerator = exponent == 0 ? fraction : fraction | (std::uint64_t(1) << 52);
    if ((bits >> 63) != 0)
    {
        numerator = -numerator;
    }
    int const shift = exponent == 0 ? -1074 : static_cast<int>(exponent) - 1023 - 52;
    ec::exact_number result(numerator);
    if (shift >= 0)
    {
        result *= cpp_int(1) << shift;
    }
    else
    {
        result /= cpp_int(1) << -shift;
    }
    return result;
}

void test_exact_value()
{
    static_assert(has_exact_value<float>::value, "float supported");
    static_assert(has_exact_value<double>::value, "double supported");
    static_assert(has_exact_value<std::uint64_t>::value, "integers supported");
    static_assert(!has_exact_value<long double>::value, "No long double promise");

    BOOST_TEST(ec::exact_value((std::numeric_limits<std::int64_t>::min)())
               == -(ec::exact_number(cpp_int(1) << 63)));
    BOOST_TEST(ec::exact_value((std::numeric_limits<std::uint64_t>::max)())
               == ec::exact_number((cpp_int(1) << 64) - 1));
    BOOST_TEST(ec::exact_value(std::uint64_t(9007199254740993ULL))
               == ec::exact_number((cpp_int(1) << 53) + 1));
    BOOST_TEST(ec::exact_value(true) == 1);
    BOOST_TEST(ec::exact_value(false) == 0);

    double const p53 = std::ldexp(1.0, 53);
    std::array<double, 16> const values = {{
        0.0, -0.0, 0.1, -0.1, p53, p53 - 1,
        std::nextafter(p53, (std::numeric_limits<double>::infinity)()),
        std::nextafter(p53, 0.0),
        (std::numeric_limits<double>::denorm_min)(),
        -(std::numeric_limits<double>::denorm_min)(),
        (std::numeric_limits<double>::min)(),
        std::nextafter((std::numeric_limits<double>::min)(), 0.0),
        (std::numeric_limits<double>::max)(),
        -(std::numeric_limits<double>::max)(), 1.0 / 3.0, -17.25
    }};
    for (double value : values)
    {
        BOOST_TEST(ec::exact_value(value) == binary_oracle(value));
    }
    BOOST_TEST(ec::exact_value(0.1) != ec::exact_number(1) / 10);

    std::array<float, 7> const floats = {{
        0.0f, -0.1f, (std::numeric_limits<float>::denorm_min)(),
        (std::numeric_limits<float>::min)(), (std::numeric_limits<float>::max)(),
        std::ldexp(1.0f, 24), std::nextafter(std::ldexp(1.0f, 24), 0.0f)
    }};
    for (float value : floats)
    {
        BOOST_TEST(ec::exact_value(value) == binary_oracle(static_cast<double>(value)));
    }
    BOOST_TEST_THROWS(ec::exact_value((std::numeric_limits<double>::infinity)()),
                      std::domain_error);
    BOOST_TEST_THROWS(ec::exact_value(-(std::numeric_limits<double>::infinity)()),
                      std::domain_error);
    BOOST_TEST_THROWS(ec::exact_value((std::numeric_limits<double>::quiet_NaN)()),
                      std::domain_error);
    BOOST_TEST_THROWS(ec::exact_value((std::numeric_limits<float>::infinity)()),
                      std::domain_error);
}

template <typename T>
int orientation_oracle(T ax, T ay, T bx, T by, T cx, T cy)
{
    cpp_int const determinant = (cpp_int(bx) - cpp_int(ax)) * (cpp_int(cy) - cpp_int(ay))
                              - (cpp_int(by) - cpp_int(ay)) * (cpp_int(cx) - cpp_int(ax));
    return determinant > 0 ? 1 : determinant < 0 ? -1 : 0;
}

template <typename T>
void check_orientation(T ax, T ay, T bx, T by, T cx, T cy)
{
    int const expected = orientation_oracle(ax, ay, bx, by, cx, cy);
    BOOST_TEST_EQ(ec::integer_orientation(ax, ay, bx, by, cx, cy), expected);
    BOOST_TEST_EQ(ec::integer_orientation(bx, by, ax, ay, cx, cy), -expected);
    BOOST_TEST_EQ(ec::integer_orientation(cx, cy, ax, ay, bx, by), expected);
}

template <typename T>
void test_integer_type(std::mt19937_64& random)
{
    T const lo = (std::numeric_limits<T>::min)(), hi = (std::numeric_limits<T>::max)();
    T const zero = 0, one = 1;
    BOOST_TEST_EQ(ec::integer_orientation(lo, lo, hi, lo, lo, hi), 1);
    check_orientation(lo, lo, hi, hi, hi, lo);
    check_orientation(lo, hi, hi, lo, lo, lo);
    check_orientation(lo, lo, hi, hi, zero, zero);
    check_orientation(lo, hi, lo, hi, hi, lo);
    check_orientation(hi, hi, T(hi - one), hi, hi, T(hi - one));
    check_orientation(lo, lo, T(lo + one), lo, lo, T(lo + one));

    // Defined arithmetic conversion of random unsigned bits to the full signed
    // range, without an out-of-range unsigned-to-signed cast.
    std::uniform_int_distribution<T> distribution(lo, hi);
    for (unsigned i = 0; i < 1000; ++i)
    {
        std::array<T, 6> v;
        for (T& coordinate : v)
        {
            coordinate = distribution(random);
        }
        check_orientation(v[0], v[1], v[2], v[3], v[4], v[5]);
    }
}

void check_contains(ec::interval const& interval, ec::exact_number const& value)
{
    BOOST_TEST(!std::isnan(interval.lower()));
    BOOST_TEST(!std::isnan(interval.upper()));
    if (interval.valid())
    {
        ec::exact_number const center = binary_oracle(interval.approximate);
        ec::exact_number const error = binary_oracle(interval.error);
        BOOST_TEST(center - error <= value);
        BOOST_TEST(center + error >= value);
    }
    if (std::isfinite(interval.lower()))
    {
        BOOST_TEST(binary_oracle(interval.lower()) <= value);
    }
    else
    {
        BOOST_TEST(interval.lower() < 0);
    }
    if (std::isfinite(interval.upper()))
    {
        BOOST_TEST(binary_oracle(interval.upper()) >= value);
    }
    else
    {
        BOOST_TEST(interval.upper() > 0);
    }
    ec::filter_result const result = ec::sign(interval);
    if (result.certain)
    {
        BOOST_TEST_EQ(result.value, value > 0 ? 1 : value < 0 ? -1 : 0);
    }
}

void check_operation_bounds(ec::interval const& a, ec::interval const& b)
{
    ec::exact_number const ac = binary_oracle(a.approximate), ae = binary_oracle(a.error);
    ec::exact_number const bc = binary_oracle(b.approximate), be = binary_oracle(b.error);
    std::array<ec::exact_number, 3> const av = {{ac - ae, ac, ac + ae}};
    std::array<ec::exact_number, 3> const bv = {{bc - be, bc, bc + be}};
    ec::interval const sum = a + b, difference = a - b, product = a * b, quotient = a / b;
    for (auto const& x : av)
    {
        check_contains(a, x);
        for (auto const& y : bv)
        {
            check_contains(sum, x + y);
            check_contains(difference, x - y);
            check_contains(product, x * y);
            if (y != 0)
            {
                check_contains(quotient, x / y);
            }
            ec::filter_result const order = ec::compare(a, b);
            if (order.certain)
            {
                BOOST_TEST_EQ(order.value, x > y ? 1 : x < y ? -1 : 0);
            }
        }
    }
    if (bv[0] <= 0 && bv[2] >= 0)
    {
        BOOST_TEST(!ec::sign(quotient).certain);
        BOOST_TEST(!quotient.valid());
    }
}

void test_intervals()
{
    ec::interval const zero{0, 0}, one{1, 0};
    BOOST_TEST(ec::sign(zero).certain);
    BOOST_TEST_EQ(ec::sign(zero).value, 0);
    BOOST_TEST(ec::sign(one).certain);
    BOOST_TEST_EQ(ec::sign(one).value, 1);
    BOOST_TEST_EQ(ec::sign(ec::interval{-1, 0}).value, -1);
    BOOST_TEST(ec::compare(one, one).certain);
    BOOST_TEST_EQ(ec::compare(one, one).value, 0);
    BOOST_TEST(!ec::sign(ec::interval{0, 1}).certain);
    BOOST_TEST(!ec::compare(ec::interval{1, 1}, one).certain);
    BOOST_TEST(ec::compare(zero, one).certain);
    BOOST_TEST_EQ(ec::compare(zero, one).value, -1);
    BOOST_TEST_EQ(ec::compare(one, zero).value, 1);

    ec::interval const cancel = (ec::interval{1e16, 0} + one) - ec::interval{1e16, 0};
    check_contains(cancel, 1);
    BOOST_TEST(!ec::sign(cancel).certain);
    check_contains(one - one, 0);
    BOOST_TEST(ec::sign(one / ec::interval{3, 0}).certain);
    BOOST_TEST_EQ(ec::sign(one / ec::interval{3, 0}).value, 1);
    BOOST_TEST(ec::sign(one / ec::interval{-3, 0}).certain);
    BOOST_TEST_EQ(ec::sign(one / ec::interval{-3, 0}).value, -1);

    double const tiny = (std::numeric_limits<double>::denorm_min)();
    double const huge = (std::numeric_limits<double>::max)();
    check_operation_bounds(ec::interval{tiny, 0}, ec::interval{0.5, 0});
    check_operation_bounds(ec::interval{tiny, tiny}, ec::interval{-tiny, tiny});
    check_operation_bounds(ec::interval{1, tiny}, ec::interval{-1, tiny});
    check_operation_bounds(ec::interval{tiny, 0}, ec::interval{tiny, 0});
    check_operation_bounds(one, ec::interval{tiny, 0});
    check_operation_bounds(ec::interval{huge, 0}, ec::interval{huge, 0});
    check_operation_bounds(ec::interval{huge, huge}, ec::interval{-huge, 0});
    check_operation_bounds(ec::interval{1, 0.25}, ec::interval{2, 0.125});
    check_operation_bounds(ec::interval{-1, 0.25}, ec::interval{-2, 0.125});
    check_operation_bounds(ec::interval{0, 0}, ec::interval{0, 1});
    check_operation_bounds(one, ec::interval{0, 0});
    check_operation_bounds(one, ec::interval{0, 1});
    check_operation_bounds(one, ec::interval{1, 1});
    check_operation_bounds(one, ec::interval{-1, 1});
    BOOST_TEST(!ec::sign(one / ec::interval{0, 1}).certain);
    BOOST_TEST(!ec::sign(one / zero).certain);
    BOOST_TEST(!ec::sign(ec::interval{huge, 0} * ec::interval{2, 0}).certain);

    std::array<ec::interval, 5> const invalid = {{
        ec::interval::uncertain(), {1, -1},
        {(std::numeric_limits<double>::infinity)(), 0},
        {(std::numeric_limits<double>::quiet_NaN)(), 0},
        {0, (std::numeric_limits<double>::quiet_NaN)()}
    }};
    for (auto const& value : invalid)
    {
        BOOST_TEST(!value.valid());
        BOOST_TEST(!ec::sign(value).certain);
        BOOST_TEST(!ec::compare(value, one).certain);
        BOOST_TEST(!ec::compare(one, value).certain);
        BOOST_TEST(!ec::sign(value + one).certain);
        BOOST_TEST(!ec::sign(one - value).certain);
        BOOST_TEST(!ec::sign(value * zero).certain);
        BOOST_TEST(!ec::sign(one / value).certain);
        BOOST_TEST(!ec::sign(value / one).certain);
    }
}

double random_finite(std::mt19937_64& random)
{
    // Sample all finite exponents and significands, not only a small real range.
    std::uint64_t const bits = (random() & ~(std::uint64_t(0x7ff) << 52))
                             | ((random() % 2047) << 52);
    double value;
    std::memcpy(&value, &bits, sizeof(value));
    return value;
}

void test_random_intervals(std::mt19937_64& random)
{
    for (unsigned i = 0; i < 400; ++i)
    {
        double const a = random_finite(random), b = random_finite(random);
        BOOST_TEST(ec::exact_value(a) == binary_oracle(a));
        BOOST_TEST(ec::exact_value(b) == binary_oracle(b));
        check_operation_bounds(ec::interval{a, 0}, ec::interval{b, 0});
        check_operation_bounds(ec::interval{a, std::abs(random_finite(random))},
                               ec::interval{b, std::abs(random_finite(random))});
    }

    std::uniform_int_distribution<int> mantissa(-1000000, 1000000), exponent(-200, 200);
    for (unsigned i = 0; i < 400; ++i)
    {
        ec::interval const a{std::ldexp(double(mantissa(random)), exponent(random)),
                             std::ldexp(1.0, exponent(random))};
        ec::interval const b{std::ldexp(double(mantissa(random)), exponent(random)),
                             std::ldexp(1.0, exponent(random))};
        check_operation_bounds(a, b);
        // Verify a chain, including the original uncertainty, against rationals.
        ec::exact_number const x = binary_oracle(a.approximate) + binary_oracle(a.error);
        ec::exact_number const y = binary_oracle(b.approximate) - binary_oracle(b.error);
        check_contains((a + b) * (a - b), (x + y) * (x - y));
    }
}

int main()
{
    test_exact_value();
    std::mt19937_64 random(0x4a6f7921);
    test_integer_type<std::int16_t>(random);
    test_integer_type<std::uint16_t>(random);
    test_integer_type<std::int32_t>(random);
    test_integer_type<std::uint32_t>(random);
    test_integer_type<std::int64_t>(random);
    test_integer_type<std::uint64_t>(random);
    check_orientation<std::int64_t>(0, 0, 3037000500LL, 3037000501LL,
                                    3037000499LL, 3037000500LL);
    check_orientation<std::int64_t>(0, 0, 9007199254740992LL, 9007199254740993LL,
                                    9007199254740991LL, 9007199254740992LL);
    check_orientation<std::int64_t>(0, 0, 9223372036854775806LL, 9223372036854775807LL,
                                    9223372036854775805LL, 9223372036854775806LL);
    check_orientation<std::int32_t>(0, 0, 2147483646, 2147483647,
                                    2147483645, 2147483646);
    check_orientation<std::uint64_t>(0, 0, 18446744073709551614ULL, 18446744073709551615ULL,
                                     18446744073709551613ULL, 18446744073709551614ULL);

    int const original_rounding = std::fegetround();
    for (int mode : {FE_TONEAREST, FE_DOWNWARD, FE_UPWARD, FE_TOWARDZERO})
    {
        if (std::fesetround(mode) == 0)
        {
            test_intervals();
            test_random_intervals(random);
        }
    }
    BOOST_TEST_EQ(std::fesetround(original_rounding), 0);
    return boost::report_errors();
}
