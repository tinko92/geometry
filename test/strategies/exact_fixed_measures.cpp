// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.

#include <boost/core/lightweight_test.hpp>
#include <boost/geometry/strategies/cartesian/exact/fixed_measures.hpp>
#include <boost/geometry/strategies/cartesian/exact/arithmetic.hpp>

#include <array>
#include <cfenv>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <new>
#include <random>
#include <type_traits>
#include <utility>

namespace ec = boost::geometry::strategies::detail::exact_cartesian;
using boost::multiprecision::cpp_int;
using rational = ec::exact_number;

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

template <typename Function>
auto without_allocations(Function function) -> decltype(function())
{
    decltype(function()) result{};
    {
        allocation_scope scope;
        result = function();
    }
    BOOST_TEST_EQ(allocations, 0u);
    return result;
}

template <std::size_t N, std::size_t D>
std::array<rational, 2> coordinates(ec::fixed_homogeneous<N, D> const& point)
{
    return {{rational(cpp_int(point.X)) / rational(cpp_int(point.W)),
             rational(cpp_int(point.Y)) / rational(cpp_int(point.W))}};
}

template <std::size_t N, std::size_t D>
rational value(ec::fixed_rational_pair<N, D> const& pair)
{
    BOOST_TEST(pair.D > 0);
    return rational(cpp_int(pair.N)) / rational(cpp_int(pair.D));
}

template <std::size_t N, std::size_t D, typename Coordinate>
void check_native(ec::fixed_rational_pair<N, D> const& pair, Coordinate native)
{
    auto const expected = value(pair);
    auto const scalar = ec::exact_value(native);
    BOOST_TEST_EQ(without_allocations([&] { return ec::fixed_pair_compare(pair, native); }),
                  (expected > scalar) - (expected < scalar));
}

template <typename A, typename B, typename P>
void check_measures(A const& a, B const& b, P const& point)
{
    auto const av = coordinates(a), bv = coordinates(b), pv = coordinates(point);
    rational const dx = bv[0] - av[0], dy = bv[1] - av[1];
    rational const expected_distance = dx * dx + dy * dy;
    rational const expected_determinant = dx * (pv[1] - av[1]) - dy * (pv[0] - av[0]);
    auto const distance = without_allocations([&] { return ec::fixed_squared_distance(a, b); });
    auto const determinant = without_allocations([&] { return ec::fixed_signed_determinant(a, b, point); });
    BOOST_TEST(value(distance) == expected_distance);
    BOOST_TEST(value(determinant) == expected_determinant);
    BOOST_TEST_EQ(without_allocations([&] { return ec::fixed_pair_compare(distance, determinant); }),
                  (expected_distance > expected_determinant) - (expected_distance < expected_determinant));
    auto reverse = without_allocations([&] { return ec::fixed_signed_determinant(b, a, point); });
    reverse.N = -reverse.N;
    BOOST_TEST_EQ(without_allocations([&] { return ec::fixed_pair_compare(determinant, reverse); }), 0);
    auto const symmetric = without_allocations([&] { return ec::fixed_squared_distance(b, a); });
    BOOST_TEST_EQ(without_allocations([&] { return ec::fixed_pair_compare(distance, symmetric); }), 0);
    auto const zero = without_allocations([&] { return ec::fixed_signed_determinant(a, b, a); });
    BOOST_TEST(zero.N == 0);
    check_native(determinant, 0);
    check_native(distance, 1);
    check_native(determinant, -1);
}

template <typename Coordinate>
void test_native_type()
{
    Coordinate const low = (std::numeric_limits<Coordinate>::lowest)();
    Coordinate const high = (std::numeric_limits<Coordinate>::max)();
    auto const a = ec::fixed_input_homogeneous(low, low);
    auto const b = ec::fixed_input_homogeneous(high, low);
    auto const c = ec::fixed_input_homogeneous(low, high);
    check_measures(a, b, c);
    check_measures(a, a, a);
    check_measures(c, b, a);
}

template <typename Coordinate>
void test_constructed()
{
    auto const a = ec::direct_native_homogeneous<Coordinate>({{0, 0, 1, 1, 0, 1, 1, -1}});
    auto const b = ec::direct_native_homogeneous<Coordinate>({{1, 0, 2, 1, 1, 1, 2, -1}});
    auto const c = ec::direct_native_homogeneous<Coordinate>({{0, 1, 1, 2, 0, 2, 1, 0}});
    check_measures(a, b, c);
    auto const leaf = ec::fixed_input_homogeneous(Coordinate(0), Coordinate(0));
    check_measures(leaf, b, c);
    check_measures(a, leaf, c);
    check_measures(a, b, leaf);
    auto const distance = ec::fixed_squared_distance(a, b);
    check_native(distance, (std::numeric_limits<double>::denorm_min)());
    check_native(distance, (std::numeric_limits<double>::max)());
    check_native(distance, (std::numeric_limits<std::int64_t>::min)());
    check_native(distance, (std::numeric_limits<std::uint64_t>::max)());
}

template <std::size_t Bits>
ec::fixed_integer<Bits> maximum_fixed()
{
    ec::fixed_integer<Bits> value(1);
    value <<= Bits - 1;
    --value;
    value += value;
    ++value;
    return value;
}

void test_capacities()
{
    using point = ec::fixed_homogeneous<195, 131>;
    using traits = ec::fixed_measure_traits<195, 131>;
    using pair = typename traits::pair_type;
    static_assert(std::is_same<decltype(ec::fixed_squared_distance(std::declval<point>(), std::declval<point>())),
                               pair>::value, "Same-format distance capacity");
    static_assert(std::is_same<decltype(ec::fixed_signed_determinant(std::declval<point>(), std::declval<point>(),
                                                                   std::declval<point>())), pair>::value,
                  "Same-format determinant capacity");
    static_assert(std::numeric_limits<pair::numerator_type>::digits == 655, "2N+2D+3 numerator");
    static_assert(std::numeric_limits<pair::denominator_type>::digits == 524, "4D denominator");
    static_assert(std::numeric_limits<traits::word_type>::digits == 1184, "2N+6D+8 word");
    point a, b, c;
    a.X = a.Y = maximum_fixed<195>();
    a.W = maximum_fixed<131>();
    b = a;
    b.X = -b.X;
    b.Y = -b.Y;
    c = a;
    c.X = -c.X;
    check_measures(a, b, c);
    pair positive, negative;
    positive.N = maximum_fixed<655>();
    positive.D = maximum_fixed<524>();
    negative = positive;
    negative.N = -negative.N;
    BOOST_TEST_EQ(without_allocations([&] { return ec::fixed_pair_compare(positive, negative); }), 1);
    BOOST_TEST_EQ(without_allocations([&] { return ec::fixed_pair_compare(negative, positive); }), -1);
    // Include both native numerator and native grid-denominator width. The
    // latter dominates for an unusually wide numerator with a tiny denominator.
    ec::fixed_rational_pair<2005, 4> large;
    large.N = 1;
    large.N <<= 2000;
    check_native(large, (std::numeric_limits<double>::denorm_min)());
    ec::fixed_homogeneous<3000, 2500> wide_a, wide_b, wide_c;
    wide_a.X = wide_a.Y = maximum_fixed<3000>();
    wide_a.W = maximum_fixed<2500>();
    wide_b = wide_a;
    wide_b.X = -wide_b.X;
    wide_b.Y = -wide_b.Y;
    wide_c = wide_a;
    wide_c.X = -wide_c.X;
    check_measures(wide_a, wide_b, wide_c);
}

void test_random()
{
    std::mt19937_64 generator(0x43c171a9);
    std::uniform_int_distribution<std::int64_t> random_integer(
        (std::numeric_limits<std::int64_t>::min)(), (std::numeric_limits<std::int64_t>::max)());
    for (int i = 0; i < 32; ++i)
    {
        auto const a = ec::fixed_input_homogeneous(random_integer(generator), random_integer(generator));
        auto const b = ec::fixed_input_homogeneous(random_integer(generator), random_integer(generator));
        auto const c = ec::fixed_input_homogeneous(random_integer(generator), random_integer(generator));
        check_measures(a, b, c);
        std::array<double, 6> values;
        for (auto& value : values)
        {
            double const mantissa = double(generator() & ((std::uint64_t(1) << 53) - 1));
            value = std::ldexp(mantissa, int(generator() % 120) - 100);
            if (generator() & 1) { value = -value; }
        }
        auto const x = ec::fixed_input_homogeneous(values[0], values[1]);
        auto const y = ec::fixed_input_homogeneous(values[2], values[3]);
        auto const z = ec::fixed_input_homogeneous(values[4], values[5]);
        check_measures(x, y, z);
    }
}

void test_invalid()
{
    ec::fixed_homogeneous<64> invalid, valid;
    invalid.W = 0;
    BOOST_TEST_THROWS(ec::fixed_squared_distance(invalid, valid), std::domain_error);
    BOOST_TEST_THROWS(ec::fixed_signed_determinant(valid, invalid, valid), std::domain_error);
    ec::fixed_rational_pair<64, 64> pair;
    pair.D = -1;
    BOOST_TEST_THROWS(ec::fixed_pair_compare(pair, pair), std::domain_error);
}

int main()
{
    test_native_type<std::int64_t>();
    test_native_type<std::uint64_t>();
    test_capacities();
    test_invalid();
    int const original = std::fegetround();
    for (int mode : {FE_TONEAREST, FE_UPWARD, FE_DOWNWARD, FE_TOWARDZERO})
    {
        BOOST_TEST_EQ(std::fesetround(mode), 0);
        test_native_type<float>();
        test_native_type<double>();
        test_constructed<double>();
        test_constructed<std::int64_t>();
    }
    BOOST_TEST_EQ(std::fesetround(original), 0);
    test_random();
    return boost::report_errors();
}
