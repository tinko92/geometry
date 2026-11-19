// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.

#include <boost/core/lightweight_test.hpp>
#include <boost/geometry/strategies/cartesian/exact/fixed_arithmetic.hpp>
#include <boost/geometry/strategies/cartesian/exact/arithmetic.hpp>

#include <array>
#include <cfenv>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
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

template <typename T, typename = void>
struct has_fixed_value : std::false_type {};
template <typename T>
struct has_fixed_value<T, decltype(void(ec::fixed_value(std::declval<T>())))> : std::true_type {};

static_assert(!has_fixed_value<rational>::value, "Arbitrary literals need a separate slow mode");
static_assert(!has_fixed_value<cpp_int>::value, "Arbitrary integers need a separate slow mode");
static_assert(!has_fixed_value<long double>::value, "Long double is not silently narrowed");
static_assert(ec::native_fixed_traits<double>::coordinate_bits == 2098, "Binary64 grid width");
static_assert(ec::native_fixed_traits<double>::intersection_numerator_bits == 6297, "Degree-three bound");
static_assert(ec::native_fixed_traits<double>::intersection_denominator_bits == 5273, "Degree-two plus grid");
static_assert(ec::native_fixed_traits<float>::coordinate_bits == 277, "Binary32 grid width");
static_assert(ec::native_fixed_traits<std::int64_t>::intersection_numerator_bits == 195, "Integer numerator");
static_assert(ec::native_fixed_traits<std::uint64_t>::intersection_denominator_bits == 131, "Integer denominator");
static_assert(std::numeric_limits<ec::fixed_integer<2098>>::digits == 2098, "No modular padding capacity");
static_assert(ec::native_fixed_polynomial<double, 8, 256>::bits == 16792, "Explicit degree-eight term bound");

template <std::size_t N, std::size_t D>
std::array<rational, 2> coordinates(ec::fixed_homogeneous<N, D> const& point)
{
    cpp_int const x(point.X), y(point.Y), w(point.W);
    BOOST_TEST(w > 0);
    cpp_int const ax = x < 0 ? cpp_int(-x) : x;
    cpp_int const ay = y < 0 ? cpp_int(-y) : y;
    BOOST_TEST(boost::multiprecision::gcd(boost::multiprecision::gcd(ax, ay), w) == 1);
    // Rational division keeps a zero numerator canonical.
    return {{rational(x) / rational(w), rational(y) / rational(w)}};
}

template <typename Coordinate>
void check_input(Coordinate x, Coordinate y)
{
    auto const value = without_allocations([&] { return ec::fixed_value(x); });
    rational scaled = ec::exact_value(x);
    scaled *= cpp_int(1) << ec::native_fixed_traits<Coordinate>::scale_bits;
    BOOST_TEST(boost::multiprecision::denominator(scaled) == 1);
    BOOST_TEST(cpp_int(value) == boost::multiprecision::numerator(scaled));
    auto const point = without_allocations([&] { return ec::fixed_input_homogeneous(x, y); });
    ec::fixed_homogeneous<128> compact;
    bool const fits = ec::fixed_fits<128>(point.X) && ec::fixed_fits<128>(point.Y) && ec::fixed_fits<128>(point.W);
    BOOST_TEST_EQ(without_allocations([&] { return ec::try_fixed_input_homogeneous(x, y, compact); }), fits);
    if (fits)
    {
        BOOST_TEST(cpp_int(point.X) == cpp_int(compact.X));
        BOOST_TEST(cpp_int(point.Y) == cpp_int(compact.Y));
        BOOST_TEST(cpp_int(point.W) == cpp_int(compact.W));
    }
    auto const actual = coordinates(point);
    BOOST_TEST(actual[0] == ec::exact_value(x));
    BOOST_TEST(actual[1] == ec::exact_value(y));
}

template <typename Coordinate>
void check_predicates(std::array<Coordinate, 6> const& values)
{
    rational const ax = ec::exact_value(values[0]), ay = ec::exact_value(values[1]);
    rational const bx = ec::exact_value(values[2]), by = ec::exact_value(values[3]);
    rational const cx = ec::exact_value(values[4]), cy = ec::exact_value(values[5]);
    rational const determinant = (bx - ax) * (cy - ay) - (by - ay) * (cx - ax);
    BOOST_TEST_EQ(without_allocations([&] { return ec::fixed_orientation(values); }), determinant.sign());
    auto const a = without_allocations([&] { return ec::fixed_input_homogeneous(values[0], values[1]); });
    auto const b = without_allocations([&] { return ec::fixed_input_homogeneous(values[2], values[3]); });
    auto const c = without_allocations([&] { return ec::fixed_input_homogeneous(values[4], values[5]); });
    BOOST_TEST_EQ(without_allocations([&] { return ec::fixed_orientation(a, b, c); }), determinant.sign());
    BOOST_TEST_EQ(without_allocations([&] { return ec::fixed_coordinate_compare(a, b, 0); }), (ax > bx) - (ax < bx));
    BOOST_TEST_EQ(without_allocations([&] { return ec::fixed_coordinate_compare(a, b, 1); }), (ay > by) - (ay < by));
    std::array<Coordinate, 4> const products{{values[0], values[1], values[2], values[3]}};
    rational const difference = ax * ay - bx * by;
    BOOST_TEST_EQ(without_allocations([&] { return ec::fixed_product_difference(products); }), difference.sign());
}

template <typename Coordinate>
void check_intersection(std::array<Coordinate, 8> values)
{
    std::array<rational, 8> q;
    for (std::size_t i = 0; i < 8; ++i) { q[i] = ec::exact_value(values[i]); }
    rational const aa = q[1] - q[3], ab = q[2] - q[0], ac = q[0] * q[3] - q[1] * q[2];
    rational const ba = q[5] - q[7], bb = q[6] - q[4], bc = q[4] * q[7] - q[5] * q[6];
    rational const w = aa * bb - ab * ba;
    if (w == 0)
    {
        BOOST_TEST_THROWS(ec::direct_native_homogeneous(values), std::domain_error);
        return;
    }
    auto const result = without_allocations([&] { return ec::direct_native_homogeneous(values); });
    auto const actual = coordinates(result);
    BOOST_TEST(actual[0] == (ab * bc - ac * bb) / w);
    BOOST_TEST(actual[1] == (ac * ba - aa * bc) / w);
    std::swap(values[0], values[2]);
    std::swap(values[1], values[3]);
    auto const reversed = without_allocations([&] { return ec::direct_native_homogeneous(values); });
    BOOST_TEST(result.X == reversed.X && result.Y == reversed.Y && result.W == reversed.W);
    auto const origin = ec::fixed_input_homogeneous(Coordinate(0), Coordinate(0));
    auto const unit = ec::fixed_input_homogeneous(Coordinate(1), Coordinate(0));
    BOOST_TEST_EQ(without_allocations([&] { return ec::fixed_orientation(origin, unit, result); }), actual[1].sign());
    BOOST_TEST_EQ(without_allocations([&] { return ec::fixed_coordinate_compare(result, origin, 0); }), actual[0].sign());
    BOOST_TEST_EQ(without_allocations([&] { return ec::fixed_coordinate_compare(result, reversed, 1); }), 0);
}

template <typename Coordinate>
void test_integer()
{
    static_assert(has_fixed_value<Coordinate>::value, "Native integers supported");
    Coordinate const low = (std::numeric_limits<Coordinate>::min)();
    Coordinate const high = (std::numeric_limits<Coordinate>::max)();
    check_input(low, high);
    check_input(high, low);
    check_input(Coordinate(0), Coordinate(0));
    check_predicates<Coordinate>({{low, low, high, low, low, high}});
    check_predicates<Coordinate>({{low, high, high, low, high, high}});
    check_intersection<Coordinate>({{0, 0, 1, 1, 0, 1, 1, 0}});
    check_intersection<Coordinate>({{low, low, high, high, low, high, high, low}});
}

template <typename Coordinate>
void test_floating()
{
    Coordinate const tiny = (std::numeric_limits<Coordinate>::denorm_min)();
    Coordinate const normal = (std::numeric_limits<Coordinate>::min)();
    Coordinate const large = (std::numeric_limits<Coordinate>::max)();
    std::array<Coordinate, 12> const values{{0, -Coordinate(0), tiny, -tiny, normal, -normal,
        large, -large, Coordinate(0.1), Coordinate(-3.5), Coordinate(1), std::nextafter(normal, Coordinate(0))}};
    for (std::size_t i = 0; i < values.size(); ++i) { check_input(values[i], values[(i + 1) % values.size()]); }
    check_predicates<Coordinate>({{-large, -large, large, -large, -large, large}});
    check_predicates<Coordinate>({{tiny, 0, 0, tiny, 0, 0}});
    check_intersection<Coordinate>({{-large, -large, large, large, -large, large, large, -large}});
    check_intersection<Coordinate>({{0, 0, large, tiny, 0, tiny, large, 0}});
    check_intersection<Coordinate>({{large, 0, large, 1, 0, tiny, 1, tiny}});
    Coordinate const half = std::ldexp(large, -1);
    check_intersection<Coordinate>({{tiny, large, large, Coordinate(2) * tiny,
                                    0, Coordinate(3) * tiny, half, half}});
    check_intersection<Coordinate>({{0, 0, 1, 1, 0, 1, 1, -1}});
    check_intersection<Coordinate>({{0, 0, 1, 1, 0, 1, 1, 2}});
    BOOST_TEST_THROWS(ec::fixed_value((std::numeric_limits<Coordinate>::infinity)()), std::domain_error);
    BOOST_TEST_THROWS(ec::fixed_value((std::numeric_limits<Coordinate>::quiet_NaN)()), std::domain_error);
}

void test_random()
{
    std::mt19937_64 generator(0x86c3d799);
    for (int i = 0; i < 64; ++i)
    {
        std::array<double, 8> values;
        for (auto& value : values)
        {
            std::uint64_t bits = generator();
            if (((bits >> 52) & 0x7ff) == 0x7ff) { bits ^= std::uint64_t(1) << 52; }
            std::memcpy(&value, &bits, sizeof(value));
        }
        check_input(values[0], values[1]);
        check_predicates<double>({{values[0], values[1], values[2], values[3], values[4], values[5]}});
        if (i < 16) { check_intersection(values); }
        std::array<float, 8> floats;
        for (auto& value : floats)
        {
            std::uint32_t bits = static_cast<std::uint32_t>(generator());
            if (((bits >> 23) & 0xff) == 0xff) { bits ^= std::uint32_t(1) << 23; }
            std::memcpy(&value, &bits, sizeof(value));
        }
        check_input(floats[0], floats[1]);
        check_predicates<float>({{floats[0], floats[1], floats[2], floats[3], floats[4], floats[5]}});
        if (i < 16) { check_intersection(floats); }
    }
    std::uniform_int_distribution<std::int64_t> integers(
        (std::numeric_limits<std::int64_t>::min)(), (std::numeric_limits<std::int64_t>::max)());
    for (int i = 0; i < 64; ++i)
    {
        std::array<std::int64_t, 8> values;
        for (auto& value : values) { value = integers(generator); }
        check_predicates<std::int64_t>({{values[0], values[1], values[2], values[3], values[4], values[5]}});
        if (i < 16) { check_intersection(values); }
        std::array<std::uint64_t, 8> unsigned_values;
        for (auto& value : unsigned_values) { value = generator(); }
        check_predicates<std::uint64_t>({{unsigned_values[0], unsigned_values[1], unsigned_values[2],
                                         unsigned_values[3], unsigned_values[4], unsigned_values[5]}});
        if (i < 16) { check_intersection(unsigned_values); }
    }
}

void test_fixed_operations()
{
    using integer = ec::fixed_integer<4096>;
    integer a(1), b(1);
    a <<= 3000;
    a -= 1;
    b <<= 2500;
    b += 17;
    auto const product = without_allocations([&] { return ec::fixed_multiply<8192>(a, -b); });
    BOOST_TEST(cpp_int(product) == -cpp_int(a) * cpp_int(b));
    auto const common = without_allocations([&] { return ec::fixed_gcd(a, b); });
    BOOST_TEST(cpp_int(common) == boost::multiprecision::gcd(cpp_int(a), cpp_int(b)));
    ec::fixed_homogeneous<4, 128> origin;
    origin.W <<= 100;
    auto const normalized = without_allocations([&] { auto point = origin; ec::fixed_normalize(point); return point; });
    BOOST_TEST(normalized.X == 0 && normalized.Y == 0 && normalized.W == 1);
    BOOST_TEST_THROWS((ec::fixed_multiply<8>(ec::fixed_integer<8>(16), ec::fixed_integer<8>(16))), std::overflow_error);
    ec::fixed_homogeneous<8> invalid;
    invalid.W = 0;
    BOOST_TEST_THROWS(ec::fixed_normalize(invalid), std::domain_error);
    BOOST_TEST_THROWS(ec::fixed_coordinate_compare(invalid, invalid, 0), std::domain_error);
}

template <std::size_t Bits, std::size_t A, std::size_t B>
void check_multiply(ec::fixed_integer<A> const& a, ec::fixed_integer<B> const& b)
{
    cpp_int const expected = cpp_int(a) * cpp_int(b);
    auto const actual = without_allocations([&] { return ec::fixed_multiply<Bits>(a, b); });
    BOOST_TEST(cpp_int(actual) == expected);
}

void test_full_width_multiply()
{
    using integer = ec::fixed_integer<32768>;
    integer const half = (integer(1) << 32767) - 1;
    integer const maximum = (half << 1) + 1;
    integer alternating = 0;
    for (unsigned i = 0; i < 1024; ++i)
    {
        alternating <<= 32;
        alternating += (i & 1) ? 0x55555555u : 0xaaaaaaaau;
    }
    for (int a_sign : {-1, 1})
    for (int b_sign : {-1, 1})
    {
        integer const a = a_sign < 0 ? integer(-maximum) : maximum;
        integer const b = b_sign < 0 ? integer(-alternating) : alternating;
        check_multiply<65536>(a, b);
        check_multiply<65536>(b, a);
    }
    check_multiply<65536>(maximum, maximum);
    check_multiply<32768>(maximum, ec::fixed_integer<1>(-1));
    check_multiply<1>(maximum, ec::fixed_integer<1>(0));
    check_multiply<1>(ec::fixed_integer<1>(0), maximum);
    check_multiply<33>(ec::fixed_integer<17>(65535), ec::fixed_integer<17>(-65535));
    // Overflow can occur during multiplication, a later digit shift/add, or
    // narrowing an operand. Positive partial sums cannot hide excess bits.
    BOOST_TEST_THROWS((ec::fixed_multiply<32768>(maximum, ec::fixed_integer<2>(2))), std::overflow_error);
    BOOST_TEST_THROWS((ec::fixed_multiply<32768>(ec::fixed_integer<2>(-2), maximum)), std::overflow_error);
    BOOST_TEST_THROWS((ec::fixed_multiply<32>(ec::fixed_integer<33>(65536),
                                           ec::fixed_integer<33>(65536))), std::overflow_error);
    BOOST_TEST_THROWS((ec::fixed_multiply<1>(maximum, ec::fixed_integer<1>(1))), std::range_error);
}

void test_degree_eight_capacity()
{
    constexpr unsigned bits = ec::native_fixed_traits<double>::coordinate_bits;
    auto const value = ec::fixed_value((std::numeric_limits<double>::max)());
    auto const square = without_allocations([&] { return ec::fixed_multiply<2 * bits>(value, value); });
    auto const fourth = without_allocations([&] { return ec::fixed_multiply<4 * bits>(square, square); });
    auto const eighth = without_allocations([&] { return ec::fixed_multiply<8 * bits>(fourth, fourth); });
    using integer = ec::native_fixed_polynomial<double, 8, 256>::type;
    auto const polynomial = without_allocations([&]
    {
        integer result(0);
        for (unsigned i = 0; i < 256; ++i) { result += integer(eighth); }
        return result;
    });
    cpp_int expected(value);
    expected *= expected;
    expected *= expected;
    expected *= expected;
    expected <<= 8;
    BOOST_TEST(cpp_int(polynomial) == expected);
    BOOST_TEST_THROWS((integer(polynomial + integer(eighth))), std::overflow_error);
}

template <std::size_t Bits>
void check_narrow_capacity()
{
    ec::fixed_integer<1024> const limit = ec::fixed_integer<1024>(1) << Bits;
    for (int sign : {-1, 1})
    {
        ec::fixed_homogeneous<1024> source;
        source.X = sign * (limit - 1);
        source.Y = -source.X;
        source.W = limit - 1;
        ec::fixed_homogeneous<Bits> small;
        BOOST_TEST(without_allocations([&] { return ec::try_fixed_narrow(source, small); }));
        BOOST_TEST(cpp_int(source.X) == cpp_int(small.X));
        BOOST_TEST(cpp_int(source.Y) == cpp_int(small.Y));
        BOOST_TEST(cpp_int(source.W) == cpp_int(small.W));
        auto const saved = small;
        for (int component = 0; component < 3; ++component)
        {
            auto oversized = source;
            if (component == 0) oversized.X = sign * limit;
            if (component == 1) oversized.Y = sign * limit;
            if (component == 2) oversized.W = limit;
            BOOST_TEST(!without_allocations([&] { return ec::try_fixed_narrow(oversized, small); }));
            BOOST_TEST(small.X == saved.X && small.Y == saved.Y && small.W == saved.W);
        }
    }
    BOOST_TEST(ec::fixed_fits<Bits>(ec::fixed_integer<1>(0)));
    BOOST_TEST(ec::fixed_fits<Bits>(ec::fixed_integer<128>(-1)));
}

int main()
{
    test_integer<bool>();
    test_integer<std::int8_t>();
    test_integer<std::uint8_t>();
    test_integer<std::int16_t>();
    test_integer<std::uint16_t>();
    test_integer<std::int32_t>();
    test_integer<std::uint32_t>();
    test_integer<std::int64_t>();
    test_integer<std::uint64_t>();
    test_fixed_operations();
    test_full_width_multiply();
    test_degree_eight_capacity();
    check_narrow_capacity<1>();
    check_narrow_capacity<63>();
    check_narrow_capacity<64>();
    check_narrow_capacity<65>();
    check_narrow_capacity<127>();
    check_narrow_capacity<128>();
    check_narrow_capacity<129>();
    check_narrow_capacity<512>();
    int const original = std::fegetround();
    for (int mode : {FE_TONEAREST, FE_UPWARD, FE_DOWNWARD, FE_TOWARDZERO})
    {
        BOOST_TEST_EQ(std::fesetround(mode), 0);
        test_floating<float>();
        test_floating<double>();
    }
    BOOST_TEST_EQ(std::fesetround(original), 0);
    test_random();
    auto const hi = (std::numeric_limits<std::uint64_t>::max)();
    check_intersection<std::uint64_t>({{0, 0, hi, hi - 1, 0, 1, hi - 1, hi - 1}});
    check_intersection<std::uint64_t>({{0, 0, hi, 1, 0, hi / 2, hi - 1, hi / 2 + 1}});
    return boost::report_errors();
}
