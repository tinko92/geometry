// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.

#include <boost/core/lightweight_test.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/strategies/cartesian/exact/native_predicates.hpp>
#include <boost/multiprecision/cpp_int.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <memory>
#include <new>
#include <random>
#include <type_traits>
#include <utility>

namespace bg = boost::geometry;
namespace ec = bg::strategies::detail::exact_cartesian;
using boost::multiprecision::cpp_int;

namespace
{
bool count_allocations = false;
std::size_t allocations = 0;
}

void* operator new(std::size_t size)
{
    if (count_allocations) { ++allocations; }
    if (void* result = std::malloc(size == 0 ? 1 : size)) { return result; }
    throw std::bad_alloc();
}
void operator delete(void* value) noexcept { std::free(value); }
void* operator new[](std::size_t size) { return ::operator new(size); }
void operator delete[](void* value) noexcept { ::operator delete(value); }
#if defined(__cpp_sized_deallocation)
void operator delete(void* value, std::size_t) noexcept { std::free(value); }
void operator delete[](void* value, std::size_t) noexcept { std::free(value); }
#endif

template <typename C>
struct observed_point
{
    C values[2];
    mutable unsigned reads[2]{};

    observed_point(C x, C y) : values{x, y} {}
    observed_point(observed_point const&) = delete;
    observed_point& operator=(observed_point const&) = delete;
};

namespace boost { namespace geometry { namespace traits
{
template <typename C> struct tag<observed_point<C>> { using type = point_tag; };
template <typename C> struct coordinate_type<observed_point<C>> { using type = C; };
template <typename C> struct coordinate_system<observed_point<C>> { using type = cs::cartesian; };
template <typename C> struct dimension<observed_point<C>> : std::integral_constant<int, 2> {};
template <typename C, std::size_t D>
struct access<observed_point<C>, D>
{
    static C get(observed_point<C> const& point)
    {
        ++point.reads[D];
        return point.values[D];
    }
};
}}}

template <typename C>
using point = bg::model::d2::point_xy<C>;

template <typename Function>
int without_allocations(Function function)
{
    allocations = 0;
    count_allocations = true;
    int const result = function();
    count_allocations = false;
    BOOST_TEST_EQ(allocations, 0u);
    return result;
}

// Scale dyadic coordinates to a common integer lattice, independently of the
// predicate's expansion arithmetic and fixed-capacity numerical code.
template <typename C>
int orientation_oracle(std::array<C, 6> const& values)
{
    std::array<cpp_int, 6> integers;
    if constexpr (std::is_integral<C>::value)
    {
        for (std::size_t i = 0; i != values.size(); ++i) { integers[i] = values[i]; }
    }
    else
    {
        std::array<int, 6> exponents{};
        std::array<std::int64_t, 6> mantissas{};
        int minimum = 0;
        for (std::size_t i = 0; i != values.size(); ++i)
        {
            double const fraction = std::frexp(double(values[i]), &exponents[i]);
            mantissas[i] = std::int64_t(std::ldexp(fraction, std::numeric_limits<C>::digits));
            exponents[i] -= std::numeric_limits<C>::digits;
            if (exponents[i] < minimum) { minimum = exponents[i]; }
        }
        for (std::size_t i = 0; i != values.size(); ++i)
        {
            integers[i] = mantissas[i];
            integers[i] <<= exponents[i] - minimum;
        }
    }
    cpp_int const determinant = (integers[2] - integers[0]) * (integers[5] - integers[1])
                              - (integers[3] - integers[1]) * (integers[4] - integers[0]);
    return (determinant > 0) - (determinant < 0);
}

template <typename C>
void check_orientation(std::array<C, 6> const& values)
{
    int const expected = orientation_oracle(values);
    point<C> const a(values[0], values[1]);
    bg::model::point<C, 2, bg::cs::cartesian> const b(values[2], values[3]);
    observed_point<C> const c(values[4], values[5]);
    BOOST_TEST_EQ(without_allocations([&] { return ec::native_orientation<C>(a, b, c); }), expected);
    BOOST_TEST_EQ(without_allocations([&] { return ec::native_orientation<C>(b, c, a); }), expected);
    BOOST_TEST_EQ(without_allocations([&] { return ec::native_orientation<C>(c, a, b); }), expected);
    BOOST_TEST_EQ(without_allocations([&] { return ec::native_orientation<C>(b, a, c); }), -expected);
    BOOST_TEST_EQ(without_allocations([&] { return ec::native_orientation<C>(a, c, b); }), -expected);
    BOOST_TEST_EQ(without_allocations([&] { return ec::native_orientation<C>(c, b, a); }), -expected);
}

template <typename C>
void test_borrowing()
{
    using native = observed_point<C>;
    using view = ec::native_point_view<C, native>;
    using points = ec::predicate_points<C, native, point<C>>;
    static_assert(sizeof(view) == sizeof(native const*), "A native view only borrows its input");
    static_assert(std::is_trivially_copy_constructible<view>::value, "Copies retain the reference");
    static_assert(!std::is_constructible<view, native&&>::value, "Do not store references to temporaries");
    static_assert(std::is_same<decltype(points::make(std::declval<native const&>())), view>::value,
                  "Native predicates do not construct owning exact points");

    native const a(0, 0), b(2, 2);
    auto const av = points::make(a);
    auto const bv = points::make(b);
    auto const copy = av;
    BOOST_TEST(std::addressof(av.point()) == std::addressof(a));
    BOOST_TEST(std::addressof(copy.point()) == std::addressof(a));
    BOOST_TEST_EQ(a.reads[0] + a.reads[1] + b.reads[0] + b.reads[1], 0u);
    BOOST_TEST_EQ(without_allocations([&] { return ec::compare_x(av, bv); }), -1);
    BOOST_TEST_EQ(a.reads[0], 1u);
    BOOST_TEST_EQ(b.reads[0], 1u);
    BOOST_TEST_EQ(a.reads[1] + b.reads[1], 0u);
    BOOST_TEST_EQ(without_allocations([&] { return ec::compare_y(av, bv); }), -1);
    BOOST_TEST_EQ(a.reads[1], 1u);
    BOOST_TEST_EQ(b.reads[1], 1u);
    point<C> const same(2, 2);
    BOOST_TEST(ec::equal_points(bv, points::make(same)));
    BOOST_TEST(!ec::equal_points(av, bv));
}

template <typename C>
void test_floating()
{
    test_borrowing<C>();
    check_orientation<C>({{0, 0, 1, 0, 0, 1}});
    check_orientation<C>({{0, 0, 1, 1, 2, 2}});
    check_orientation<C>({{0, -C(0), -C(0), 0, 1, 1}});
    check_orientation<C>({{0, 0, 1, 1, 1, std::nextafter(C(1), C(0))}});
    C const large = std::ldexp(C(1), std::numeric_limits<C>::digits - 1);
    for (int exponent : {-70, -10, 0, 10, 70})
    {
        C const scale = std::ldexp(C(1), exponent);
        check_orientation<C>({{0, 0, large * scale, (large - 1) * scale,
                               (-large + 1) * scale, (-large + 2) * scale}});
        check_orientation<C>({{0, 0, scale, scale, 2 * scale, 2 * scale}});
    }
    std::mt19937_64 random(1729);
    for (int i = 0; i != 1024; ++i)
    {
        std::array<C, 6> values;
        for (C& value : values)
        {
            C const mantissa = C(int(random() % 65537) - 32768);
            value = std::ldexp(mantissa, int(random() % 141) - 70);
        }
        check_orientation(values);
    }
}

void test_adaptive_orientation()
{
    double const large = std::ldexp(1.0, 52);
    point<double> const a(0, 0), b(large, large - 1), c(-large + 1, -large + 2);
    using naive = bg::strategy::side::side_robust<double, bg::strategy::side::fp_equals_policy, 0>;
    BOOST_TEST_EQ(naive::apply(a, b, c), 0);
    BOOST_TEST_EQ(without_allocations([&] { return ec::native_orientation<double>(a, b, c); }), 1);
    // An exact nonzero determinant must not be classified using the default
    // side_robust wrapper's tolerance-based zero policy.
    BOOST_TEST_EQ(bg::strategy::side::side_robust<double>::apply(a, b, c), 0);
}

void test_lossless_promotions()
{
    using f = observed_point<float>;
    using i = point<std::int32_t>;
    using d = point<double>;
    using points = ec::predicate_points<double, f, i, d>;
    static_assert(points::use_native::value, "Lossless native promotions stay borrowed");
    static_assert(std::is_same<decltype(points::make(std::declval<f const&>())),
                              ec::native_point_view<double, f>>::value, "Retain the native point type");
    static_assert(!ec::native_predicates_detail::native_inputs<double, point<std::int64_t>>::value,
                  "Do not round int64 inputs to double");
    static_assert(!ec::native_predicates_detail::native_inputs<float, d>::value,
                  "Do not narrow double inputs");
    f const a(0, 0);
    i const b(2147483647, 2147483647);
    d const c(2147483647.0, 2147483646.0);
    auto const av = points::make(a);
    auto const bv = points::make(b);
    BOOST_TEST_EQ(a.reads[0] + a.reads[1], 0u);
    BOOST_TEST_EQ(without_allocations([&] { return ec::native_orientation<double>(a, b, c); }), -1);
    BOOST_TEST_EQ(without_allocations([&] { return ec::compare_x(av, bv); }), -1);
    i const negative(-1, 0);
    point<std::uint32_t> const positive(4294967295U, 0);
    using integer_points = ec::predicate_points<std::int64_t, i, decltype(positive)>;
    static_assert(integer_points::use_native::value, "Promote signed/unsigned inputs losslessly");
    BOOST_TEST_EQ(ec::compare_x(integer_points::make(negative), integer_points::make(positive)), -1);
    BOOST_TEST_EQ(without_allocations([&]
    {
        return ec::native_orientation<std::int64_t>(negative, positive, point<std::int16_t>(0, 1));
    }), 1);
}

template <typename C>
void test_integral()
{
    test_borrowing<C>();
    C const low = (std::numeric_limits<C>::min)(), high = (std::numeric_limits<C>::max)();
    check_orientation<C>({{low, low, high, low, low, high}});
    check_orientation<C>({{low, low, high, high, 0, 0}});
    check_orientation<C>({{low, low, high, C(high - 1), C(high - 1), C(high - 2)}});
    check_orientation<C>({{high, high, high, high, low, low}});
    std::mt19937_64 random(5489);
    for (int i = 0; i != 128; ++i)
    {
        std::array<C, 6> values;
        for (C& value : values) { value = static_cast<C>(random()); }
        check_orientation(values);
    }
}

int main()
{
    test_floating<float>();
    test_floating<double>();
    test_adaptive_orientation();
    test_lossless_promotions();
    test_integral<std::int8_t>();
    test_integral<std::uint8_t>();
    test_integral<std::int16_t>();
    test_integral<std::uint16_t>();
    test_integral<std::int32_t>();
    test_integral<std::uint32_t>();
    test_integral<std::int64_t>();
    test_integral<std::uint64_t>();
    return boost::report_errors();
}
