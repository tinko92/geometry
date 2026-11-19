// Boost.Geometry
// Unit Test
// Copyright (c) 2026 Tinko Bartels.
// Licensed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt

#include <boost/core/lightweight_test.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/strategies/relate/cartesian_exact.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <new>
#include <random>
#include <type_traits>

namespace bg = boost::geometry;
namespace ec = bg::strategies::detail::exact_cartesian;
using rational = ec::exact_number;

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

void operator delete(void* pointer) noexcept { std::free(pointer); }
void* operator new[](std::size_t size) { return ::operator new(size); }
void operator delete[](void* pointer) noexcept { ::operator delete(pointer); }

#if defined(__cpp_sized_deallocation)
void operator delete(void* pointer, std::size_t) noexcept { std::free(pointer); }
void operator delete[](void* pointer, std::size_t) noexcept { std::free(pointer); }
#endif

template <typename Coordinate>
using point = bg::model::d2::point_xy<Coordinate>;

template <std::size_t Dimension, typename Point>
rational value(Point const& p)
{
    return ec::point_detail::input_exact_value(bg::get<Dimension>(p));
}

int compare(rational const& a, rational const& b)
{
    return (a > b) - (a < b);
}

template <typename A, typename B, typename P>
int orientation_oracle(A const& a, B const& b, P const& p)
{
    rational const determinant = (value<0>(b) - value<0>(a)) * (value<1>(p) - value<1>(a))
                               - (value<1>(b) - value<1>(a)) * (value<0>(p) - value<0>(a));
    return determinant.sign();
}

template <typename Coordinate, typename A, typename B>
void check_pair(A const& a, B const& b)
{
    int const x = compare(value<0>(a), value<0>(b));
    int const y = compare(value<1>(a), value<1>(b));
    int const order = x == 0 ? y : x;
    auto const ea = ec::make_point<Coordinate>(a);
    auto const eb = ec::make_point<Coordinate>(b);
    BOOST_TEST_EQ(ec::point_equal<Coordinate>::apply(a, b), order == 0);
    BOOST_TEST_EQ(ec::point_equal<Coordinate>::apply(a, b), ec::equal_points(ea, eb));
    BOOST_TEST_EQ((ec::point_compare<Coordinate, bg::strategy::compare::less>::apply(a, b)), order < 0);
    BOOST_TEST_EQ((ec::point_compare<Coordinate, bg::strategy::compare::greater>::apply(a, b)), order > 0);
    BOOST_TEST_EQ((ec::point_compare<Coordinate, bg::strategy::compare::equal_to>::apply(a, b)), order == 0);
}

template <typename Coordinate, typename A, typename B, typename P>
void check_orientation(A const& a, B const& b, P const& p)
{
    int const expected = orientation_oracle(a, b, p);
    bool const native = ec::native_predicates_detail::native_inputs<Coordinate, A, B, P>::value;
    allocations = 0;
    count_allocations = native;
    int const actual = ec::side<Coordinate>::apply(a, b, p);
    int const direct = ec::native_orientation<Coordinate>(a, b, p);
    std::array<int, 5> const permutations{{ec::side<Coordinate>::apply(b, p, a),
        ec::side<Coordinate>::apply(p, a, b), ec::side<Coordinate>::apply(b, a, p),
        ec::side<Coordinate>::apply(a, p, b), ec::side<Coordinate>::apply(p, b, a)}};
    count_allocations = false;
    if (native) { BOOST_TEST_EQ(allocations, 0u); }
    BOOST_TEST_EQ(actual, expected);
    BOOST_TEST_EQ(direct, expected);
    BOOST_TEST_EQ(permutations[0], expected);
    BOOST_TEST_EQ(permutations[1], expected);
    BOOST_TEST_EQ(permutations[2], -expected);
    BOOST_TEST_EQ(permutations[3], -expected);
    BOOST_TEST_EQ(permutations[4], -expected);
    BOOST_TEST_EQ(ec::side<Coordinate>::apply(a, b, p),
                  ec::orientation(ec::make_point<Coordinate>(a), ec::make_point<Coordinate>(b),
                                  ec::make_point<Coordinate>(p)));
}

template <typename Coordinate, typename P, typename SegmentPoint>
void check_winding(P const& p, SegmentPoint const& a, SegmentPoint const& b)
{
    using strategy = ec::winding<Coordinate>;
    int const side = orientation_oracle(a, b, p);
    int const ax = compare(value<0>(a), value<0>(p));
    int const bx = compare(value<0>(b), value<0>(p));
    int const ay = compare(value<1>(a), value<1>(p));
    int const by = compare(value<1>(b), value<1>(p));
    bool const boundary = side == 0 && ax * bx <= 0 && ay * by <= 0;
    int const increment = !boundary && ay <= 0 && by > 0 && side > 0 ? 1
                        : !boundary && by <= 0 && ay > 0 && side < 0 ? -1 : 0;
    for (bool initial_boundary : {false, true})
    {
        for (int initial : {-2, 0, 3})
        {
            typename strategy::state_type state{initial, initial_boundary};
            typename strategy::state_type reference{initial, initial_boundary};
            bool const keep_going = strategy::apply(p, a, b, state);
            BOOST_TEST_EQ(keep_going, !boundary);
            BOOST_TEST_EQ(state.boundary, initial_boundary || boundary);
            BOOST_TEST_EQ(state.count, initial + increment);
            BOOST_TEST_EQ(keep_going, strategy::apply(ec::make_point<Coordinate>(p),
                ec::make_point<Coordinate>(a), ec::make_point<Coordinate>(b), reference));
            BOOST_TEST_EQ(state.count, reference.count);
            BOOST_TEST_EQ(state.boundary, reference.boundary);
            BOOST_TEST_EQ(strategy::result(state), strategy::result(reference));
        }
    }
}

template <typename Coordinate>
void check_case(point<Coordinate> const& a, point<Coordinate> const& b, point<Coordinate> const& p)
{
    check_orientation<Coordinate>(a, b, p);
    check_pair<Coordinate>(a, b);
    check_pair<Coordinate>(b, p);
    check_pair<Coordinate>(p, a);
    check_pair<Coordinate>(p, p);
    check_winding<Coordinate>(p, a, b);
    check_winding<Coordinate>(p, b, a);
}

template <typename Function>
void check_no_allocations(Function const& function)
{
    allocations = 0;
    count_allocations = true;
    function();
    count_allocations = false;
    BOOST_TEST_EQ(allocations, 0u);
}

void test_orientation_scales()
{
    ec::point_cache_scope<double> operation;
    auto const check = [](std::array<double, 6> const& v)
    {
        check_orientation<double>(point<double>(v[0], v[1]), point<double>(v[2], v[3]),
                                  point<double>(v[4], v[5]));
    };
    check({{0, 0, 1, 0, 0, 1}});
    check({{0, 0, 0, 1, 1, 0}});
    check({{0, 0, 1, 1, 2, 2}});
    check({{0, 0, 0, 0, 1, 1}});

    double const large = std::ldexp(1.0, 52);
    check({{0, 0, large, large - 1, large - 1, large - 2}});
    double const minimum = std::ldexp(1.0, -200);
    double const maximum = std::ldexp(1.0, 200);
    // Keep finite values around the former custom filter's artificial guard;
    // the adaptive predicate must handle both sides of that guard exactly.
    for (double outside : {std::nextafter(minimum, 0.0),
                           std::nextafter(maximum, (std::numeric_limits<double>::infinity)())})
    {
        for (std::size_t i = 0; i != 6; ++i)
        {
            std::array<double, 6> values{{0, 0, 1, 0, 0, 1}};
            values[i] = outside;
            check(values);
            values[i] = -outside;
            check(values);
        }
    }
    // Include the lattice's smallest nonzero differences at the lower guard,
    // upper-guard products, and products spanning the full guarded range.
    double const next_minimum = std::nextafter(minimum, maximum);
    check({{minimum, minimum, next_minimum, minimum, minimum, next_minimum}});
    check({{-maximum, -maximum, maximum, -maximum, -maximum, maximum}});
    check({{0, 0, minimum, 0, 0, maximum}});

    for (int exponent = -200; exponent <= 200; ++exponent)
    {
        double const scale = std::ldexp(1.0, exponent);
        check({{0, 0, scale, 0, 0, scale}});
        check({{0, minimum, scale, scale, scale, std::nextafter(scale, 0.0)}});
    }
    // Random dyadic inputs inside the guard exercise certification, while
    // nearly collinear integer triples exercise cancellation at many scales.
    std::mt19937_64 generator(314159);
    for (int i = 0; i != 1024; ++i)
    {
        std::array<double, 6> values;
        for (auto& v : values)
        {
            int const exponent = static_cast<int>(generator() % 389) - 194;
            int const significand = static_cast<int>(generator() % 65) - 32;
            v = std::ldexp(double(significand), exponent);
        }
        check(values);
    }
    for (int exponent : {-200, -100, 0, 100, 147})
    {
        double const scale = std::ldexp(1.0, exponent);
        for (int offset = -8; offset <= 8; ++offset)
        {
            double const x = large + offset;
            std::array<double, 6> const values{{0, 0, x * scale, (x + 1) * scale,
                                               (x + 1) * scale, (x + 2) * scale}};
            check(values);
        }
    }
    check_no_allocations([&]
    {
        BOOST_TEST_EQ(ec::native_orientation<double>(point<double>(0, 0),
            point<double>(large, large - 1), point<double>(-large + 1, -large + 2)), 1);
    });
}

template <typename Coordinate, typename QueryRange, typename SegmentRange>
void check_winding_grid(QueryRange const& queries, SegmentRange const& segments)
{
    for (auto const& p : queries)
    {
        for (auto const& a : segments)
        {
            for (auto const& b : segments)
            {
                check_winding<Coordinate>(p, a, b);
            }
        }
    }
}

template <typename Coordinate>
void test_winding_grid(Coordinate offset = 0)
{
    ec::point_cache_scope<Coordinate> operation;
    std::array<point<Coordinate>, 9> grid;
    std::size_t index = 0;
    for (int x = -1; x <= 1; ++x)
    {
        for (int y = -1; y <= 1; ++y)
        {
            grid[index++] = point<Coordinate>(offset + Coordinate(x), offset + Coordinate(y));
        }
    }
    check_winding_grid<Coordinate>(grid, grid);
}

void test_constructed_winding_grid()
{
    ec::point_cache_scope<double> operation;
    using implicit = ec::exact_point<double>;
    std::array<std::array<implicit, 4>, 9> producers;
    std::array<implicit, 9> constructed;
    std::array<point<double>, 9> native;
    std::size_t index = 0;
    for (int x = -1; x <= 1; ++x)
    {
        for (int y = -1; y <= 1; ++y)
        {
            // The two integer-input lines meet at (x + 1/3, y + 1/3).
            auto& support = producers[index];
            support = {{implicit(x, y), implicit(x + 1, y + 1),
                        implicit(x, y + 1), implicit(x + 1, y - 1)}};
            constructed[index] = ec::construct_intersection(support[0], support[1], support[2], support[3]);
            native[index++] = point<double>(x, y);
        }
    }
    check_winding_grid<double>(constructed, constructed);
    check_winding_grid<double>(native, constructed);
    check_winding_grid<double>(constructed, native);
}

template <typename Coordinate>
void test_no_allocations()
{
    point<Coordinate> const a(0, 0), b(2, 2), p(0, 1), on(1, 1);
    bg::model::point<Coordinate, 2, bg::cs::cartesian> const other(2, 2);
    using inputs = ec::predicate_points<Coordinate, decltype(a), decltype(other)>;
    static_assert(inputs::use_native::value, "Different point models may share native coordinates");
    check_no_allocations([&]
    {
        BOOST_TEST_EQ(ec::side<Coordinate>::apply(a, b, p), 1);
        BOOST_TEST_EQ(ec::native_orientation<Coordinate>(a, other, p), 1);
        BOOST_TEST(ec::point_equal<Coordinate>::apply(b, other));
        BOOST_TEST(!ec::point_equal<Coordinate>::apply(a, p));
        BOOST_TEST((ec::point_compare<Coordinate, bg::strategy::compare::less>::apply(a, b)));
        BOOST_TEST((ec::point_compare<Coordinate, bg::strategy::compare::greater>::apply(b, a)));
        BOOST_TEST((ec::point_compare<Coordinate, bg::strategy::compare::equal_to>::apply(b, other)));
        typename ec::winding<Coordinate>::state_type state;
        BOOST_TEST(ec::winding<Coordinate>::apply(p, a, b, state));
        BOOST_TEST_EQ(state.count, 1);
    });
    // Collinearity and the corresponding boundary winding case stay exact.
    check_no_allocations([&]
    {
        BOOST_TEST_EQ(ec::side<Coordinate>::apply(a, b, on), 0);
        typename ec::winding<Coordinate>::state_type state;
        BOOST_TEST(!ec::winding<Coordinate>::apply(on, a, b, state));
        BOOST_TEST(state.boundary);
    });
}

template <typename Coordinate>
void test_integers()
{
    ec::point_cache_scope<Coordinate> operation;
    using p = point<Coordinate>;
    Coordinate const lo = (std::numeric_limits<Coordinate>::min)();
    Coordinate const hi = (std::numeric_limits<Coordinate>::max)();
    check_case(p(lo, lo), p(hi, lo), p(lo, hi));
    check_case(p(lo, hi), p(hi, lo), p(hi, hi));
    check_case(p(lo, lo), p(hi, hi), p(0, 0));
    // Full-span differences exceed the native product guard; cancellation
    // leaves determinant -1 or zero rather than a large well-separated sign.
    check_case(p(lo, lo), p(hi, Coordinate(hi - 1)), p(Coordinate(hi - 1), Coordinate(hi - 2)));
    check_case(p(lo, lo), p(hi, hi), p(Coordinate(hi - 1), Coordinate(hi - 1)));
    check_case(p(hi, hi), p(Coordinate(hi - 1), hi), p(hi, Coordinate(hi - 1)));
    check_case(p(lo, lo), p(Coordinate(lo + 1), lo), p(lo, Coordinate(lo + 1)));
    check_case(p(hi, lo), p(hi, lo), p(lo, hi));
    std::mt19937_64 generator(5489);
    for (int i = 0; i != 32; ++i)
    {
        std::array<Coordinate, 6> values;
        for (auto& v : values) { v = static_cast<Coordinate>(generator()); }
        check_case(p(values[0], values[1]), p(values[2], values[3]), p(values[4], values[5]));
    }
    test_no_allocations<Coordinate>();
}

template <typename Coordinate>
void test_floating()
{
    ec::point_cache_scope<Coordinate> operation;
    using p = point<Coordinate>;
    Coordinate const large = std::ldexp(Coordinate(1), std::numeric_limits<Coordinate>::digits - 1);
    // All binary32 inputs fit the binary64 arithmetic range. Binary64 fixtures
    // avoid arithmetic underflow/overflow, which are outside the contract.
    Coordinate const tiny = sizeof(Coordinate) == sizeof(float)
        ? (std::numeric_limits<Coordinate>::denorm_min)() : std::ldexp(Coordinate(1), -200);
    Coordinate const normal = sizeof(Coordinate) == sizeof(float)
        ? (std::numeric_limits<Coordinate>::min)() : std::ldexp(Coordinate(1), -180);
    Coordinate const maximum = sizeof(Coordinate) == sizeof(float)
        ? (std::numeric_limits<Coordinate>::max)() : std::ldexp(Coordinate(1), 200);
    Coordinate const below_one = std::nextafter(Coordinate(1), Coordinate(0));
    check_case(p(0, 0), p(large, large - 1), p(large - 1, large - 2));
    check_case(p(0, 0), p(1, 1), p(1, below_one));
    check_case(p(0, 0), p(1, 1), p(2, 2));
    check_case(p(0, -Coordinate(0)), p(-Coordinate(0), 0), p(1, 1));
    check_case(p(0, 0), p(tiny, 0), p(0, tiny));
    check_case(p(0, 0), p(normal, normal), p(normal, std::nextafter(normal, maximum)));
    check_case(p(-maximum, -maximum), p(maximum, maximum), p(maximum, -maximum));
    check_case(p(0, 0), p(maximum, maximum), p(std::nextafter(maximum, Coordinate(0)), maximum));
    check_case(p(0, 0), p(maximum, tiny), p(maximum, normal));
    check_case(p(tiny, tiny), p(tiny, tiny), p(maximum, -maximum));
    for (int exponent : {-100, -10, 0, 10, 100})
    {
        Coordinate const scale = std::ldexp(Coordinate(1), exponent);
        check_case(p(0, 0), p(scale, scale), p(2 * scale, 2 * scale));
    }
    // Retain diverse significands and exponents without binary64 range errors.
    using bits_type = typename std::conditional<sizeof(Coordinate) == 4,
                                               std::uint32_t, std::uint64_t>::type;
    std::mt19937_64 generator(1729);
    for (int i = 0; i != 32; ++i)
    {
        std::array<Coordinate, 6> values;
        for (auto& v : values)
        {
            do
            {
                bits_type const bits = static_cast<bits_type>(generator());
                std::memcpy(&v, &bits, sizeof(v));
            }
            while (!std::isfinite(v) || (v != 0 && (std::abs(v) < tiny || std::abs(v) > maximum)));
        }
        check_case(p(values[0], values[1]), p(values[2], values[3]), p(values[4], values[5]));
    }
    test_no_allocations<Coordinate>();
    check_no_allocations([&]
    {
        BOOST_TEST_EQ(ec::side<Coordinate>::apply(p(0, 0), p(large, large - 1),
                                                 p(large - 1, large - 2)), -1);
    });
}

void test_fallbacks()
{
    ec::point_cache_scope<double> operation;
    using p = point<double>;
    using i = point<std::int64_t>;
    using f = point<float>;
    using implicit = ec::exact_point<double>;
    static_assert(!ec::predicate_points<double, p, i>::use_native::value, "Reject lossy int64/double promotion");
    static_assert(!ec::predicate_points<double, i, i>::use_native::value, "Do not round int64 to the strategy type");
    static_assert(!ec::predicate_points<double, implicit, p>::use_native::value, "Retain implicit exact data");
    static_assert(!ec::predicate_points<float, f, p>::use_native::value, "Do not narrow double to float");
    static_assert(ec::predicate_points<double, f, p>::use_native::value, "Lossless promotions stay native");
    using small = point<std::int32_t>;
    static_assert(ec::predicate_points<double, small, p>::use_native::value, "Int32 fits double exactly");
    check_pair<double>(small(2147483647, 0), p(2147483647.0, 0));
    check_orientation<double>(small(2147483647, 0), p(2147483647.0, 1), f(0, 0));
    check_winding<double>(small(0, 1), p(0, 0), p(2, 2));
    double const next = std::nextafter(1.0, 2.0);
    check_pair<double>(f(1, 0), p(next, 0));
    check_orientation<double>(f(0, 0), f(1, 1), p(1, next));
    check_winding<double>(p(1, next), f(0, 0), f(2, 2));

    // Retain unit differences beyond double's consecutive-integer range under
    // the integer strategy rather than requesting a now-rejected lossy mix.
    {
        ec::point_cache_scope<std::int64_t> integer_operation;
        std::int64_t const n = 9007199254740992LL;
        check_pair<std::int64_t>(i(n, 0), i(n + 1, 0));
        check_orientation<std::int64_t>(i(n, n), i(n + 1, n), i(n, n + 1));
        check_winding<std::int64_t>(i(n + 1, n + 1), i(n, n), i(n + 2, n + 2));
    }

    implicit const a(0, 0), b(1, 1), c(0, 1), d(1, -1), e(1, 0), f0(0, -1);
    auto const third = ec::construct_intersection(a, b, c, d);
    auto const nested = ec::construct_intersection(third, e, f0, c);
    check_pair<double>(third, p(1.0 / 3.0, 1.0 / 3.0));
    check_pair<double>(third, third);
    check_pair<double>(third, nested);
    check_orientation<double>(p(0, 0), p(1, 1), third);
    check_orientation<double>(third, p(1, 0), nested);
    check_winding<double>(third, p(0, 0), p(1, 1));
    check_winding<double>(p(0, 0), third, nested);
}

int main()
{
    test_integers<std::int8_t>();
    test_integers<std::uint8_t>();
    test_integers<std::int16_t>();
    test_integers<std::uint16_t>();
    test_integers<std::int32_t>();
    test_integers<std::uint32_t>();
    test_integers<std::int64_t>();
    test_integers<std::uint64_t>();
    test_winding_grid<std::int16_t>();
    test_winding_grid<std::int64_t>();
    test_winding_grid<std::int64_t>(9007199254740992LL);
    test_orientation_scales();
    test_floating<float>();
    test_floating<double>();
    test_winding_grid<double>();
    test_winding_grid<double>(std::ldexp(1.0, 52));
    test_constructed_winding_grid();
    test_fallbacks();
    return boost::report_errors();
}
