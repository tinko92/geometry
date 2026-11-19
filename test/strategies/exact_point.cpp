// Boost.Geometry
// Unit Test
// Copyright (c) 2026
// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include <boost/core/lightweight_test.hpp>

#include <atomic>
#include <cfenv>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <new>
#include <thread>
#include <type_traits>
#include <vector>

#include <boost/geometry/geometries/concepts/point_concept.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/strategies/cartesian/exact/point.hpp>
#include "exact_test_context.hpp"

namespace
{
thread_local bool count_allocations = false;
thread_local std::size_t allocations = 0;
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

namespace bg = boost::geometry;
namespace ec = bg::strategies::detail::exact_cartesian;
using rational = ec::exact_number;
using integer = boost::multiprecision::cpp_int;
using approximate_output = bg::strategies::relate::approximate_output;
using nearest_output = bg::strategies::relate::nearest_output;

template <typename Coordinate>
using input_point = bg::model::d2::point_xy<Coordinate>;

template <typename Coordinate>
using point = ec::exact_point<Coordinate>;

template <typename Coordinate>
point<Coordinate> third()
{
    return ec::construct_intersection(
        ec::make_point(input_point<Coordinate>(0, 0)),
        ec::make_point<Coordinate>(input_point<Coordinate>(1, 1)),
        ec::make_point(input_point<Coordinate>(0, 1)),
        ec::make_point(input_point<Coordinate>(1, -1)));
}

template <typename Coordinate>
void check_bounds(point<Coordinate> const& p)
{
    for (std::size_t dimension = 0; dimension != 2; ++dimension)
    {
        rational const value = p.exact_coordinate(dimension);
        for (auto const& bounds : {p.approximation(dimension), p.refined_approximation(dimension)})
        {
            if (std::isfinite(bounds.lower()))
            {
                BOOST_TEST(ec::exact_value(bounds.lower()) <= value);
            }
            if (std::isfinite(bounds.upper()))
            {
                BOOST_TEST(value <= ec::exact_value(bounds.upper()));
            }
        }
    }
}

template <typename Coordinate>
void test_constructed_points()
{
    ec::point_cache_scope<Coordinate> operation;
    BOOST_CONCEPT_ASSERT((bg::concepts::Point<point<Coordinate>>));
    BOOST_CONCEPT_ASSERT((bg::concepts::ConstPoint<point<Coordinate>>));
    static_assert(std::is_same<typename bg::coordinate_type<point<Coordinate>>::type,
                               rational>::value, "Traits must expose true rational coordinates");

    point<Coordinate> const a(0, 0), b(1, 1), c(0, 1), d(1, -1);
    auto const p = third<Coordinate>();
    rational const one_third = rational(1) / 3;
    BOOST_TEST(p.exact_coordinate(0) == one_third);
    BOOST_TEST(p.exact_coordinate(1) == one_third);
    BOOST_TEST(bg::get<0>(p) == one_third);
    BOOST_TEST(bg::get<1>(p) == one_third);
    BOOST_TEST(p.homogeneous().W > 0);
    BOOST_TEST(p.homogeneous().X == p.homogeneous().Y);
    BOOST_TEST_EQ(ec::orientation(a, b, p), 0);
    BOOST_TEST_EQ(ec::orientation(c, d, p), 0);
    BOOST_TEST(ec::on_segment(a, b, p));
    BOOST_TEST(ec::on_segment(d, c, p));
    BOOST_TEST(ec::on_segment(p, p, p));
    BOOST_TEST(!ec::on_segment(a, p, b));

    auto const other = ec::construct_intersection(point<Coordinate>(0, 0), point<Coordinate>(2, 2),
                                                 point<Coordinate>(0, 2), point<Coordinate>(1, -3));
    BOOST_TEST(ec::equal_points(p, other));
    BOOST_TEST_EQ(ec::compare_x(p, other), 0);
    BOOST_TEST_EQ(ec::compare_y(p, other), 0);
    BOOST_TEST_EQ(p.refined_approximation(0).approximate, 1.0 / 3.0);
    auto const shared = p;
    BOOST_TEST_EQ(p.refined_approximation(0).approximate, shared.refined_approximation(0).approximate);
    BOOST_TEST_EQ(p.refined_approximation(0).error, shared.refined_approximation(0).error);
    BOOST_TEST(ec::equal_points(p, ec::construct_intersection(b, a, c, d)));
    BOOST_TEST(ec::equal_points(p, ec::construct_intersection(a, b, d, c)));
    BOOST_TEST(ec::equal_points(p, ec::construct_intersection(c, d, a, b)));
    BOOST_TEST(ec::equal_points(p, ec::construct_intersection(d, c, b, a)));

    auto const nested = ec::construct_intersection(p, point<Coordinate>(1, 0),
                                                  point<Coordinate>(0, -1), point<Coordinate>(0, 1));
    BOOST_TEST(nested.exact_coordinate(0) == 0);
    BOOST_TEST(nested.exact_coordinate(1) == rational(1) / 2);
    BOOST_TEST(!nested.has_fixed_representation());
    BOOST_TEST_EQ(ec::orientation(p, point<Coordinate>(1, 0), nested), 0);
    check_bounds(p);
    check_bounds(other);
    check_bounds(nested);

    auto copy = p;
    auto const immutable_bounds = p.approximation(0);
    double const approximate = immutable_bounds.approximate;
    double const error = immutable_bounds.error;
    BOOST_TEST(ec::equal_points(copy, p));
    auto const identity = ec::make_point<Coordinate>(p);
    BOOST_TEST(ec::equal_points(identity, p));
    auto const inferred_identity = ec::make_point(p);
    BOOST_TEST(ec::equal_points(inferred_identity, p));
    bg::set<0>(copy, rational(7) / 11);
    BOOST_TEST(bg::get<0>(copy) == rational(7) / 11);
    BOOST_TEST(bg::get<1>(copy) == one_third);
    BOOST_TEST(bg::get<0>(p) == one_third);
    BOOST_TEST_EQ(p.approximation(0).approximate, approximate);
    BOOST_TEST_EQ(p.approximation(0).error, error);
    check_bounds(copy);

    auto const shared_construction = ec::construct_intersection(a, b, c, d);
    auto const copied_construction = shared_construction;
    BOOST_TEST(ec::equal_points(copied_construction, p));
    BOOST_TEST_EQ(ec::orientation(a, b, copied_construction), 0);

    // Initial zero-denominator refinement rejects parallel supporting lines:
    // a unique finite intersection is required at construction.
    auto parallel = [&]
    {
        return ec::construct_intersection(a, b, point<Coordinate>(0, 1), point<Coordinate>(1, 2));
    };
    BOOST_TEST_THROWS(parallel(), std::domain_error);
    BOOST_TEST_THROWS(parallel(), std::domain_error);
}

template <typename Coordinate>
void test_unrepresentable_third()
{
    ec::point_cache_scope<Coordinate> operation;
    auto const p = third<Coordinate>();
    input_point<Coordinate> rounded;
    ec::assign_output<nearest_output>(p, rounded);
    auto const q = ec::make_point(rounded);
    rational const discrepancy = 3 * ec::exact_value(bg::get<0>(rounded)) - 1;
    BOOST_TEST(discrepancy != 0);
    int const expected = discrepancy > 0 ? 1 : -1;
    BOOST_TEST_EQ(ec::orientation(point<Coordinate>(0, 1), point<Coordinate>(1, -1), q), expected);
    BOOST_TEST_EQ(ec::orientation(point<Coordinate>(0, 1), point<Coordinate>(1, -1), p), 0);
    BOOST_TEST(!ec::on_segment(point<Coordinate>(0, 1), point<Coordinate>(1, -1), q));
    BOOST_TEST(!ec::equal_points(p, q));
    BOOST_TEST_EQ(ec::compare_x(q, p), expected);
    BOOST_TEST_EQ(ec::compare_y(p, q), -expected);
}

void test_cancellation()
{
    ec::point_cache_scope<double> operation;
    double const large = std::ldexp(1.0, 52);
    point<double> const a(large, large), b(large + 4, large + 4);
    point<double> const c(large, large + 3), d(large + 4, large);
    auto const p = ec::construct_intersection(a, b, c, d);
    rational const expected = ec::exact_value(large) + rational(12) / 7;
    BOOST_TEST(p.exact_coordinate(0) == expected);
    BOOST_TEST(p.exact_coordinate(1) == expected);
    BOOST_TEST_EQ(ec::compare_x(p, point<double>(large + 2, large + 2)), -1);
    BOOST_TEST_EQ(ec::orientation(c, d, p), 0);
    BOOST_TEST(ec::on_segment(a, b, p));
    check_bounds(p);

    // Almost parallel lines: the tiny crossing cannot survive rounded line data.
    double const tiny = std::ldexp(1.0, -100);
    double const below_one = std::nextafter(1.0, 0.0);
    auto const close = ec::construct_intersection(point<double>(0, 0), point<double>(1, 1),
                                                 point<double>(0, tiny), point<double>(1, below_one));
    rational const parameter = ec::exact_value(tiny)
                             / (ec::exact_value(tiny) + 1 - ec::exact_value(below_one));
    BOOST_TEST(close.exact_coordinate(0) == parameter);
    BOOST_TEST(close.exact_coordinate(1) == parameter);
    BOOST_TEST_EQ(ec::compare_x(close, point<double>(0, 0)), 1);
    BOOST_TEST_EQ(ec::orientation(point<double>(0, tiny), point<double>(1, below_one), close), 0);
    BOOST_TEST(ec::on_segment(point<double>(0, 0), point<double>(1, 1), close));
    check_bounds(close);
}

template <typename Coordinate>
void check_leaf_orientation(point<Coordinate> const& a, point<Coordinate> const& b,
                            point<Coordinate> const& c)
{
    rational const determinant = (b.exact_coordinate(0) - a.exact_coordinate(0))
                               * (c.exact_coordinate(1) - a.exact_coordinate(1))
                               - (b.exact_coordinate(1) - a.exact_coordinate(1))
                               * (c.exact_coordinate(0) - a.exact_coordinate(0));
    int const expected = determinant > 0 ? 1 : determinant < 0 ? -1 : 0;
    BOOST_TEST_EQ(ec::orientation(a, b, c), expected);
    BOOST_TEST_EQ(ec::orientation(b, c, a), expected);
    BOOST_TEST_EQ(ec::orientation(c, a, b), expected);
    BOOST_TEST_EQ(ec::orientation(b, a, c), -expected);
    BOOST_TEST_EQ(ec::orientation(a, c, b), -expected);
    BOOST_TEST_EQ(ec::orientation(c, b, a), -expected);
}

void test_floating_leaf_fallbacks()
{
    ec::point_cache_scope<double> operation;
    ec::point_cache_scope<float> float_operation;
    double const large = std::ldexp(1.0, 52);
    check_leaf_orientation(point<double>(0, 0), point<double>(large, large - 1),
                           point<double>(large - 1, large - 2));
    float const large_float = std::ldexp(1.0f, 23);
    check_leaf_orientation(point<float>(0, 0), point<float>(large_float, large_float - 1),
                           point<float>(large_float - 1, large_float - 2));
    int const exponents[] = {0, 100, -100, 190, -190, 450, -450};
    for (int exponent : exponents)
    {
        double const scale = std::ldexp(1.0, exponent);
        point<double> const a(scale, scale), b(3 * scale, 3 * scale);
        point<double> const c(2 * scale, std::nextafter(2 * scale,
                                                       (std::numeric_limits<double>::infinity)()));
        check_leaf_orientation(a, b, c);
        check_leaf_orientation(a, b, point<double>(2 * scale, 2 * scale));
        check_bounds(c);
    }
}

void test_large_integers()
{
    ec::point_cache_scope<std::int64_t> operation;
    ec::point_cache_scope<std::uint64_t> unsigned_operation;
    using coordinate = std::int64_t;
    coordinate const low = (std::numeric_limits<coordinate>::min)();
    coordinate const high = (std::numeric_limits<coordinate>::max)();
    point<coordinate> const a(low, low), b(high, high - 1), c(high - 1, high - 2);
    BOOST_TEST_EQ(ec::orientation(a, b, c), -1);
    BOOST_TEST_EQ(ec::orientation(a, c, b), 1);
    BOOST_TEST_EQ(ec::orientation(a, b, b), 0);
    BOOST_TEST_EQ(ec::compare_x(b, c), 1);
    BOOST_TEST_EQ(ec::compare_y(b, c), 1);
    BOOST_TEST(!ec::equal_points(b, c));
    check_bounds(a);
    check_bounds(b);
    check_bounds(c);

    auto const p = ec::construct_intersection(point<coordinate>(low, 0), point<coordinate>(high, 0),
                                             point<coordinate>(0, low), point<coordinate>(0, high));
    BOOST_TEST(ec::equal_points(p, point<coordinate>(0, 0)));
    check_bounds(p);

    using unsigned_coordinate = std::uint64_t;
    unsigned_coordinate const maximum = (std::numeric_limits<unsigned_coordinate>::max)();
    point<unsigned_coordinate> const u(0, 0), v(maximum, maximum - 1), w(maximum - 1, maximum - 2);
    BOOST_TEST_EQ(ec::orientation(u, v, w), -1);
    BOOST_TEST_EQ(ec::compare_x(v, w), 1);
    check_bounds(v);
}

template <typename Coordinate>
void test_integral_leaf_no_allocations()
{
    ec::point_cache_scope<Coordinate> operation;
    Coordinate const low = (std::numeric_limits<Coordinate>::lowest)();
    Coordinate const high = (std::numeric_limits<Coordinate>::max)();
    std::array<std::array<Coordinate, 6>, 3> const inputs{{
        {{low, low, high, Coordinate(high - 1), Coordinate(high - 1), Coordinate(high - 2)}},
        {{low, high, high, low, high, high}},
        {{0, 0, 3, 2, 2, 3}}}};
    for (auto const& input : inputs)
    {
        point<Coordinate> const a(input[0], input[1]), b(input[2], input[3]), c(input[4], input[5]);
        rational const determinant = (ec::exact_value(input[2]) - ec::exact_value(input[0]))
                                   * (ec::exact_value(input[5]) - ec::exact_value(input[1]))
                                   - (ec::exact_value(input[3]) - ec::exact_value(input[1]))
                                   * (ec::exact_value(input[4]) - ec::exact_value(input[0]));
        int const expected = determinant.sign();
        allocations = 0;
        count_allocations = true;
        int const abc = ec::orientation(a, b, c), acb = ec::orientation(a, c, b);
        int const bca = ec::orientation(b, c, a);
        count_allocations = false;
        BOOST_TEST_EQ(allocations, 0u);
        BOOST_TEST_EQ(abc, expected);
        BOOST_TEST_EQ(acb, -expected);
        BOOST_TEST_EQ(bca, expected);
    }
}

void test_constructed_forward_predicates()
{
    ec::point_cache_scope<double> operation;
    auto const a = third<double>();
    auto const b = ec::construct_intersection(point<double>(0, 0), point<double>(2, 2),
        point<double>(0, 1), point<double>(1, 0));
    auto const c = ec::construct_intersection(point<double>(0, 1), point<double>(1, 0),
        point<double>(0, 0), point<double>(1, 2));
    allocations = 0;
    count_allocations = true;
    int const orientation = ec::orientation(a, b, c);
    int const x = ec::compare_x(a, b), y = ec::compare_y(a, b);
    bool const equal = ec::equal_points(a, c);
    count_allocations = false;
    BOOST_TEST_EQ(allocations, 0u);
    BOOST_TEST_EQ(orientation, 1);
    BOOST_TEST_EQ(x, -1);
    BOOST_TEST_EQ(y, -1);
    BOOST_TEST(!equal);
}

point<double> literal(rational const& value)
{
    point<double> result;
    bg::set<0>(result, value);
    bg::set<1>(result, value);
    return result;
}

template <typename Output>
Output round_output(rational const& value)
{
    input_point<Output> output;
    ec::assign_output<nearest_output>(literal(value), output);
    BOOST_TEST_EQ(bg::get<0>(output), bg::get<1>(output));
    return bg::get<0>(output);
}

template <typename Output>
void test_binary_rounding()
{
    ec::point_cache_scope<double> operation;
    constexpr int digits = std::numeric_limits<Output>::digits;
    rational const ulp = rational(1) / (integer(1) << (digits - 1));
    rational const half = ulp / 2;
    rational const small = ulp / 1024;
    Output const next = std::nextafter(Output(1), Output(2));
    Output const next_next = std::nextafter(next, Output(2));
    BOOST_TEST_EQ(round_output<Output>(1 + half), Output(1));
    BOOST_TEST_EQ(round_output<Output>(1 + half + small), next);
    BOOST_TEST_EQ(round_output<Output>(1 + half - small), Output(1));
    BOOST_TEST_EQ(round_output<Output>(1 + ulp + half), next_next);
    BOOST_TEST_EQ(round_output<Output>(-1 - half), Output(-1));
    BOOST_TEST_EQ(round_output<Output>(-1 - ulp - half), -next_next);

    // Exercise binade boundaries without floating underflow or overflow.
    rational const binade = ec::exact_value(std::ldexp(Output(1), -32));
    rational const binade_step = binade * ulp / 2;
    Output const lower = std::nextafter(std::ldexp(Output(1), -32), Output(0));
    BOOST_TEST_EQ(round_output<Output>(binade - binade_step / 2), std::ldexp(Output(1), -32));
    BOOST_TEST_EQ(round_output<Output>(binade - binade_step * 3 / 4), lower);
    BOOST_TEST_EQ(round_output<Output>(-binade + binade_step * 3 / 4), -lower);
}

template <typename Output, typename OutputPolicy, typename Coordinate>
void check_integral_output(point<Coordinate> const& p, rational const& x, rational const& y)
{
    BOOST_TEST(p.exact_coordinate(0) == x);
    BOOST_TEST(p.exact_coordinate(1) == y);
    input_point<Output> output(17, 23);
    auto convert = [&](std::size_t dimension)
    {
        return p.is_input() ? bg::util::numeric_cast<Output>(p.input_coordinate(dimension))
                            : bg::util::numeric_cast<Output>(p.output_approximation(dimension));
    };
    Output expected_x, expected_y;
    try
    {
        expected_x = convert(0);
        expected_y = convert(1);
    }
    catch (boost::numeric::positive_overflow const&)
    {
        BOOST_TEST_THROWS(ec::assign_output<OutputPolicy>(p, output), boost::numeric::positive_overflow);
        BOOST_TEST_EQ(bg::get<0>(output), Output(17));
        BOOST_TEST_EQ(bg::get<1>(output), Output(23));
        return;
    }
    catch (boost::numeric::negative_overflow const&)
    {
        BOOST_TEST_THROWS(ec::assign_output<OutputPolicy>(p, output), boost::numeric::negative_overflow);
        BOOST_TEST_EQ(bg::get<0>(output), Output(17));
        BOOST_TEST_EQ(bg::get<1>(output), Output(23));
        return;
    }
    auto const before = ec::filtering_statistics();
    allocations = 0;
    count_allocations = true;
    ec::assign_output<OutputPolicy>(p, output);
    count_allocations = false;
    BOOST_TEST_EQ(allocations, 0u);
    BOOST_TEST_EQ(bg::get<0>(output), expected_x);
    BOOST_TEST_EQ(bg::get<1>(output), expected_y);
    BOOST_TEST_EQ(ec::filtering_statistics().construction_filters, before.construction_filters);
    BOOST_TEST_EQ(ec::filtering_statistics().point_coordinate_values, before.point_coordinate_values);
    BOOST_TEST_EQ(ec::filtering_statistics().homogeneous_values, before.homogeneous_values);
}

template <typename OutputPolicy>
void test_integral_conversion()
{
    ec::point_cache_scope<double> operation;
    std::array<rational, 6> const values{{rational(1) / 3, rational(1) / 2, rational(3) / 2,
                                       rational(5) / 2, rational(-3) / 2, rational(-5) / 2}};
    for (auto const& value : values)
        check_integral_output<int, OutputPolicy>(literal(value), value, value);
    auto const high = (std::numeric_limits<std::int64_t>::max)();
    auto const low = (std::numeric_limits<std::int64_t>::min)();
    // Preserve the rational fixtures. The stored binary64 approximation of
    // INT64_MAX is out of range; this is an ordinary numeric_cast failure.
    check_integral_output<std::int64_t, OutputPolicy>(literal(rational(high)), rational(high), rational(high));
    check_integral_output<std::int64_t, OutputPolicy>(literal(rational(low)), rational(low), rational(low));
    rational const above_int = rational((std::numeric_limits<int>::max)()) + 1;
    check_integral_output<int, OutputPolicy>(literal(above_int), above_int, above_int);
    check_integral_output<unsigned, OutputPolicy>(literal(rational(-1)), -1, -1);
    auto p = literal(rational(1));
    bg::set<1>(p, rational(integer(1) << 100));
    check_integral_output<int, OutputPolicy>(p, rational(1), rational(integer(1) << 100));

    // Native leaves retain their original full-width coordinates, rather than
    // being converted through the approximation used for constructed points.
    point<std::int64_t> const native(high, low);
    check_integral_output<std::int64_t, OutputPolicy>(native, rational(high), rational(low));
    auto const maximum = (std::numeric_limits<std::uint64_t>::max)();
    point<std::uint64_t> const unsigned_native(maximum, 0);
    check_integral_output<std::uint64_t, OutputPolicy>(unsigned_native, rational(maximum), rational(0));
}

void test_concurrent_fallback()
{
    // Each thread owns a complete operation, including all working points.
    std::atomic<bool> ready{false};
    std::atomic<bool> okay{true};
    std::vector<std::thread> workers;
    for (int worker = 0; worker < 8; ++worker)
    {
        workers.emplace_back([&]
        {
            while (!ready.load(std::memory_order_acquire))
            {
                std::this_thread::yield();
            }
            try
            {
                ec::point_cache_scope<double> operation;
                auto const p = third<double>();
                auto const q = ec::construct_intersection(point<double>(0, 0), point<double>(2, 2),
                    point<double>(0, 2), point<double>(1, -3));
                for (int iteration = 0; iteration < 40; ++iteration)
                {
                    if (!ec::equal_points(p, q)
                        || ec::orientation(point<double>(0, 1), point<double>(1, -1), p) != 0
                        || p.exact_coordinate(0) != rational(1) / 3)
                    {
                        okay.store(false, std::memory_order_relaxed);
                    }
                }
            }
            catch (...)
            {
                okay.store(false, std::memory_order_relaxed);
            }
        });
    }
    ready.store(true, std::memory_order_release);
    for (auto& worker : workers)
    {
        worker.join();
    }
    BOOST_TEST(okay.load());
}

template <typename Coordinate>
void test_inline_native_points()
{
    ec::point_cache_scope<Coordinate> operation;
    ec::operation_detail::require<Coordinate>();
    prepare_exact_cache_repeat<Coordinate>();
    allocations = 0;
    count_allocations = true;
    point<Coordinate> const a(0, 0), b(3, 3), c(0, 1), d(2, 0);
    auto const p = ec::construct_intersection(a, b, c, d);
    // Both nested intersections reuse a native supporting line, not a rounded
    // coordinate of an earlier intersection.
    auto const q = ec::construct_intersection(p, b, point<Coordinate>(0, 1), point<Coordinate>(2, 1));
    auto const r = ec::construct_intersection(q, p, point<Coordinate>(2, 0), point<Coordinate>(2, 3));
    auto copy = r;
    auto const h = copy.fixed_homogeneous();
    int sign = 9;
    bool const fixed = point<Coordinate>::try_fixed_orientation(a, p, point<Coordinate>(0, 1), sign);
    bool const equal = ec::equal_points(copy, point<Coordinate>(2, 2));
    std::array<Coordinate, 8> producers;
    bool const retained_native = r.native_producers(producers, std::true_type());
    count_allocations = false;
    BOOST_TEST_EQ(allocations, 0u);
    BOOST_TEST(fixed);
    BOOST_TEST_EQ(sign, 1);
    BOOST_TEST(equal);
    BOOST_TEST(retained_native);
    BOOST_TEST(p.has_fixed_representation());
    BOOST_TEST(q.has_fixed_representation());
    BOOST_TEST(r.has_fixed_representation());
    BOOST_TEST(h.X == 2 * h.W && h.Y == 2 * h.W);
    BOOST_TEST(!copy.shares_node(p));
    BOOST_TEST(copy.shares_node(r));
}

void test_indexed_point_cache()
{
    using cache = ec::point_cache_detail::native_cache<double>;
    BOOST_TEST(ec::operation_detail::current<double>() == nullptr);
    {
        ec::point_cache_scope<double> outer;
        auto* const root = ec::operation_detail::active<double>();
        BOOST_TEST(ec::operation_detail::current<double>() == nullptr);
        {
            ec::point_cache_scope<double> inner;
            BOOST_TEST(ec::operation_detail::active<double>() == root);
            BOOST_TEST(ec::operation_detail::current<double>() == nullptr);
            third<double>().fixed_homogeneous();
        }
        BOOST_TEST(ec::operation_detail::active<double>() == root);
        BOOST_TEST(ec::operation_detail::current<double>() != nullptr);
        BOOST_TEST_EQ(cache::size(), 1u);
    }
    BOOST_TEST(ec::operation_detail::active<double>() == nullptr);
    {
        ec::point_cache_scope<double> operation;
        BOOST_TEST(ec::operation_detail::current<double>() == nullptr);
        auto const p = third<double>();
        auto* const owner = &ec::operation_detail::require<double>();
        BOOST_TEST_EQ(cache::size(), 0u);
        auto const first = p.fixed_homogeneous();
        BOOST_TEST_EQ(cache::size(), 1u);
        auto const copy = p;
        copy.fixed_homogeneous();
        allocations = 0;
        count_allocations = true;
        auto const cached = copy.fixed_homogeneous();
        count_allocations = false;
        BOOST_TEST_EQ(allocations, 0u);
        BOOST_TEST(cached.X == first.X && cached.Y == first.Y && cached.W == first.W);
        BOOST_TEST_EQ(cache::size(), 1u);
        {
            ec::point_cache_scope<double> nested;
            BOOST_TEST(ec::operation_detail::current<double>() == owner);
            std::array<point<double>, 40> points;
            for (std::size_t i = 0; i != points.size(); ++i)
            {
                double const x = static_cast<double>(i + 2);
                points[i] = ec::construct_intersection(point<double>(0, 0), point<double>(x, x),
                    point<double>(0, x), point<double>(x, 0));
                auto const value = points[i].fixed_homogeneous();
                BOOST_TEST(value.X * 2 == value.W * (i + 2));
            }
            BOOST_TEST_EQ(cache::size(), 41u);
        }
        BOOST_TEST(ec::operation_detail::current<double>() == owner);
        auto const after_growth = copy.fixed_homogeneous();
        BOOST_TEST(after_growth.X == first.X && after_growth.W == first.W);
        BOOST_TEST_EQ(cache::size(), 41u);
    }
    BOOST_TEST(ec::operation_detail::current<double>() == nullptr);
    {
        ec::point_cache_scope<double> next_operation;
        BOOST_TEST(ec::operation_detail::current<double>() == nullptr);
        BOOST_TEST_EQ(cache::size(), 0u);
        auto const p = third<double>();
        auto const value = p.fixed_homogeneous();
        BOOST_TEST(value.X * 3 == value.W);
        BOOST_TEST_EQ(cache::size(), 1u);
    }
    BOOST_TEST(ec::operation_detail::current<double>() == nullptr);
}

template <typename Coordinate>
void test_native_fixed_extremes()
{
    ec::point_cache_scope<Coordinate> operation;
    ec::operation_detail::require<Coordinate>();
    prepare_exact_cache_repeat<Coordinate>();
    // Full binary64 extrema overflow floating differences/products, outside D6.
    Coordinate const low = std::is_same<Coordinate, double>::value
        ? Coordinate(-std::ldexp(1.0, 200)) : (std::numeric_limits<Coordinate>::lowest)();
    Coordinate const high = std::is_same<Coordinate, double>::value
        ? Coordinate(std::ldexp(1.0, 200)) : (std::numeric_limits<Coordinate>::max)();
    allocations = 0;
    count_allocations = true;
    point<Coordinate> const a(low, low), b(high, high), c(low, high), d(high, low);
    auto const p = ec::construct_intersection(a, b, c, d);
    auto const h = p.fixed_homogeneous();
    int sign = 9;
    bool const fixed = point<Coordinate>::try_fixed_orientation(a, b, p, sign);
    int const from_low = ec::compare_x(p, a), from_high = ec::compare_y(p, b);
    count_allocations = false;
    BOOST_TEST_EQ(allocations, 0u);
    BOOST_TEST(h.X == h.Y && h.W > 0);
    BOOST_TEST(fixed);
    BOOST_TEST_EQ(sign, 0);
    BOOST_TEST_EQ(from_low, 1);
    BOOST_TEST_EQ(from_high, -1);
    rational const expected = (ec::exact_value(low) + ec::exact_value(high)) / 2;
    BOOST_TEST(p.exact_coordinate(0) == expected);
    BOOST_TEST(p.exact_coordinate(1) == expected);
}

void test_dense_fixed_orientation_no_allocations()
{
    ec::point_cache_scope<double> operation;
    // Dense bounded diagnostic literals stress fixed arithmetic without relying
    // on the former fixture's floating underflow/overflow during construction.
    integer const numerator = integer(1) << 4096;
    integer const denominator = integer(1) << 4095;
    auto const a = point<double>::literal(rational(numerator - 1) / (denominator + 17),
        rational(numerator - 5) / (denominator + 17));
    auto const b = point<double>::literal(rational(numerator - 7) / (denominator + 19),
        rational(numerator - 11) / (denominator + 19));
    auto const c = point<double>::literal(rational(numerator - 13) / (denominator + 23),
        rational(numerator - 17) / (denominator + 23));
    for (auto const* p : {&a, &b, &c})
    {
        auto const fixed = p->precision_homogeneous();
        std::size_t maximum_width = 0, maximum_density = 0;
        for (integer value : {integer(fixed.X), integer(fixed.Y)})
        {
            if (value < 0) { value = -value; }
            if (value == 0) { continue; }
            std::size_t const width = boost::multiprecision::msb(value) + 1;
            std::size_t density = 0;
            for (std::size_t bit = 0; bit != width; ++bit)
            {
                density += boost::multiprecision::bit_test(value, static_cast<unsigned>(bit));
            }
            maximum_width = (std::max)(maximum_width, width);
            maximum_density = (std::max)(maximum_density, density);
        }
        // These operands exceed Karatsuba's scratch-allocation threshold;
        // unlike small dyadics, their dense limbs exercise the full backend.
        BOOST_TEST(maximum_width > 3000);
        BOOST_TEST(maximum_density > 1000);
    }
    auto const ah = a.homogeneous(), bh = b.homogeneous(), ch = c.homogeneous();
    auto const line = ec::point_detail::cross(ah, bh);
    integer const determinant = line.X * ch.X + line.Y * ch.Y + line.W * ch.W;
    int const expected = determinant.sign();
    BOOST_TEST(expected != 0);
    int abc = 7, acb = 7, bca = 7;
    allocations = 0;
    count_allocations = true;
    bool const first = point<double>::try_fixed_orientation(a, b, c, abc);
    bool const second = point<double>::try_fixed_orientation(a, c, b, acb);
    bool const third = point<double>::try_fixed_orientation(b, c, a, bca);
    count_allocations = false;
    BOOST_TEST_EQ(allocations, 0u);
    BOOST_TEST(first && second && third);
    BOOST_TEST_EQ(abc, expected);
    BOOST_TEST_EQ(acb, -expected);
    BOOST_TEST_EQ(bca, expected);
}

void test_expansion_refinement_cache()
{
    using cache = ec::point_cache_detail::native_cache<double>;
    ec::point_cache_scope<double> scope;
    auto const p = third<double>();
    auto const initial = p.approximation(0);
    BOOST_TEST_EQ(cache::size(), 0u);
#ifdef BOOST_GEOMETRY_EXACT_FILTER_STATISTICS
    auto const filters = ec::filtering_statistics().construction_filters;
    auto const coordinates = ec::filtering_statistics().point_coordinate_values;
    auto const homogeneous = ec::filtering_statistics().homogeneous_values;
#endif
    auto const x = p.refined_approximation(0);
    auto const y = p.refined_approximation(1);
    BOOST_TEST(x.error < initial.error);
    BOOST_TEST_EQ(cache::size(), 1u);
#ifdef BOOST_GEOMETRY_EXACT_FILTER_STATISTICS
    BOOST_TEST_EQ(ec::filtering_statistics().construction_filters, filters + 2);
#endif
    auto const copy = p;
    copy.refined_approximation(0);
    copy.refined_approximation(1);
    allocations = 0;
    count_allocations = true;
    auto const cached_x = copy.refined_approximation(0);
    auto const cached_y = copy.refined_approximation(1);
    auto const comparison = point<double>::compare_native(copy, 0, 1.0 / 3.0);
    count_allocations = false;
    BOOST_TEST_EQ(allocations, 0u);
    BOOST_TEST_EQ(cached_x.approximate, x.approximate);
    BOOST_TEST_EQ(cached_x.error, x.error);
    BOOST_TEST_EQ(cached_y.approximate, y.approximate);
    BOOST_TEST_EQ(cached_y.error, y.error);
    BOOST_TEST(comparison.certain);
    BOOST_TEST_EQ(comparison.value, 1);
    BOOST_TEST_EQ(cache::size(), 1u);
    BOOST_TEST_EQ(p.approximation(0).approximate, x.approximate);
    BOOST_TEST_EQ(p.approximation(0).error, x.error);
#ifdef BOOST_GEOMETRY_EXACT_FILTER_STATISTICS
    BOOST_TEST_EQ(ec::filtering_statistics().construction_filters, filters + 2);
    BOOST_TEST_EQ(ec::filtering_statistics().point_coordinate_values, coordinates);
    BOOST_TEST_EQ(ec::filtering_statistics().homogeneous_values, homogeneous);
#endif
    auto const other = ec::construct_intersection(point<double>(0, 0), point<double>(2, 2),
        point<double>(0, 2), point<double>(1, -3));
    auto const equality = point<double>::compare_points(p, other, 0);
    BOOST_TEST(equality.certain);
    BOOST_TEST_EQ(equality.value, 0);
    BOOST_TEST_EQ(cache::size(), 2u);
    auto const fixed = p.fixed_homogeneous();
    BOOST_TEST(fixed.X * 3 == fixed.W);
    BOOST_TEST_EQ(cache::size(), 2u);
}

void test_cached_approximation_copies()
{
    using cache = ec::point_cache_detail::native_cache<double>;
    ec::point_cache_scope<double> scope;
    auto const p = third<double>();
    auto const original = p.approximation(0);
    BOOST_TEST_EQ(cache::size(), 0u);
    auto const refined = p.refined_approximation(0);
    BOOST_TEST(refined.error < original.error);
    auto const copy = p;
    {
        ec::point_cache_scope<double> nested;
        allocations = 0;
        count_allocations = true;
        auto const reused = p.approximation(0);
        auto const copied = copy.approximation(0);
        count_allocations = false;
        BOOST_TEST_EQ(allocations, 0u);
        BOOST_TEST_EQ(reused.approximate, refined.approximate);
        BOOST_TEST_EQ(reused.error, refined.error);
        BOOST_TEST_EQ(copied.approximate, refined.approximate);
        BOOST_TEST_EQ(copied.error, refined.error);
        BOOST_TEST_EQ(cache::size(), 1u);
    }
    BOOST_TEST_EQ(p.approximation(0).approximate, refined.approximate);
    BOOST_TEST_EQ(p.approximation(0).error, refined.error);
    BOOST_TEST_EQ(copy.approximation(0).error, refined.error);
    BOOST_TEST_EQ(cache::size(), 1u);
}

void test_fixed_cache_bounds()
{
    ec::point_cache_scope<std::int64_t> integral_scope;
    ec::point_cache_scope<double> floating_scope;
    auto const high = (std::numeric_limits<std::int64_t>::max)();
    auto const p = ec::construct_intersection(point<std::int64_t>(high - 4, 0),
        point<std::int64_t>(high, 3), point<std::int64_t>(high - 4, 3), point<std::int64_t>(high, 0));
    std::array<double, 8> values;
    BOOST_TEST(!p.native_producers(values));
    auto const original_y = p.approximation(1);
    p.fixed_homogeneous();
    auto const x = p.approximation(0), y = p.approximation(1);
    BOOST_TEST(x.valid() && y.valid());
    // Unrepresentable native int64 endpoints may trigger refinement already
    // during construction, so this lookup need not improve the bound again.
    BOOST_TEST(!original_y.valid() || y.error <= original_y.error);
    BOOST_TEST_EQ(y.approximate, 1.5);
    BOOST_TEST_EQ(p.refined_approximation(1).error, y.error);
    BOOST_TEST(ec::exact_value(x.lower()) <= rational(high - 2));
    BOOST_TEST(rational(high - 2) <= ec::exact_value(x.upper()));
    BOOST_TEST(ec::exact_value(y.lower()) <= rational(3) / 2);
    BOOST_TEST(rational(3) / 2 <= ec::exact_value(y.upper()));
    auto const f = third<double>();
    auto const expanded = f.refined_approximation(0);
    f.fixed_homogeneous();
    BOOST_TEST_EQ(f.approximation(0).error, expanded.error);
    double const small = std::ldexp(1.0, -100);
    double const large = std::ldexp(1.0, 100);
    auto const cancellation = ec::construct_intersection(point<double>(small, large), point<double>(large, -small),
        point<double>(-large, small), point<double>(small, -large));
    cancellation.fixed_homogeneous();
    BOOST_TEST(cancellation.approximation(0).valid());
    BOOST_TEST(cancellation.approximation(1).valid());
    check_bounds(cancellation);
}

template <typename Output>
void check_retained_quotient(std::array<double, 8> const& values, std::size_t dimension)
{
    ec::construction_detail::expansion<32> numerator;
    ec::construction_detail::expansion<16> denominator;
    BOOST_TEST(ec::construction_detail::homogeneous_coordinate(values, dimension, numerator, denominator));
    Output retained = Output(17), rebuilt = Output(17);
    allocations = 0;
    count_allocations = true;
    bool const certified = ec::try_round_quotient(numerator, denominator, retained);
    bool const original = ec::try_round_intersection(values, dimension, rebuilt);
    count_allocations = false;
    BOOST_TEST_EQ(allocations, 0u);
    BOOST_TEST(certified);
    BOOST_TEST_EQ(certified, original);
    BOOST_TEST(retained == rebuilt);
    auto const p = ec::construct_intersection(point<double>(values[0], values[1]),
        point<double>(values[2], values[3]), point<double>(values[4], values[5]),
        point<double>(values[6], values[7]));
    BOOST_TEST(retained == ec::point_detail::output_rounder<Output>::apply(p.exact_coordinate(dimension)));
}

void test_retained_rounding_quotient()
{
    ec::point_cache_scope<double> scope;
    std::array<std::array<double, 8>, 4> const values{{
        {{0, 0, 1, 1, 0, 1, 1, -1}},
        {{0, 0, -1, -1, 0, -1, -1, 1}},
        {{0, 0, 0, 1, -1, 0, 2, 1}},
        {{0, 0, 1, 1, 0, 1, 1, std::nextafter(2.0, 3.0)}}}};
    for (auto const& input : values)
    {
        for (std::size_t dimension = 0; dimension != 2; ++dimension)
        {
            check_retained_quotient<double>(input, dimension);
            check_retained_quotient<float>(input, dimension);
            check_retained_quotient<std::int64_t>(input, dimension);
        }
    }
    ec::construction_detail::expansion<32> numerator;
    ec::construction_detail::expansion<16> denominator;
    BOOST_TEST(ec::construction_detail::homogeneous_coordinate(values[0], 0, numerator, denominator));
    double control = 17;
    BOOST_TEST(ec::try_round_quotient(numerator, denominator, control));
    BOOST_TEST_EQ(control, 1.0 / 3.0);
    double rounded = 17;
    ec::construction_detail::expansion<16> const zero_denominator;
    BOOST_TEST(!ec::try_round_quotient(numerator, zero_denominator, rounded));
    BOOST_TEST_EQ(rounded, 17);
    numerator.valid = false;
    BOOST_TEST(!ec::try_round_quotient(numerator, denominator, rounded));
    BOOST_TEST_EQ(rounded, 17);
    long double unsupported = 17;
    BOOST_TEST(!ec::try_round_quotient(numerator, denominator, unsupported));
    BOOST_TEST(unsupported == 17);
}

void test_uncertain_initial_denominator()
{
    using cache = ec::point_cache_detail::native_cache<double>;
    ec::point_cache_scope<double> scope;
    double const above_two = std::nextafter(2.0, (std::numeric_limits<double>::infinity)());
    std::array<double, 8> const values{{0, 0, 1, 1, 0, 1, 1, above_two}};
    ec::forward_error::point_approximation approximation;
    BOOST_TEST(!ec::forward_error::try_intersection(values, approximation));
#ifdef BOOST_GEOMETRY_EXACT_FILTER_STATISTICS
    auto const coordinates = ec::filtering_statistics().point_coordinate_values;
    auto const homogeneous = ec::filtering_statistics().homogeneous_values;
#endif
    auto const p = ec::construct_intersection(point<double>(0, 0), point<double>(1, 1),
        point<double>(0, 1), point<double>(1, above_two));
    BOOST_TEST(p.approximation(0).valid());
    BOOST_TEST(p.approximation(1).valid());
    BOOST_TEST_EQ(p.approximation(0).approximate, -std::ldexp(1.0, 51));
    BOOST_TEST_EQ(cache::size(), 1u);
    auto const copy = p;
    copy.refined_approximation(0);
    allocations = 0;
    count_allocations = true;
    auto const refined = copy.refined_approximation(0);
    count_allocations = false;
    BOOST_TEST_EQ(allocations, 0u);
    BOOST_TEST_EQ(refined.approximate, p.approximation(0).approximate);
    BOOST_TEST_EQ(refined.error, p.approximation(0).error);
#ifdef BOOST_GEOMETRY_EXACT_FILTER_STATISTICS
    BOOST_TEST_EQ(ec::filtering_statistics().point_coordinate_values, coordinates);
    BOOST_TEST_EQ(ec::filtering_statistics().homogeneous_values, homogeneous);
#endif
    check_bounds(p);
    auto parallel = []
    {
        return ec::construct_intersection(point<double>(0, 0), point<double>(1, 1),
            point<double>(0, 1), point<double>(1, 2));
    };
    BOOST_TEST_THROWS(parallel(), std::domain_error);
    BOOST_TEST_THROWS(parallel(), std::domain_error);
}

void test_same_operation_native_producers()
{
    ec::point_cache_scope<double> operation;
    point<double> const a(0, 0), b(3, 3), c(0, 1), d(1, -1);
    auto const p = ec::construct_intersection(a, b, c, d);
    p.fixed_homogeneous();
    auto const copy = p;
    std::array<double, 8> producers;
    BOOST_TEST(copy.native_producers(producers));
    BOOST_TEST((producers == std::array<double, 8>{{0, 0, 3, 3, 0, 1, 1, -1}}));
    auto const before = ec::point_cache_detail::native_cache<double>::size();
    {
        ec::point_cache_scope<double> nested;
        BOOST_TEST(copy.exact_coordinate(0) == rational(1) / 3);
        BOOST_TEST(copy.exact_coordinate(1) == rational(1) / 3);
        BOOST_TEST_EQ(ec::orientation(a, b, copy), 0);
        BOOST_TEST_EQ(ec::orientation(c, d, copy), 0);
        check_bounds(copy);
    }
    BOOST_TEST_EQ(ec::point_cache_detail::native_cache<double>::size(), before);
    BOOST_TEST(ec::equal_points(p, copy));
}

void test_independent_point_cache_concurrency()
{
    std::atomic<bool> ready{false}, okay{true};
    std::vector<std::thread> workers;
    for (int worker = 0; worker != 8; ++worker)
    {
        workers.emplace_back([&]
        {
            while (!ready.load(std::memory_order_acquire)) { std::this_thread::yield(); }
            try
            {
                for (int iteration = 0; iteration != 20; ++iteration)
                {
                    ec::point_cache_scope<double> scope;
                    auto const p = third<double>();
                    auto const early_copy = p;
                    auto const h = p.fixed_homogeneous();
                    auto const copy = p;
                    auto const other = copy.fixed_homogeneous();
                    auto const early = early_copy.fixed_homogeneous();
                    auto const refined = copy.refined_approximation(0);
                    auto const current = p.approximation(0);
                    if (h.X * 3 != h.W || other.X != h.X || other.W != h.W
                        || early.X != h.X || early.W != h.W || !refined.valid() || !current.valid()
                        || p.exact_coordinate(1) != rational(1) / 3)
                    {
                        okay.store(false, std::memory_order_relaxed);
                    }
                }
            }
            catch (...) { okay.store(false, std::memory_order_relaxed); }
        });
    }
    ready.store(true, std::memory_order_release);
    for (auto& worker : workers) { worker.join(); }
    BOOST_TEST(okay.load());
}

void test_constructed_round_to_nearest()
{
    ec::point_cache_scope<double> double_scope;
    ec::point_cache_scope<float> float_scope;
    BOOST_TEST_EQ(std::fegetround(), FE_TONEAREST);
    auto const p = third<double>();
    auto const f = third<float>();
    check_bounds(p);
    check_bounds(f);
    BOOST_TEST_EQ(ec::compare_x(p, point<double>(0, 0)), 1);
    BOOST_TEST_EQ(ec::compare_x(p, point<double>(1, 1)), -1);
    BOOST_TEST_EQ(ec::orientation(point<double>(0, 1), point<double>(1, -1), p), 0);
    BOOST_TEST_EQ(ec::orientation(point<float>(0, 1), point<float>(1, -1), f), 0);
}

template <typename Output, typename Coordinate>
void check_native_output(point<Coordinate> const& p, rational const& x, rational const& y)
{
    ec::point_cache_scope<Coordinate> scope;
    if constexpr (std::is_integral<Output>::value)
    {
        check_integral_output<Output, approximate_output>(p, x, y);
        check_integral_output<Output, nearest_output>(p, x, y);
        return;
    }
    Output const expected_x = ec::point_detail::output_rounder<Output>::apply(x);
    Output const expected_y = ec::point_detail::output_rounder<Output>::apply(y);
    input_point<Output> output;
#ifdef BOOST_GEOMETRY_EXACT_FILTER_STATISTICS
    auto const coordinates = ec::filtering_statistics().point_coordinate_values;
    auto const homogeneous = ec::filtering_statistics().homogeneous_values;
    auto const native = ec::filtering_statistics().native_exact_values;
#endif
    prepare_exact_cache_repeat<Coordinate>();
    allocations = 0;
    count_allocations = true;
    ec::assign_output<nearest_output>(p, output);
    count_allocations = false;
    BOOST_TEST_EQ(allocations, 0u);
    BOOST_TEST(bg::get<0>(output) == expected_x);
    BOOST_TEST(bg::get<1>(output) == expected_y);
    if (expected_x == Output(0))
    {
        BOOST_TEST_EQ(std::signbit(double(bg::get<0>(output))), std::signbit(double(expected_x)));
    }
    if (expected_y == Output(0))
    {
        BOOST_TEST_EQ(std::signbit(double(bg::get<1>(output))), std::signbit(double(expected_y)));
    }
#ifdef BOOST_GEOMETRY_EXACT_FILTER_STATISTICS
    BOOST_TEST_EQ(ec::filtering_statistics().point_coordinate_values, coordinates);
    BOOST_TEST_EQ(ec::filtering_statistics().homogeneous_values, homogeneous);
    BOOST_TEST_EQ(ec::filtering_statistics().native_exact_values, native);
#endif
}

void test_fixed_native_output()
{
    ec::point_cache_scope<double> double_scope;
    ec::point_cache_scope<std::int64_t> integral_scope;
    ec::point_cache_scope<std::uint64_t> unsigned_scope;
    auto const p = third<double>();
    check_native_output<double>(p, rational(1) / 3, rational(1) / 3);
    check_native_output<float>(p, rational(1) / 3, rational(1) / 3);
    check_native_output<int>(p, rational(1) / 3, rational(1) / 3);
    auto const high = (std::numeric_limits<std::int64_t>::max)();
    auto const low = (std::numeric_limits<std::int64_t>::min)();
    check_native_output<double>(point<std::int64_t>(9007199254740993LL, low),
        rational(9007199254740993LL), rational(low));
    auto const maximum = (std::numeric_limits<std::uint64_t>::max)();
    check_native_output<double>(point<std::uint64_t>(maximum, 0), rational(maximum), rational(0));
    auto const crossing = ec::construct_intersection(point<std::int64_t>(high - 2, 0),
        point<std::int64_t>(high, 2), point<std::int64_t>(high - 2, 2), point<std::int64_t>(high, 0));
    check_native_output<double>(crossing, rational(high - 1), rational(1));
    check_native_output<std::int64_t>(crossing, rational(high - 1), rational(1));
    for (int exponent : {100, -100})
    {
        double const scale = std::ldexp(1.0, exponent);
        auto const large = ec::construct_intersection(point<double>(0, 0), point<double>(scale, scale),
            point<double>(0, scale), point<double>(scale, -scale));
        rational const expected = ec::exact_value(scale) / 3;
        check_native_output<double>(large, expected, expected);
    }
    input_point<int> unchanged(17, 23);
    BOOST_TEST_THROWS(ec::assign_output(point<std::int64_t>(7, high), unchanged), boost::numeric::positive_overflow);
    BOOST_TEST_EQ(bg::get<0>(unchanged), 17);
    BOOST_TEST_EQ(bg::get<1>(unchanged), 23);
    input_point<rational> exact;
    ec::assign_output(crossing, exact);
    BOOST_TEST(bg::get<0>(exact) == rational(high - 1));
    BOOST_TEST(bg::get<1>(exact) == 1);
}

template <typename Coordinate>
void test_materialized_precision_cache()
{
    using p = point<Coordinate>;
    ec::point_cache_scope<Coordinate> owner;
    // This diagnostic construction also splices a derived rational point.
    // Keep that expression within this coordinate type's auxiliary capacity.
    constexpr std::size_t bits = (std::min)(std::size_t(4096), ec::precision_denominator_bits<Coordinate> / 2 - 8);
    rational const dense = rational((integer(1) << bits) - 1) / ((integer(1) << (bits - 1)) + 17);
    rational const third_value = rational(1) / 3;
    {
        // Reserve reusable precision slots; all later arithmetic must allocate
        // neither temporary Q values nor separate ownership records.
        std::array<p, 8> warm;
        for (auto& value : warm) { value = p::literal(dense, third_value); }
    }
    input_point<double> output;
    prepare_exact_cache_repeat<Coordinate>();
    allocations = 0;
    count_allocations = true;
    auto const original = p::literal(dense, third_value);
    auto changed = original;
    changed.set_exact_coordinate(1, dense);
    auto const derived = ec::construct_intersection(original, p(1, 0), p(0, -1), p(0, 1));
    auto const spliced = p::splice_coordinates(original, 0, derived, 1);
    bool const different = !ec::equal_points(original, changed);
    int const collinear = ec::orientation(original, p(1, 0), derived);
    ec::assign_output<nearest_output>(spliced, output);
    count_allocations = false;
    BOOST_TEST_EQ(allocations, 0u);
    BOOST_TEST(different);
    BOOST_TEST_EQ(collinear, 0);
    BOOST_TEST(original.exact_coordinate(1) == third_value);
    BOOST_TEST(changed.exact_coordinate(1) == dense);
    BOOST_TEST(derived.exact_coordinate(0) == 0);
    BOOST_TEST(spliced.exact_coordinate(0) == dense);
    BOOST_TEST(spliced.exact_coordinate(1) == derived.exact_coordinate(1));
    BOOST_TEST_EQ(bg::get<0>(output), ec::point_detail::output_rounder<double>::apply(dense));

    p shared;
    {
        ec::point_cache_scope<Coordinate> scope;
        auto temporary = p::literal(dense, third_value);
        shared = temporary;
        temporary.set_exact_coordinate(0, third_value);
        BOOST_TEST(shared.exact_coordinate(0) == dense);
        BOOST_TEST(temporary.exact_coordinate(0) == third_value);
    }
    BOOST_TEST(shared.exact_coordinate(0) == dense);
    auto const identity = ec::make_point<Coordinate>(shared);
    BOOST_TEST(identity.shares_node(shared));
    BOOST_TEST(identity.exact_coordinate(0) == dense);
    std::atomic<bool> okay{true};
    std::vector<std::thread> workers;
    for (int i = 0; i != 4; ++i)
    {
        workers.emplace_back([&]
        {
            try
            {
                ec::point_cache_scope<Coordinate> operation;
                auto const local = p::literal(dense, third_value);
                for (int iteration = 0; iteration != 8; ++iteration)
                {
                    auto const copy = local;
                    auto const restored = p::from_homogeneous(copy.precision_homogeneous());
                    if (ec::compare_x(copy, restored) != 0 || ec::compare_y(copy, restored) != 0)
                    {
                        okay.store(false, std::memory_order_relaxed);
                    }
                }
            }
            catch (...) { okay.store(false, std::memory_order_relaxed); }
        });
    }
    for (auto& worker : workers) { worker.join(); }
    BOOST_TEST(okay.load());
    auto const native_corner = point<std::int64_t>::splice_coordinates(
        point<std::int64_t>(9007199254740993LL, 0), 0, point<std::int64_t>(0, -7), 1);
    BOOST_TEST(native_corner.is_input());
    BOOST_TEST_EQ(native_corner.input_coordinate(0), 9007199254740993LL);
    BOOST_TEST_EQ(native_corner.input_coordinate(1), -7);

    rational const excessive(integer(1) << (ec::precision_numerator_bits<Coordinate> + 1));
    BOOST_TEST_THROWS(p::literal(excessive, third_value), std::overflow_error);
    ec::fixed_homogeneous<ec::precision_numerator_bits<Coordinate> + 2> wide;
    wide.X = 1;
    wide.X <<= ec::precision_numerator_bits<Coordinate> + 1;
    BOOST_TEST_THROWS(p::from_homogeneous(wide), std::overflow_error);
    if constexpr (std::is_same<Coordinate, float>::value)
    {
        // The old cross-format float copy retained a binary64-sized payload.
        // That 4096-bit fixture is not a binary32 auxiliary expression.
        rational const binary64_fixture = rational((integer(1) << 4096) - 1)
            / ((integer(1) << 4095) + 17);
        BOOST_TEST_THROWS(p::literal(binary64_fixture, third_value), std::overflow_error);
    }
}

void test_derived_cache_capacities()
{
    static_assert(sizeof(ec::point_cache_detail::native_cache_handle) == sizeof(std::size_t),
                  "Native cache handles are one plain index");
    static_assert(std::is_same<ec::auxiliary_point_type<std::int16_t>, ec::fixed_homogeneous<87, 71>>::value,
                  "Int16 auxiliary capacity follows the expression degree");
    static_assert(std::is_same<ec::auxiliary_point_type<std::int64_t>, ec::fixed_homogeneous<327, 263>>::value,
                  "Int64 auxiliary capacity follows the expression degree");
    static_assert(std::is_same<ec::auxiliary_point_type<double>, ec::fixed_homogeneous<11571, 10547>>::value,
                  "Binary64 auxiliary capacity follows the expression degree");
    // Legacy native box-corner setters must not create operation-owned values.
    point<double> corner;
    bg::set<0>(corner, rational(std::nextafter(1.0, 2.0)));
    bg::set<1>(corner, rational(-7));
    BOOST_TEST(corner.is_input());
    BOOST_TEST(bg::get<0>(corner) == rational(std::nextafter(1.0, 2.0)));
    BOOST_TEST_EQ(corner.input_coordinate(1), -7.0);
    ec::point_cache_scope<std::int64_t> operation;
    using p = point<std::int64_t>;
    auto const lo = (std::numeric_limits<std::int64_t>::min)();
    auto const hi = (std::numeric_limits<std::int64_t>::max)();
    auto const a = ec::construct_intersection(p(lo, lo), p(hi, hi - 1), p(lo, hi), p(hi - 1, lo));
    auto const b = ec::construct_intersection(p(lo, lo + 1), p(hi - 1, hi), p(lo + 1, hi), p(hi, lo));
    auto const splice = p::splice_coordinates(a, 0, b, 1);
    auto const middle = p::from_homogeneous(ec::fixed_midpoint(a.fixed_homogeneous(), b.fixed_homogeneous()));
    BOOST_TEST(splice.exact_coordinate(0) == a.exact_coordinate(0));
    BOOST_TEST(splice.exact_coordinate(1) == b.exact_coordinate(1));
    BOOST_TEST(middle.exact_coordinate(0) * 2 == a.exact_coordinate(0) + b.exact_coordinate(0));
    BOOST_TEST(middle.exact_coordinate(1) * 2 == a.exact_coordinate(1) + b.exact_coordinate(1));
}

void test_point_copies()
{
    ec::point_cache_scope<double> double_scope;
    ec::point_cache_scope<std::int64_t> integral_scope;
    input_point<std::int64_t> const native(9007199254740993LL,
        (std::numeric_limits<std::int64_t>::min)());
    point<std::int64_t> copy;
    bg::detail::conversion::point_to_point<decltype(native), decltype(copy), 0, 2>::apply(native, copy);
    BOOST_TEST(copy.is_input());
    input_point<std::int64_t> output;
    ec::assign_output(copy, output);
    BOOST_TEST_EQ(bg::get<0>(output), bg::get<0>(native));
    BOOST_TEST_EQ(bg::get<1>(output), bg::get<1>(native));

    static_assert(!ec::lossless_coordinate_conversion<std::int64_t, double>::value,
                  "Int64 to double loses native inputs and is not supported");
    static_assert(!ec::lossless_coordinate_conversion<double, float>::value,
                  "Narrowing floating input is not supported");
    static_assert(ec::lossless_coordinate_conversion<std::int32_t, double>::value,
                  "Int32 to double is lossless");
    static_assert(ec::lossless_coordinate_conversion<float, double>::value,
                  "Float to double is lossless");
    input_point<std::int32_t> const narrower((std::numeric_limits<std::int32_t>::max)(),
        (std::numeric_limits<std::int32_t>::min)());
    point<double> promoted;
    bg::detail::conversion::point_to_point<decltype(narrower), decltype(promoted), 0, 2>::apply(narrower, promoted);
    BOOST_TEST(promoted.is_input());
    BOOST_TEST(bg::get<0>(promoted) == rational(bg::get<0>(narrower)));
    BOOST_TEST(bg::get<1>(promoted) == rational(bg::get<1>(narrower)));
    input_point<float> const single(std::nextafter(1.0f, 2.0f), -7.0f);
    auto const widened = ec::make_point<double>(single);
    BOOST_TEST(bg::get<0>(widened) == ec::exact_value(bg::get<0>(single)));
    BOOST_TEST(bg::get<1>(widened) == -7);

    // Preserve the large-coordinate geometry as an explicit bounded diagnostic
    // value, not as an allegedly lossless native int64-to-double promotion.
    auto const diagnostic = point<double>::from_homogeneous(
        ec::fixed_input_homogeneous(bg::get<0>(native), bg::get<1>(native)));
    BOOST_TEST(bg::get<0>(diagnostic) == rational(9007199254740993LL));
    BOOST_TEST(bg::get<1>(diagnostic) == rational(bg::get<1>(native)));

    auto const constructed = third<double>();
    point<double> retained;
    bg::detail::conversion::point_to_point<decltype(constructed), decltype(retained), 0, 2>
        ::apply(constructed, retained);
    BOOST_TEST(retained.shares_node(constructed));
    input_point<double> rounded;
    ec::assign_output<nearest_output>(retained, rounded);
    BOOST_TEST_EQ(bg::get<0>(rounded), 1.0 / 3.0);

    input_point<double> const partial(17, 23);
    bg::detail::conversion::point_to_point<decltype(partial), decltype(retained), 1, 2>
        ::apply(partial, retained);
    BOOST_TEST(bg::get<0>(retained) == rational(1) / 3);
    BOOST_TEST(bg::get<1>(retained) == 23);
    BOOST_TEST(bg::get<1>(constructed) == rational(1) / 3);

    input_point<double> const floating(-0.0, std::ldexp(1.0, -100));
    bg::detail::conversion::point_to_point<decltype(floating), decltype(retained), 0, 2>
        ::apply(floating, retained);
    BOOST_TEST(retained.is_input());
    ec::assign_output(retained, rounded);
    BOOST_TEST(std::signbit(bg::get<0>(rounded)));
    BOOST_TEST_EQ(bg::get<1>(rounded), bg::get<1>(floating));
}

int main()
{
    BOOST_TEST_EQ(std::fegetround(), FE_TONEAREST);
    {
        ec::point_cache_scope<double> operation;
        auto const p = ec::construct_intersection(point<double>(0, 1), point<double>(11, -13),
            point<double>(-17, 6), point<double>(5, -1));
        BOOST_TEST(p.exact_coordinate(0) == rational(3) / 7);
        BOOST_TEST(p.exact_coordinate(1) == rational(5) / 11);
    }
    test_constructed_points<float>();
    test_constructed_points<double>();
    test_constructed_points<int>();
    test_unrepresentable_third<float>();
    test_unrepresentable_third<double>();
    test_cancellation();
    test_floating_leaf_fallbacks();
    test_large_integers();
    test_integral_leaf_no_allocations<std::int16_t>();
    test_integral_leaf_no_allocations<std::int32_t>();
    test_integral_leaf_no_allocations<std::int64_t>();
    test_integral_leaf_no_allocations<std::uint64_t>();
    test_constructed_forward_predicates();
    test_binary_rounding<float>();
    test_binary_rounding<double>();
    test_integral_conversion<approximate_output>();
    test_integral_conversion<nearest_output>();
    test_concurrent_fallback();
    test_inline_native_points<float>();
    test_inline_native_points<double>();
    test_inline_native_points<std::int64_t>();
    test_inline_native_points<std::uint64_t>();
    test_native_fixed_extremes<float>();
    test_native_fixed_extremes<double>();
    test_native_fixed_extremes<std::int64_t>();
    test_native_fixed_extremes<std::uint64_t>();
    test_dense_fixed_orientation_no_allocations();
    test_expansion_refinement_cache();
    test_cached_approximation_copies();
    test_fixed_cache_bounds();
    test_retained_rounding_quotient();
    test_uncertain_initial_denominator();
    test_indexed_point_cache();
    test_same_operation_native_producers();
    test_independent_point_cache_concurrency();
    test_constructed_round_to_nearest();
    test_fixed_native_output();
    test_materialized_precision_cache<float>();
    test_materialized_precision_cache<double>();
    test_materialized_precision_cache<std::int16_t>();
    test_materialized_precision_cache<std::int64_t>();
    test_derived_cache_capacities();
    test_point_copies();
    BOOST_TEST_THROWS((point<double>((std::numeric_limits<double>::infinity)(), 0)), std::domain_error);
    BOOST_TEST_THROWS((point<double>((std::numeric_limits<double>::quiet_NaN)(), 0)), std::domain_error);
    return boost::report_errors();
}
