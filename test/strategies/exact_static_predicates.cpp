// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.

#define BOOST_GEOMETRY_EXACT_FILTER_STATISTICS
#include <boost/core/lightweight_test.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/strategies/cartesian/exact/static_predicates.hpp>
#include <boost/geometry/strategies/detail/distance_measure.hpp>
#include "exact_test_context.hpp"

#include <array>
#include <atomic>
#include <cfenv>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <new>
#include <random>
#include <type_traits>
#include <thread>
#include <utility>

namespace bg = boost::geometry;
namespace ec = bg::strategies::detail::exact_cartesian;
using rational = ec::exact_number;
using boost::multiprecision::cpp_int;
template <typename C> using point = bg::model::d2::point_xy<C>;
using exact_coordinates = std::array<rational, 2>;

namespace
{
bool counting = false;
std::size_t allocations = 0;
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
auto no_allocations(Function function) -> decltype(function())
{
    function(); // Cache growth is allowed; numerical objects have no own pool.
    prepare_exact_cache_repeat();
    allocations = 0;
    counting = true;
    auto result = function();
    counting = false;
    BOOST_TEST_EQ(allocations, 0u);
    return result;
}

template <typename C>
exact_coordinates coordinates(ec::exact_point<C> const& point)
{
    auto convert = [](auto const& h)
    {
        return exact_coordinates{{rational(cpp_int(h.X)) / rational(cpp_int(h.W)),
                                  rational(cpp_int(h.Y)) / rational(cpp_int(h.W))}};
    };
    return point.has_fixed_representation() ? convert(point.fixed_homogeneous()) : convert(point.precision_homogeneous());
}
template <typename Point>
exact_coordinates coordinates(Point const& point)
{
    return {{ec::exact_value(bg::get<0>(point)), ec::exact_value(bg::get<1>(point))}};
}
exact_coordinates difference(exact_coordinates const& a, exact_coordinates const& b)
{ return {{a[0] - b[0], a[1] - b[1]}}; }
rational dot(exact_coordinates const& a, exact_coordinates const& b) { return a[0] * b[0] + a[1] * b[1]; }
rational cross(exact_coordinates const& a, exact_coordinates const& b) { return a[0] * b[1] - a[1] * b[0]; }
int compare(rational const& a, rational const& b) { return (a > b) - (a < b); }

template <typename Measure>
void check_bounds(Measure const& value, rational const& expected)
{
    if (value.valid)
    {
        rational error = rational(value.approximate) - expected;
        if (error < 0) { error = -error; }
        BOOST_TEST(error <= rational(value.error));
    }
}

template <typename A, typename B, typename P, typename Q>
void check_case(A const& a, B const& b, P const& p, Q const& q)
{
    auto const x = coordinates(a), y = coordinates(b), z = coordinates(p), w = coordinates(q);
    auto const ab = difference(y, x), ap = difference(z, x), bp = difference(z, y), pq = difference(w, z);
    rational const direction = dot(ab, bp), determinant = cross(ab, ap);
    rational const distance = dot(ab, ab), other_distance = dot(pq, pq);
    ec::filtering_statistics() = {};
    BOOST_TEST_EQ(no_allocations([&] { return ec::static_direction(a, b, p); }), direction.sign());
    BOOST_TEST_EQ(no_allocations([&] { return ec::static_orientation(a, b, p); }), determinant.sign());
    BOOST_TEST_EQ(no_allocations([&] { return ec::static_orientation(b, a, p); }), -determinant.sign());
    for (std::size_t d = 0; d < 2; ++d)
    {
        BOOST_TEST_EQ(no_allocations([&] { return ec::static_compare_coordinate(a, b, d); }), compare(x[d], y[d]));
    }
    BOOST_TEST_EQ(no_allocations([&] { return ec::static_compare_coordinate<0,double>(a,b); }), compare(x[0],y[0]));
    BOOST_TEST_EQ(no_allocations([&] { return ec::static_compare_coordinate<1,double>(a,b); }), compare(x[1],y[1]));
    auto first = no_allocations([&] { return ec::static_squared_distance(a, b); });
    auto second = no_allocations([&] { return ec::static_squared_distance(p, q); });
    auto measure = no_allocations([&] { return ec::static_signed_determinant(a, b, p); });
    check_bounds(first, distance);
    check_bounds(second, other_distance);
    check_bounds(measure, determinant);
    BOOST_TEST_EQ(no_allocations([&] { return ec::compare_distance(first, second); }), compare(distance, other_distance));
    BOOST_TEST_EQ(no_allocations([&] { return first < second; }), distance < other_distance);
    BOOST_TEST_EQ(no_allocations([&] { return first == second; }), distance == other_distance);
    BOOST_TEST_EQ(no_allocations([&] { return first >= second; }), distance >= other_distance);
    BOOST_TEST_EQ(no_allocations([&] { return first == 0; }), distance == 0);
    BOOST_TEST_EQ(no_allocations([&] { return first.is_zero(); }), distance == 0);
    BOOST_TEST_EQ(no_allocations([&] { return 0 < first; }), distance > 0);
    BOOST_TEST_EQ(no_allocations([&] { return measure > 0; }), determinant > 0);
    BOOST_TEST_EQ(no_allocations([&] { return measure < 0; }), determinant < 0);
    BOOST_TEST_EQ(no_allocations([&] { return -measure > 0; }), determinant < 0);
    BOOST_TEST(no_allocations([&] { return -(-measure) == measure; }));
    BOOST_TEST_EQ(no_allocations([&] { return ec::compare_distance(first, measure); }), compare(distance, determinant));
    auto const wrapped = no_allocations([&] { return ec::static_distance_measure(a,b,p); });
    BOOST_TEST_EQ(no_allocations([&] { return wrapped.is_zero(); }), determinant == 0);
    BOOST_TEST_EQ(no_allocations([&] { return wrapped.is_positive(); }), determinant > 0);
    BOOST_TEST_EQ(no_allocations([&] { return wrapped.is_negative(); }), determinant < 0);
    BOOST_TEST_EQ(no_allocations([&] { return ec::static_coordinate<0>(a) < ec::static_coordinate<1>(b); }), x[0] < y[1]);
    BOOST_TEST_EQ(no_allocations([&] { return ec::static_coordinate<1>(a) == ec::static_coordinate<0>(b); }), x[1] == y[0]);
    auto const statistics = ec::filtering_statistics();
    BOOST_TEST_EQ(statistics.native_rationals, 0u);
    BOOST_TEST_EQ(statistics.point_coordinate_values, 0u);
    BOOST_TEST_EQ(statistics.homogeneous_values, 0u);
}

void segment_order()
{
    using p = ec::exact_point<double>;
    p const a(0, 0), b(1, 0);
    double const tiny = std::ldexp(1.0, -150);
    p const near_start = ec::construct_intersection(a, b, p(0, -tiny), p(1, 1));
    p const near_end = ec::construct_intersection(a, b, p(1, -tiny), p(2, 1));
    std::array<p, 7> const points{{p(-tiny, 0), a, p(tiny, 0), near_start, b, near_end, p(0, tiny)}};
    std::array<int, 7> const expected{{0, 1, 2, 2, 3, 4, 7}};
    for (std::size_t i = 0; i < points.size(); ++i)
    {
        BOOST_TEST_EQ(no_allocations([&] { return ec::classify_point_on_segment(a, b, points[i]); }), expected[i]);
        int const reverse = expected[i] <= 4 ? 4 - expected[i] : expected[i];
        BOOST_TEST_EQ(no_allocations([&] { return ec::classify_point_on_segment(b, a, points[i]); }), reverse);
    }
    for (auto const& first : points)
    for (auto const& second : points)
    {
        auto const x = coordinates(first), y = coordinates(second);
        BOOST_TEST_EQ(no_allocations([&] { return ec::compare_along_segment(a, b, first, second); }), compare(x[0], y[0]));
        BOOST_TEST_EQ(no_allocations([&] { return ec::compare_along_segment(b, a, first, second); }), compare(y[0], x[0]));
    }
    BOOST_TEST_EQ(ec::classify_point_on_segment(a, a, a), 5);
    BOOST_TEST_EQ(ec::classify_point_on_segment(a, a, b), 6);
    BOOST_TEST_EQ(ec::compare_along_segment(a, a, a, a), 0);
    BOOST_TEST_THROWS(ec::compare_along_segment(a, a, a, b), std::domain_error);
    BOOST_TEST_EQ(ec::classify_point_on_segment(p(0, 1), a, p(0, tiny)), 2);
    BOOST_TEST_EQ(ec::compare_along_segment(p(0, 1), a, p(0, tiny), a), -1);
    check_case(a, b, near_start, near_end);
    check_case(near_end, near_start, p(0, 1), p(0, -1));
}

void native_extremes()
{
    auto const lo = (std::numeric_limits<std::int64_t>::min)();
    auto const hi = (std::numeric_limits<std::int64_t>::max)();
    auto const top = (std::numeric_limits<std::uint64_t>::max)();
    using i = point<std::int64_t>;
    using u = point<std::uint64_t>;
    check_case(i(lo, lo), i(hi, hi - 1), i(hi - 1, hi - 2), i(lo, hi));
    check_case(u(0, 0), u(top, top - 1), u(top - 1, top - 2), u(top, 0));
    check_case(i(-1, lo), u(top, top), point<double>(0.1, -0.2), point<float>(1, 2));
    std::int64_t const large = std::int64_t(1) << 53;
    check_case(i(large + 1, large), point<double>(double(large), 0), i(large + 2, large + 3), u(top, 1));
    auto const mixed = ec::static_squared_distance(i(large + 1, 0), point<double>(double(large), 0));
    static_assert(std::is_same<decltype(mixed.a), ec::exact_point<std::int64_t>>::value, "Retain native integer precision");
    static_assert(std::is_same<decltype(mixed.b), ec::exact_point<double>>::value, "Retain native binary precision");
    BOOST_TEST(mixed > 0);
    static_assert(!ec::lossless_coordinate_conversion<std::int64_t, double>::value,
                  "The strategy rejects lossy mixed-type promotion");
    auto const homogeneous = ec::static_squared_distance<double>(point<std::int32_t>(18, 0), point<double>(17, 0));
    static_assert(std::is_same<typename std::decay<decltype(homogeneous)>::type,
        ec::squared_distance_value<double>>::value, "Explicit strategy coordinate gives a homogeneous typed result");
    BOOST_TEST(no_allocations([&] { return homogeneous == ec::static_squared_distance(point<double>(0,0),point<double>(1,0)); }));
    BOOST_TEST_EQ(ec::static_compare_coordinate(i(large + 1, 0), point<double>(double(large), 0), 0), 1);
    BOOST_TEST_EQ(ec::classify_point_on_segment(u(0,0),u(top,0),u(top-1,0)),2);
    BOOST_TEST_EQ(ec::compare_along_segment(i(hi,0),i(lo,0),i(hi-1,0),i(lo+1,0)),-1);
    double const maximum = std::ldexp(1.0, 200);
    double const tiny = std::ldexp(1.0, -200);
    check_case(point<double>(-maximum, -tiny), point<double>(maximum, tiny),
        point<double>(maximum, -tiny), point<double>(0, 0));
}

void randomized()
{
    std::mt19937_64 random(0x535441544943);
    for (unsigned n = 0; n < 64; ++n)
    {
        std::array<std::int64_t, 8> integers;
        for (auto& value : integers) { value = std::int64_t(random()); }
        check_case(point<std::int64_t>(integers[0], integers[1]), point<std::int64_t>(integers[2], integers[3]),
            point<std::int64_t>(integers[4], integers[5]), point<std::int64_t>(integers[6], integers[7]));
        std::array<double, 8> binary;
        for (auto& value : binary)
        {
            value = std::ldexp(double(std::int64_t(random() % 2001) - 1000), int(random() % 201) - 100);
        }
        check_case(point<double>(binary[0], binary[1]), point<double>(binary[2], binary[3]),
            point<double>(binary[4], binary[5]), point<double>(binary[6], binary[7]));
    }
}

void lifetime_and_preconditions()
{
    auto retained = []
    {
        using p = ec::exact_point<double>;
        auto third = ec::construct_intersection(p(0, 0), p(1, 1), p(0, 1), p(1, -1));
        return ec::static_signed_determinant(p(0, 0), p(1, 0), third);
    }();
    BOOST_TEST(retained > 0);
    BOOST_TEST(-retained < 0);
    BOOST_TEST_THROWS(decltype(retained)(1), std::invalid_argument);
    ec::squared_distance_value<double> const unset = -1, zero = 0;
    BOOST_TEST(unset < 0);
    BOOST_TEST(!unset.is_zero());
    BOOST_TEST(zero.is_zero());
    BOOST_TEST(unset < zero);
    BOOST_TEST(zero > unset);
    BOOST_TEST(unset == unset);
    BOOST_TEST_THROWS(ec::squared_distance_value<double>(2), std::invalid_argument);
    auto literal = ec::exact_point<double>::literal(rational(1) / 3, 0);
    BOOST_TEST(ec::static_coordinate<0>(literal) > ec::static_coordinate<0>(point<double>(0, 0)));
    check_case(literal, point<double>(0, 0), point<std::int64_t>(1, 0), point<double>(1, 1));
    static_assert(!std::is_constructible<ec::static_coordinate_view<double>, int>::value, "Coordinates are views, not scalar expressions");
}

void retained_rational_inputs()
{
    using rp = point<rational>;
    rational const tiny = rational(1) / rational(cpp_int(1) << 150);
    rp const a(0,0), b(1,0), p(tiny,0), after(1+tiny,0), off(0,tiny);
    auto const working = ec::make_point<double>(p);
    ec::filtering_statistics() = {};
    BOOST_TEST_EQ((ec::static_compare_coordinate<0,double>(p,a)),1);
    BOOST_TEST_EQ((ec::static_compare_coordinate<0,double>(a,p)),-1);
    BOOST_TEST_EQ((ec::static_compare_coordinate<0,double>(p,point<double>(0,0))),1);
    BOOST_TEST_EQ((ec::static_compare_coordinate<0,double>(point<double>(0,0),p)),-1);
    BOOST_TEST_EQ((ec::static_compare_coordinate<0,double>(working,p)),0);
    BOOST_TEST_EQ((ec::static_compare_coordinate<1,double>(off,a)),1);
    BOOST_TEST(ec::static_equal_points<double>(working,p));
    BOOST_TEST(!ec::static_equal_points<double>(a,p));
    BOOST_TEST_EQ(ec::static_orientation<double>(a,b,off),1);
    BOOST_TEST_EQ(ec::static_orientation<double>(a,b,p),0);
    BOOST_TEST_EQ(ec::static_direction<double>(a,b,p),-1);
    BOOST_TEST_EQ(ec::static_direction_code<double>(a,b,after),1);
    BOOST_TEST_EQ(ec::classify_point_on_segment<double>(a,b,p),2);
    BOOST_TEST_EQ(ec::classify_point_on_segment<double>(a,b,after),4);
    BOOST_TEST_EQ(ec::compare_along_segment<double>(a,b,p,after),-1);
    BOOST_TEST((ec::static_coordinate<0,double>(p)>ec::static_coordinate<0,double>(a)));
    BOOST_TEST(ec::static_distance_measure<double>(a,b,off).is_positive());
    auto const distance = ec::static_squared_distance<double>(a,p);
    BOOST_TEST(!distance.is_zero());
    BOOST_TEST(distance > 0);
    BOOST_TEST_EQ(ec::filtering_statistics().point_coordinate_values,0u);
    BOOST_TEST_EQ(ec::filtering_statistics().homogeneous_values,0u);
}

void compact_comparisons()
{
    ec::point_cache_scope operation;
    using p = ec::exact_point<double>;
    p const origin(0, 0);
    auto const first = ec::static_squared_distance(origin, p(3, 4));
    auto const tied = ec::static_squared_distance(origin, p(-4, 3));
    decltype(first)::compact_pair_type compact;
    BOOST_TEST(no_allocations([&] { return first.compact_pair(compact); }));
    BOOST_TEST_EQ(compact.N, 25);
    BOOST_TEST_EQ(compact.D, 1);
    BOOST_TEST_EQ(no_allocations([&] { return ec::compare_distance(first, tied); }), 0);

    std::int64_t const limit = std::int64_t(1) << 30;
    point<std::int64_t> const low(-limit, -limit), high(limit, limit);
    std::uint64_t native = 0;
    BOOST_TEST(no_allocations([&] { return ec::static_predicates_detail::small_integer_squared_distance(
        ec::make_point<std::int64_t>(low), ec::make_point<std::int64_t>(high), native); }));
    BOOST_TEST_EQ(native, std::uint64_t(1) << 63);
    check_case(low, high, high, low);
    for (auto const& outside : {p(double(limit + 1), 0), p(double(-limit - 1), 0), p(0.5, 0)})
    {
        native = 17;
        BOOST_TEST(!no_allocations([&]
            { return ec::static_predicates_detail::small_integer_squared_distance(origin, outside, native); }));
        BOOST_TEST_EQ(native, 17u);
        check_case(origin, outside, outside, origin);
    }
    check_case(point<std::uint64_t>(0, 0), point<std::uint64_t>(limit, limit), low, point<std::int64_t>(0, 0));

    p const intersection = ec::construct_intersection(p(0, 0), p(2, 2), p(0, 2), p(1, 0));
    ec::fixed_homogeneous<128> coordinate;
    BOOST_TEST(no_allocations([&] { return intersection.copy_compact_homogeneous(coordinate); }));
    BOOST_TEST_EQ(coordinate.X, 2);
    BOOST_TEST_EQ(coordinate.Y, 2);
    BOOST_TEST_EQ(coordinate.W, 3);
    auto const rational_intersection = p::literal(rational(2) / 3, rational(2) / 3);
    BOOST_TEST(no_allocations([&] { return rational_intersection.copy_compact_homogeneous(coordinate); }));
    BOOST_TEST_EQ(no_allocations([&]
        { return ec::static_squared_distance(origin, intersection)
            == ec::static_squared_distance(origin, rational_intersection); }), true);

    for (unsigned bits : {63, 64, 65, 126, 127, 128, 129, 255, 512, 2048})
    {
        rational const tiny = rational(1) / rational(cpp_int(1) << bits);
        auto const low = p::literal(1 - tiny, 0);
        auto const high = p::literal(1 + tiny, 0);
        auto const below = ec::static_squared_distance(origin, low);
        auto const above = ec::static_squared_distance(origin, high);
        bool const fits = bits < 128; // 1+2^-bits has a bits+1-bit numerator.
        BOOST_TEST_EQ(no_allocations([&] { return above.compact_pair(compact); }), fits);
        BOOST_TEST_EQ(no_allocations([&] { return ec::compare_distance(below, above); }), -1);
        BOOST_TEST_EQ(no_allocations([&] { return ec::compare_distance(above, below); }), 1);
        BOOST_TEST_EQ(no_allocations([&] { return ec::compare_distance(above, above); }), 0);
        check_case(origin, high, origin, low);
    }

    ec::fixed_homogeneous<128> maximum;
    maximum.X = (ec::fixed_integer<128>(1) << 127) - 1;
    maximum.X = (maximum.X << 1) + 1;
    maximum.Y = -maximum.X + 1;
    maximum.W = maximum.X - 2;
    auto const full = p::from_homogeneous(maximum);
    check_case(full, p(2, -2), p(-1, 1), p(-2, 2));

    std::atomic<unsigned> errors{0};
    std::array<std::thread, 4> workers;
    for (auto& worker : workers) worker = std::thread([&]
    {
        ec::point_cache_scope thread_operation;
        auto const local_intersection = ec::construct_intersection(p(0, 0), p(2, 2), p(0, 2), p(1, 0));
        auto const expected = ec::static_squared_distance(origin, p::literal(rational(2) / 3, rational(2) / 3));
        for (unsigned i = 0; i < 100; ++i)
            if (ec::compare_distance(ec::static_squared_distance(origin, local_intersection), expected) != 0) ++errors;
    });
    for (auto& worker : workers) worker.join();
    BOOST_TEST_EQ(errors.load(), 0u);
}

int main()
{
    ec::point_cache_scope operation;
    int const original = std::fegetround();
    for (int mode : {FE_TONEAREST})
    {
        if (std::fesetround(mode) != 0) { continue; }
        native_extremes();
        segment_order();
        compact_comparisons();
        check_case(point<double>(0, 0), point<double>(1, 1), point<double>(0, 2), point<double>(1, 0));
    }
    BOOST_TEST_EQ(std::fesetround(original), 0);
    randomized();
    lifetime_and_preconditions();
    retained_rational_inputs();
    return boost::report_errors();
}
