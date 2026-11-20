// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.

#include <boost/core/lightweight_test.hpp>
#include <boost/geometry/strategies/cartesian/exact/borrowed_query.hpp>
#include <boost/geometry/strategies/cartesian/exact/static_predicates.hpp>
#include <boost/geometry/geometries/point_xy.hpp>

#include <array>
#include <cfenv>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <new>
#include <type_traits>
#include <utility>

namespace bg = boost::geometry;
namespace ec = bg::strategies::detail::exact_cartesian;

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
    if (counting) ++allocations;
    if (void* result = std::malloc(size ? size : 1)) return result;
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
auto without_allocations(Function function)
{
    auto result = [&]
    {
        allocation_scope scope;
        return function();
    }();
    BOOST_TEST_EQ(allocations, 0u);
    return result;
}

template <typename Root, typename A, typename B, typename = void>
struct has_midpoint : std::false_type {};
template <typename Root, typename A, typename B>
struct has_midpoint<Root, A, B, std::void_t<decltype(ec::query_point<Root>::midpoint(
    std::declval<A const&>(), std::declval<B const&>()))>> : std::true_type {};
template <typename Root, typename A, typename B, typename = void>
struct has_splice : std::false_type {};
template <typename Root, typename A, typename B>
struct has_splice<Root, A, B, std::void_t<decltype(ec::query_point<Root>::splice_coordinates(
    std::declval<A const&>(), 0, std::declval<B const&>(), 1))>> : std::true_type {};

template <typename Query>
void check_oracle(Query const& query, ec::exact_number const& x, ec::exact_number const& y)
{
    auto const value = without_allocations([&] { return query.homogeneous(); });
    // Independent unbounded rational arithmetic is confined to this test oracle.
    ec::exact_number actual_x(value.X), actual_y(value.Y);
    actual_x /= value.W;
    actual_y /= value.W;
    BOOST_TEST(actual_x == x);
    BOOST_TEST(actual_y == y);

    ec::fixed_homogeneous<128, 128> compact;
    compact.X = 7; compact.Y = -5; compact.W = 3;
    bool const copied = without_allocations([&] { return query.copy_compact_homogeneous(compact); });
    bool const fits = ec::fixed_fits<128>(value.X) && ec::fixed_fits<128>(value.Y)
        && ec::fixed_fits<128>(value.W);
    BOOST_TEST_EQ(copied, fits);
    if (copied)
    {
        ec::exact_number compact_x(compact.X), compact_y(compact.Y);
        compact_x /= compact.W;
        compact_y /= compact.W;
        BOOST_TEST(compact_x == x);
        BOOST_TEST(compact_y == y);
    }
    else BOOST_TEST(compact.X == 7 && compact.Y == -5 && compact.W == 3);
}

template <typename Root, typename A, typename B>
void check_native_pair(A const& a, B const& b)
{
    using query = ec::query_point<Root>;
    using C = typename Root::coordinate_type;
    auto& storage = ec::operation_detail::require<C>().storage();
    std::array<ec::exact_number, 2> const first{{ec::exact_number(bg::get<0>(a)), ec::exact_number(bg::get<1>(a))}};
    std::array<ec::exact_number, 2> const second{{ec::exact_number(bg::get<0>(b)), ec::exact_number(bg::get<1>(b))}};
    auto const before = storage.auxiliary.size();
    auto const middle = [&]
    {
        allocation_scope scope;
        return query::midpoint(a, b);
    }();
    ec::exact_number const x = (first[0] + second[0]) / 2;
    ec::exact_number const y = (first[1] + second[1]) / 2;
    if (storage.auxiliary.size() == before)
    {
        BOOST_TEST_EQ(allocations, 0u);
        BOOST_TEST(ec::exact_number(middle.approximate(0)) == x);
        BOOST_TEST(ec::exact_number(middle.approximate(1)) == y);
        for (std::size_t d = 0; d != 2; ++d)
        {
            BOOST_TEST_EQ(middle.approximation(d).error, 0.0);
            BOOST_TEST_EQ(middle.refine(d).error, 0.0);
            BOOST_TEST(middle.approximation(d).valid());
        }
    }
    else BOOST_TEST_EQ(storage.auxiliary.size(), before + 1);
    check_oracle(middle, x, y);

    for (std::size_t i = 0; i != 2; ++i)
        for (std::size_t j = 0; j != 2; ++j)
        {
            double const ax = i == 0 ? double(bg::get<0>(a)) : double(bg::get<1>(a));
            double const by = j == 0 ? double(bg::get<0>(b)) : double(bg::get<1>(b));
            bool const exact = ec::exact_number(ax) == first[i] && ec::exact_number(by) == second[j];
            auto const previous = storage.auxiliary.size();
            auto const corner = [&]
            {
                allocation_scope scope;
                return query::splice_coordinates(a, i, b, j);
            }();
            BOOST_TEST_EQ(storage.auxiliary.size(), previous + (exact ? 0 : 1));
            if (exact)
            {
                BOOST_TEST_EQ(allocations, 0u);
                for (std::size_t d = 0; d != 2; ++d)
                {
                    BOOST_TEST_EQ(corner.approximation(d).error, 0.0);
                    BOOST_TEST_EQ(corner.refine(d).error, 0.0);
                    BOOST_TEST(corner.approximation(d).valid());
                }
            }
            check_oracle(corner, first[i], second[j]);
        }
}

template <typename C, typename Policy>
void test_native_queries()
{
    using native = bg::model::d2::point_xy<C>;
    using root = ec::borrowed_root<C, Policy, native, native>;
    using point = ec::borrowed_point<root>;
    using query = ec::query_point<root>;
    static_assert(sizeof(std::size_t) != 8 || sizeof(query) == 24, "No extra exactness flag or error fields");
    ec::borrowed_context<root> context;
    auto& storage = ec::operation_detail::require<C>().storage();
    native const a(2, 4), b(4, 6);
    auto const first = point::borrow(a), second = point::borrow(b);
    auto const middle = without_allocations([&] { return query::midpoint(first, second); });
    auto const corner = without_allocations([&] { return query::splice_coordinates(first, 1, b, 0); });
    check_oracle(middle, ec::exact_number(3), ec::exact_number(5));
    check_oracle(corner, ec::exact_number(4), ec::exact_number(4));
    ec::fixed_homogeneous<3, 1> small;
    BOOST_TEST(without_allocations([&] { return middle.copy_compact_homogeneous(small); }));
    BOOST_TEST(small.X == 3 && small.Y == 5 && small.W == 1);
    ec::fixed_homogeneous<2, 1> too_small;
    too_small.X = 2; too_small.Y = -1; too_small.W = 1;
    BOOST_TEST(!without_allocations([&] { return corner.copy_compact_homogeneous(too_small); }));
    BOOST_TEST(too_small.X == 2 && too_small.Y == -1 && too_small.W == 1);
    C const lowest = (std::numeric_limits<C>::lowest)();
    native const low_a(lowest, 0), low_b(lowest, 1);
    auto const mixed_denominator = without_allocations([&] { return query::midpoint(low_a, low_b); });
    check_oracle(mixed_denominator, ec::exact_number(lowest), ec::exact_number(1) / 2);
    without_allocations([&] { return middle.refine(0); });
    BOOST_TEST(storage.native.empty());
    BOOST_TEST(storage.auxiliary.empty());
    BOOST_TEST_EQ(storage.auxiliary.capacity(), 0u);

    std::array<C, 10> values{{C(0), C(1), C(2), C(3), C(4), C(5), C(6), C(7),
        (std::numeric_limits<C>::lowest)(), (std::numeric_limits<C>::max)()}};
    if constexpr (std::is_floating_point<C>::value)
    {
        values[2] = (std::numeric_limits<C>::denorm_min)();
        values[3] = C(2) * values[2];
        values[4] = (std::numeric_limits<C>::min)();
        values[5] = std::nextafter(C(1), C(2));
        values[6] = std::nextafter(values[4], C(1));
        values[7] = -(std::numeric_limits<C>::denorm_min)();
    }
    else if constexpr (std::numeric_limits<C>::digits > 53)
    {
        values[4] = C(std::uint64_t(1) << 53);
        values[5] = values[4] + 1;
        values[6] = values[4] + 2;
        values[7] = (std::numeric_limits<C>::max)() - C(2047);
    }
    for (std::size_t i = 0; i != values.size(); ++i)
        for (std::size_t j = 0; j != values.size(); ++j)
            check_native_pair<root>(native(values[i], values[j]), native(values[j], values[i]));

    if constexpr (std::is_floating_point<C>::value)
    {
        C const maximum = (std::numeric_limits<C>::max)();
        C const tiny = (std::numeric_limits<C>::denorm_min)();
        native const high(maximum, -maximum), low(tiny, -tiny), zero(0, 0);
        auto const before = storage.auxiliary.size();
        std::feclearexcept(FE_OVERFLOW | FE_INVALID);
        auto const same = without_allocations([&] { return query::midpoint(high, high); });
        BOOST_TEST_EQ(std::fetestexcept(FE_OVERFLOW | FE_INVALID), 0);
        check_oracle(same, ec::exact_number(maximum), ec::exact_number(-maximum));
        auto const subnormal = without_allocations([&] { return query::midpoint(low, low); });
        check_oracle(subnormal, ec::exact_number(tiny), ec::exact_number(-tiny));
        BOOST_TEST_EQ(storage.auxiliary.size(), before);
        auto const half = query::midpoint(low, zero);
        check_oracle(half, ec::exact_number(tiny) / 2, -ec::exact_number(tiny) / 2);
        // A binary32 half-subnormal fits binary64; a binary64 half-subnormal does not.
        BOOST_TEST_EQ(storage.auxiliary.size(), before + (std::is_same<C, double>::value ? 1 : 0));
        ec::fixed_homogeneous<128, 128> rejected;
        BOOST_TEST(!ec::borrowed_query_detail::copy_compact_geometry<root>(low, rejected));
        BOOST_TEST(!without_allocations([&] { return half.copy_compact_homogeneous(rejected); }));
    }
}

template <typename C, typename Policy>
void test_compact_rational_queries()
{
    using native = bg::model::d2::point_xy<C>;
    using polygon = bg::model::polygon<native>;
    using root = ec::borrowed_root<C, Policy, polygon, polygon>;
    using point = ec::borrowed_point<root>;
    using query = ec::query_point<root>;
    using compact = ec::fixed_homogeneous<128, 128>;
    static_assert(std::is_same<decltype(ec::fixed_midpoint(std::declval<compact>(), std::declval<compact>())),
        ec::fixed_homogeneous<257, 257>>::value, "Compact midpoint scratch includes both carry bits");
    static_assert(std::is_same<decltype(ec::fixed_splice(std::declval<compact>(), 0, std::declval<compact>(), 1)),
        ec::fixed_homogeneous<256, 256>>::value, "Compact splice uses its actual expression format");
    polygon a, b;
    a.outer() = {{0, 0}, {3, 3}, {3, 0}, {0, 0}};
    b.outer() = {{0, 1}, {2, 0}, {0, 0}, {0, 1}};
    ec::borrowed_context<root> context;
    ec::borrowed_stage<root> stage(context, a, b);
    typename root::operation_type operations[2];
    operations[0].seg_id = bg::segment_identifier(0, -1, -1, 0);
    operations[1].seg_id = bg::segment_identifier(1, -1, -1, 0);
    ec::borrowed_record<root> record;
    record.stage = &stage;
    record.operations = {{&operations[0], &operations[1]}};
    point crossing(record);
    record.numeric.init(crossing.producers());
    compact value;
    BOOST_TEST(crossing.copy_compact_homogeneous(value));
    BOOST_TEST(value.X == 2 && value.Y == 2 && value.W == 3);
    auto& storage = ec::operation_detail::require<C>().storage();
    storage.auxiliary.reserve(2);
    auto const native_before = storage.native.size();
    auto const before = storage.auxiliary.size();
    auto const middle = without_allocations([&] { return query::midpoint(crossing, point::borrow(a.outer()[0])); });
    BOOST_TEST_EQ(storage.auxiliary.size(), before + 1);
    auto const& cached_middle = storage.auxiliary.back();
    BOOST_TEST(cached_middle.X == 1 && cached_middle.Y == 1 && cached_middle.W == 3);
    check_oracle(middle, ec::exact_number(1) / 3, ec::exact_number(1) / 3);
    BOOST_TEST(middle.approximation(0).error > 0);
    native const target(1, 1);
    auto const corner = without_allocations([&] { return query::splice_coordinates(crossing, 0, target, 1); });
    BOOST_TEST_EQ(storage.auxiliary.size(), before + 2);
    auto const& cached_corner = storage.auxiliary.back();
    BOOST_TEST(cached_corner.X == 2 && cached_corner.Y == 3 && cached_corner.W == 3);
    check_oracle(corner, ec::exact_number(2) / 3, ec::exact_number(1));
    BOOST_TEST_EQ(storage.native.size(), native_before);
}

template <typename Policy>
void test_axis_origin()
{
    using native = bg::model::d2::point_xy<double>;
    using polygon = bg::model::polygon<native>;
    using root = ec::borrowed_root<double, Policy, polygon, polygon>;
    using point = ec::borrowed_point<root>;
    using query = ec::query_point<root>;
    for (bool horizontal : {false, true})
    {
        polygon a, b;
        a.outer() = {{2, -4}, {2, 4}, {0, 0}, {2, -4}};
        b.outer() = {{0, 0}, {4, horizontal ? 0.0 : 2.0}, {0, 4}, {0, 0}};
        ec::borrowed_context<root> context;
        ec::borrowed_stage<root> stage(context, a, b);
        typename root::operation_type operations[2];
        operations[0].seg_id = bg::segment_identifier(0, -1, -1, 0);
        operations[1].seg_id = bg::segment_identifier(1, -1, -1, 0);
        ec::borrowed_record<root> record;
        record.stage = &stage;
        record.operations = {{&operations[0], &operations[1]}};
        point crossing(record);
        record.numeric.init(crossing.producers());
        auto& storage = ec::operation_detail::require<double>().storage();
        auto const before = storage.auxiliary.size(), native_before = storage.native.size();
        native const target(4, 6);
        auto const corner = without_allocations([&] { return query::splice_coordinates(crossing, 0, target, 1); });
        check_oracle(corner, ec::exact_number(2), ec::exact_number(6));
        BOOST_TEST_EQ(storage.auxiliary.size(), before);
        BOOST_TEST_EQ(storage.native.size(), native_before);
        auto const middle = query::midpoint(crossing, target);
        check_oracle(middle, ec::exact_number(3), horizontal ? ec::exact_number(3) : ec::exact_number(7) / 2);
        BOOST_TEST_EQ(storage.auxiliary.size(), before + (horizontal ? 0 : 1));
    }
}

template <typename C, typename Policy>
void test_query()
{
    using native = bg::model::d2::point_xy<C>;
    using polygon = bg::model::polygon<native>;
    using root = ec::borrowed_root<C, Policy, polygon, polygon>;
    using other_root = ec::borrowed_root<C, Policy, polygon, bg::model::ring<native>>;
    using point = ec::borrowed_point<root>;
    using query = ec::query_point<root>;
    static_assert(has_midpoint<root, point, native>::value && has_midpoint<root, native, point>::value,
                  "Geometry/native operands are supported in either order");
    static_assert(has_midpoint<root, native, native>::value, "An explicitly bound query can use two native points");
    static_assert(!has_midpoint<root, query, point>::value && !has_midpoint<root, point, query>::value,
                  "A midpoint cannot recursively consume a query");
    static_assert(!has_splice<root, query, native>::value && !has_splice<root, native, query>::value,
                  "A splice cannot recursively consume a query");
    static_assert(!has_midpoint<root, ec::borrowed_point<other_root>, native>::value,
                  "A different Root is not a geometry operand of this query");
    static_assert(!std::is_default_constructible<query>::value
        && !std::is_constructible<query, typename query::homogeneous_type>::value,
        "Only the proved construction factories create queries");
    static_assert(!std::is_constructible<point, query>::value && !std::is_convertible<query, point>::value,
                  "Queries are not retained geometry references");
    struct expected_layout
    {
        ec::point_precision_detail::point_handle<C> handle;
        std::array<double, 2> approximate;
    };
    static_assert(sizeof(query) == sizeof(expected_layout), "Only a typed cache handle and two doubles");

    polygon a, b;
    a.outer() = {{0, 0}, {3, 3}, {3, 0}, {0, 0}};
    b.outer() = {{0, 3}, {3, 0}, {0, 0}, {0, 3}};
    ec::borrowed_context<root> context;
    ec::borrowed_stage<root> stage(context, a, b);
    typename root::operation_type operations[2];
    operations[0].seg_id = bg::segment_identifier(0, -1, -1, 0);
    operations[1].seg_id = bg::segment_identifier(1, -1, -1, 0);
    ec::borrowed_record<root> record;
    record.stage = &stage;
    record.operations = {{&operations[0], &operations[1]}};
    point crossing(record);
    record.numeric.init(crossing.producers());
    crossing.homogeneous();
    auto& storage = ec::operation_detail::require<C>().storage();
    storage.auxiliary.reserve(8);

    auto const origin = point::from_iterator(a.outer().cbegin());
    auto const middle = without_allocations([&] { return query::midpoint(crossing, origin); });
    auto const exact_middle = middle.homogeneous();
    BOOST_TEST(exact_middle.X == 3 && exact_middle.Y == 3 && exact_middle.W == 4);
    ec::fixed_homogeneous<2, 3> compact;
    BOOST_TEST(without_allocations([&] { return middle.copy_compact_homogeneous(compact); }));
    BOOST_TEST(compact.X == 3 && compact.Y == 3 && compact.W == 4);
    ec::fixed_homogeneous<2, 2> too_small;
    too_small.X = 2; too_small.Y = -1; too_small.W = 1;
    BOOST_TEST(!without_allocations([&] { return middle.copy_compact_homogeneous(too_small); }));
    BOOST_TEST(too_small.X == 2 && too_small.Y == -1 && too_small.W == 1);
    BOOST_TEST_EQ(middle.approximate(0), 0.75);
    BOOST_TEST_EQ(middle.approximation(0).approximate, middle.refine(0).approximate);
    BOOST_TEST(middle.approximation(0).error > 0);
    auto const count = storage.auxiliary.size();
    without_allocations([&] { return middle.refine(1); });
    BOOST_TEST_EQ(storage.auxiliary.size(), count);

    native const target(2, 6);
    auto const corner = without_allocations([&] { return query::splice_coordinates(crossing, 1, target, 0); });
    auto const exact_corner = corner.homogeneous();
    BOOST_TEST(exact_corner.X == 3 && exact_corner.Y == 4 && exact_corner.W == 2);
    auto const mixed = without_allocations([&] { return query::midpoint(target, crossing); });
    auto const exact_mixed = mixed.homogeneous();
    BOOST_TEST(exact_mixed.X == 7 && exact_mixed.Y == 15 && exact_mixed.W == 4);
    BOOST_TEST_THROWS(query::splice_coordinates(crossing, 2, target, 0), std::out_of_range);

#ifdef BOOST_GEOMETRY_TEST_REJECT_QUERY_BORROW
    auto const forbidden = point::borrow(middle);
    (void)forbidden;
#endif
#ifdef BOOST_GEOMETRY_TEST_REJECT_MIXED_ROOT_PREDICATE
    auto const other = ec::borrowed_point<other_root>::borrow(a.outer()[0]);
    (void)ec::static_compare_coordinate<0>(origin, other);
#endif
}

int main()
{
    test_query<std::int64_t, bg::strategies::relate::approximate_output>();
    test_query<std::int64_t, bg::strategies::relate::nearest_output>();
    test_query<double, bg::strategies::relate::approximate_output>();
    test_query<double, bg::strategies::relate::nearest_output>();
    test_native_queries<std::int8_t, bg::strategies::relate::approximate_output>();
    test_native_queries<std::int64_t, bg::strategies::relate::approximate_output>();
    test_native_queries<std::uint64_t, bg::strategies::relate::approximate_output>();
    test_native_queries<float, bg::strategies::relate::approximate_output>();
    test_native_queries<double, bg::strategies::relate::approximate_output>();
    test_native_queries<std::int64_t, bg::strategies::relate::nearest_output>();
    test_native_queries<float, bg::strategies::relate::nearest_output>();
    test_native_queries<double, bg::strategies::relate::nearest_output>();
    test_axis_origin<bg::strategies::relate::approximate_output>();
    test_axis_origin<bg::strategies::relate::nearest_output>();
    test_compact_rational_queries<std::int64_t, bg::strategies::relate::approximate_output>();
    test_compact_rational_queries<float, bg::strategies::relate::approximate_output>();
    test_compact_rational_queries<double, bg::strategies::relate::approximate_output>();
    test_compact_rational_queries<std::int64_t, bg::strategies::relate::nearest_output>();
    test_compact_rational_queries<float, bg::strategies::relate::nearest_output>();
    test_compact_rational_queries<double, bg::strategies::relate::nearest_output>();
    return boost::report_errors();
}
