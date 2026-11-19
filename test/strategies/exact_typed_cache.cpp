// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.

#include <array>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <new>
#include <thread>
#include <type_traits>

#include <boost/core/lightweight_test.hpp>
#include <boost/geometry/strategies/cartesian/exact/point_cache.hpp>
#include <boost/geometry/strategies/cartesian/exact/point_precision_cache.hpp>
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

namespace ec = boost::geometry::strategies::detail::exact_cartesian;
namespace op = ec::operation_detail;
namespace pc = ec::point_cache_detail;
namespace pp = ec::point_precision_detail;

template <typename C>
std::array<C, 8> crossing()
{
    return {{C(0), C(0), C(1), C(1), C(0), C(1), C(1), C(0)}};
}

template <typename C>
void test_scopes()
{
    BOOST_TEST(op::active<C>() == nullptr);
    BOOST_TEST(op::current<C>() == nullptr);
    allocations = 0;
    count_allocations = true;
    {
        ec::point_cache_scope<C> root;
        auto* const owner = op::active<C>();
        BOOST_TEST(owner != nullptr);
        BOOST_TEST(op::current<C>() == nullptr);
        pc::native_cache_handle empty;
        ec::interval bounds;
        BOOST_TEST(!pc::native_cache<C>::cached_bounds(empty, 0, bounds));
        prepare_exact_cache_repeat<C>();
        BOOST_TEST(op::current<C>() == nullptr);
        {
            ec::point_cache_scope<C> nested;
            BOOST_TEST(op::active<C>() == owner);
            BOOST_TEST(op::current<C>() == nullptr);
        }
        auto& context = op::require<C>();
        auto const& const_context = context;
        BOOST_TEST(&const_context.storage() == &context.template storage<C>());
        BOOST_TEST(context.storage().native.empty());
        BOOST_TEST(context.storage().auxiliary.empty());
        {
            ec::point_cache_scope<C> nested;
            BOOST_TEST(op::active<C>() == owner);
            BOOST_TEST(op::current<C>() == &context);
        }
        BOOST_TEST(op::current<C>() == &context);
    }
    count_allocations = false;
    BOOST_TEST_EQ(allocations, 0u);
    BOOST_TEST(op::active<C>() == nullptr);
    BOOST_TEST(op::current<C>() == nullptr);
    {
        ec::point_cache_scope<C> next;
        BOOST_TEST(op::require<C>().storage().native.empty());
        BOOST_TEST(op::require<C>().storage().auxiliary.empty());
    }
}

template <typename C>
void test_cache_entries()
{
    using cache = pc::native_cache<C>;
    using auxiliary_handle = pp::point_handle<C>;
    static_assert(sizeof(pc::native_cache_handle) == sizeof(std::size_t), "Only a local index");
    static_assert(sizeof(auxiliary_handle) == sizeof(std::size_t), "Only a typed local index");
    static_assert(std::is_trivially_copyable<pc::native_cache_handle>::value, "Copy preserves cache hint");
    static_assert(std::is_trivially_copyable<auxiliary_handle>::value, "Copy preserves cache hint");
    static_assert(std::is_same<decltype(std::declval<auxiliary_handle>().get()),
                               ec::auxiliary_point_type<C>>::value, "No double-sized promotion");

    ec::point_cache_scope<C> scope;
    auto& storage = op::require<C>().storage();
    pc::native_cache_handle native;
    auto const coordinates = crossing<C>();
    auto const value = cache::get(native, true, coordinates);
    auto auxiliary = auxiliary_handle::create(value);
    auto native_copy = native;
    pc::native_cache_handle native_assigned;
    native_assigned = native;
    auto auxiliary_copy = auxiliary;
    auxiliary_handle auxiliary_assigned;
    auxiliary_assigned = auxiliary;
    BOOST_TEST(!bool(auxiliary_handle()));
    BOOST_TEST(auxiliary.same(auxiliary_copy));
    BOOST_TEST(auxiliary.same(auxiliary_assigned));

    prepare_exact_cache_repeat<C>();
    auto const native_capacity = storage.native.capacity();
    auto const auxiliary_capacity = storage.auxiliary.capacity();
    ec::fixed_homogeneous<32, 32> compact;
    ec::interval bounds;
    allocations = 0;
    count_allocations = true;
    for (unsigned i = 0; i != 8; ++i)
    {
        auto const cached = cache::get(native_copy, true, coordinates);
        BOOST_TEST(cached.X * 2 == cached.W && cached.Y * 2 == cached.W);
        BOOST_TEST(cache::copy_compact(native_assigned, true, coordinates, compact));
        BOOST_TEST(cache::cached_bounds(native_copy, 0, bounds));
        BOOST_TEST_EQ(bounds.approximate, 0.5);
        BOOST_TEST(bounds.lower() <= 0.5 && bounds.upper() >= 0.5);
        BOOST_TEST(auxiliary_copy.copy_compact(compact));
        BOOST_TEST(compact.X * 2 == compact.W && compact.Y * 2 == compact.W);
        BOOST_TEST(auxiliary_assigned.get().X * 2 == auxiliary_assigned.get().W);
    }
    BOOST_TEST_EQ(cache::size(), 1u);
    // Add previously uncached values without permitting vector growth.
    for (unsigned i = 0; i != 8; ++i)
    {
        pc::native_cache_handle next;
        cache::get(next, true, coordinates);
        auxiliary_handle::create(value);
    }
    count_allocations = false;
    BOOST_TEST_EQ(allocations, 0u);
    BOOST_TEST_EQ(storage.native.capacity(), native_capacity);
    BOOST_TEST_EQ(storage.auxiliary.capacity(), auxiliary_capacity);

    // Growth must not invalidate either original or copied index handles.
    for (unsigned i = 0; i != 32; ++i)
    {
        pc::native_cache_handle next;
        cache::get(next, true, coordinates);
        auxiliary_handle::create(value);
    }
    BOOST_TEST(storage.native.capacity() > native_capacity);
    BOOST_TEST(storage.auxiliary.capacity() > auxiliary_capacity);
    BOOST_TEST(auxiliary_copy.copy_compact(compact));
    BOOST_TEST(compact.X * 2 == compact.W);
    BOOST_TEST(cache::copy_compact(native_copy, true, coordinates, compact));
    BOOST_TEST(compact.X * 2 == compact.W);
    BOOST_TEST_EQ(cache::size(), 41u);
}

template <typename C>
void test_capacity()
{
    constexpr auto n = ec::native_fixed_traits<C>::intersection_numerator_bits;
    constexpr auto d = ec::native_fixed_traits<C>::intersection_denominator_bits;
    static_assert(ec::precision_numerator_bits<C> == n + d + 1, "Splice/midpoint numerator");
    static_assert(ec::precision_denominator_bits<C> == 2 * d + 1, "Splice/midpoint denominator");
    static_assert(std::is_same<ec::precision_homogeneous_point<C>,
        ec::fixed_homogeneous<n + d + 1, 2 * d + 1>>::value, "Coordinate/degree-derived payload");

    ec::point_cache_scope<C> scope;
    ec::fixed_homogeneous<n + d + 2, 2 * d + 2> wide;
    wide.X = 1;
    wide.X <<= n + d + 1;
    bool rejected = false;
    try { pp::point_handle<C>::create(wide); }
    catch (std::overflow_error const&) { rejected = true; }
    BOOST_TEST(rejected);
    BOOST_TEST(op::require<C>().storage().auxiliary.empty());
    wide.X >>= 1;
    auto const handle = pp::point_handle<C>::create(wide);
    BOOST_TEST(handle.get().X == wide.X);
    ec::fixed_homogeneous<n, d> compact;
    BOOST_TEST(!handle.copy_compact(compact));

    wide.X = 0;
    wide.W = 1;
    wide.W <<= 2 * d + 1;
    rejected = false;
    try { pp::point_handle<C>::create(wide); }
    catch (std::overflow_error const&) { rejected = true; }
    BOOST_TEST(rejected);
    BOOST_TEST_EQ(op::require<C>().storage().auxiliary.size(), 1u);
}

void test_separate_coordinates()
{
    ec::point_cache_scope<> floating;
    auto const floating_value = ec::fixed_input_homogeneous(1.0, 2.0);
    auto const floating_handle = pp::point_handle<double>::create(floating_value);
    auto* const floating_owner = op::current<double>();
    {
        ec::point_cache_scope<std::int16_t> integral;
        BOOST_TEST(op::current<std::int16_t>() == nullptr);
        prepare_exact_cache_repeat<std::int16_t>();
        BOOST_TEST(op::current<std::int16_t>() == nullptr);
        auto const integral_handle = pp::point_handle<std::int16_t>::create(ec::fixed_input_homogeneous(3, 4));
        BOOST_TEST_EQ(integral_handle.get().X, 3);
        BOOST_TEST_EQ(floating_handle.get().X, 1);
        BOOST_TEST(op::current<double>() == floating_owner);
        BOOST_TEST_EQ(op::require<double>().storage().auxiliary.size(), 1u);
        BOOST_TEST_EQ(op::require<std::int16_t>().storage().auxiliary.size(), 1u);
    }
    BOOST_TEST(op::current<std::int16_t>() == nullptr);
    BOOST_TEST_EQ(floating_handle.get().Y, 2);
}

void test_expansion_cache()
{
    using cache = pc::native_cache<double>;
    ec::point_cache_scope<> scope;
    op::require<double>().storage().native.reserve(1);
    auto const first = crossing<double>();
    std::array<double, 8> const second{{0, 0, 1, 1, 0, 2, 2, 0}};
    pc::native_cache_handle a, b;
    int result = 0;
    // Obtaining b grows the vector after a's expansion has been evaluated.
    BOOST_TEST(cache::compare(a, first, b, second, 0, result));
    BOOST_TEST_EQ(result, -1);
    auto copy = a;
    prepare_exact_cache_repeat();
    ec::interval bounds;
    allocations = 0;
    count_allocations = true;
    BOOST_TEST(cache::refine(copy, first, 0, bounds));
    BOOST_TEST(cache::cached_bounds(a, 0, bounds));
    BOOST_TEST_EQ(bounds.approximate, 0.5);
    BOOST_TEST(cache::compare(copy, first, 0, 0.5, result));
    BOOST_TEST_EQ(result, 0);
    BOOST_TEST(cache::compare(copy, first, b, second, 0, result));
    BOOST_TEST_EQ(result, -1);
    count_allocations = false;
    BOOST_TEST_EQ(allocations, 0u);
    BOOST_TEST_EQ(cache::size(), 2u);
}

void test_thread_isolation()
{
    ec::point_cache_scope<double> scope;
    auto* const owner = &op::require<double>();
    bool isolated = false;
    std::thread worker([&]
    {
        bool const starts_empty = op::active<double>() == nullptr;
        ec::point_cache_scope<double> independent;
        auto& context = op::require<double>();
        isolated = starts_empty && &context != owner && context.storage().native.empty();
        pc::native_cache_handle handle;
        auto const value = pc::native_cache<double>::get(handle, true, crossing<double>());
        isolated = isolated && value.X * 2 == value.W && context.storage().native.size() == 1;
    });
    worker.join();
    BOOST_TEST(isolated);
    BOOST_TEST(op::current<double>() == owner);
    BOOST_TEST(owner->storage().native.empty());
}

template <typename C>
void test_coordinate()
{
    test_scopes<C>();
    test_cache_entries<C>();
    test_capacity<C>();
}

int main()
{
    test_coordinate<bool>();
    test_coordinate<std::int8_t>();
    test_coordinate<std::uint8_t>();
    test_coordinate<std::int16_t>();
    test_coordinate<std::uint16_t>();
    test_coordinate<std::int32_t>();
    test_coordinate<std::uint32_t>();
    test_coordinate<std::int64_t>();
    test_coordinate<std::uint64_t>();
    test_coordinate<float>();
    test_coordinate<double>();
    test_separate_coordinates();
    test_expansion_cache();
    test_thread_isolation();
    return boost::report_errors();
}
