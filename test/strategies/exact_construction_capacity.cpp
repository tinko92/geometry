// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.

#include <boost/core/lightweight_test.hpp>
#include <boost/geometry/strategies/cartesian/exact/construction_capacity.hpp>
#include <boost/geometry/strategies/cartesian/exact/fixed_measures.hpp>

#include <array>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <new>
#include <type_traits>
#include <utility>

namespace ec = boost::geometry::strategies::detail::exact_cartesian;
using family = ec::construction_family;
template <typename C, family F> using capacity = ec::construction_capacity<C, F>;

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

template <typename C, typename = void> struct has_capacity : std::false_type {};
template <typename C>
struct has_capacity<C, std::void_t<typename capacity<C, family::native_vertex>::homogeneous_type>>
    : std::true_type {};
static_assert(!has_capacity<long double>::value, "No unproved extended format");
static_assert(!has_capacity<ec::fixed_integer<128>>::value, "No arbitrary integer leaf");

template <typename C, std::size_t B, std::size_t S>
void check_formats()
{
    using native = capacity<C, family::native_vertex>;
    using root = capacity<C, family::original_support_intersection>;
    using midpoint = capacity<C, family::midpoint_query>;
    using splice = capacity<C, family::splice_query>;
    using clip = capacity<C, family::native_axis_clip_query>;
    using native_traits = ec::native_fixed_traits<C>;
    static_assert(native::numerator_bits == B && native::denominator_bits == S + 1, "Native grid");
    static_assert(root::numerator_bits == 3 * B + 3, "Original support numerator");
    static_assert(root::denominator_bits == 2 * B + S + 3, "Original support denominator");
    static_assert(midpoint::numerator_bits == 5 * B + S + 7, "One midpoint numerator");
    static_assert(midpoint::denominator_bits == 4 * B + 2 * S + 7, "One midpoint denominator");
    static_assert(splice::numerator_bits == 5 * B + S + 6, "One splice numerator");
    static_assert(splice::denominator_bits == 4 * B + 2 * S + 6, "One splice denominator");
    static_assert(clip::numerator_bits == 2 * B + 2, "Grid-cancelled native clip numerator");
    static_assert(clip::denominator_bits == B + S + 1, "Grid-cancelled native clip denominator");
    static_assert(native::is_geometry_point && root::is_geometry_point, "Retained geometry families");
    static_assert(!midpoint::is_geometry_point && !splice::is_geometry_point && !clip::is_geometry_point,
                  "Queries do not certify original-edge provenance");
    static_assert(std::is_same<typename native::homogeneous_type, typename native_traits::input_type>::value,
                  "Reuse native storage");
    using point = typename root::homogeneous_type;
    static_assert(std::is_same<point, typename native_traits::intersection_type>::value, "Reuse root storage");
    static_assert(std::is_same<decltype(ec::fixed_midpoint(std::declval<point>(), std::declval<point>())),
                               typename midpoint::homogeneous_type>::value, "Actual midpoint return format");
    static_assert(std::is_same<decltype(ec::fixed_splice(std::declval<point>(), 0, std::declval<point>(), 1)),
                               typename splice::homogeneous_type>::value, "Actual splice return format");
    using scratch = ec::construction_predicate_capacity<root, root, midpoint>;
    static_assert(scratch::coordinate_compare_bits == root::numerator_bits + root::denominator_bits + 1,
                  "Comparison uses its two operands");
    static_assert(scratch::orientation_bits == 2 * root::numerator_bits + 2 * root::denominator_bits + 4,
                  "Six-term root/root/midpoint determinant");
    static_assert(scratch::direction_bits == 2 * root::numerator_bits + 3 * root::denominator_bits + 4,
                  "Two affine differences and a dot-product carry");
    static_assert(std::numeric_limits<typename scratch::orientation_type>::digits == scratch::orientation_bits,
                  "Scratch is checked fixed capacity, not modular limb padding");
    using measures = ec::fixed_measure_traits<root::numerator_bits, root::denominator_bits>;
    static_assert(std::is_same<decltype(ec::fixed_squared_distance(std::declval<point>(), std::declval<point>())),
                               typename measures::pair_type>::value, "Reuse existing measure capacity");
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

// Full-format values test arithmetic bounds, not geometric reachability.
template <typename C>
void test_saturated_scratch()
{
    using root = capacity<C, family::original_support_intersection>;
    using midpoint = capacity<C, family::midpoint_query>;
    using splice = capacity<C, family::splice_query>;
    typename root::homogeneous_type a;
    a.X = a.Y = maximum_fixed<root::numerator_bits>();
    a.W = maximum_fixed<root::denominator_bits>();
    auto b = a, c = a;
    b.X = -b.X;
    c.Y = -c.Y;
    auto middle = without_allocations([&] { return ec::fixed_midpoint(a, b); });
    auto corner = without_allocations([&] { return ec::fixed_splice(a, 0, c, 1); });
    BOOST_TEST(middle.X == 0 && middle.W > 0 && middle.Y > 0);
    BOOST_TEST_EQ(without_allocations([&] { return ec::fixed_coordinate_compare(a, b, 0); }), 1);
    BOOST_TEST_EQ(without_allocations([&] { return ec::fixed_orientation(a, b, c); }), 1);
    BOOST_TEST_EQ(without_allocations([&] { return ec::fixed_orientation(a, middle, corner); }), 1);
    BOOST_TEST_EQ(without_allocations([&]
    {
        return ec::fixed_orientation(ec::fixed_midpoint(a, a),
            ec::fixed_midpoint(b, b), ec::fixed_midpoint(c, c));
    }), 1);
    BOOST_TEST_EQ(without_allocations([&] { return ec::fixed_coordinate_compare(c, corner, 1); }), 0);
    auto const distance = without_allocations([&] { return ec::fixed_squared_distance(a, b); });
    auto const determinant = without_allocations([&] { return ec::fixed_signed_determinant(a, b, c); });
    BOOST_TEST_EQ(without_allocations([&] { return ec::fixed_pair_compare(distance, determinant); }), 0);
    BOOST_TEST(without_allocations([&]
    {
        ec::fixed_normalize(middle);
        ec::fixed_normalize(corner);
        typename midpoint::homogeneous_type m;
        typename splice::homogeneous_type s;
        return ec::try_fixed_narrow(middle, m) && ec::try_fixed_narrow(corner, s);
    }));
}

// This direct formula is the precondition of native_axis_clip_query. It is not
// a claim that the old generic four-point intersection uses this raw format.
template <std::size_t Axis, typename C>
auto native_clip(std::array<C, 4> const& support, C bound)
    -> typename capacity<C, family::native_axis_clip_query>::homogeneous_type
{
    using traits = capacity<C, family::native_axis_clip_query>;
    auto const line = ec::fixed_native_line(support);
    auto const p = ec::fixed_value(bound);
    typename traits::homogeneous_type result;
    if constexpr (Axis == 0)
    {
        result.X = ec::fixed_multiply<traits::numerator_bits>(line.B, p);
        result.Y = -(ec::fixed_multiply<traits::numerator_bits>(line.A, p) + line.C);
        result.W = line.B;
    }
    else
    {
        result.X = -(ec::fixed_multiply<traits::numerator_bits>(line.B, p) + line.C);
        result.Y = ec::fixed_multiply<traits::numerator_bits>(line.A, p);
        result.W = line.A;
    }
    if (result.W < 0)
    {
        result.X = -result.X;
        result.Y = -result.Y;
        result.W = -result.W;
    }
    result.W <<= ec::native_fixed_traits<C>::scale_bits;
    ec::fixed_normalize(result);
    return result;
}

template <std::size_t Axis, typename C>
void check_clip(std::array<C, 4> support, C bound)
{
    auto const result = without_allocations([&] { return native_clip<Axis>(support, bound); });
    std::array<C, 8> values{{support[0], support[1], support[2], support[3],
        Axis == 0 ? bound : C(0), Axis == 1 ? bound : C(0),
        Axis == 0 ? bound : C(1), Axis == 1 ? bound : C(1)}};
    auto const expected = without_allocations([&] { return ec::direct_native_homogeneous(values); });
    BOOST_TEST_EQ(without_allocations([&] { return ec::fixed_coordinate_compare(result, expected, 0); }), 0);
    BOOST_TEST_EQ(without_allocations([&] { return ec::fixed_coordinate_compare(result, expected, 1); }), 0);
    // Positive-W primitive triples are unique, not just coordinate-equivalent.
    BOOST_TEST(result.X == expected.X && result.Y == expected.Y && result.W == expected.W);
    std::swap(support[0], support[2]);
    std::swap(support[1], support[3]);
    auto const reversed = without_allocations([&] { return native_clip<Axis>(support, bound); });
    BOOST_TEST(result.X == reversed.X && result.Y == reversed.Y && result.W == reversed.W);
}

template <typename C>
void test_native_clips()
{
    C const low = (std::numeric_limits<C>::lowest)();
    C const high = (std::numeric_limits<C>::max)();
    check_clip<0, C>({{low, low, high, high}}, C(0));
    check_clip<1, C>({{low, high, high, low}}, high / C(2));
    if constexpr (std::is_floating_point<C>::value)
    {
        C const tiny = (std::numeric_limits<C>::denorm_min)();
        check_clip<0, C>({{0, 0, high, tiny}}, tiny);
        check_clip<1, C>({{0, 0, tiny, high}}, high);
    }
}

template <typename C>
void test_normalized_splices()
{
    auto const a = ec::direct_native_homogeneous<C>({{0, 0, 1, 1, 0, 1, 1, 0}});
    auto const b = ec::direct_native_homogeneous<C>({{0, 0, 1, 1, 0, 1, 2, 0}});
    auto current = ec::fixed_splice(a, 0, b, 1);
    BOOST_TEST(without_allocations([&]
    {
        for (unsigned i = 0; i != 16; ++i)
        {
            auto next = ec::fixed_splice(current, 0, b, 1);
            ec::fixed_normalize(next);
            if (!ec::try_fixed_narrow(next, current)) return false;
        }
        return current.X == 3 && current.Y == 4 && current.W == 6;
    }));
}

void test_recursive_midpoint_nonclosure()
{
    using query = capacity<std::int64_t, family::midpoint_query>;
    ec::fixed_homogeneous<2, query::denominator_bits + 2> current, zero;
    current.X = 1;
    BOOST_TEST(without_allocations([&]
    {
        for (std::size_t i = 1; i <= query::denominator_bits; ++i)
        {
            auto next = ec::fixed_midpoint(current, zero);
            ec::fixed_normalize(next);
            if (!ec::try_fixed_narrow(next, current)) return false;
            typename query::homogeneous_type query_storage;
            bool const fits = ec::try_fixed_narrow(current, query_storage);
            if (fits != (i < query::denominator_bits)) return false;
        }
        return current.X == 1 && current.Y == 0;
    }));
}

int main()
{
    check_formats<std::int8_t, 8, 0>();
    check_formats<std::uint8_t, 8, 0>();
    check_formats<std::int16_t, 16, 0>();
    check_formats<std::uint16_t, 16, 0>();
    check_formats<std::int32_t, 32, 0>();
    check_formats<std::uint32_t, 32, 0>();
    check_formats<std::int64_t, 64, 0>();
    check_formats<std::uint64_t, 64, 0>();
    check_formats<float, 277, 149>();
    check_formats<double, 2098, 1074>();
    test_saturated_scratch<std::int64_t>();
    test_saturated_scratch<std::uint64_t>();
    test_saturated_scratch<float>();
    test_saturated_scratch<double>();
    test_native_clips<std::int64_t>();
    test_native_clips<std::uint64_t>();
    test_native_clips<float>();
    test_native_clips<double>();
    test_normalized_splices<std::int64_t>();
    test_normalized_splices<float>();
    test_normalized_splices<double>();
    test_recursive_midpoint_nonclosure();
    return boost::report_errors();
}
