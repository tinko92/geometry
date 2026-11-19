// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.

#include <boost/core/lightweight_test.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/strategies/cartesian/exact/borrowed_clipping.hpp>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <new>
#include <random>
#include <vector>

namespace bg = boost::geometry;
namespace ec = bg::strategies::detail::exact_cartesian;
namespace policy = bg::strategies::relate;
using rational = boost::multiprecision::cpp_rational;

namespace
{
bool count_allocations = false;
std::size_t allocations = 0;
}

void* operator new(std::size_t size)
{
    if (count_allocations) ++allocations;
    if (void* value = std::malloc(size == 0 ? 1 : size)) return value;
    throw std::bad_alloc();
}
void operator delete(void* value) noexcept { std::free(value); }
void* operator new[](std::size_t size) { return ::operator new(size); }
void operator delete[](void* value) noexcept { std::free(value); }
#if defined(__cpp_sized_deallocation)
void operator delete(void* value, std::size_t) noexcept { std::free(value); }
void operator delete[](void* value, std::size_t) noexcept { std::free(value); }
#endif

template <typename Root>
struct clipping_adapter : bg::strategies::relate::cartesian<>
{
    ec::borrowed_context<Root>& context;
    explicit clipping_adapter(ec::borrowed_context<Root>& context) : context(context) {}
    template <typename Point, typename Box>
    auto clipping(Box const& box) const { return ec::borrowed_clipping<Root>(context, box); }
};

template <typename C, typename Policy>
struct types
{
    using native = bg::model::d2::point_xy<C>;
    using line = bg::model::linestring<native>;
    using box = bg::model::box<native>;
    using root = ec::borrowed_root<C, Policy, line, box>;
    using point = ec::borrowed_point<root>;
    using working_line = typename ec::borrowed_geometry_types<root>::line;
};

template <typename Point>
void check_point(Point const& point, rational const& x, rational const& y)
{
    BOOST_TEST(point.exact_coordinate(0) == x);
    BOOST_TEST(point.exact_coordinate(1) == y);
}

template <typename Root>
void check_definition(ec::borrowed_point<Root> const& point,
                       ec::original_edge input, std::size_t box, std::size_t face)
{
    BOOST_TEST(!point.is_native());
    if (point.is_native()) return;
    auto const& record = point.record();
    BOOST_TEST(record.stage->finalized);
    BOOST_TEST(!record.is_seed);
    BOOST_TEST(record.stage->model[0] == nullptr && record.stage->model[1] == nullptr);
    BOOST_TEST(record.stage->support(record.operations[0]->seg_id) == input);
    BOOST_TEST((record.stage->support(record.operations[1]->seg_id) == ec::original_edge{box, face}));
}

template <typename C, typename Policy>
void test_lifetime_and_original_supports()
{
    using T = types<C, Policy>;
    using root = typename T::root;
    typename T::line const original{{-3, -4}, {4, 3}};
    typename T::box const box({0, 0}, {2, 2});
    typename T::box const unrelated_box1({10, 10}, {11, 11});
    typename T::box const unrelated_box2({20, 20}, {21, 21});
    ec::borrowed_context<root> context;
    context.register_range(unrelated_box1);
    context.register_range(unrelated_box2);
    auto const input_id = context.register_range(original);
    auto const box_id = context.register_range(box);
    BOOST_TEST(input_id > 1 && box_id > 1);
    ec::original_edge const input{input_id, 0};
    std::vector<typename T::working_line> output;
    {
        ec::borrowed_stage<root> unrelated(context, original, box);
        context.active_stage = &unrelated;
        bg::intersection(original, box, output, clipping_adapter<root>(context));
        BOOST_TEST(context.active_stage == &unrelated);
        context.active_stage = nullptr;
    }
    BOOST_TEST_EQ(output.size(), 1u);
    if (output.empty()) return;
    auto& line = output.front();
    BOOST_TEST_EQ(line.size(), 2u);
    if (line.size() != 2) return;
    // Left clipping first creates (0,-1), then bottom clipping replaces it.
    // The replacement still uses the native input edge, never the derived line.
    check_point(line.front(), 1, 0);
    check_point(line.back(), 2, 1);
    check_definition(line.front(), input, box_id, 3);
    check_definition(line.back(), input, box_id, 2);
    BOOST_TEST(!line.incoming.front().valid());
    BOOST_TEST(line.incoming.back() == input);
    BOOST_TEST(line.pending == input);
    auto const producers = line.front().producers();
    BOOST_TEST_EQ(producers[0], C(-3)); BOOST_TEST_EQ(producers[1], C(-4));
    BOOST_TEST_EQ(producers[2], C(4)); BOOST_TEST_EQ(producers[3], C(3));

    if constexpr (std::is_floating_point<C>::value)
    {
        typename T::box const narrower({1, 0}, {1.5, 2});
        std::vector<typename T::working_line> second;
        {
            ec::borrowed_clipping<root> clipping(context, narrower);
            bg::detail::intersection::clip_range_with_box<typename T::working_line>(
                narrower, line, std::back_inserter(second), clipping);
        }
        line.clear();
        BOOST_TEST_EQ(second.size(), 1u);
        if (second.empty()) return;
        auto const& recut = second.front();
        check_point(recut.front(), 1, 0);
        check_point(recut.back(), rational(3) / 2, rational(1) / 2);
        check_definition(recut.back(), input, context.register_range(narrower), 2);
        BOOST_TEST(recut.incoming.back() == input);
        auto const original_producers = recut.back().producers();
        for (std::size_t i = 0; i < 4; ++i) BOOST_TEST_EQ(original_producers[i], producers[i]);
    }
}

template <typename C, typename Policy>
void test_metadata_and_no_allocation()
{
    using T = types<C, Policy>;
    using root = typename T::root;
    typename T::line const input{{-2, 0}, {1, 0}, {4, 0}, {4, 0}};
    typename T::box const box({0, -1}, {2, 1});
    ec::borrowed_context<root> context;
    auto const range = context.register_range(input);
    std::vector<typename T::working_line> output;
    bg::intersection(input, box, output, clipping_adapter<root>(context));
    BOOST_TEST_EQ(output.size(), 1u);
    if (output.empty()) return;
    auto const& line = output.front();
    BOOST_TEST_EQ(line.size(), 3u);
    if (line.size() != 3) return;
    BOOST_TEST((line.incoming[1] == ec::original_edge{range, 0}));
    BOOST_TEST((line.incoming[2] == ec::original_edge{range, 1}));
    check_point(line[0], 0, 0); check_point(line[1], 1, 0); check_point(line[2], 2, 0);
    BOOST_TEST(line[1].is_native());

    // Registration is outside the measured primitive. Native predicates and an
    // unchanged clip must not introduce an owning point or a record allocation.
    ec::borrowed_clipping<root> clipping(context, box);
    clipping.prepare_segment(input, 1);
    auto a = context.native_vertex(input, 1), b = a;
    bg::model::referring_segment<typename T::point> segment(a, b);
    bool first = true, second = true;
    allocations = 0;
    count_allocations = true;
    bool const visible = clipping.clip_segment(box, segment, first, second);
    count_allocations = false;
    BOOST_TEST(visible && !first && !second);
    BOOST_TEST_EQ(allocations, 0u);
    BOOST_TEST(a.is_native() && b.is_native());

    // Construct=false preserves the original endpoints and flags.
    a = context.native_vertex(input, 0); b = context.native_vertex(input, 1);
    clipping.prepare_segment(input, 0);
    first = false; second = true;
    BOOST_TEST(clipping.template clip_segment<false>(box, segment, first, second));
    BOOST_TEST(!first && second && a.is_native() && b.is_native());
    check_point(a, -2, 0); check_point(b, 1, 0);
}

template <typename C, typename Policy>
void test_closed_support()
{
    using native = bg::model::d2::point_xy<C>;
    using ring = bg::model::ring<native>;
    using box_type = bg::model::box<native>;
    using root = ec::borrowed_root<C, Policy, ring, box_type>;
    using point = ec::borrowed_point<root>;
    using working_ring = typename ec::borrowed_geometry_types<root>::template ring<true, true>;
    using open_ring = typename ec::borrowed_geometry_types<root>::template ring<true, false>;
    ring const input{{0, 0}, {0, 4}, {4, 4}, {4, 0}, {0, 0}};
    box_type const box({1, -1}, {3, 1});
    ec::borrowed_context<root> context;
    ec::original_edge const original{context.register_range(input), 3};
    working_ring working;
    open_ring open;
    ec::borrowed_geometry_converter<root> converter(context);
    converter.apply(input, working);
    converter.apply(input, open);
    working.incoming.front() = {}; // Only the appended closing label is authoritative.
    for (bool closed : {true, false})
    {
        ec::borrowed_clipping<root> clipping(context, box);
        point a, b;
        if (closed)
        {
            clipping.prepare_segment(working, 3);
            a = working[3]; b = working[4];
        }
        else
        {
            clipping.prepare_segment(open, 3);
            a = open[3]; b = open[0];
        }
        bg::model::referring_segment<point> segment(a, b);
        bool first = false, second = false;
        BOOST_TEST(clipping.clip_segment(box, segment, first, second));
        BOOST_TEST(first && second);
        check_point(a, 3, 0); check_point(b, 1, 0);
        check_definition(a, original, context.register_range(box), 2);
        check_definition(b, original, context.register_range(box), 0);
    }
}

template <typename C, typename Policy>
void compare_native(std::array<C, 4> const& coordinates, C width = C(4), C height = C(4))
{
    using T = types<C, Policy>;
    using root = typename T::root;
    typename T::line const input{{coordinates[0], coordinates[1]}, {coordinates[2], coordinates[3]}};
    typename T::box const box({0, 0}, {width, height});
    rational lo = 0, hi = 1;
    bool visible = true;
    for (std::size_t d = 0; d < 2; ++d)
    {
        rational const origin(coordinates[d]), delta = rational(coordinates[d + 2]) - origin;
        rational const bound = d == 0 ? rational(width) : rational(height);
        if (delta == 0) { if (origin < 0 || origin > bound) visible = false; }
        else
        {
            rational a = -origin / delta, b = (bound - origin) / delta;
            if (a > b) std::swap(a, b);
            lo = (std::max)(lo, a); hi = (std::min)(hi, b);
        }
    }
    visible = visible && lo <= hi;
    bool const isolated = coordinates[0] == coordinates[2] && coordinates[1] == coordinates[3];
    bool const line_visible = visible && lo < hi && !isolated;
    ec::borrowed_context<root> context;
    std::vector<typename T::working_line> output;
    bg::intersection(input, box, output, clipping_adapter<root>(context));
    BOOST_TEST_EQ(output.size(), line_visible ? 1u : 0u);
    if (!line_visible || output.empty()) return;
    auto const& line = output.front();
    BOOST_TEST_EQ(line.size(), 2u);
    if (line.size() != 2) return;
    rational const dx = rational(coordinates[2]) - coordinates[0];
    rational const dy = rational(coordinates[3]) - coordinates[1];
    check_point(line.front(), coordinates[0] + lo * dx, coordinates[1] + lo * dy);
    check_point(line.back(), coordinates[0] + hi * dx, coordinates[1] + hi * dy);
    BOOST_TEST((line.incoming.back() == ec::original_edge{context.register_range(input), 0}));
}

template <typename C, typename Policy>
void test_scaled(C origin, C step)
{
    using T = types<C, Policy>;
    using root = typename T::root;
    typename T::line const input{{origin - 3 * step, origin - 4 * step},
                                 {origin + 4 * step, origin + 3 * step}};
    typename T::box const box({origin, origin}, {origin + 2 * step, origin + 2 * step});
    ec::borrowed_context<root> context;
    std::vector<typename T::working_line> output;
    bg::intersection(input, box, output, clipping_adapter<root>(context));
    BOOST_TEST_EQ(output.size(), 1u);
    if (output.empty()) return;
    rational const base(origin), unit(step);
    check_point(output.front().front(), base + unit, base);
    check_point(output.front().back(), base + 2 * unit, base + unit);
}

template <typename C, typename Policy>
void test_cases()
{
    test_lifetime_and_original_supports<C, Policy>();
    test_metadata_and_no_allocation<C, Policy>();
    test_closed_support<C, Policy>();
    if constexpr (std::is_integral<C>::value) test_scaled<C, Policy>(C(1) << 60, 1);
    else
    {
        test_scaled<C, Policy>(std::ldexp(C(1), 500), std::ldexp(C(1), 450));
        test_scaled<C, Policy>(0, (std::numeric_limits<C>::denorm_min)());
    }
    for (auto const& value : {std::array<C, 4>{{-3, 1, 7, 2}}, {{7, 2, -3, 1}},
            {{1, -3, 2, 7}}, {{2, 7, 1, -3}}, {{-3, -4, 7, 6}}, {{-2, 0, 0, 0}},
            {{0, 0, 4, 0}}, {{-2, 0, 6, 0}}, {{0, -2, 0, 6}},
            {{1, 1, 1, 1}}, {{-1, -1, -1, -1}}, {{0, 1, -1, 1}}})
    {
        compare_native<C, Policy>(value);
        compare_native<C, Policy>(value, 0, 4);
        compare_native<C, Policy>(value, 4, 0);
        compare_native<C, Policy>(value, 0, 0);
    }
    std::mt19937 random(31829);
    for (int i = 0; i < 200; ++i)
    {
        std::array<C, 4> value;
        for (auto& coordinate : value) coordinate = C(int(random() % 17) - 6);
        compare_native<C, Policy>(value);
        compare_native<C, Policy>(value, 0, 4);
        compare_native<C, Policy>(value, 4, 0);
        compare_native<C, Policy>(value, 0, 0);
    }
}

int main()
{
    test_cases<double, policy::approximate_output>();
    test_cases<double, policy::eager_nearest_output>();
    test_cases<double, policy::deferred_nearest_output>();
    test_cases<std::int64_t, policy::approximate_output>();
    test_cases<std::int64_t, policy::eager_nearest_output>();
    test_cases<std::int64_t, policy::deferred_nearest_output>();
    return boost::report_errors();
}
