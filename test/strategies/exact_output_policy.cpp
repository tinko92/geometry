// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.

#ifndef BOOST_GEOMETRY_EXACT_FILTER_STATISTICS
#define BOOST_GEOMETRY_EXACT_FILTER_STATISTICS
#endif
#include <boost/core/lightweight_test.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/adapted/std_array.hpp>
#include <boost/geometry/strategies/relate/cartesian_exact.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>

#include <array>
#include <cfenv>
#include <cstdlib>
#include <new>
#include <tuple>
#include <vector>

BOOST_GEOMETRY_REGISTER_STD_ARRAY_CS(::boost::geometry::cs::cartesian)

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
void operator delete(void* p) noexcept { std::free(p); }
void* operator new[](std::size_t n) { return ::operator new(n); }
void operator delete[](void* p) noexcept { ::operator delete(p); }
#if defined(__cpp_sized_deallocation)
void operator delete(void* p, std::size_t) noexcept { std::free(p); }
void operator delete[](void* p, std::size_t) noexcept { std::free(p); }
#endif

namespace bg = boost::geometry;
namespace ec = bg::strategies::detail::exact_cartesian;
namespace sr = bg::strategies::relate;
using exact_point = ec::exact_point<double>;
template <typename C> using point = bg::model::d2::point_xy<C>;
using default_strategy = sr::cartesian_exact<double>;
using nearest_strategy = sr::cartesian_exact<double, sr::nearest_output>;

static_assert(std::is_same<default_strategy::output_policy, sr::approximate_output>::value,
              "Approximate output is the public default");
static_assert(std::is_same<default_strategy::intersection_point_type<point<double>>,
                          nearest_strategy::intersection_point_type<point<double>>>::value,
              "This stage does not change the point representation");

exact_point crossing()
{
    return ec::construct_intersection(exact_point(43, 7), exact_point(-24, 81),
        exact_point(-78, 22), exact_point(63, 95));
}

template <typename Output>
point<Output> reference()
{
    using high = boost::multiprecision::cpp_bin_float_100;
    // Cramer's rule gives (-74541, 917427)/15325. This independent high-
    // precision division has ample precision to resolve either target format.
    return {(high(-74541) / 15325).template convert_to<Output>(),
            (high(917427) / 15325).template convert_to<Output>()};
}

void same_statistics(ec::filter_statistics const& a, ec::filter_statistics const& b)
{
    BOOST_TEST_EQ(a.native_exact_values, b.native_exact_values);
    BOOST_TEST_EQ(a.native_rationals, b.native_rationals);
    BOOST_TEST_EQ(a.lazy_exact_values, b.lazy_exact_values);
    BOOST_TEST_EQ(a.point_exact_values, b.point_exact_values);
    BOOST_TEST_EQ(a.point_coordinate_values, b.point_coordinate_values);
    BOOST_TEST_EQ(a.homogeneous_values, b.homogeneous_values);
    BOOST_TEST_EQ(a.ratio_values, b.ratio_values);
    BOOST_TEST_EQ(a.construction_filters, b.construction_filters);
    BOOST_TEST_EQ(a.native_point_constructions, b.native_point_constructions);
}

void test_approximate_emission()
{
    ec::point_cache_scope<double> scope;
    auto const p = crossing();
    auto const expected = reference<double>();
    double const x = p.output_approximation(0), y = p.output_approximation(1);
    BOOST_TEST(x != bg::get<0>(expected) || y != bg::get<1>(expected));
    BOOST_TEST(ec::operation_detail::current<double>() == nullptr);
    point<double> output;
    point<float> single;
    exact_point internal;
    auto const before = ec::filtering_statistics();
    allocations = 0;
    count_allocations = true;
    ec::assign_output(p, output);
    BOOST_TEST_EQ(bg::get<0>(output), x);
    BOOST_TEST_EQ(bg::get<1>(output), y);
    default_strategy::convert_geometry(p, single);
    default_strategy::convert_geometry(p, internal);
    count_allocations = false;
    BOOST_TEST_EQ(allocations, 0u);
    BOOST_TEST_EQ(bg::get<0>(single), static_cast<float>(x));
    BOOST_TEST_EQ(bg::get<1>(single), static_cast<float>(y));
    BOOST_TEST(internal.shares_node(p));
    BOOST_TEST(ec::operation_detail::current<double>() == nullptr);
    same_statistics(before, ec::filtering_statistics());

    // Refinement for a predicate updates the stored approximation. Emission
    // must use that updated value, without initiating any further refinement.
    p.refined_approximation(0);
    p.refined_approximation(1);
    double const refined_x = p.output_approximation(0), refined_y = p.output_approximation(1);
    auto const refined_before = ec::filtering_statistics();
    auto const entries = ec::point_cache_detail::native_cache<double>::size();
    allocations = 0;
    count_allocations = true;
    ec::assign_output(p, output);
    count_allocations = false;
    BOOST_TEST_EQ(allocations, 0u);
    BOOST_TEST_EQ(bg::get<0>(output), refined_x);
    BOOST_TEST_EQ(bg::get<1>(output), refined_y);
    BOOST_TEST(refined_x != x || refined_y != y);
    BOOST_TEST_EQ(ec::point_cache_detail::native_cache<double>::size(), entries);
    same_statistics(refined_before, ec::filtering_statistics());
}

void test_nearest_and_exact_destinations()
{
    ec::point_cache_scope<double> scope;
    auto const p = crossing();
    point<double> output;
    point<float> single;
    nearest_strategy::convert_geometry(p, output);
    nearest_strategy::convert_geometry(p, single);
    auto const expected = reference<double>();
    auto const expected_single = reference<float>();
    BOOST_TEST_EQ(bg::get<0>(output), bg::get<0>(expected));
    BOOST_TEST_EQ(bg::get<1>(output), bg::get<1>(expected));
    BOOST_TEST_EQ(bg::get<0>(single), bg::get<0>(expected_single));
    BOOST_TEST_EQ(bg::get<1>(single), bg::get<1>(expected_single));

    exact_point internal;
    auto const before = ec::filtering_statistics();
    ec::assign_output<sr::nearest_output>(p, internal);
    BOOST_TEST(internal.shares_node(p));
    same_statistics(before, ec::filtering_statistics());
    BOOST_TEST(internal.exact_coordinate(0) == ec::exact_number(-74541) / 15325);
    BOOST_TEST(internal.exact_coordinate(1) == ec::exact_number(917427) / 15325);
    point<ec::exact_number> diagnostic;
    default_strategy::convert_geometry(p, diagnostic);
    BOOST_TEST(bg::get<0>(diagnostic) == internal.exact_coordinate(0));
    BOOST_TEST(bg::get<1>(diagnostic) == internal.exact_coordinate(1));

    ec::exact_number const ulp = ec::exact_number(1) / (ec::point_detail::output_integer(1) << 52);
    auto const tie = exact_point::literal(1 + ulp / 2, 1 + 3 * ulp / 2);
    ec::assign_output<sr::nearest_output>(tie, output);
    BOOST_TEST_EQ(bg::get<0>(output), 1.0);
    BOOST_TEST_EQ(bg::get<1>(output), 1.0 + 2 * std::ldexp(1.0, -52));
}

void test_integer_conversion()
{
    ec::point_cache_scope<double> scope;
    for (double x : {-2.5, -1.5, -0.5, -0.3, 0.3, 0.5, 1.5, 2.5})
    {
        auto const p = exact_point::literal(ec::exact_value(x), ec::exact_value(-x));
        point<int> ordinary, approximate, nearest;
        bg::convert(point<double>(x, -x), ordinary);
        default_strategy::convert_geometry(p, approximate);
        nearest_strategy::convert_geometry(p, nearest);
        BOOST_TEST_EQ(bg::get<0>(ordinary), bg::util::numeric_cast<int>(x));
        BOOST_TEST_EQ(bg::get<0>(approximate), bg::get<0>(ordinary));
        BOOST_TEST_EQ(bg::get<1>(approximate), bg::get<1>(ordinary));
        BOOST_TEST_EQ(bg::get<0>(nearest), bg::get<0>(ordinary));
        BOOST_TEST_EQ(bg::get<1>(nearest), bg::get<1>(ordinary));
    }
    using line = bg::model::linestring<point<double>>;
    for (double origin : {-2.0, -1.0, 0.0, 1.0, 2.0})
    {
        line const a{{origin, origin}, {origin + 1, origin + 1}};
        line const b{{origin, origin + 1}, {origin + 1, origin}};
        std::vector<point<int>> ordinary, exact;
        bg::intersection(a, b, ordinary);
        bg::intersection(a, b, exact, default_strategy());
        BOOST_TEST_EQ(ordinary.size(), 1u);
        BOOST_TEST_EQ(exact.size(), ordinary.size());
        BOOST_TEST_EQ(bg::get<0>(exact.front()), bg::get<0>(ordinary.front()));
        BOOST_TEST_EQ(bg::get<1>(exact.front()), bg::get<1>(ordinary.front()));
    }
    // Native integer construction is not a scalar conversion: it rounds a
    // displacement from a selected endpoint, even differently in x and y.
    using int_line = bg::model::linestring<point<int>>;
    std::vector<point<int>> integer_construction;
    bg::intersection(int_line{{0, 0}, {1, 1}}, int_line{{0, 1}, {1, 0}}, integer_construction);
    BOOST_TEST_EQ(integer_construction.size(), 1u);
    BOOST_TEST_EQ(bg::get<0>(integer_construction.front()), 1);
    BOOST_TEST_EQ(bg::get<1>(integer_construction.front()), 0);
    BOOST_TEST_EQ(bg::math::divide(1, 2), 1);
    BOOST_TEST_EQ(bg::math::divide(-1, 2), -1);

    point<int> unchanged(17, 23);
    auto const out_of_range = exact_point::literal(7, ec::exact_number(ec::point_detail::output_integer(1) << 100));
    BOOST_TEST_THROWS(default_strategy::convert_geometry(out_of_range, unchanged), boost::numeric::bad_numeric_cast);
    BOOST_TEST_EQ(bg::get<0>(unchanged), 17);
    BOOST_TEST_EQ(bg::get<1>(unchanged), 23);
    point<std::int64_t> large;
    ec::exact_point<std::int64_t> const native(9007199254740993LL, -7);
    ec::assign_output(native, large);
    BOOST_TEST_EQ(bg::get<0>(large), 9007199254740993LL);
}

template <typename Native>
struct native_range
{
    using point_type = Native;
    std::array<Native, 2> points;
    Native const& at(std::size_t i) const { return points.at(i); }
};

template <typename Native>
Native native_point(int x, int y)
{
    Native p;
    bg::set<0>(p, x);
    bg::set<1>(p, y);
    return p;
}

template <typename C, typename Native>
void test_direct_native_policy()
{
    ec::point_cache_scope<C> scope;
    native_range<Native> const a{{native_point<Native>(43, 7), native_point<Native>(-24, 81)}};
    native_range<Native> const b{{native_point<Native>(-78, 22), native_point<Native>(63, 95)}};
    using result = ec::segment_intersection_result<point<double>>;
    using policy = bg::policies::relate::segments_intersection_points<result>;
    using combined_policy = bg::policies::relate::segments_intersection_policy<result>;
    using approximate = sr::cartesian_exact<C>;
    using nearest = sr::cartesian_exact<C, sr::nearest_output>;
    auto const raw = approximate::relate().apply(a, b, policy());
    auto const rounded = nearest::relate().apply(a, b, policy());
    auto const combined = nearest::relate().apply(a, b, combined_policy());
    BOOST_TEST_EQ(raw.count, 1u);
    BOOST_TEST_EQ(rounded.count, 1u);
    BOOST_TEST_EQ(combined.intersection_points.count, 1u);
    auto const expected = reference<double>();
    BOOST_TEST_EQ(bg::get<0>(rounded.intersections[0]), bg::get<0>(expected));
    BOOST_TEST_EQ(bg::get<1>(rounded.intersections[0]), bg::get<1>(expected));
    BOOST_TEST_EQ(bg::get<0>(combined.intersection_points.intersections[0]), bg::get<0>(expected));
    BOOST_TEST(bg::get<0>(raw.intersections[0]) != bg::get<0>(expected)
        || bg::get<1>(raw.intersections[0]) != bg::get<1>(expected));
    point<double> converted;
    auto const before = ec::filtering_statistics();
    approximate::convert_geometry(a.at(0), converted);
    BOOST_TEST_EQ(bg::get<0>(converted), 43);
    BOOST_TEST_EQ(bg::get<1>(converted), 7);
    same_statistics(before, ec::filtering_statistics());

    native_range<Native> const common{{a.at(1), native_point<Native>(10, 70)}};
    auto const shared = nearest::relate().apply(a, common, policy());
    BOOST_TEST_EQ(shared.count, 1u);
    BOOST_TEST_EQ(bg::get<0>(shared.intersections[0]), -24);
    BOOST_TEST_EQ(bg::get<1>(shared.intersections[0]), 81);
}

void test_bg_integer_construction_inventory()
{
    using int_line = bg::model::linestring<point<int>>;
    using float_line = bg::model::linestring<point<double>>;
    for (int origin : {-2, -1, 0, 1, 2})
    for (int sign : {-1, 1})
    for (unsigned reverse = 0; reverse != 4; ++reverse)
    {
        int_line a{{origin, origin}, {origin + sign, origin + sign}};
        int_line b{{origin, origin + sign}, {origin + sign, origin}};
        if (reverse & 1) std::reverse(a.begin(), a.end());
        if (reverse & 2) std::reverse(b.begin(), b.end());
        float_line af, bf;
        bg::convert(a, af);
        bg::convert(b, bf);
        std::vector<point<int>> native_result, floating_result, approximate, nearest;
        bg::intersection(a, b, native_result);
        bg::intersection(af, bf, floating_result);
        bg::intersection(a, b, approximate, sr::cartesian_exact<int>());
        bg::intersection(a, b, nearest, sr::cartesian_exact<int, sr::nearest_output>());
        BOOST_TEST_EQ(native_result.size(), 1u);
        BOOST_TEST_EQ(floating_result.size(), 1u);
        BOOST_TEST_EQ(approximate.size(), 1u);
        BOOST_TEST_EQ(nearest.size(), 1u);
        // Native integer construction is a different BG operation. Both exact
        // strategies intentionally use ordinary conversion of their stored
        // approximation, not endpoint-relative displacement rounding.
        BOOST_TEST_EQ(bg::get<0>(native_result[0]), (reverse & 2) ? origin : origin + sign);
        BOOST_TEST_EQ(bg::get<1>(native_result[0]), (reverse & 2) ? origin + sign : origin);
        int const truncated = bg::util::numeric_cast<int>(origin + sign * 0.5);
        BOOST_TEST_EQ(bg::get<0>(floating_result[0]), truncated);
        BOOST_TEST_EQ(bg::get<1>(floating_result[0]), truncated);
        BOOST_TEST_EQ(bg::get<0>(approximate[0]), truncated);
        BOOST_TEST_EQ(bg::get<1>(approximate[0]), truncated);
        BOOST_TEST_EQ(bg::get<0>(nearest[0]), truncated);
        BOOST_TEST_EQ(bg::get<1>(nearest[0]), truncated);
    }
}

template <typename Strategy, bool SourceClockwise, bool SourceClosed, bool DestinationClosed, bool Reverse>
void test_ring_closure()
{
    ec::point_cache_scope<double> scope;
    using source_ring = bg::model::ring<exact_point, SourceClockwise, SourceClosed>;
    using output_ring = bg::model::ring<point<double>, SourceClockwise != Reverse, DestinationClosed>;
    source_ring source{crossing(), exact_point(0, 70), exact_point(10, 60)};
    if constexpr (!SourceClockwise) { std::swap(source[1], source[2]); }
    double const x = source.front().output_approximation(0);
    double const y = source.front().output_approximation(1);
    if constexpr (SourceClosed)
    {
        source.push_back(source.front());
        source.back().refined_approximation(0);
        source.back().refined_approximation(1);
        BOOST_TEST(source.back().output_approximation(0) != x
            || source.back().output_approximation(1) != y);
        BOOST_TEST_EQ(source.front().output_approximation(0), x);
        BOOST_TEST_EQ(source.front().output_approximation(1), y);
    }

    point<double> expected;
    Strategy::convert_geometry(source.front(), expected);
    output_ring output;
    output.reserve(4);
    auto const before = ec::filtering_statistics();
    auto const entries = ec::point_cache_detail::native_cache<double>::size();
    allocations = 0;
    count_allocations = true;
    Strategy::convert_geometry(source, output);
    count_allocations = false;
    BOOST_TEST_EQ(allocations, 0u);
    if constexpr (std::is_same<Strategy, default_strategy>::value)
    {
        same_statistics(before, ec::filtering_statistics());
        BOOST_TEST_EQ(ec::point_cache_detail::native_cache<double>::size(), entries);
        BOOST_TEST_EQ(bg::get<0>(output.front()), x);
        BOOST_TEST_EQ(bg::get<1>(output.front()), y);
    }
    BOOST_TEST_EQ(output.size(), DestinationClosed ? 4u : 3u);
    BOOST_TEST_EQ(bg::get<0>(output.front()), bg::get<0>(expected));
    BOOST_TEST_EQ(bg::get<1>(output.front()), bg::get<1>(expected));
    for (std::size_t i = 1; i < 3; ++i)
    {
        point<double> vertex;
        Strategy::convert_geometry(source[Reverse ? 3 - i : i], vertex);
        BOOST_TEST_EQ(bg::get<0>(output[i]), bg::get<0>(vertex));
        BOOST_TEST_EQ(bg::get<1>(output[i]), bg::get<1>(vertex));
    }
    if constexpr (DestinationClosed)
    {
        BOOST_TEST_EQ(bg::get<0>(output.back()), bg::get<0>(output.front()));
        BOOST_TEST_EQ(bg::get<1>(output.back()), bg::get<1>(output.front()));
    }

    using polygon = bg::model::polygon<point<double>, SourceClockwise != Reverse, DestinationClosed>;
    polygon polygon_output;
    polygon_output.inners().resize(1);
    Strategy::convert_geometry(source, polygon_output);
    BOOST_TEST(polygon_output.inners().empty());
    BOOST_TEST_EQ(polygon_output.outer().size(), output.size());
    for (std::size_t i = 0; i < output.size(); ++i)
    {
        BOOST_TEST_EQ(bg::get<0>(polygon_output.outer()[i]), bg::get<0>(output[i]));
        BOOST_TEST_EQ(bg::get<1>(polygon_output.outer()[i]), bg::get<1>(output[i]));
    }

    source.clear();
    Strategy::convert_geometry(source, output);
    Strategy::convert_geometry(source, polygon_output);
    BOOST_TEST(output.empty());
    BOOST_TEST(polygon_output.outer().empty());
    BOOST_TEST(polygon_output.inners().empty());
}

template <typename Strategy, bool SourceClockwise>
void test_ring_closure_traits()
{
    test_ring_closure<Strategy, SourceClockwise, false, false, false>();
    test_ring_closure<Strategy, SourceClockwise, false, false, true>();
    test_ring_closure<Strategy, SourceClockwise, false, true, false>();
    test_ring_closure<Strategy, SourceClockwise, false, true, true>();
    test_ring_closure<Strategy, SourceClockwise, true, false, false>();
    test_ring_closure<Strategy, SourceClockwise, true, false, true>();
    test_ring_closure<Strategy, SourceClockwise, true, true, false>();
    test_ring_closure<Strategy, SourceClockwise, true, true, true>();
}

template <typename Strategy, bool Reverse>
void test_exact_ring_closure_provenance()
{
    ec::point_cache_scope<double> scope;
    auto const first = ec::construct_intersection(exact_point(0, 0), exact_point(2, 2),
        exact_point(0, 2), exact_point(2, 0));
    auto const closing = ec::construct_intersection(exact_point(1, 0), exact_point(1, 2),
        exact_point(0, 1), exact_point(2, 1));
    BOOST_TEST(!first.shares_node(closing));
    std::array<double, 8> first_values, closing_values;
    BOOST_TEST(first.native_producers(first_values));
    BOOST_TEST(closing.native_producers(closing_values));
    BOOST_TEST(first_values != closing_values);
    using source_ring = bg::model::ring<exact_point>;
    using output_ring = bg::model::ring<exact_point, !Reverse>;
    source_ring const source{first, exact_point(0, 3), exact_point(1, 3), closing};
    output_ring output;
    output.reserve(4);
    auto const before = ec::filtering_statistics();
    Strategy::convert_geometry(source, output);
    same_statistics(before, ec::filtering_statistics());
    BOOST_TEST_EQ(output.size(), 4u);
    BOOST_TEST(output.front().shares_node(Reverse ? closing : first));
    BOOST_TEST(output.back().shares_node(Reverse ? first : closing));
    BOOST_TEST(!output.front().shares_node(output.back()));
    auto const& retained_closing = Reverse ? output.front() : output.back();
    std::array<double, 8> retained_values;
    BOOST_TEST(retained_closing.native_producers(retained_values));
    BOOST_TEST(retained_values == closing_values);
    BOOST_TEST(ec::equal_points(output.front(), output.back()));
    auto const next = ec::construct_intersection(retained_closing, exact_point(1, 3),
        exact_point(0, 2), exact_point(2, 2));
    BOOST_TEST(next.native_producers(retained_values));
    BOOST_TEST(ec::equal_points(next, exact_point(1, 2)));
}

template <typename Strategy>
void test_geometry_propagation()
{
    ec::point_cache_scope<double> scope;
    auto const p = crossing();
    point<double> expected;
    Strategy::convert_geometry(p, expected);
    using working_polygon = bg::model::polygon<exact_point>;
    using polygon = bg::model::polygon<point<double>>;
    working_polygon source;
    source.outer() = {p, exact_point(0, 70), exact_point(10, 60), p};
    source.inners().resize(1);
    source.inners()[0] = {exact_point(0, 61), exact_point(1, 61), exact_point(0, 62), exact_point(0, 61)};
    polygon converted;
    Strategy::convert_geometry(source, converted);
    BOOST_TEST_EQ(converted.outer().size(), source.outer().size());
    BOOST_TEST_EQ(converted.inners().size(), 1u);
    BOOST_TEST_EQ(bg::get<0>(converted.outer().front()), bg::get<0>(expected));
    working_polygon internal;
    Strategy::convert_geometry(source, internal);
    BOOST_TEST(internal.outer().front().shares_node(p));
    using working_points = bg::model::multi_point<exact_point>;
    using points = bg::model::multi_point<point<double>>;
    std::tuple<working_points, std::vector<working_polygon>> const tuple_source{{p}, {source}};
    std::tuple<points, std::vector<polygon>> tuple_destination;
    Strategy::append_output(tuple_source, tuple_destination);
    BOOST_TEST_EQ(bg::get<0>(std::get<0>(tuple_destination).front()), bg::get<0>(expected));
    BOOST_TEST_EQ(bg::get<1>(std::get<1>(tuple_destination).front().outer().front()), bg::get<1>(expected));
    using working_variant = boost::variant<exact_point, working_polygon>;
    using output_variant = boost::variant<point<double>, polygon>;
    bg::model::geometry_collection<working_variant> collection_source;
    collection_source.push_back(p);
    collection_source.push_back(source);
    bg::model::geometry_collection<output_variant> collection_output;
    Strategy::convert_geometry(collection_source, collection_output);
    BOOST_TEST_EQ(collection_output.size(), 2u);
    BOOST_TEST_EQ(bg::get<0>(boost::get<point<double>>(collection_output[0])), bg::get<0>(expected));
    BOOST_TEST_EQ(bg::get<1>(boost::get<polygon>(collection_output[1]).outer().front()), bg::get<1>(expected));
    bg::model::box<exact_point> const box(p, exact_point(50, 100));
    bg::model::box<point<double>> box_output;
    Strategy::convert_geometry(box, box_output);
    BOOST_TEST_EQ((bg::get<bg::min_corner, 0>(box_output)), bg::get<0>(expected));
    BOOST_TEST_EQ((bg::get<bg::min_corner, 1>(box_output)), bg::get<1>(expected));
    bg::model::box<exact_point> const containing(exact_point(-100, -100), exact_point(200, 200));
    BOOST_TEST(Strategy::box_intersection().apply(box, containing, box_output));
    BOOST_TEST_EQ((bg::get<bg::min_corner, 0>(box_output)), bg::get<0>(expected));
    BOOST_TEST_EQ((bg::get<bg::min_corner, 1>(box_output)), bg::get<1>(expected));
}

void test_public_set_operations()
{
    using line = bg::model::linestring<point<double>>;
    line const a{{43, 7}, {-24, 81}}, b{{-78, 22}, {63, 95}};
    std::vector<point<double>> approximate, nearest;
    bg::intersection(a, b, approximate, default_strategy());
    bg::intersection(a, b, nearest, nearest_strategy());
    BOOST_TEST_EQ(approximate.size(), 1u);
    BOOST_TEST_EQ(nearest.size(), 1u);
    ec::point_cache_scope<double> scope;
    auto const p = crossing();
    BOOST_TEST_EQ(bg::get<0>(approximate.front()), p.output_approximation(0));
    BOOST_TEST_EQ(bg::get<1>(approximate.front()), p.output_approximation(1));
    auto const expected = reference<double>();
    BOOST_TEST_EQ(bg::get<0>(nearest.front()), bg::get<0>(expected));
    BOOST_TEST_EQ(bg::get<1>(nearest.front()), bg::get<1>(expected));
    BOOST_TEST_EQ(bg::relation(a, b, default_strategy()).str(), bg::relation(a, b, nearest_strategy()).str());
}

void test_polygon_operation_emission()
{
    ec::point_cache_scope<double> scope;
    using polygon = bg::model::polygon<point<double>>;
    polygon a, b;
    a.outer() = {{43, 7}, {-24, 81}, {100, 100}, {43, 7}};
    b.outer() = {{-78, 22}, {63, 95}, {-100, -100}, {-78, 22}};
    auto bits = [](double value) { return boost::core::bit_cast<std::uint64_t>(value); };
    auto current = [](auto const& p, std::size_t dimension)
    {
        return p.is_native() ? p.native_coordinate(dimension)
                             : p.record().numeric.approximate(dimension);
    };
    auto cache_snapshot = [&]
    {
        auto const* owner = ec::operation_detail::current<double>();
        void const* native_data = nullptr;
        void const* auxiliary_data = nullptr;
        std::vector<std::uint64_t> state;
        auto integer = [&](auto const& value)
        {
            auto const& backend = value.backend();
            state.push_back(backend.size());
            state.push_back(backend.sign());
            for (std::size_t i = 0; i < backend.size(); ++i) state.push_back(backend.limbs()[i]);
        };
        auto homogeneous = [&](auto const& p) { integer(p.X); integer(p.Y); integer(p.W); };
        auto expansion = [&](auto const& value)
        {
            state.push_back(value.valid);
            state.push_back(value.size);
            for (std::size_t i = 0; i < value.size; ++i) state.push_back(bits(value.components[i]));
        };
        if (owner)
        {
            auto const& cache = owner->storage();
            native_data = cache.native.data();
            auxiliary_data = cache.auxiliary.data();
            state = {cache.native.size(), cache.native.capacity(), cache.auxiliary.size(), cache.auxiliary.capacity()};
            for (auto const& entry : cache.native)
            {
                state.push_back(entry.exact_ready);
                homogeneous(entry.value);
                for (auto const& coordinate : entry.expanded)
                {
                    state.push_back(coordinate.attempted);
                    state.push_back(coordinate.valid);
                    state.push_back(coordinate.bounded);
                    expansion(coordinate.numerator);
                    expansion(coordinate.denominator);
                    state.push_back(bits(coordinate.bounds.approximate));
                    state.push_back(bits(coordinate.bounds.error));
                }
            }
            for (auto const& p : cache.auxiliary) homogeneous(p);
        }
        return std::make_tuple(owner, native_data, auxiliary_data, std::move(state));
    };
    auto check = [&](auto operation)
    {
        auto check_policy = [&](auto const& strategy)
        {
            strategy.with_sources(a, b, [&](auto const& bound)
            {
                using bound_type = std::decay_t<decltype(bound)>;
                using working_polygon = typename bound_type::template overlay_geometry_type<polygon>;
                using output_policy = typename bound_type::output_policy;
                using traits = ec::construction_rounding::policy_traits<output_policy>;
                std::vector<working_polygon> working;
                operation(a, b, working, bound);
                BOOST_TEST(!working.empty());
                std::vector<polygon> expected(working.size()), output;
                std::size_t constructed = 0;
                auto capture_ring = [&](auto const& source, auto& destination)
                {
                    for (auto const& p : source)
                    {
                        if (!p.is_native())
                        {
                            ++constructed;
                            // Deferred emission deliberately requests nearest.
                            // Warm it explicitly BEFORE the export snapshot;
                            // ordinary/eager need no such output refinement.
                            if constexpr (traits::nearest_output && !traits::eager)
                            { p.refine(0); p.refine(1); }
                        }
                        destination.emplace_back(current(p, 0), current(p, 1));
                    }
                };
                for (std::size_t i = 0; i < working.size(); ++i)
                {
                    capture_ring(working[i].outer(), expected[i].outer());
                    expected[i].inners().resize(working[i].inners().size());
                    for (std::size_t j = 0; j < working[i].inners().size(); ++j)
                        capture_ring(working[i].inners()[j], expected[i].inners()[j]);
                }
                BOOST_TEST(constructed > 0);
                // Never materialize legacy exact_point auxiliaries: doing so
                // would replace the stored approximation we intend to test.
                auto const cache_before = cache_snapshot();
                auto const statistics_before = ec::filtering_statistics();
                bound.convert_geometry(working, output);
                same_statistics(statistics_before, ec::filtering_statistics());
                BOOST_TEST(cache_snapshot() == cache_before);
                BOOST_TEST_EQ(output.size(), expected.size());
                if (output.size() != expected.size()) return;
                auto compare_ring = [&](auto const& source, auto const& actual, auto const& reference)
                {
                    BOOST_TEST_EQ(actual.size(), reference.size());
                    if (actual.size() != reference.size()) return;
                    for (std::size_t i = 0; i < actual.size(); ++i)
                    {
                        BOOST_TEST_EQ(bits(bg::get<0>(actual[i])), bits(bg::get<0>(reference[i])));
                        BOOST_TEST_EQ(bits(bg::get<1>(actual[i])), bits(bg::get<1>(reference[i])));
                        BOOST_TEST_EQ(bits(current(source[i], 0)), bits(bg::get<0>(reference[i])));
                        BOOST_TEST_EQ(bits(current(source[i], 1)), bits(bg::get<1>(reference[i])));
                    }
                };
                for (std::size_t i = 0; i < working.size(); ++i)
                {
                    compare_ring(working[i].outer(), output[i].outer(), expected[i].outer());
                    BOOST_TEST_EQ(output[i].inners().size(), expected[i].inners().size());
                    if (output[i].inners().size() != expected[i].inners().size()) continue;
                    for (std::size_t j = 0; j < working[i].inners().size(); ++j)
                        compare_ring(working[i].inners()[j], output[i].inners()[j], expected[i].inners()[j]);
                }
            });
        };
        check_policy(default_strategy());
        check_policy(sr::cartesian_exact<double, sr::eager_nearest_output>());
        check_policy(sr::cartesian_exact<double, sr::deferred_nearest_output>());
        check_policy(nearest_strategy());
    };
    check([](auto const& a, auto const& b, auto& out, auto const& strategy) { bg::intersection(a, b, out, strategy); });
    check([](auto const& a, auto const& b, auto& out, auto const& strategy) { bg::union_(a, b, out, strategy); });
    check([](auto const& a, auto const& b, auto& out, auto const& strategy) { bg::difference(a, b, out, strategy); });
    check([](auto const& a, auto const& b, auto& out, auto const& strategy) { bg::sym_difference(a, b, out, strategy); });
}

int main()
{
    BOOST_TEST_EQ(std::fegetround(), FE_TONEAREST);
    test_approximate_emission();
    test_nearest_and_exact_destinations();
    test_integer_conversion();
    test_bg_integer_construction_inventory();
    test_direct_native_policy<double, point<double>>();
    test_direct_native_policy<double, std::array<double, 2>>();
    test_direct_native_policy<double, point<float>>();
    test_direct_native_policy<double, std::array<float, 2>>();
    test_direct_native_policy<int, point<int>>();
    test_direct_native_policy<int, std::array<int, 2>>();
    test_direct_native_policy<int, std::array<std::int16_t, 2>>();
    test_ring_closure_traits<default_strategy, false>();
    test_ring_closure_traits<default_strategy, true>();
    test_ring_closure_traits<nearest_strategy, false>();
    test_ring_closure_traits<nearest_strategy, true>();
    test_exact_ring_closure_provenance<default_strategy, false>();
    test_exact_ring_closure_provenance<default_strategy, true>();
    test_exact_ring_closure_provenance<nearest_strategy, false>();
    test_exact_ring_closure_provenance<nearest_strategy, true>();
    test_geometry_propagation<default_strategy>();
    test_geometry_propagation<nearest_strategy>();
    test_public_set_operations();
    test_polygon_operation_emission();
    return boost::report_errors();
}
