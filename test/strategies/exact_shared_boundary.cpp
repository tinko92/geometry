// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.

#define BOOST_GEOMETRY_EXACT_FILTER_STATISTICS
#include <boost/core/lightweight_test.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/strategies/relate/cartesian_exact.hpp>

#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <string>
#include <type_traits>
#include <vector>

namespace bg = boost::geometry;
namespace ec = bg::strategies::detail::exact_cartesian;
using rational_point = bg::model::d2::point_xy<ec::exact_number>;
using rational_polygon = bg::model::polygon<rational_point>;
using rational_multi = bg::model::multi_polygon<rational_polygon>;

enum class operation { unite, intersect, subtract, symmetric };
constexpr std::array<operation, 4> operations{{operation::unite, operation::intersect,
                                             operation::subtract, operation::symmetric}};

struct rectangle { int left, bottom, right, top; };
struct fixture
{
    char const* name;
    rectangle first, second;
    std::array<std::vector<rectangle>, 4> expected;
};

std::array<fixture, 4> fixtures()
{
    return {{
        {"overlap", {0, 0, 32, 16}, {8, 0, 40, 16}, {{
            {{0, 0, 40, 16}}, {{8, 0, 32, 16}}, {{0, 0, 8, 16}},
            {{0, 0, 8, 16}, {32, 0, 40, 16}}}}},
        {"adjacent", {0, 0, 32, 16}, {32, 0, 64, 16}, {{
            {{0, 0, 64, 16}}, {}, {{0, 0, 32, 16}}, {{0, 0, 64, 16}}}}},
        {"identical", {0, 0, 32, 16}, {0, 0, 32, 16}, {{
            {{0, 0, 32, 16}}, {{0, 0, 32, 16}}, {}, {}}}},
        {"contained", {0, 0, 32, 16}, {0, 0, 24, 16}, {{
            {{0, 0, 32, 16}}, {{0, 0, 24, 16}}, {{24, 0, 32, 16}}, {{24, 0, 32, 16}}}}}
    }};
}

template <typename Polygon>
Polygon make_rectangle(rectangle const& box, int stride, int transform)
{
    using coordinate = bg::coordinate_type_t<Polygon>;
    Polygon result;
    auto& ring = result.outer();
    std::array<std::array<int, 2>, 4> const corners{{
        {{box.left, box.bottom}}, {{box.left, box.top}},
        {{box.right, box.top}}, {{box.right, box.bottom}}}};
    for (std::size_t side = 0; side != corners.size(); ++side)
    {
        auto const& first = corners[side];
        auto const& second = corners[(side + 1) % corners.size()];
        int const dx = (second[0] > first[0]) - (second[0] < first[0]);
        int const dy = (second[1] > first[1]) - (second[1] < first[1]);
        int const length = dx ? (second[0] - first[0]) * dx : (second[1] - first[1]) * dy;
        for (int offset = 0; offset < length; offset += stride)
        {
            int x = first[0] + offset * dx, y = first[1] + offset * dy;
            if (transform == 1) { int const old_x = x; x = -y; y = old_x; }
            if (transform == 2) { x += y; }
            ring.emplace_back(coordinate(x), coordinate(y));
        }
    }
    if constexpr (bg::point_order<Polygon>::value == bg::counterclockwise)
        std::reverse(ring.begin(), ring.end());
    if constexpr (bg::closure<Polygon>::value == bg::closed)
        ring.push_back(ring.front());
    return result;
}

template <typename A, typename B, typename Result, typename Strategy>
void overlay(operation op, A const& a, B const& b, Result& result, Strategy const& strategy)
{
    switch (op)
    {
    case operation::unite: bg::union_(a, b, result, strategy); break;
    case operation::intersect: bg::intersection(a, b, result, strategy); break;
    case operation::subtract: bg::difference(a, b, result, strategy); break;
    case operation::symmetric: bg::sym_difference(a, b, result, strategy); break;
    }
}

void check_set(rational_multi const& actual, rational_multi const& expected, std::string const& label)
{
    bg::strategies::relate::cartesian<ec::exact_number> const strategy;
    bool const equal = actual.empty() || expected.empty()
        ? actual.empty() && expected.empty()
        : bg::relate(actual, expected, bg::de9im::mask("T*F**FFF*"), strategy);
    if (!equal)
        std::cerr << label << " actual=" << bg::wkt(actual) << " expected=" << bg::wkt(expected) << '\n';
    BOOST_TEST(equal);
    std::string reason;
    bool const valid = bg::is_valid(actual, reason, strategy);
    if (!valid) std::cerr << label << " invalid exact output: " << reason << '\n';
    BOOST_TEST(valid);
}

template <typename C, typename Policy>
void check(char const* policy_name)
{
    using native_point = bg::model::d2::point_xy<C>;
    using closed_polygon = bg::model::polygon<native_point, true, true>;
    using open_reversed_polygon = bg::model::polygon<native_point, false, false>;
    using native_multi = bg::model::multi_polygon<closed_polygon>;
    bg::strategies::relate::cartesian_exact<C, Policy> const strategy;
    bg::strategies::relate::cartesian<ec::exact_number> const exact_control;
    std::string const prefix = std::string(std::is_integral<C>::value ? "int64/" : "double/") + policy_name;
    std::uint64_t small = 0, compact = 0, full = 0;
    for (auto const& item : fixtures())
    for (int transform = 0; transform != 3; ++transform)
    {
        auto const a = make_rectangle<closed_polygon>(item.first, 1, transform);
        auto const b = make_rectangle<open_reversed_polygon>(item.second, 3, transform);
        BOOST_TEST(bg::is_valid(a));
        BOOST_TEST(bg::is_valid(b));
        rational_polygon qa, qb;
        bg::convert(a, qa);
        bg::convert(b, qb);
        for (std::size_t index = 0; index != operations.size(); ++index)
        {
            std::string const label = prefix + "/" + item.name + "/" + std::to_string(transform)
                + "/" + std::to_string(index);
            rational_multi expected, control;
            for (auto const& box : item.expected[index])
                expected.push_back(make_rectangle<rational_polygon>(box, 128, transform));
            overlay(operations[index], qa, qb, control, exact_control);
            check_set(control, expected, label + "/control");

            auto collect_statistics = [&]
            {
                auto const counts = ec::filtering_statistics();
                if (counts.full_measure_comparisons)
                    std::cerr << label << " full measure comparisons=" << counts.full_measure_comparisons << '\n';
                BOOST_TEST_EQ(counts.full_measure_comparisons, 0u);
                small += counts.small_integer_measure_comparisons;
                compact += counts.compact_measure_comparisons;
                full += counts.full_measure_comparisons;
            };

            rational_multi exact;
            ec::filtering_statistics() = {};
            overlay(operations[index], a, b, exact, strategy);
            collect_statistics();
            check_set(exact, expected, label + "/exact");
            check_set(exact, control, label + "/control-match");

            // Every analytic result lies on the small integer lattice, so all
            // output policies must also export the same native set.
            native_multi native;
            ec::filtering_statistics() = {};
            overlay(operations[index], a, b, native, strategy);
            collect_statistics();
            std::string reason;
            bool const valid = bg::is_valid(native, reason);
            if (!valid) std::cerr << label << " invalid native output: " << reason << '\n';
            BOOST_TEST(valid);
            rational_multi lifted;
            bg::convert(native, lifted);
            check_set(lifted, expected, label + "/native");
        }
    }
    BOOST_TEST(small > 0);
    BOOST_TEST_EQ(full, 0u);
    std::cout << prefix << " small_integer_measure_comparisons=" << small
              << " compact_measure_comparisons=" << compact << " full_measure_comparisons=" << full << '\n';
}

int main()
{
    check<double, bg::strategies::relate::approximate_output>("approximate");
    check<double, bg::strategies::relate::eager_nearest_output>("eager");
    check<double, bg::strategies::relate::deferred_nearest_output>("deferred");
    check<double, bg::strategies::relate::lean_bounded_ulp_output<100>>("lean100");
    check<std::int64_t, bg::strategies::relate::approximate_output>("approximate");
    check<std::int64_t, bg::strategies::relate::eager_nearest_output>("eager");
    check<std::int64_t, bg::strategies::relate::deferred_nearest_output>("deferred");
    return boost::report_errors();
}
