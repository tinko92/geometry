// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.

#include <boost/core/lightweight_test.hpp>
#include <boost/geometry/algorithms/area.hpp>
#include <boost/geometry/algorithms/convert.hpp>
#include <boost/geometry/algorithms/correct.hpp>
#include <boost/geometry/algorithms/difference.hpp>
#include <boost/geometry/algorithms/intersection.hpp>
#include <boost/geometry/algorithms/is_valid.hpp>
#include <boost/geometry/algorithms/num_interior_rings.hpp>
#include <boost/geometry/algorithms/prepare_filtered_decimal_polygon.hpp>
#include <boost/geometry/algorithms/prepare_rational_polygon.hpp>
#include <boost/geometry/algorithms/relation.hpp>
#include <boost/geometry/io/wkt/read.hpp>
#include <boost/geometry/strategies/relate/cartesian_decimal.hpp>
#include <boost/geometry/strategies/relate/cartesian_exact.hpp>

#include "exact_additional_issue_cases.hpp"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <string>
#include <utility>

namespace bg = boost::geometry;
using rational = boost::multiprecision::cpp_rational;
using decimal = bg::exact_coordinate<std::string>;

template <typename C>
using polygon = bg::model::polygon<bg::model::point<C, 2, bg::cs::cartesian>, false>;
template <typename C>
using polygons = bg::model::multi_polygon<polygon<C>>;
using output = polygons<rational>;

bg::strategies::relate::cartesian<rational> const plain;
bg::strategies::relate::cartesian_decimal const filtered;
bg::strategies::relate::cartesian_exact<double> const binary;

template <typename G>
G read(char const* wkt)
{
    G result;
    bg::read_wkt(wkt, result);
    return result;
}

template <typename A, typename B, typename Strategy>
output common(A const& a, B const& b, Strategy const& strategy)
{
    output result;
    BOOST_TEST(bg::intersection(a, b, result, strategy));
    BOOST_TEST_EQ(result.size(), 1u);
    BOOST_TEST(bg::is_valid(result, plain));
    return result;
}

template <typename G, typename Strategy>
G normalized(G geometry, Strategy const& strategy)
{
    // Match the historical replay's winding/closure normalization, without
    // rounding the exact decimal ordinates or rewriting the source tokens.
    bg::correct(geometry, strategy);
    return geometry;
}

struct intersection_observation
{
    rational decimal_area;
    rational binary_area;
};

intersection_observation replay_intersection(exact_issue_cases::wkt_pair const& pair)
{
    auto const sa = read<polygon<decimal>>(pair.first);
    auto const sb = read<polygon<decimal>>(pair.second);
    auto const a = normalized(bg::prepare_rational_polygon(sa), plain);
    auto const b = normalized(bg::prepare_rational_polygon(sb), plain);
    auto const fa = normalized(bg::prepare_filtered_decimal_polygon(sa), plain);
    auto const fb = normalized(bg::prepare_filtered_decimal_polygon(sb), plain);
    auto const ba = normalized(read<polygon<double>>(pair.first), binary);
    auto const bb = normalized(read<polygon<double>>(pair.second), binary);

    auto const reference = common(a, b, plain);
    auto const decimal_result = common(fa, fb, filtered);
    auto const binary_result = common(ba, bb, binary);
    BOOST_TEST_EQ(bg::relation(a, b, plain).str(), "212101212");
    BOOST_TEST_EQ(bg::relation(fa, fb, filtered).str(), "212101212");
    BOOST_TEST_EQ(bg::relation(ba, bb, binary).str(), "212101212");
    BOOST_TEST_EQ(bg::relation(decimal_result, reference, plain).str(), "2FFF1FFF2");

    return {bg::area(decimal_result, plain), bg::area(binary_result, plain)};
}

void report(unsigned issue, char const* matrix, intersection_observation const& observation)
{
    std::cout << '#' << issue << " matrix (both)=" << matrix << " overlap_area: report decimals="
              << std::setprecision(17) << observation.decimal_area.convert_to<double>()
              << ", double-preparsed=" << observation.binary_area.convert_to<double>()
              << ", decimal minus binary="
              << rational(observation.decimal_area - observation.binary_area).convert_to<double>() << '\n';
}

void issue_1201()
{
    auto const& pairs = exact_issue_cases::cases();
    auto const found = std::find_if(pairs.begin(), pairs.end(), [](auto const& p) { return p.issue == 1201; });
    BOOST_TEST(found != pairs.end());
    if (found == pairs.end())
        return;

    auto const observation = replay_intersection(*found);
    // For the corrected comment pair, x=0 clips a triangle of width e.
    // Its vertical side is the difference of the two edge intercepts.
    auto const expected = [](rational const& e)
    {
        return rational(e * (50 * e / (10 + e) - 50 * e / (20 + e)) / 2);
    };
    BOOST_TEST_EQ(observation.decimal_area, expected(rational(1, 100000000000000LL)));
    BOOST_TEST_EQ(observation.binary_area, expected(rational(1e-14)));
    BOOST_TEST(observation.decimal_area != observation.binary_area);
    report(1201, "212101212", observation);
}

void issue_704()
{
    auto const& pairs = exact_additional_issue_cases::pairs();
    auto const found = std::find_if(pairs.begin(), pairs.end(), [](auto const& p) { return p.data.issue == 704; });
    BOOST_TEST(found != pairs.end());
    if (found == pairs.end())
        return;

    auto const observation = replay_intersection(found->data);
    BOOST_TEST(observation.decimal_area > 0);
    BOOST_TEST(observation.binary_area > 0);
    BOOST_TEST(observation.decimal_area != observation.binary_area);
    report(704, "212101212", observation);
}

void issue_669()
{
    namespace data = exact_additional_issue_cases::data;
    auto const sa = read<polygons<decimal>>(data::issue_669_0);
    auto const sb = read<polygons<decimal>>(data::issue_669_1);
    auto const a = normalized(bg::prepare_rational_polygon(sa), plain);
    auto const b = normalized(bg::prepare_rational_polygon(sb), plain);
    auto const fa = normalized(bg::prepare_filtered_decimal_polygon(sa), plain);
    auto const fb = normalized(bg::prepare_filtered_decimal_polygon(sb), plain);
    auto const ba = normalized(read<polygons<double>>(data::issue_669_0), binary);
    auto const bb = normalized(read<polygons<double>>(data::issue_669_1), binary);

    auto const check = [](auto const& lhs, auto const& rhs, auto const& strategy, rational const& input_area)
    {
        output result, overlap;
        bg::difference(lhs, rhs, result, strategy);
        bg::intersection(lhs, rhs, overlap, strategy);
        BOOST_TEST(bg::is_valid(result, plain));
        BOOST_TEST(bg::is_valid(overlap, plain));
        BOOST_TEST_EQ(result.size(), 1u);
        BOOST_TEST_EQ(bg::num_interior_rings(result), 1u);
        BOOST_TEST_EQ(bg::area(result, plain) + bg::area(overlap, plain), input_area);
        return std::make_pair(std::move(result), std::move(overlap));
    };
    rational const input_area = bg::area(a, plain);
    auto const reference = check(a, b, plain, input_area);
    auto const decimal_result = check(fa, fb, filtered, input_area);
    output binary_input;
    // cpp_rational's floating constructor preserves the exact binary value.
    bg::convert(ba, binary_input);
    auto const binary_result = check(ba, bb, binary, bg::area(binary_input, plain));

    BOOST_TEST_EQ(bg::relation(decimal_result.first, reference.first, plain).str(), "2FFF1FFF2");
    BOOST_TEST_EQ(bg::relation(decimal_result.second, reference.second, plain).str(), "2FFF1FFF2");
    BOOST_TEST_EQ(bg::relation(a, b, plain).str(), "212101212");
    BOOST_TEST_EQ(bg::relation(fa, fb, filtered).str(), "212101212");
    BOOST_TEST_EQ(bg::relation(ba, bb, binary).str(), "212101212");
    rational const decimal_area = bg::area(decimal_result.second, plain);
    rational const binary_area = bg::area(binary_result.second, plain);
    BOOST_TEST(decimal_area > 0);
    BOOST_TEST(binary_area > 0);
    BOOST_TEST(decimal_area != binary_area);
    report(669, "212101212", {decimal_area, binary_area});
}

int main()
{
    issue_1201();
    issue_704();
    issue_669();
    return boost::report_errors();
}
