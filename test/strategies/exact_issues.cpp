// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Licensed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt

#include <boost/core/lightweight_test.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/strategies/relate/cartesian.hpp>
#include <boost/geometry/strategies/relate/cartesian_exact.hpp>

#include "exact_issue_cases.hpp"

#include <cstddef>
#include <cstdint>
#include <string>

namespace bg = boost::geometry;
using exact_number = bg::strategies::detail::exact_cartesian::exact_number;

template <typename Geometry>
void test_exact_validity(Geometry const& result, exact_issue_cases::wkt_pair const& item)
{
    bg::strategies::detail::exact_cartesian::point_cache_scope operation;
    // is_valid does not yet support cartesian_exact's turn-ratio types.
    bg::strategies::relate::cartesian<exact_number> validation;
    std::string reason;
    if (!bg::is_valid(result, reason, validation))
    {
        BOOST_ERROR((std::string(item.name) + ": invalid exact result: " + reason).c_str());
    }
}

void test_predicate(bool result, exact_issue_cases::wkt_pair const& item)
{
    if (item.expected_predicate >= 0)
    {
        BOOST_TEST_EQ(result, item.expected_predicate != 0);
    }
}

template <typename Coordinate>
void test_issue(exact_issue_cases::wkt_pair const& item)
{
    bg::strategies::detail::exact_cartesian::point_cache_scope<Coordinate> operation_scope;
    using point = bg::model::d2::point_xy<Coordinate>;
    using polygon = bg::model::polygon<point>;
    using multi = bg::model::multi_polygon<polygon>;
    using implicit_point = bg::strategies::detail::exact_cartesian::exact_point<Coordinate>;
    using exact_polygon = bg::model::polygon<implicit_point>;
    using exact_multi = bg::model::multi_polygon<exact_polygon>;

    bg::strategies::relate::cartesian_exact<Coordinate> strategy;
    bg::strategies::relate::cartesian<exact_number> validation;
    std::string const operation(item.operation);
    if (operation == "within")
    {
        using line = bg::model::linestring<point>;
        bg::model::multi_linestring<line> a, b;
        bg::read_wkt(item.first, a);
        bg::read_wkt(item.second, b);
        test_predicate(bg::within(a, b, strategy), item);
        return;
    }

    auto read = [&](char const* wkt)
    {
        multi result;
        if (std::string(wkt).find("MULTIPOLYGON") == 0)
        {
            bg::read_wkt(wkt, result);
        }
        else
        {
            result.emplace_back();
            bg::read_wkt(wkt, result.back());
        }
        if (item.correct_inputs)
        {
            bg::correct(result, strategy);
        }
        return result;
    };
    auto const a = read(item.first);
    auto const b = read(item.second);
    if (operation == "covered_by")
    {
        test_predicate(bg::covered_by(a, b, strategy), item);
        return;
    }

    exact_multi result;
    if (operation == "intersection")
    {
        bg::intersection(a, b, result, strategy);
    }
    else if (operation == "union")
    {
        bg::union_(a, b, result, strategy);
    }
    else if (operation == "difference")
    {
        bg::difference(a, b, result, strategy);
    }
    else if (operation == "sym_difference")
    {
        bg::sym_difference(a, b, result, strategy);
    }
    else
    {
        BOOST_ERROR((std::string(item.name) + ": unknown operation: " + operation).c_str());
        return;
    }

    if (item.expected_components >= 0)
    {
        BOOST_TEST_EQ(result.size(), static_cast<std::size_t>(item.expected_components));
    }
    test_exact_validity(result, item);

    if (operation == "union")
    {
        exact_multi reverse;
        bg::union_(b, a, reverse, strategy);
        test_exact_validity(reverse, item);
        BOOST_TEST_EQ(reverse.size(), result.size());
        exact_number const area = bg::area(result, validation);
        BOOST_TEST_EQ(bg::area(reverse, validation), area);
        if (item.issue == 1506)
        {
            exact_number const input_area = bg::area(a, validation) + bg::area(b, validation);
            BOOST_TEST_EQ(area, input_area);
        }
    }
}

int main()
{
    for (auto const& item : exact_issue_cases::cases())
    {
        std::string const coordinate(item.coordinate_type);
        if (coordinate == "double")
        {
            test_issue<double>(item);
        }
        else if (coordinate == "float")
        {
            test_issue<float>(item);
        }
        else if (coordinate == "int32_t")
        {
            test_issue<std::int32_t>(item);
        }
        else if (coordinate == "int64_t")
        {
            test_issue<std::int64_t>(item);
        }
        else
        {
            BOOST_ERROR((std::string(item.name) + ": unknown coordinate type: " + coordinate).c_str());
        }
    }
    return boost::report_errors();
}
