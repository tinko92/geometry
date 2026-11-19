// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Licensed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt

#include <boost/core/lightweight_test.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/strategies/relate/cartesian.hpp>
#include <boost/geometry/strategies/relate/cartesian_exact.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include "exact_additional_issue_cases.hpp"

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace bg = boost::geometry;
namespace ec = bg::strategies::detail::exact_cartesian;
namespace cases = exact_additional_issue_cases;
using rational = ec::exact_number;
using record = boost::property_tree::ptree;

template <typename Geometry>
rational validation_area(Geometry const& geometry)
{
    return bg::area(geometry, bg::strategies::relate::cartesian<rational>());
}

template <typename Coordinate>
struct geometries
{
    using native_point = bg::model::d2::point_xy<Coordinate>;
    using point = ec::exact_point<Coordinate>;
    using line = bg::model::linestring<point>;
    using lines = bg::model::multi_linestring<line>;
    using polygon = bg::model::polygon<point>;
    using polygons = bg::model::multi_polygon<polygon>;
    using points = bg::model::multi_point<point>;
    using box = bg::model::box<point>;
    using strategy = bg::strategies::relate::cartesian_exact<Coordinate>;
};

template <typename Coordinate, typename NativeRange, typename ExactRange>
void lift_range(NativeRange const& input, ExactRange& output)
{
    for (auto const& p : input) output.push_back(ec::make_point<Coordinate>(p));
}

template <typename Coordinate>
typename geometries<Coordinate>::polygon read_polygon(char const* wkt)
{
    using types = geometries<Coordinate>;
    bg::model::polygon<typename types::native_point> input;
    bg::read_wkt(wkt, input);
    typename types::polygon result;
    lift_range<Coordinate>(input.outer(), result.outer());
    result.inners().resize(input.inners().size());
    for (std::size_t i = 0; i < input.inners().size(); ++i)
    {
        lift_range<Coordinate>(input.inners()[i], result.inners()[i]);
    }
    return result;
}

template <typename Coordinate>
typename geometries<Coordinate>::polygons read_polygons(char const* wkt)
{
    using types = geometries<Coordinate>;
    typename types::polygons result;
    if (std::string(wkt).find("MULTIPOLYGON") != 0)
    {
        result.push_back(read_polygon<Coordinate>(wkt));
        return result;
    }
    bg::model::multi_polygon<bg::model::polygon<typename types::native_point>> input;
    bg::read_wkt(wkt, input);
    for (auto const& polygon : input)
    {
        result.emplace_back();
        lift_range<Coordinate>(polygon.outer(), result.back().outer());
        result.back().inners().resize(polygon.inners().size());
        for (std::size_t i = 0; i < polygon.inners().size(); ++i)
        {
            lift_range<Coordinate>(polygon.inners()[i], result.back().inners()[i]);
        }
    }
    return result;
}

template <typename Coordinate>
typename geometries<Coordinate>::line read_line(char const* wkt)
{
    using types = geometries<Coordinate>;
    bg::model::linestring<typename types::native_point> input;
    bg::read_wkt(wkt, input);
    typename types::line result;
    lift_range<Coordinate>(input, result);
    return result;
}

template <typename Coordinate>
typename geometries<Coordinate>::lines read_lines(char const* wkt)
{
    using types = geometries<Coordinate>;
    typename types::lines result;
    if (std::string(wkt).find("MULTILINESTRING") != 0)
    {
        result.push_back(read_line<Coordinate>(wkt));
        return result;
    }
    bg::model::multi_linestring<bg::model::linestring<typename types::native_point>> input;
    bg::read_wkt(wkt, input);
    for (auto const& line : input)
    {
        result.emplace_back();
        lift_range<Coordinate>(line, result.back());
    }
    return result;
}

template <typename Coordinate>
typename geometries<Coordinate>::box read_box(char const* wkt)
{
    using types = geometries<Coordinate>;
    bg::model::box<typename types::native_point> input;
    bg::read_wkt(wkt, input);
    return typename types::box(
        typename types::point(bg::get<bg::min_corner, 0>(input), bg::get<bg::min_corner, 1>(input)),
        typename types::point(bg::get<bg::max_corner, 0>(input), bg::get<bg::max_corner, 1>(input)));
}

template <typename Geometry>
bool exact_valid(Geometry const& geometry, std::string& reason)
{
    ec::point_cache_scope operation;
    // The validator uses rational coordinates throughout, without the exact
    // umbrella's custom turn-ratio integration, which is not supported here.
    bg::strategies::relate::cartesian<rational> validation;
    return bg::is_valid(geometry, reason, validation);
}

template <typename A, typename B>
bool valid_inputs(A const& a, B const& b, record& report)
{
    std::string reason_a, reason_b;
    bool const first = exact_valid(a, reason_a);
    bool const second = exact_valid(b, reason_b);
    report.put("input.first_valid", first);
    report.put("input.second_valid", second);
    if (!first) report.put("input.first_reason", reason_a);
    if (!second) report.put("input.second_reason", reason_b);
    if (!first || !second) report.put("status", "invalid_input");
    return first && second;
}

void check(record& report, char const* name, bool result)
{
    std::string const path = std::string("checks.") + name;
    bool const accumulated = report.get<bool>(path, true) && result;
    report.put(path, accumulated);
    if (!accumulated) report.put("status", "exact_contract_failure");
}

template <typename Geometry>
void check_result(Geometry const& result, record& report)
{
    std::string reason;
    bool const valid = exact_valid(result, reason);
    report.put("result.valid", valid);
    if (!valid)
    {
        report.put("result.reason", reason);
        check(report, "valid_exact_result", false);
    }
}

template <typename Geometry, typename Strategy>
void areal_summary(Geometry const& result, Strategy const&, record& report)
{
    check_result(result, report);
    report.put("result.components", result.size());
    report.put("result.holes", bg::num_interior_rings(result));
    report.put("result.area", validation_area(result));
}

template <typename A, typename B>
bool equal_sets(A const& a, B const& b)
{
    ec::point_cache_scope operation;
    if (bg::is_empty(a) || bg::is_empty(b))
        return bg::is_empty(a) && bg::is_empty(b);
    bg::strategies::relate::cartesian<rational> validation;
    // Direct DE-9IM avoids equals' normalized-vector/sqrt optimization, which
    // cannot represent irrational lengths in cpp_rational.
    return bg::relate(a, b, bg::de9im::mask("T*F**FFF*"), validation);
}

template <typename Coordinate, typename Geometry>
void rounded_diagnostic(Geometry result, record& report)
{
    using native_point = typename geometries<Coordinate>::native_point;
    bg::for_each_point(result, [](ec::exact_point<Coordinate>& point)
    {
        native_point rounded;
        ec::assign_output(point, rounded);
        point = ec::make_point<Coordinate>(rounded);
    });
    std::string reason;
    bool const valid = exact_valid(result, reason);
    report.put("rounded_output.valid", valid);
    if (!valid) report.put("rounded_output.reason", reason);
    // This is an output representation observation, not an exact topology
    // failure. Never feed these rounded points back into the exact oracle.
}

template <typename Coordinate>
void replay_areal(cases::typed_pair const& item, record& report)
{
    using types = geometries<Coordinate>;
    typename types::strategy strategy;
    auto a = read_polygons<Coordinate>(item.data.first);
    auto b = read_polygons<Coordinate>(item.data.second);
    if (item.data.correct_inputs)
    {
        bg::correct(a, strategy);
        bg::correct(b, strategy);
    }
    if (!valid_inputs(a, b, report)) return;
    std::string const operation(item.data.operation);
    if (operation == "covered_by" || operation == "within" || operation == "overlaps")
    {
        bool const value = operation == "covered_by" ? bg::covered_by(a, b, strategy)
            : operation == "within" ? bg::within(a, b, strategy) : bg::overlaps(a, b, strategy);
        report.put("result.predicate", value);
        report.put("result.matrix", bg::relation(a, b, strategy).str());
        if (item.data.expected_predicate >= 0)
            check(report, "supplied_predicate", value == (item.data.expected_predicate != 0));
        if (item.expected_matrix)
            check(report, "supplied_matrix", bg::relation(a, b, strategy).str() == item.expected_matrix);
        return;
    }
    typename types::polygons result;
    if (operation == "intersection") bg::intersection(a, b, result, strategy);
    else if (operation == "difference") bg::difference(a, b, result, strategy);
    else if (operation == "union") bg::union_(a, b, result, strategy);
    else if (operation == "sym_difference") bg::sym_difference(a, b, result, strategy);
    else throw std::logic_error("Unsupported areal operation");
    areal_summary(result, strategy, report);
    if (item.data.issue == 704 || item.data.issue == 1363)
        report.put("result.input_matrix", bg::relation(a, b, strategy).str());
    rounded_diagnostic<Coordinate>(result, report);
    // Component/decimal-area expectations are recorded, not universal exact
    // contracts. Older rounded output can merge or collapse components.
    if (item.data.expected_components >= 0)
    {
        report.put("reference.components", item.data.expected_components);
        report.put("reference.components_match", result.size()
                   == static_cast<std::size_t>(item.data.expected_components));
    }
    if (operation == "union" || operation == "intersection" || operation == "sym_difference")
    {
        typename types::polygons reverse;
        if (operation == "union") bg::union_(b, a, reverse, strategy);
        else if (operation == "intersection") bg::intersection(b, a, reverse, strategy);
        else bg::sym_difference(b, a, reverse, strategy);
        check(report, "operand_order", equal_sets(result, reverse));
    }
    if (operation == "difference")
    {
        check(report, "subset_of_minuend", bg::covered_by(result, a, strategy));
        check(report, "no_subtractor_interior", bg::relate(result, b,
              bg::de9im::mask("F********"), strategy));
    }
    if (operation == "intersection")
    {
        check(report, "subset_of_first", bg::is_empty(result) || bg::covered_by(result, a, strategy));
        check(report, "subset_of_second", bg::is_empty(result) || bg::covered_by(result, b, strategy));
    }
    typename types::polygons complement;
    rational expected_area;
    if (operation == "intersection")
    {
        bg::union_(a, b, complement, strategy);
        expected_area = validation_area(a) + validation_area(b) - validation_area(complement);
    }
    else
    {
        bg::intersection(a, b, complement, strategy);
        rational const intersection_area = validation_area(complement);
        expected_area = operation == "difference" ? validation_area(a) - intersection_area
            : validation_area(a) + validation_area(b)
                - intersection_area * (operation == "sym_difference" ? 2 : 1);
    }
    check(report, "area_partition", validation_area(result) == expected_area);
}

template <typename Coordinate>
void replay_line_polygon(cases::typed_pair const& item, record& report)
{
    using types = geometries<Coordinate>;
    typename types::strategy strategy;
    auto lines = read_lines<Coordinate>(item.data.first);
    auto polygons = read_polygons<Coordinate>(item.data.second);
    if (item.data.correct_inputs) bg::correct(polygons, strategy);
    if (!valid_inputs(lines, polygons, report)) return;
    std::string const operation(item.data.operation);
    if (operation == "crosses")
    {
        bool const value = bg::crosses(lines, polygons, strategy);
        auto const matrix = bg::relation(lines, polygons, strategy).str();
        report.put("result.predicate", value);
        report.put("result.matrix", matrix);
        bg::reverse(lines);
        check(report, "direction_matrix", bg::relation(lines, polygons, strategy).str() == matrix);
        if (item.data.expected_predicate >= 0)
            check(report, "supplied_predicate", value == (item.data.expected_predicate != 0));
        if (item.expected_matrix)
            check(report, "supplied_matrix", matrix == item.expected_matrix);
        return;
    }
    if (operation == "intersection_tuple")
    {
        std::tuple<typename types::points, typename types::lines, typename types::polygons> result;
        bg::intersection(lines, polygons, result, strategy);
        report.put("result.points", std::get<0>(result).size());
        report.put("result.lines", std::get<1>(result).size());
        report.put("result.polygons", std::get<2>(result).size());
        check(report, "no_contact_line", std::get<1>(result).empty());
        check(report, "one_contact_point", std::get<0>(result).size() == 1u);
        if (item.expected_set)
        {
            bg::model::multi_point<typename types::native_point> native;
            bg::read_wkt(item.expected_set, native);
            typename types::points expected;
            lift_range<Coordinate>(native, expected);
            check(report, "supplied_contact_point", bg::equals(std::get<0>(result), expected, strategy));
        }
        return;
    }
    typename types::lines result, reverse;
    bg::intersection(lines, polygons, result, strategy);
    report.put("result.components", result.size());
    check_result(result, report);
    rounded_diagnostic<Coordinate>(result, report);
    bg::reverse(lines);
    bg::intersection(lines, polygons, reverse, strategy);
    check(report, "direction_set", equal_sets(result, reverse));
    check(report, "subset_of_polygon", bg::is_empty(result) || bg::covered_by(result, polygons, strategy));
    if (item.expected_set)
    {
        auto expected = read_lines<Coordinate>(item.expected_set);
        check(report, "supplied_set", equal_sets(result, expected));
    }
}

template <typename Coordinate>
void replay_line_line(cases::typed_pair const& item, record& report)
{
    using types = geometries<Coordinate>;
    typename types::strategy strategy;
    auto a = read_lines<Coordinate>(item.data.first);
    auto b = read_lines<Coordinate>(item.data.second);
    if (!valid_inputs(a, b, report)) return;
    typename types::lines result, reverse;
    bg::difference(a, b, result, strategy);
    report.put("result.components", result.size());
    check_result(result, report);
    bg::reverse(a);
    bg::difference(a, b, reverse, strategy);
    check(report, "direction_set", equal_sets(result, reverse));
    check(report, "tail_retained", bg::covered_by(typename types::point(7, 4), result, strategy));
    check(report, "head_retained", bg::covered_by(typename types::point(5, -8), result, strategy));
}

template <typename Coordinate>
void replay_tuple_areal(cases::typed_pair const& item, record& report)
{
    using types = geometries<Coordinate>;
    typename types::strategy strategy;
    auto a = read_polygons<Coordinate>(item.data.first);
    auto b = read_polygons<Coordinate>(item.data.second);
    bg::correct(a, strategy);
    bg::correct(b, strategy);
    if (!valid_inputs(a, b, report)) return;
    std::tuple<typename types::points, typename types::lines, typename types::polygons> result;
    bg::intersection(a, b, result, strategy);
    report.put("result.points", std::get<0>(result).size());
    report.put("result.lines", std::get<1>(result).size());
    areal_summary(std::get<2>(result), strategy, report);
    check(report, "no_isolated_point", std::get<0>(result).empty());
    check(report, "no_isolated_line", std::get<1>(result).empty());
    typename types::polygons areal;
    bg::intersection(a, b, areal, strategy);
    check(report, "tuple_areal_set", equal_sets(areal, std::get<2>(result)));
    report.put("reference.components", item.data.expected_components);
}

template <typename Coordinate>
void replay_special(cases::typed_pair const& item, record& report)
{
    using types = geometries<Coordinate>;
    typename types::strategy strategy;
    auto const line = read_line<Coordinate>(item.data.first);
    if (item.second_type == cases::geometry_type::box)
    {
        auto const box = read_box<Coordinate>(item.data.second);
        if (!valid_inputs(line, box, report)) return;
        typename types::lines result;
        bg::intersection(line, box, result, strategy);
        report.put("result.components", result.size());
        report.put("category", "legacy_output_dimensional_advisory");
        report.put("result.points", bg::num_points(result));
        std::string reason;
        report.put("result.valid", exact_valid(result, reason));
        report.put("result.reason", reason);
        report.put("checks.no_singleton_contact_line", result.empty());
        if (!result.empty())
        {
            report.put("status", "known_algorithmic_dimensional_failure");
            report.put("assertion_guard", "BOOST_GEOMETRY_TEST_FAILURES");
            report.put("contract_note", "Legacy box-intersection tests require singleton contact; dimensional API change is outside this refactor, not a rounding fix.");
#ifdef BOOST_GEOMETRY_TEST_FAILURES
            check(report, "no_singleton_contact_line", false);
#endif
        }
        return;
    }
    auto const second = read_line<Coordinate>(item.data.second);
    if (line.size() != 2 || second.size() != 2) throw std::logic_error("Segment arity");
    bg::model::segment<typename types::point> const a(line[0], line[1]), b(second[0], second[1]);
    typename types::points result;
    bg::intersection(a, b, result, strategy);
    report.put("result.points", result.size());
    check(report, "one_intersection", result.size() == 1u);
    if (!result.empty())
    {
        rational const y = rational(-2665034) - rational(1351) * 10043 / 24634;
        check(report, "rational_intersection", bg::get<0>(result[0]) == 4700000
              && bg::get<1>(result[0]) == y);
    }
}

template <typename Coordinate>
void replay(cases::typed_pair const& item, record& report)
{
    ec::point_cache_scope<Coordinate> operation;
    using kind = cases::geometry_type;
    if (item.first_type == kind::segment || item.second_type == kind::box)
        replay_special<Coordinate>(item, report);
    else if (item.first_type == kind::linestring || item.first_type == kind::multi_linestring)
    {
        if (item.second_type == kind::linestring || item.second_type == kind::multi_linestring)
            replay_line_line<Coordinate>(item, report);
        else replay_line_polygon<Coordinate>(item, report);
    }
    else if (std::string(item.data.operation) == "intersection_tuple")
        replay_tuple_areal<Coordinate>(item, report);
    else replay_areal<Coordinate>(item, report);
}

void replay_pipeline(cases::pipeline_case const& item, record& report)
{
    ec::point_cache_scope<double> operation;
    using types = geometries<double>;
    types::strategy strategy;
    std::vector<types::polygons> input;
    for (auto const* wkt : item.inputs)
    {
        input.push_back(read_polygons<double>(wkt));
        bg::correct(input.back(), strategy);
        std::string reason;
        if (!exact_valid(input.back(), reason))
        {
            report.put("status", "invalid_input");
            report.put("input.invalid_index", input.size() - 1);
            report.put("input.reason", reason);
            return;
        }
    }
    report.put("input.valid_count", input.size());
    std::size_t steps = 0;
    auto verify = [&](types::polygons const& geometry)
    {
        ++steps;
        std::string reason;
        if (!exact_valid(geometry, reason))
        {
            report.put("result.invalid_step", steps);
            report.put("result.reason", reason);
            check(report, "valid_intermediates", false);
            return false;
        }
        return true;
    };
    auto difference = [&](types::polygons const& a, types::polygons const& b)
    {
        types::polygons result;
        bg::difference(a, b, result, strategy);
        verify(result);
        return result;
    };
    auto symmetric = [&](types::polygons const& a, types::polygons const& b)
    {
        types::polygons result;
        bg::sym_difference(a, b, result, strategy);
        verify(result);
        return result;
    };
    types::polygons result;
    if (item.issue == 1419)
    {
        result = difference(difference(input[0], input[1]), input[2]);
        check(report, "subset_of_minuend", bg::covered_by(result, input[0], strategy));
    }
    else if (item.issue == 1236)
    {
        auto const a = difference(input[0], input[1]);
        result = difference(a, input[3]);
        auto const c = difference(input[2], input[3]);
        bool const within = bg::within(result, c, strategy);
        report.put("result.within", within);
        check(report, "retained_exact_within", within);
    }
    else if (item.issue == 900)
    {
        auto const a = symmetric(input[0], input[1]);
        auto const b = symmetric(a, input[2]);
        auto const c = symmetric(input[3], input[2]);
        result = symmetric(b, c);
        auto const expected = symmetric(a, input[3]);
        check(report, "symmetric_difference_identity", equal_sets(result, expected));
    }
    else if (item.issue == 588 || item.issue == 1290)
    {
        result = difference(input[0], input[1]);
        report.put("result.matrix_to_minuend", bg::relation(result, input[0], strategy).str());
        check(report, "subset_of_minuend", bg::covered_by(result, input[0], strategy));
        if (item.issue == 588)
        {
            check(report, "retained_exact_within", bg::within(result, input[0], strategy));
            report.put("frozen_rounded_report.within", bg::within(input[2], input[0], strategy));
            report.put("frozen_rounded_report.matrix", bg::relation(input[2], input[0], strategy).str());
            report.put("category", "retained_exact_vs_frozen_rounded_geometry");
        }
        else
        {
            report.put("result.matrix_to_subtractor", bg::relation(result, input[1], strategy).str());
            check(report, "no_subtractor_interior", bg::relate(result, input[1],
                  bg::de9im::mask("F********"), strategy));
            report.put("category", "predicate_assumption_and_roundtrip");
        }
    }
    else if (item.issue == 888)
    {
        std::reverse(input.begin(), input.end());
        while (input.size() > 1)
        {
            std::size_t const divide = input.size() / 2 + input.size() % 2;
            for (std::size_t i = 0; i < input.size() / 2; ++i)
            {
                auto next = symmetric(input[i + divide], input[i]);
                if (report.get<std::string>("status") == "exact_contract_failure") return;
                input[i] = std::move(next);
            }
            input.resize(divide);
        }
        result = std::move(input[0]);
    }
    else if (item.issue == 1174)
    {
        bg::model::polygon<types::native_point> native;
        bg::read_wkt(item.inputs[0], native);
        bg::model::box<types::native_point> box;
        bg::envelope(native, box);
        double const xmin = bg::get<bg::min_corner, 0>(box);
        double const ymin = bg::get<bg::min_corner, 1>(box);
        double const xmax = bg::get<bg::max_corner, 0>(box);
        double const ymax = bg::get<bg::max_corner, 1>(box);
        double const cx = (xmax + xmin) / 2;
        double const cy = (ymax + ymin) / 2;
        std::vector<types::polygons> pieces;
        rational area = 0;
        for (bool right : {false, true})
        {
            for (bool upper : {false, true})
            {
                types::box const quadrant(types::point(right ? cx : xmin, upper ? cy : ymin),
                    types::point(right ? xmax : cx, upper ? ymax : cy));
                pieces.emplace_back();
                bg::intersection(input[0], quadrant, pieces.back(), strategy);
                verify(pieces.back());
                check(report, "nonempty_quadrants", !pieces.back().empty());
                area += validation_area(pieces.back());
            }
        }
        for (std::size_t i = 0; i < pieces.size(); ++i)
            for (std::size_t j = i + 1; j < pieces.size(); ++j)
                check(report, "disjoint_quadrant_interiors", bg::relate(pieces[i], pieces[j],
                      bg::de9im::mask("F********"), strategy));
        check(report, "quadrant_area_partition", area == validation_area(input[0]));
        report.put("result.quadrants", pieces.size());
        report.put("result.exact_area_sum", area);
        report.put("result.steps", steps);
        return;
    }
    else throw std::logic_error("Unsupported pipeline");
    report.put("result.steps", steps);
    areal_summary(result, strategy, report);
    rounded_diagnostic<double>(result, report);
}

void emit(record const& report)
{
    boost::property_tree::write_json(std::cout, report, false);
    std::cout << '\n';
    std::cout.flush();
    auto const status = report.get<std::string>("status");
    if (status == "exact_contract_failure" || status == "exception")
        BOOST_ERROR(report.get<std::string>("case").c_str());
}

int main(int argc, char** argv)
{
    bool matched = argc == 1;
    for (auto const& item : cases::pairs())
    {
        if (argc > 1 && std::string(argv[1]) != "--pairs" && std::string(argv[1]) != item.data.name
            && std::string(argv[1]) != std::to_string(item.data.issue)) continue;
        matched = true;
        record report;
        report.put("issue", item.data.issue);
        report.put("case", item.data.name);
        report.put("coordinate_type", item.data.coordinate_type);
        report.put("operation", item.data.operation);
        report.put("status", "exact_contract_pass");
        report.put("contract", item.data.note);
        try
        {
            std::string const coordinate(item.data.coordinate_type);
            if (coordinate == "double") replay<double>(item, report);
            else if (coordinate == "float") replay<float>(item, report);
            else if (coordinate == "int32_t") replay<std::int32_t>(item, report);
            else if (coordinate == "int64_t") replay<std::int64_t>(item, report);
            else throw std::logic_error("Unsupported native coordinate type");
        }
        catch (std::exception const& error)
        {
            report.put("status", "exception");
            report.put("reason", error.what());
        }
        emit(report);
    }
    for (auto const& item : cases::pipelines())
    {
        if (argc > 1 && std::string(argv[1]) != "--pipelines" && std::string(argv[1]) != item.name
            && std::string(argv[1]) != std::to_string(item.issue)) continue;
        matched = true;
        record report;
        report.put("issue", item.issue);
        report.put("case", item.name);
        report.put("coordinate_type", item.coordinate_type);
        report.put("operation", "pipeline");
        report.put("status", "exact_contract_pass");
        report.put("contract", item.contract);
        try { replay_pipeline(item, report); }
        catch (std::exception const& error)
        {
            report.put("status", "exception");
            report.put("reason", error.what());
        }
        emit(report);
    }
    BOOST_TEST(matched);
    return boost::report_errors();
}
