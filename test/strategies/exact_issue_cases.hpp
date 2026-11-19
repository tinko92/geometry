// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Licensed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_GEOMETRY_TEST_STRATEGIES_EXACT_ISSUE_CASES_HPP
#define BOOST_GEOMETRY_TEST_STRATEGIES_EXACT_ISSUE_CASES_HPP

#include <array>

namespace exact_issue_cases
{

// Data only, not evidence that any issue is fixed. Expectations come from the
// linked reports/comments or existing tests; the runner must validate inputs and
// compare exact internal geometry before interpreting rounded public output.
// Parse using coordinate_type (in particular float), then correct orientation
// and closure if correct_inputs is true. -1 means no supplied expectation.
struct wkt_pair
{
    unsigned issue;
    char const* name;
    char const* operation;
    char const* coordinate_type;
    char const* first;
    char const* second;
    bool correct_inputs;
    int expected_components;
    int expected_predicate;
    double reference_area; // Reported/test value, not an exact area oracle.
    char const* note;
};

inline std::array<wkt_pair, 16> const& cases()
{
    static std::array<wkt_pair, 16> const result = {{
        {1487, "issue_1487", "intersection", "double",
         "POLYGON((5613.4732841326168 -2310.7179505603272,5514.0601555991116 -2350.5086913635655,"
         "5508.3308534808984 -2352.8018811960578,5507.275791306095 -2346.7215468157037,"
         "5501.2357570754521 -2311.9127661882035,5484.909411045599 -2280.4588267460131,"
         "5484.1182457186051 -2278.9345869535609,5484.1182457171381 -2495.9701154741069,"
         "5616.9999440144747 -2442.8237527817514,5616.9999440144647 -2309.3063823768307,"
         "5613.4732841326168 -2310.7179505603272))",
         "POLYGON((5110.0842867836363 -2645.5659190779725,5768.7833420500028 -2382.1176086992587,"
         "5611.4496151045878 -2258.7436467911998,5616.9999440144638 -2309.3063823768302,"
         "5508.3308534808984 -2352.8018811960578,5501.2357570754521 -2311.9127661882053,"
         "5484.1182457186051 -2278.9345869535618,5508.7291498683508 -2145.887437971141,"
         "5110.0842867836363 -2645.5659190779725))",
         true, 1, -1, 18983.884971862572,
         "Priority: exact epsilon clustering disabled; PR 1508 is not assumed present or unnecessary."},
        {1506, "issue_1506", "union", "double",
         "MULTIPOLYGON(((7.5766965679702691 55.1955384390566266,7.5769710638206007 55.1955418668735760,"
         "7.5776907371234223 55.1955508539257096,7.5777041477885838 55.1955510213938609,"
         "7.5724506758568699 55.1833987553893692,7.5766965679702691 55.1955384390566266)))",
         "MULTIPOLYGON(((7.5776907371234223 55.1955508539257096,7.5769710638206007 55.1955418668735760,"
         "7.5786734102424926 55.1979584470965179,7.5776907371234223 55.1955508539257096)))",
         true, 1, -1, 6.951063106e-6,
         "Run both argument orders; compare union area with the exact input-area sum."},
        {1490, "issue_1490", "difference", "int64_t",
         "MULTIPOLYGON(((0 0,20 0,20 20,0 20,0 0),(5 15,14 14,9 6,5 15)))",
         "MULTIPOLYGON(((5 15,5 5,15 5,15 15,5 15)))",
         true, 1, -1, 300,
         "Existing difference_multi test; small-integer graph/topology control, not a pure FP claim."},
        {1420, "issue_1420", "union", "int32_t",
         "MULTIPOLYGON(((4 6,2 6,2 4,0 4,0 8,4 8,6 8,6 6,6 4,4 4,4 6)),"
         "((2 2,2 4,4 4,4 2,2 2)),((2 2,2 0,0 0,0 2,2 2)),((4 2,6 2,6 0,4 0,4 2)))",
         "MULTIPOLYGON(((2 2,2 4,4 4,4 2,2 2)))",
         false, 4, -1, -1,
         "Reporter explicitly identifies a small-integer traversal regression; no FP-only diagnosis."},
        {1421, "issue_1421_difference", "difference", "int32_t",
         "MULTIPOLYGON(((8 8,8 6,8 4,6 4,6 2,6 0,0 0,0 6,2 6,2 8,8 8),"
         "(2 4,2 2,4 2,4 4,2 4),(6 6,4 6,4 4,6 4,6 6)))",
         "MULTIPOLYGON(((6 6,6 4,4 4,4 6,6 6)),((6 6,6 8,8 8,8 6,6 6)),"
         "((4 6,2 6,2 8,4 8,4 6)),((4 4,4 2,2 2,2 4,4 4)),((4 2,6 2,6 0,4 0,4 2)),"
         "((2 6,2 4,0 4,0 6,2 6)),((2 2,2 0,0 0,0 2,2 2)),((6 2,6 4,8 4,8 2,6 2)))",
         false, 6, -1, -1, "Also run sym_difference: the report expects one component."},
        {1395, "issue_1395", "intersection", "double",
         "MULTIPOLYGON(((2601.6180555555557 201.34722222222223,2601.6180555555557 582.8472222222223,"
         "372.41666666666663 1065.3993055555554,372.41666666666663 201.34722222222223,"
         "2601.6180555555557 201.34722222222223)))",
         "MULTIPOLYGON(((2482.777777777778 201.34722222222223,2400.270833333333 1421.5416666666665,"
         "1983.1944444444443 1421.5416666666665,1900.6875 201.34722222222223,"
         "2482.777777777778 201.34722222222223)))",
         true, -1, -1, -1,
         "Original CCW/open custom model and /fp:fast; normalize model and rerun strict FP."},
        {1241, "issue_1241", "difference", "double",
         "POLYGON((-124.19999999845255 51.901455507812500,-124.19999999460376 51.935823093966235,"
         "-123.99999999648789 51.935823093966235,-123.99999999317222 51.902564290309876,"
         "-124.19999999845255 51.901455507812500))",
         "POLYGON((-123.99999999367975 51.907655109375000,-123.99999999291659 51.900000006653443,"
         "-124.19999999861555 51.900000005468293,-124.19999999792353 51.906179355468751,"
         "-123.99999999367975 51.907655109375000))",
         true, 1, -1, -1, "Missing near-collinear turns; keep every input digit."},
        {1211, "issue_1211", "intersection", "float",
         "POLYGON((0.875 0.4375,0.875 0.375,0.9375 0.375))",
         "POLYGON((1 0,1 0.581662834,0.824054539 0))",
         true, 1, -1, -1,
         "Original CCW/open float model; report bounds area below 0.01, not an exact value."},
        {1201, "issue_1201", "intersection", "double",
         "POLYGON((1e-14 0,-10 50,-20 50,1e-14 0))",
         "POLYGON((0 1000,0 0,1000 0,1000 1000,0 1000))",
         true, 1, -1, -1,
         "Comment corrects the original expectation: the exact result is a tiny nonempty triangle."},
        {1178, "issue_1178", "difference", "int32_t",
         "POLYGON((-1 -1,-1 9,9 9,9 -1,-1 -1))",
         "MULTIPOLYGON(((0 1,0 2,1 2,1 1,0 1)),((5 0,1 0,1 1,4 1,4 2,1 2,1 3,5 3,5 0)))",
         false, 2, -1, -1, "Touching subtractor components enclose an island; graph/topology control."},
        {631, "issue_631", "sym_difference", "int32_t",
         "POLYGON((0 0,0 100,100 100,100 0,0 0))",
         "POLYGON((0 -2147483648,0 100,100 100,100 -2147483648,0 -2147483648))",
         false, 1, -1, -1, "INT_MIN input; expected rectangle spans y=[INT_MIN,0], x=[0,100]."},
        {852, "issue_852", "intersection", "double",
         "POLYGON((0.33333327000590729 0.33333344706803125,0.50967067548140821 0,0 0,"
         "0 0.93199284523078707,0.33333327000590729 0.33333344706803125))",
         "POLYGON((0 0.8554426406197928,0 1,1 0,0.54614530500375813 0,0 0.8554426406197928))",
         false, 1, -1, -1, "Near-duplicate output vertices; check exact result separately from rounded output."},
        {1053, "issue_1053_self", "covered_by", "double",
         "POLYGON((0.3 0.3,0.5 0.75,0.6 0.4,0.3 0.3))",
         "POLYGON((0.3 0.3,0.5 0.75,0.6 0.4,0.3 0.3))",
         true, -1, 1, -1, "Self-relation should hold; also test within and the exact DE-9IM matrix."},
        {1238, "issue_1238", "within", "double",
         "MULTILINESTRING((10 10,20 20),(60 -60,30 15))",
         "MULTILINESTRING((10 10,20 20),(15 15,30 15))",
         false, -1, 0, -1, "Body calls within(second named geometry, first); this ordering preserves that call."},
        {1212, "issue_1212_simplified", "difference", "int64_t",
         "POLYGON((850000 0,1350000 500000,790000 1280000,30000 770000,850000 0),"
         "(510711 675561,620711 945561,680711 735561,510711 675561),"
         "(1035911 461661,715911 591661,755911 691661,1035911 461661))",
         "POLYGON((711861 626536,741861 656536,651861 836536,711861 626536))",
         true, -1, -1, -1,
         "Simplified comment data; original uses mpq_rational, and comments identify a combinatorial cause."},
        {494, "issue_494_representation_control", "union", "int32_t",
         "POLYGON((-8 0,-10 0,-10 4,10 4,10 0,8 0,4 0,0 0,-8 0))",
         "POLYGON((0 2,4 2,4 0,4 -2,0 -2,0 0,0 2))",
         false, 1, -1, -1,
         "Representation-only control: removing a collinear vertex does not change the geometric set."}
    }};
    return result;
}

} // namespace exact_issue_cases

#endif // BOOST_GEOMETRY_TEST_STRATEGIES_EXACT_ISSUE_CASES_HPP
