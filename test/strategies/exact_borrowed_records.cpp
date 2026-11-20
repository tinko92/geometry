// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.

#include <boost/core/lightweight_test.hpp>
#include <boost/geometry/strategies/cartesian/exact/borrowed_records.hpp>
#include <boost/geometry/geometries/point_xy.hpp>

namespace bg = boost::geometry;
namespace ec = bg::strategies::detail::exact_cartesian;

template <typename Policy>
void test_records()
{
    using coordinate = std::int64_t;
    using native = bg::model::d2::point_xy<coordinate>;
    using polygon = bg::model::polygon<native>;
    using root = ec::borrowed_root<coordinate, Policy, polygon, polygon>;
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
    ec::borrowed_point<root> crossing(record);
    record.numeric.init(crossing.producers());
    BOOST_TEST(crossing.exact_coordinate(0) == ec::exact_number(3) / 2);
    BOOST_TEST(crossing.exact_coordinate(1) == ec::exact_number(3) / 2);
    BOOST_TEST_EQ(crossing.template emit<int>(0), 1);
    BOOST_TEST_EQ(crossing.template emit<double>(1), 1.5);

    using line = typename ec::borrowed_geometry_types<root>::line;
    line intermediate;
    intermediate.push_back(stage.vertex(operations[0].seg_id, 0));
    intermediate.pending = stage.support(operations[0].seg_id);
    intermediate.push_back(crossing);
    ec::borrowed_stage<root> next(context, intermediate, b);
    auto const original = next.support(bg::segment_identifier(0, -1, -1, 0));
    BOOST_TEST(original == intermediate.incoming[1]);
    auto const original_coordinates = context.support_coordinates(original);
    BOOST_TEST_EQ(original_coordinates[0], 0);
    BOOST_TEST_EQ(original_coordinates[2], 3);

    // Finalized definitions use the existing operation IDs, not the transient
    // intermediate geometry. Retain native source references through export.
    for (auto& operation : operations)
    {
        auto const edge = stage.support(operation.seg_id);
        operation.seg_id = bg::segment_identifier(edge.range, -1, -1, edge.segment);
    }
    stage.finalized = true;
    stage.sources = {};
    intermediate.clear();
    BOOST_TEST(crossing.exact_coordinate(0) == ec::exact_number(3) / 2);
    BOOST_TEST_EQ(crossing.refine(0).approximate, 1.5);
}

int main()
{
    test_records<bg::strategies::relate::approximate_output>();
    test_records<bg::strategies::relate::nearest_output>();
    return boost::report_errors();
}
