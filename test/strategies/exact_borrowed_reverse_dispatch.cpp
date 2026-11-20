// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.

#include <boost/core/lightweight_test.hpp>
#include <boost/geometry/algorithms/intersection.hpp>
#include <boost/geometry/strategies/relate/cartesian_exact.hpp>

#include <array>
#include <cstdint>
#include <vector>

namespace bg = boost::geometry;
namespace ec = bg::strategies::detail::exact_cartesian;

template <typename Policy>
void test_reverse_dispatch()
{
    using C = std::int64_t;
    using native = bg::model::d2::point_xy<C>;
    using ccw_polygon = bg::model::polygon<native, false, true>;
    using cw_polygon = bg::model::polygon<native, true, true>;
    using native_multi = bg::model::multi_polygon<cw_polygon>;
    using root = ec::borrowed_root<C, Policy, ccw_polygon, native_multi>;
    using point = ec::borrowed_point<root>;
    using working_multi = typename ec::borrowed_working<root, native_multi>::type;
    using working_polygon = typename ec::borrowed_geometry_types<root>::template polygon<true, true>;
    ccw_polygon a;
    a.outer() = {{0, 0}, {4, 0}, {4, 4}, {0, 4}, {0, 0}};
    native_multi original;
    original.resize(1);
    original[0].outer() = {{2, -1}, {2, 5}, {6, 5}, {6, -1}, {2, -1}};
    ec::borrowed_context<root> context;
    ec::borrowed_strategy<root> strategy(context);
    working_multi b;
    strategy.convert_geometry(original, b);
    {
        ec::borrowed_stage<root> stage(context, a, b);
        struct model_turn { std::array<typename root::operation_type, 2> operations; } model;
        model.operations[0].seg_id = bg::segment_identifier(1, 0, -1, 0);
        model.operations[1].seg_id = bg::segment_identifier(0, -1, -1, 0);
        auto* previous = context.active_stage;
        context.active_stage = &stage;
        {
            auto const scope = strategy.template model_turn_scope<false, true>(model);
            BOOST_TEST(stage.reverse[0]);
            BOOST_TEST(!stage.reverse[1]);
            (void)scope;
        }
        context.active_stage = previous;
    }
    std::vector<working_polygon> result;
    BOOST_TEST(bg::intersection(a, b, result, strategy));
    BOOST_TEST_EQ(result.size(), 1u);
    b.clear(); // Retained records must now resolve original supports, not this frame.
    if (result.empty()) return;
    auto const& ring = result.front().outer();
    std::array<native, 4> const corners{{{2, 0}, {2, 4}, {4, 4}, {4, 0}}};
    std::array<bool, 4> found{{false, false, false, false}};
    auto const count = ec::borrowed_detail::range_size(ring);
    BOOST_TEST_EQ(count, 4u);
    for (std::size_t i = 0; i != count; ++i)
    {
        for (std::size_t j = 0; j != corners.size(); ++j)
            if (ec::root_equal<root>(ring[i], corners[j])) found[j] = true;
        auto const support = ec::borrowed_detail::incoming_at(ring, (i + 1) % count);
        BOOST_TEST(support.valid());
        if (!support.valid()) continue;
        auto const values = context.support_coordinates(support);
        native const first(values[0], values[1]), second(values[2], values[3]);
        BOOST_TEST_EQ(ec::root_orientation<root>(first, second, ring[i]), 0);
        BOOST_TEST_EQ(ec::root_orientation<root>(first, second, ring[(i + 1) % count]), 0);
    }
    for (bool value : found) BOOST_TEST(value);
}

int main()
{
    test_reverse_dispatch<bg::strategies::relate::approximate_output>();
    test_reverse_dispatch<bg::strategies::relate::nearest_output>();
    return boost::report_errors();
}
