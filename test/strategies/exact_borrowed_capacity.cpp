// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.

#include <boost/core/lightweight_test.hpp>
#include <boost/geometry/strategies/relate/cartesian_exact.hpp>
#include <boost/geometry/algorithms/detail/overlay/append_no_dups_or_spikes.hpp>
#include <boost/geometry/views/detail/boundary_view.hpp>
#ifdef BOOST_GEOMETRY_TEST_PUBLIC_BOX_OUTPUT
#include <boost/geometry/algorithms/intersection.hpp>
#endif

#include <array>
#include <cstdint>
#include <type_traits>

namespace bg = boost::geometry;
namespace ec = bg::strategies::detail::exact_cartesian;

template <typename Policy>
void test_boundary_mapping()
{
    using C = std::int64_t;
    using native = bg::model::d2::point_xy<C>;
    using polygon = bg::model::polygon<native>;
    using multi = bg::model::multi_polygon<polygon>;
    using root = ec::borrowed_root<C, Policy, multi, polygon>;
    using point = ec::borrowed_point<root>;
    using query = ec::query_point<root>;
    using bounds = ec::root_bounds<root>;
    using root_capacity = ec::construction_capacity<C, ec::construction_family::original_support_intersection>;
    using intersection_result = ec::borrowed_intersection_result<root>;
    using native_box = bg::model::box<native>;
    static_assert(std::is_same<typename ec::borrowed_output_working<root, native_box>::type,
        native_box>::value, "Standalone native box output bypasses generic working-geometry selection");
    static_assert(std::is_same<typename ec::borrowed_strategy<root>::template geometry_box_type<point>, bounds>::value,
                  "Geometry envelope metadata uses four coordinate origins");
    static_assert(std::is_same<decltype(std::declval<point>().homogeneous()),
        typename root_capacity::homogeneous_type>::value, "Root values have the original-support capacity");
    static_assert(std::is_same<decltype(std::declval<query>().homogeneous()),
        ec::auxiliary_point_type<C>>::value, "Query predicates use the actual auxiliary format");
    static_assert(std::tuple_size<decltype(std::declval<intersection_result>().definitions)>::value == 1,
                  "Straight segments need one proper-crossing numeric seed, not two");
    static_assert(std::tuple_size<decltype(std::declval<intersection_result>().intersections)>::value == 2,
                  "Collinear overlap can still return two existing point references");
    static_assert(!ec::borrowed_detail::is_borrowed_point<query>::value,
                  "Query diagnostics do not make it a geometry record");
    static_assert(std::is_same<typename bg::helper_geometry<point>::type, point>::value,
                  "Point-on-border helpers copy Root references, not rational coordinate snapshots");
    static_assert(!ec::borrowed_query_detail::geometry_operand<root, bounds>::value,
                  "Coordinate-origin metadata is not a query construction leaf");

    auto square = [](C x, C y, C size)
    {
        return typename polygon::ring_type{{x, y}, {x, y + size}, {x + size, y + size},
                                           {x + size, y}, {x, y}};
    };
    polygon a, b, other;
    a.outer() = square(0, 0, 100);
    a.inners().push_back(square(10, 10, 10));
    a.inners().push_back(square(30, 30, 10));
    b.outer() = square(200, 200, 100);
    b.inners().push_back(square(210, 210, 10));
    other.outer() = square(400, 400, 100);
    multi source{polygon{}, a, b};
    std::array<std::pair<int, int>, 5> const selected{{{1, -1}, {1, 0}, {1, 1}, {2, -1}, {2, 0}}};
    ec::borrowed_context<root> context;
    ec::borrowed_stage<root>* retained;
    {
        bg::detail::boundary_view<multi const> view(source);
        BOOST_TEST_EQ(boost::size(view), 6u);
        retained = &context.template retain<ec::borrowed_stage<root>>(context, view, other);
    }
    // The temporary boundary owner has gone; all locators resolve the original
    // native multipolygon, not its allocated array of ring views.
    for (std::size_t i = 0; i != selected.size(); ++i)
    {
        bg::segment_identifier const flat(0, i + 1, -1, 0);
        auto const mapped = retained->source_id(source, flat);
        BOOST_TEST_EQ(mapped.multi_index, selected[i].first);
        BOOST_TEST_EQ(mapped.ring_index, selected[i].second);
        auto const& polygon = source[selected[i].first];
        auto const& ring = selected[i].second < 0 ? polygon.outer() : polygon.inners()[selected[i].second];
        auto const vertex = retained->vertex(flat, 0);
        BOOST_TEST_EQ(vertex.native_coordinate(0), bg::get<0>(ring[0]));
        auto const support = context.support_coordinates(retained->support(flat));
        BOOST_TEST_EQ(support[0], bg::get<0>(ring[0]));
        BOOST_TEST_EQ(support[1], bg::get<1>(ring[0]));
        BOOST_TEST_EQ(support[2], bg::get<0>(ring[1]));
        retained->reverse[0] = true;
        auto const reversed = retained->vertex(flat, 1);
        BOOST_TEST_EQ(reversed.native_coordinate(0), bg::get<0>(ring[3]));
        BOOST_TEST_EQ(reversed.native_coordinate(1), bg::get<1>(ring[3]));
        auto const reverse_support = context.support_coordinates(retained->support(flat));
        BOOST_TEST_EQ(reverse_support[0], bg::get<0>(ring[3]));
        BOOST_TEST_EQ(reverse_support[1], bg::get<1>(ring[3]));
        retained->reverse[0] = false;
    }
    bg::detail::boundary_view<polygon const> single(source[1]);
    ec::borrowed_stage<root> polygon_stage(context, single, other);
    bg::segment_identifier const hole(0, 2, -1, 0);
    auto const mapped = polygon_stage.source_id(source[1], hole);
    BOOST_TEST_EQ(mapped.multi_index, -1);
    BOOST_TEST_EQ(mapped.ring_index, 1);
    BOOST_TEST_EQ(polygon_stage.vertex(hole, 0).native_coordinate(0), 30);

    ec::borrowed_stage<root> ordinary(context, source, other);
    bg::segment_identifier const ordinary_id(0, 2, 0, 0);
    BOOST_TEST_EQ(ordinary.source_id(source, ordinary_id).multi_index, 2);
    BOOST_TEST_EQ(ordinary.vertex(ordinary_id, 0).native_coordinate(0), 210);

    using working = typename ec::borrowed_working<root, multi>::type;
    working converted;
    ec::borrowed_geometry_converter<root>(context).apply(source, converted);
    bg::detail::boundary_view<working const> working_view(converted);
    ec::borrowed_stage<root> working_stage(context, working_view, other);
    auto const working_edge = working_stage.support(bg::segment_identifier(0, 3, -1, 0));
    auto const working_support = context.support_coordinates(working_edge);
    BOOST_TEST_EQ(working_support[0], 30);
    BOOST_TEST_EQ(working_stage.vertex(bg::segment_identifier(0, 3, -1, 0), 0).native_coordinate(0), 30);

    using box_root = ec::borrowed_root<C, Policy, multi, native_box>;
    ec::borrowed_context<box_root> box_context;
    native_box const box(native(0, 0), native(4, 4));
    ec::borrowed_point<box_root> sample;
    auto const box_status = ec::ring_topology<C>().analyze(box, sample);
    BOOST_TEST(box_status.valid && box_status.clockwise);
    BOOST_TEST_EQ(sample.native_coordinate(0), 0);
    BOOST_TEST_EQ(sample.native_coordinate(1), 0);
    auto const boundary = ec::ring_topology_detail::make_boundary<C>(box);
    static_assert(std::is_same<typename decltype(boundary)::point_type,
        ec::borrowed_detail::indexed_vertex<native_box>>::value, "Native box corners retain typed accessors");
    BOOST_TEST_EQ(ec::ring_topology_detail::locate(native(2, 2), boundary), 1);
    BOOST_TEST(ec::ring_topology<C>().representative(box, sample));

#ifdef BOOST_GEOMETRY_TEST_GENERATED_OPEN_SUPPORT
    using open_ring = typename ec::borrowed_geometry_types<root>::template ring<true, false>;
    open_ring generated;
    auto const range_id = context.register_range(a.outer());
    generated.push_back(point::borrow(a.outer()[0]));
    for (std::size_t i = 1; i != a.outer().size(); ++i)
    {
        generated.pending = {range_id, i - 1};
        generated.push_back(point::borrow(a.outer()[i]));
    }
    ec::borrowed_strategy<root> strategy(context);
    bg::detail::overlay::fix_closure(generated, strategy);
    BOOST_TEST(generated.incoming.front().valid());
#endif

#ifdef BOOST_GEOMETRY_TEST_REJECT_BOUNDS_GEOMETRY
    bounds metadata;
    ec::root_envelope<root>::apply(source[1], metadata);
    typename ec::borrowed_geometry_types<root>::template ring<true, true> forbidden;
    ec::borrowed_geometry_converter<root>(context).apply(metadata, forbidden);
#endif
#ifdef BOOST_GEOMETRY_TEST_REJECT_WORKING_BOX
    using forbidden_box = typename ec::borrowed_working<root, native_box>::type;
    (void)sizeof(forbidden_box);
#endif
#ifdef BOOST_GEOMETRY_TEST_BOUND_BOX_ALIAS
    static_assert(std::is_same<typename ec::borrowed_strategy<root>::template overlay_output_type<native_box>,
        native_box>::value, "RootBound uses the lazy standalone box output selector");
#endif
#ifdef BOOST_GEOMETRY_TEST_PUBLIC_BOX_OUTPUT
    native_box const first(native(0, 0), native(3, 3)), second(native(1, 1), native(4, 4));
    native_box box_output;
    bg::strategies::relate::cartesian_exact<C, Policy> const public_strategy;
    BOOST_TEST(bg::intersection(first, second, box_output, public_strategy));
    BOOST_TEST_EQ((bg::get<bg::min_corner, 0>(box_output)), 1);
    BOOST_TEST_EQ((bg::get<bg::max_corner, 1>(box_output)), 3);
    bg::model::linestring<native> const input_line{{-2, 0}, {2, 0}};
    native_box const clip_bounds(native(-1, -1), native(1, 1));
    std::vector<bg::model::linestring<native>> clipped;
    BOOST_TEST(bg::intersection(input_line, clip_bounds, clipped, public_strategy));
    BOOST_TEST_EQ(clipped.size(), 1u);
    if (!clipped.empty())
    {
        BOOST_TEST_EQ(bg::get<0>(clipped.front().front()), -1);
        BOOST_TEST_EQ(bg::get<0>(clipped.front().back()), 1);
    }
#endif
#ifdef BOOST_GEOMETRY_TEST_REJECT_BOUNDS_TOPOLOGY
    bounds metadata;
    ec::root_envelope<root>::apply(source[1], metadata);
    (void)ec::ring_topology<C>().analyze(metadata);
#endif
}

int main()
{
    test_boundary_mapping<bg::strategies::relate::approximate_output>();
    test_boundary_mapping<bg::strategies::relate::eager_nearest_output>();
    test_boundary_mapping<bg::strategies::relate::deferred_nearest_output>();
    return boost::report_errors();
}
