// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.

#define BOOST_GEOMETRY_EXACT_FILTER_STATISTICS
#include <boost/core/lightweight_test.hpp>
#include <boost/geometry/algorithms/covered_by.hpp>
#include <boost/geometry/algorithms/envelope.hpp>
#include <boost/geometry/strategies/cartesian/exact/borrowed_bounds.hpp>
#include <boost/geometry/geometries/pointing_segment.hpp>

#include <cstdint>
#include <limits>
#include <type_traits>

namespace bg = boost::geometry;
namespace ec = bg::strategies::detail::exact_cartesian;

template <typename Root>
struct bounds_strategy : bg::strategies::relate::cartesian_exact<typename Root::coordinate_type>
{
    ec::borrowed_context<Root>& context;
    explicit bounds_strategy(ec::borrowed_context<Root>& owner) : context(owner) {}
    template <typename Point> using geometry_box_type = ec::root_bounds<Root>;
    template <typename G, typename Box> auto envelope(G const&, Box const&) const
    { return ec::root_envelope<Root>{context}; }
    template <typename Box, typename G> auto expand(Box const&, G const&) const
    { return ec::root_expand<Root>{context}; }
    template <typename A, typename B> auto covered_by(A const&, B const&) const
    { return ec::root_box_predicate<Root, ec::box_relation::covered_by>{}; }
    template <typename A, typename B, typename Apply>
    decltype(auto) with_sources(A const&, B const&, Apply&& apply) const
    { return std::forward<Apply>(apply)(*this); }
};

template <ec::box_relation Relation, typename Root, typename A, typename B>
bool compare(A const& a, B const& b) { return ec::root_box_predicate<Root, Relation>::apply(a, b); }

template <typename Policy>
void test_bounds()
{
    using C = std::int64_t;
    using native = bg::model::d2::point_xy<C>;
    using polygon = bg::model::polygon<native>;
    using multi = bg::model::multi_polygon<polygon>;
    using box = bg::model::box<native>;
    using root = ec::borrowed_root<C, Policy, multi, box>;
    using point = ec::borrowed_point<root>;
    using bounds = ec::root_bounds<root>;
    using ring = typename ec::borrowed_geometry_types<root>::template ring<true, true>;
    static_assert(bg::dimension<bounds>::value == 2, "Root bounds expose Box metadata traits");
    static_assert(!std::is_constructible<point, bounds>::value, "Bounds are not geometry vertices");

    polygon a, b;
    a.outer() = {{0, 0}, {0, 4}, {4, 4}, {4, 0}, {0, 0}};
    b.outer() = {{-2, 1}, {-2, 6}, {3, 6}, {3, 1}, {-2, 1}};
    multi source{a, b};
    box const containing(native(-10, -10), native(10, 10));
    ec::borrowed_context<root> context;
    ec::root_envelope<root> envelope{context};
    ec::root_expand<root> expand{context};
    bounds_strategy<root> strategy(context);
    ec::filtering_statistics() = {};
    bounds exact;
    bg::envelope(source, exact, strategy);
    BOOST_TEST(exact.initialized);
    BOOST_TEST_EQ(exact.origins[0].native_coordinate(0), -2);
    BOOST_TEST_EQ(exact.origins[1].native_coordinate(1), 0);
    BOOST_TEST_EQ(exact.origins[2].native_coordinate(0), 4);
    BOOST_TEST_EQ(exact.origins[3].native_coordinate(1), 6);
    BOOST_TEST((compare<ec::box_relation::covered_by, root>(exact, containing)));
    BOOST_TEST((compare<ec::box_relation::within, root>(exact, containing)));
    BOOST_TEST(!(compare<ec::box_relation::disjoint, root>(exact, containing)));
    box const adjacent(native(4, 1), native(7, 3));
    BOOST_TEST((compare<ec::box_relation::touches, root>(exact, adjacent)));
    BOOST_TEST(!(compare<ec::box_relation::disjoint, root>(adjacent, exact)));
    bounds native_box;
    envelope.apply(containing, native_box);
    BOOST_TEST((compare<ec::box_relation::covered_by, root>(exact, native_box)));
    BOOST_TEST((compare<ec::box_relation::within, root>(native_box, native_box)));
    bounds degenerate;
    envelope.apply(a.outer().front(), degenerate);
    BOOST_TEST(!(compare<ec::box_relation::within, root>(degenerate, native_box)));
    BOOST_TEST((compare<ec::box_relation::covered_by, root>(a.outer().front(), exact)));

    ec::search_box conservative;
    envelope.apply(source, conservative);
    BOOST_TEST((compare<ec::box_relation::covered_by, root>(exact, conservative)));
    ec::search_box const fractional(ec::search_point(-2.5, -0.5), ec::search_point(4.5, 6.5));
    BOOST_TEST((compare<ec::box_relation::covered_by, root>(exact, fractional)));
    BOOST_TEST(!(compare<ec::box_relation::covered_by, root>(fractional, exact)));
    double const infinity = (std::numeric_limits<double>::infinity)();
    ec::search_box const unbounded(ec::search_point(-infinity, -infinity), ec::search_point(infinity, infinity));
    BOOST_TEST((compare<ec::box_relation::covered_by, root>(exact, unbounded)));
    bounds empty;
    BOOST_TEST(!(compare<ec::box_relation::covered_by, root>(empty, containing)));
    BOOST_TEST((compare<ec::box_relation::disjoint, root>(empty, exact)));

    ring generated;
    for (auto it = a.outer().cbegin(); it != a.outer().cend(); ++it) generated.push_back(point::from_iterator(it));
    bounds root_geometry;
    envelope.apply(generated, root_geometry);
    expand.apply(root_geometry, b);
    BOOST_TEST((compare<ec::box_relation::covered_by, root>(exact, root_geometry)));
    BOOST_TEST((compare<ec::box_relation::covered_by, root>(root_geometry, exact)));
    BOOST_TEST(bg::covered_by(generated, containing, strategy));
    BOOST_TEST(bg::covered_by(source, containing, strategy));
    using metadata_point = bg::model::point<double, 2, bg::cs::cartesian>;
    using metadata_box = bg::model::box<metadata_point>;
    metadata_box rtree;
    envelope.apply(generated, rtree);
    BOOST_TEST_EQ((bg::get<bg::max_corner, 0>(rtree)), 4.0);
    metadata_box const extra(metadata_point(-3.5, -1.5), metadata_point(8.5, 7.5));
    expand.apply(rtree, extra);
    BOOST_TEST_EQ((bg::get<bg::min_corner, 0>(rtree)), -3.5);
    envelope.apply(extra, conservative);
    BOOST_TEST_EQ((bg::get<bg::max_corner, 1>(conservative)), 7.5);
    BOOST_TEST((compare<ec::box_relation::covered_by, root>(exact, rtree)));
    using reverse_iterator = boost::reverse_iterator<typename polygon::ring_type::const_iterator>;
    ec::iterator_point_view<reverse_iterator> const reversed(reverse_iterator(a.outer().cend()));
    bounds reversed_origin;
    envelope.apply(reversed, reversed_origin);
    BOOST_TEST_EQ(reversed_origin.origins[0].native_coordinate(0), 0);
    typename ec::root_envelope<root>::template state<bounds> state;
    bounds first, second, merged;
    envelope.apply(a, first); envelope.apply(b, second);
    ec::root_envelope<root>::apply(state, first);
    ec::root_envelope<root>::apply(state, second);
    ec::root_envelope<root>::result(state, merged);
    BOOST_TEST((compare<ec::box_relation::covered_by, root>(exact, merged)));
    BOOST_TEST((compare<ec::box_relation::covered_by, root>(merged, exact)));
    multi const no_geometry;
    envelope.apply(no_geometry, merged);
    BOOST_TEST(!merged.initialized);
    BOOST_TEST_EQ(ec::filtering_statistics().native_point_constructions, 0u);
    BOOST_TEST_EQ(ec::filtering_statistics().homogeneous_values, 0u);
    BOOST_TEST_EQ(ec::filtering_statistics().point_coordinate_values, 0u);
    BOOST_TEST(ec::operation_detail::require<C>().storage().auxiliary.empty());

#ifdef BOOST_GEOMETRY_TEST_REJECT_QUERY_BOUNDS
    auto const query = ec::query_point<root>::midpoint(generated[0], generated[1]);
    envelope.apply(query, exact);
#endif
}

template <typename C, typename Policy>
void test_segments_and_lines()
{
    using native = bg::model::d2::point_xy<C>;
    using segment = bg::model::segment<native>;
    using line = bg::model::linestring<native>;
    using root = ec::borrowed_root<C, Policy, segment, line>;
    using point = ec::borrowed_point<root>;
    ec::borrowed_context<root> context;
    ec::root_envelope<root> envelope{context};
    ec::borrowed_strategy<root> strategy(context);
    ec::root_bounds<root> bounds;
    segment const native_segment(native(3, -2), native(-1, 5));
    ec::filtering_statistics() = {};
    envelope.apply(native_segment, bounds);
    BOOST_TEST_EQ(bounds.origins[0].native_coordinate(0), -1.0);
    BOOST_TEST_EQ(bounds.origins[3].native_coordinate(1), 5.0);
    bg::envelope(native_segment, bounds, strategy);
    BOOST_TEST_EQ(bounds.origins[0].native_coordinate(0), -1.0);
    line const range{{-3, 1}, {2, 7}};
    bg::model::pointing_segment<native const> const native_pointed(range[0], range[1]);
    envelope.apply(native_pointed, bounds);
    BOOST_TEST_EQ(bounds.origins[0].native_coordinate(0), -3.0);
    BOOST_TEST_EQ(bounds.origins[3].native_coordinate(1), 7.0);
    bg::envelope(native_pointed, bounds, strategy);
    BOOST_TEST_EQ(bounds.origins[0].native_coordinate(0), -3.0);
    BOOST_TEST_EQ(bounds.origins[3].native_coordinate(1), 7.0);
    bg::model::referring_segment<native const> const native_referred(range[0], range[1]);
    envelope.apply(native_referred, bounds);
    BOOST_TEST_EQ(bounds.origins[2].native_coordinate(0), 2.0);
    bg::envelope(native_referred, bounds, strategy);
    BOOST_TEST_EQ(bounds.origins[2].native_coordinate(0), 2.0);
    auto const first = point::from_iterator(range.cbegin());
    auto const second = point::from_iterator(range.cbegin() + 1);
    bg::model::referring_segment<point const> const borrowed(first, second);
    envelope.apply(borrowed, bounds);
    BOOST_TEST_EQ(bounds.origins[0].native_coordinate(0), -3.0);
    BOOST_TEST_EQ(bounds.origins[3].native_coordinate(1), 7.0);
    bg::envelope(borrowed, bounds, strategy);
    BOOST_TEST_EQ(bounds.origins[0].native_coordinate(0), -3.0);
    BOOST_TEST_EQ(bounds.origins[3].native_coordinate(1), 7.0);
    bg::model::pointing_segment<point const> const pointed(first, second);
    envelope.apply(pointed, bounds);
    BOOST_TEST_EQ(bounds.origins[0].native_coordinate(0), -3.0);
    BOOST_TEST_EQ(bounds.origins[3].native_coordinate(1), 7.0);
    bg::envelope(pointed, bounds, strategy);
    BOOST_TEST_EQ(bounds.origins[0].native_coordinate(0), -3.0);
    BOOST_TEST_EQ(bounds.origins[3].native_coordinate(1), 7.0);
    ec::search_box search;
    ec::root_envelope<root>::apply(pointed, search);
    BOOST_TEST_EQ((bg::get<bg::min_corner, 0>(search)), -3.0);
    BOOST_TEST_EQ((bg::get<bg::max_corner, 1>(search)), 7.0);
    bg::envelope(native_pointed, search, strategy);
    BOOST_TEST_EQ((bg::get<bg::min_corner, 0>(search)), -3.0);
    BOOST_TEST_EQ((bg::get<bg::max_corner, 1>(search)), 7.0);
    bg::envelope(pointed, search, strategy);
    BOOST_TEST_EQ((bg::get<bg::min_corner, 0>(search)), -3.0);
    BOOST_TEST_EQ((bg::get<bg::max_corner, 1>(search)), 7.0);
    envelope.apply(range, bounds);
    BOOST_TEST_EQ(bounds.origins[2].native_coordinate(0), 2.0);
    BOOST_TEST_EQ(ec::filtering_statistics().native_point_constructions, 0u);
    BOOST_TEST_EQ(ec::filtering_statistics().point_coordinate_values, 0u);
    BOOST_TEST(ec::operation_detail::require<C>().storage().auxiliary.empty());
}

template <typename Policy>
void test_implicit_origin()
{
    using C = double;
    using native = bg::model::d2::point_xy<C>;
    using polygon = bg::model::polygon<native>;
    using root = ec::borrowed_root<C, Policy, polygon, polygon>;
    using point = ec::borrowed_point<root>;
    C const delta = std::ldexp(1.0, -54);
    polygon a, b;
    a.outer() = {{0, -1}, {2, 1}, {1, 0}, {0, -1}};
    b.outer() = {{0, delta}, {2, delta}, {0, 0}, {0, delta}};
    ec::borrowed_context<root> context;
    ec::borrowed_stage<root> stage(context, a, b);
    typename root::operation_type operations[2];
    operations[0].seg_id = bg::segment_identifier(0, -1, -1, 0);
    operations[1].seg_id = bg::segment_identifier(1, -1, -1, 0);
    ec::borrowed_record<root> record;
    record.stage = &stage;
    record.operations = {{&operations[0], &operations[1]}};
    point crossing(record);
    record.numeric.init(crossing.producers());
    typename ec::borrowed_geometry_types<root>::line line;
    line.push_back(point::borrow(a.outer()[2]));
    line.push_back(crossing);
    ec::root_envelope<root> envelope{context};
    ec::root_bounds<root> exact;
    ec::filtering_statistics() = {};
    envelope.apply(line, exact);
    BOOST_TEST(!exact.origins[2].is_native());
    BOOST_TEST_EQ((ec::root_compare_coordinate<0, root>(exact.origins[2], a.outer()[2])), 1);
    ec::search_box const rounded_edge(ec::search_point(0, -1), ec::search_point(1, 1));
    BOOST_TEST(!(compare<ec::box_relation::covered_by, root>(exact, rounded_edge)));
    ec::search_box conservative;
    envelope.apply(line, conservative);
    BOOST_TEST((compare<ec::box_relation::covered_by, root>(exact, conservative)));
    BOOST_TEST(ec::operation_detail::require<C>().storage().auxiliary.empty());
    BOOST_TEST_EQ(ec::filtering_statistics().native_point_constructions, 0u);
    BOOST_TEST_EQ(ec::filtering_statistics().point_coordinate_values, 0u);
}

int main()
{
    test_bounds<bg::strategies::relate::approximate_output>();
    test_bounds<bg::strategies::relate::nearest_output>();
    test_segments_and_lines<std::int64_t, bg::strategies::relate::approximate_output>();
    test_segments_and_lines<std::int64_t, bg::strategies::relate::nearest_output>();
    test_segments_and_lines<double, bg::strategies::relate::approximate_output>();
    test_segments_and_lines<double, bg::strategies::relate::nearest_output>();
    test_implicit_origin<bg::strategies::relate::approximate_output>();
    test_implicit_origin<bg::strategies::relate::nearest_output>();
    return boost::report_errors();
}
