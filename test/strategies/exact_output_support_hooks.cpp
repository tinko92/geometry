// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.

#include <boost/core/lightweight_test.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/algorithms/detail/overlay/append_no_duplicates.hpp>
#include <boost/geometry/algorithms/detail/overlay/append_no_dups_or_spikes.hpp>
#include <boost/geometry/algorithms/detail/overlay/clip_linestring.hpp>
#include <boost/geometry/algorithms/detail/overlay/convert_ring.hpp>
#include <boost/geometry/algorithms/detail/overlay/copy_segment_point.hpp>
#include <boost/geometry/algorithms/detail/overlay/copy_segments.hpp>
#include <boost/geometry/algorithms/detail/overlay/graph/traverse_graph.hpp>
#include <boost/geometry/algorithms/detail/overlay/traversal_info.hpp>
#include <boost/geometry/geometries/register/point.hpp>
#include <boost/geometry/strategies/detail/output_support.hpp>

namespace bg = boost::geometry;
namespace ov = bg::detail::overlay;
using point = bg::model::d2::point_xy<double>;
using line = bg::model::linestring<point>;
using ring = bg::model::ring<point>;
using polygon = bg::model::polygon<point>;
using box = bg::model::box<point>;
using strategy = bg::strategies::relate::cartesian<>;
namespace support = bg::strategies::detail;

struct supported_range : ring
{
    bg::signed_size_type pending_support = -1;
    std::vector<bg::signed_size_type> incoming_support;
};

struct mock_strategy : strategy
{
    mutable int duplicates = 0, spikes = 0, collinear = 0;
    mutable int copies = 0, preparations = 0, reversals = 0;
    mutable bool reverse = false, remove_spikes = false;
    mutable void const* source = nullptr;
    mutable void const* other_source = nullptr;
    mutable bg::segment_identifier segment;
    mutable bg::signed_size_type to_index = -1;

    template <typename Range>
    void append_no_duplicates(Range& range, point const& p) const
    {
        ++duplicates;
        range.push_back(p);
    }

    template <typename Range>
    void append_no_dups_or_spikes(Range& range, point const& p) const
    {
        ++spikes;
        range.push_back(p);
    }

    template <typename Range>
    void append_no_collinear(Range& range, point const& p) const
    {
        ++collinear;
        range.push_back(p);
    }

    template <bool Reverse, bool RemoveSpikes, typename Geometry>
    void copy_segments(Geometry const& geometry, bg::segment_identifier const& id,
                       bg::signed_size_type index, supported_range& range) const
    {
        ++copies;
        source = &geometry;
        segment = id;
        to_index = index;
        reverse = Reverse;
        remove_spikes = RemoveSpikes;
        range.incoming_support.push_back(range.pending_support);
        range.push_back(point(9, 9));
        range.pending_support = 42;
    }

    template <bool Reverse1, bool Reverse2, typename Geometry1, typename Geometry2>
    void prepare_output_edge(supported_range& range, Geometry1 const& geometry1,
                             Geometry2 const& geometry2, bg::segment_identifier const& id) const
    {
        ++preparations;
        source = &geometry1;
        other_source = &geometry2;
        segment = id;
        reverse = id.source_index == 0 ? Reverse1 : Reverse2;
        range.pending_support = id.segment_index;
    }

    template <typename Range>
    void reverse_range(Range& range) const
    {
        ++reversals;
        boost::reverse(range);
    }
};

struct source_point
{
    double x = -1, y = -1;
    void const* source = nullptr;
    bg::segment_identifier segment;
    bg::signed_size_type offset = 0;
    bool reverse = false;

    template <bool Reverse, typename Geometry>
    static void copy_source_vertex(Geometry const& geometry, bg::segment_identifier const& id,
                                   bg::signed_size_type offset, source_point& point)
    {
        point.source = &geometry;
        point.segment = id;
        point.offset = offset;
        point.reverse = Reverse;
    }
};

BOOST_GEOMETRY_REGISTER_POINT_2D(source_point, double, bg::cs::cartesian, x, y)

void check_point(point const& actual, point const& expected)
{
    BOOST_TEST_EQ(bg::get<0>(actual), bg::get<0>(expected));
    BOOST_TEST_EQ(bg::get<1>(actual), bg::get<1>(expected));
}

template <typename Range1, typename Range2>
void check_range(Range1 const& actual, Range2 const& expected)
{
    BOOST_TEST_EQ(actual.size(), expected.size());
    for (std::size_t i = 0; i < actual.size() && i < expected.size(); ++i)
    {
        check_point(actual[i], expected[i]);
    }
}

void test_defaults()
{
    strategy const s;
    line output;
    ov::append_no_duplicates(output, point(0, 0), s);
    ov::append_no_duplicates(output, point(0, 0), s);
    check_range(output, line{{0, 0}});
    ov::append_no_dups_or_spikes(output, point(2, 0), s);
    ov::append_no_dups_or_spikes(output, point(1, 0), s);
    check_range(output, line{{0, 0}, {1, 0}});
    ov::append_no_collinear(output, point(2, 0), s);
    check_range(output, line{{0, 0}, {2, 0}});
    BOOST_TEST(!support::try_append_no_duplicates(s, output, point(3, 0)));
    BOOST_TEST(!support::try_append_no_dups_or_spikes(s, output, point(3, 0)));
    BOOST_TEST(!support::try_append_no_collinear(s, output, point(3, 0)));
    check_range(output, line{{0, 0}, {2, 0}});

    ring const input{{0, 0}, {0, 4}, {4, 4}, {4, 0}, {0, 0}};
    box const bounds(point(0, 0), point(4, 4));
    bg::segment_identifier id(0, -1, -1, 0);
    BOOST_TEST((!support::try_copy_segments<false, true>(s, input, id, 3, output)));
    support::prepare_output_edge<false, true>(s, output, input, bounds, id);
    output.clear();
    bg::copy_segments<false>(input, id, 3, s, output);
    check_range(output, line{{0, 4}, {4, 4}, {4, 0}});
    output.clear();
    bg::copy_segments<false>(bounds, id, 3, s, output);
    check_range(output, line{{0, 4}, {4, 4}, {4, 0}});
    output.clear();
    bg::copy_segments<true>(input, id, 1, s, output);
    check_range(output, line{{4, 0}});
    output.clear();
    bg::copy_segments<true>(bounds, id, 1, s, output);
    check_range(output, line{{4, 0}});

    line const source{{0, 0}, {2, 0}, {1, 0}, {1, 1}};
    id.segment_index = -1;
    output.clear();
    bg::copy_segments<false>(source, id, 3, s, output);
    check_range(output, line{{0, 0}, {1, 0}, {1, 1}});
    output.clear();
    bg::detail::copy_segments::copy_segments_linestring<false, false>::apply(
        source, id, 3, s, output);
    check_range(output, source);
    bg::copy_segments<false>(source, id, -1, s, output);
    check_range(output, source);

    point vertex;
    id.segment_index = 0;
    BOOST_TEST(!support::try_copy_source_vertex<false>(input, id, 0, vertex));
    BOOST_TEST(bg::copy_segment_point<false>(input, id, -1, vertex));
    check_point(vertex, point(4, 0));
    BOOST_TEST(bg::copy_segment_point<true>(input, id, 1, vertex));
    check_point(vertex, point(4, 0));
    BOOST_TEST(bg::copy_segment_point<false>(bounds, id, -1, vertex));
    check_point(vertex, point(4, 0));
    BOOST_TEST(bg::copy_segment_point<true>(bounds, id, 1, vertex));
    check_point(vertex, point(4, 0));

    ring const reversed(input.rbegin(), input.rend());
    ring converted;
    ov::convert_ring<bg::ring_tag>::apply(converted, input, false, true, s);
    check_range(converted, reversed);
    polygon poly;
    ov::convert_ring<bg::polygon_tag>::apply(poly, input, false, true, s);
    check_range(bg::exterior_ring(poly), reversed);
    ov::convert_ring<bg::polygon_tag>::apply(poly, input, true, true, s);
    BOOST_TEST_EQ(bg::interior_rings(poly).size(), 1u);
    check_range(bg::interior_rings(poly).front(), reversed);
}

void test_strategy_hooks()
{
    mock_strategy const s;
    line output;
    ov::append_no_duplicates(output, point(0, 0), s);
    ov::append_no_duplicates(output, point(0, 0), s);
    check_range(output, line{{0, 0}, {0, 0}});
    BOOST_TEST_EQ(s.duplicates, 2);
    output = {{0, 0}, {2, 0}};
    ov::append_no_dups_or_spikes(output, point(1, 0), s);
    check_range(output, line{{0, 0}, {2, 0}, {1, 0}});
    BOOST_TEST_EQ(s.spikes, 1);
    output = {{0, 0}, {1, 0}};
    ov::append_no_collinear(output, point(2, 0), s);
    check_range(output, line{{0, 0}, {1, 0}, {2, 0}});
    BOOST_TEST_EQ(s.collinear, 1);

    ring const input{{0, 0}, {0, 4}, {4, 4}, {4, 0}, {0, 0}};
    box const bounds(point(0, 0), point(4, 4));
    bg::segment_identifier id(0, -1, -1, 0);
    supported_range out;
    bg::copy_segments<true>(input, id, 2, s, out);
    BOOST_TEST(s.source == &input);
    BOOST_TEST(s.reverse && s.remove_spikes);
    BOOST_TEST(s.segment == id);
    BOOST_TEST_EQ(s.to_index, 2);
    bg::copy_segments<false>(bounds, id, 3, s, out);
    BOOST_TEST(s.source == &bounds);
    BOOST_TEST(!s.reverse && s.remove_spikes);
    bg::copy_segments<false>(output, id, 1, s, out);
    BOOST_TEST(s.source == &output);
    BOOST_TEST(!s.reverse && s.remove_spikes);
    // The hook runs even when the legacy linestring path would return immediately.
    bg::detail::copy_segments::copy_segments_linestring<true, false>::apply(
        output, id, -1, s, out);
    BOOST_TEST(s.reverse && !s.remove_spikes);
    BOOST_TEST_EQ(s.copies, 4);
    check_range(out, line{{9, 9}, {9, 9}, {9, 9}, {9, 9}});

    ring converted;
    polygon poly;
    ov::convert_ring<bg::ring_tag>::apply(converted, input, false, true, s);
    ov::convert_ring<bg::polygon_tag>::apply(poly, input, false, true, s);
    ov::convert_ring<bg::polygon_tag>::apply(poly, input, true, true, s);
    BOOST_TEST_EQ(s.reversals, 3);
    ring const reversed(input.rbegin(), input.rend());
    check_range(converted, reversed);
    check_range(bg::exterior_ring(poly), reversed);
    check_range(bg::interior_rings(poly).front(), reversed);
    ov::convert_ring<bg::ring_tag>::apply(converted, input, false, false, s);
    check_range(converted, input);
    BOOST_TEST_EQ(s.reversals, 3);
}

void test_source_vertex_hook()
{
    ring const input{{0, 0}, {0, 4}, {4, 4}, {4, 0}, {0, 0}};
    line const ls{{1, 2}, {3, 4}};
    box const bounds(point(0, 0), point(4, 4));
    bg::segment_identifier id(0, -1, -1, 0);
    source_point vertex;
    BOOST_TEST(bg::copy_segment_point<true>(input, id, -1, vertex));
    BOOST_TEST(vertex.source == &input);
    BOOST_TEST(vertex.segment == id);
    BOOST_TEST_EQ(vertex.offset, -1);
    BOOST_TEST(vertex.reverse);
    BOOST_TEST(bg::copy_segment_point<false>(ls, id, 1, vertex));
    BOOST_TEST(vertex.source == &ls);
    BOOST_TEST(!vertex.reverse);
    BOOST_TEST(bg::copy_segment_point<false>(bounds, id, 2, vertex));
    BOOST_TEST(vertex.source == &bounds);
    BOOST_TEST_EQ(vertex.offset, 2);
    BOOST_TEST(!vertex.reverse);
    BOOST_TEST_EQ(vertex.x, -1);
    BOOST_TEST_EQ(vertex.y, -1);
}

void test_graph_preparation()
{
    using turn = ov::traversal_turn_info<point, bg::segment_ratio<double>>;
    using turns = std::vector<turn>;
    using clusters = std::map<bg::signed_size_type, ov::cluster_info>;
    ring const input{{0, 0}, {0, 4}, {4, 4}, {4, 0}, {0, 0}};
    box const bounds(point(0, 0), point(4, 4));
    turns ts(2);
    clusters const cs;
    auto& op = ts[0].operations[0];
    op.seg_id = bg::segment_identifier(0, -1, -1, 0);
    op.enriched.travels_to_vertex_index = -1;
    op.enriched.travels_to_ip_index = 1;
    mock_strategy const s;
    ov::traverse_graph<false, true, bg::overlay_union, ring, box, turns, clusters,
                       mock_strategy> graph(input, bounds, ts, cs, s);
    supported_range out;
    ov::turn_operation_id const toi{0, 0};
    graph.use_vertices(out, toi);
    BOOST_TEST_EQ(s.preparations, 1);
    BOOST_TEST_EQ(s.copies, 0);
    BOOST_TEST_EQ(out.pending_support, 0);
    BOOST_TEST(s.source == &input && s.other_source == &bounds);
    BOOST_TEST(!s.reverse);
    // No matching next operation: the same-segment early return must also prepare.
    out.pending_support = -1;
    op.enriched.travels_to_vertex_index = 0;
    graph.use_vertices(out, toi);
    BOOST_TEST_EQ(s.preparations, 2);
    BOOST_TEST_EQ(s.copies, 0);
    BOOST_TEST_EQ(out.pending_support, 0);
    op.enriched.travels_to_vertex_index = 1;
    graph.use_vertices(out, toi);
    BOOST_TEST_EQ(s.preparations, 3);
    BOOST_TEST_EQ(s.copies, 1);
    BOOST_TEST_EQ(out.incoming_support.back(), 0);
    BOOST_TEST_EQ(out.pending_support, 42);
    op.seg_id.source_index = 1;
    graph.use_vertices(out, toi);
    BOOST_TEST_EQ(s.preparations, 4);
    BOOST_TEST_EQ(s.copies, 2);
    BOOST_TEST(s.source == &bounds && s.reverse);

    strategy const default_strategy;
    ov::traverse_graph<false, true, bg::overlay_union, ring, box, turns, clusters,
                       strategy> default_graph(input, bounds, ts, cs, default_strategy);
    ring ordinary;
    op.enriched.travels_to_vertex_index = -1;
    default_graph.use_vertices(ordinary, toi);
    BOOST_TEST(ordinary.empty());
    op.enriched.travels_to_vertex_index = 1;
    default_graph.use_vertices(ordinary, toi);
    check_range(ordinary, ring{{4, 0}});
}

struct clipping_strategy : bg::strategy::intersection::liang_barsky<box, point>
{
    using base = bg::strategy::intersection::liang_barsky<box, point>;
    explicit clipping_strategy(int source_index) : source_index(source_index) {}
    int source_index;
    mutable line const* source = nullptr;
    mutable bg::segment_identifier active;
    mutable std::vector<bg::segment_identifier> prepared;
    mutable std::size_t clip_calls = 0;

    void prepare_segment(line const& input, bg::signed_size_type index) const
    {
        source = &input;
        active = bg::segment_identifier(source_index, -1, -1, index);
        prepared.push_back(active);
    }

    bool clip_segment(box const& bounds, bg::model::referring_segment<point>& segment,
                       bool& first_clipped, bool& second_clipped) const
    {
        BOOST_TEST_EQ(prepared.size(), clip_calls + 1);
        ++clip_calls;
        BOOST_TEST(source != nullptr);
        if (source)
        {
            check_point(point(bg::get<0, 0>(segment), bg::get<0, 1>(segment)),
                        (*source)[active.segment_index]);
            check_point(point(bg::get<1, 0>(segment), bg::get<1, 1>(segment)),
                        (*source)[active.segment_index + 1]);
        }
        return base::clip_segment(bounds, segment, first_clipped, second_clipped);
    }
};

void test_clipping_preparation()
{
    box const bounds(point(0, 0), point(2, 2));
    for (line const& input : {line{}, line{{1, 1}}, line{{1, 1}, {1, 1}},
            line{{-1, 1}, {3, 1}}, line{{-1, 1}, {3, 1}, {3, 3}, {1, 1}, {1, 1}}})
    {
        clipping_strategy const s(7);
        std::vector<line> actual, expected;
        bg::detail::intersection::clip_range_with_box<line>(
            bounds, input, std::back_inserter(actual), s);
        bg::detail::intersection::clip_range_with_box<line>(
            bounds, input, std::back_inserter(expected), clipping_strategy::base{});
        auto const count = input.empty() ? 0 : input.size() - 1;
        BOOST_TEST_EQ(s.prepared.size(), count);
        BOOST_TEST_EQ(s.clip_calls, count);
        if (count != 0) BOOST_TEST(s.source == &input);
        for (std::size_t i = 0; i < s.prepared.size(); ++i)
        {
            BOOST_TEST_EQ(s.prepared[i].segment_index, bg::signed_size_type(i));
            BOOST_TEST_EQ(s.prepared[i].source_index, 7);
        }
        BOOST_TEST_EQ(actual.size(), expected.size());
        for (std::size_t i = 0; i < actual.size() && i < expected.size(); ++i)
        {
            check_range(actual[i], expected[i]);
        }
    }
}

struct conversion_strategy : strategy
{
    mutable std::vector<line const*> copied;
    void convert_geometry(line const& source, line& destination) const
    {
        copied.push_back(&source);
        bg::convert(source, destination);
    }
};

template <typename Geometry>
void check_difference_conversion(Geometry const& input, line const& cutter,
                                 std::vector<line const*> const& expected_sources)
{
    conversion_strategy const s;
    std::vector<line> actual, expected;
    bg::difference(input, cutter, actual, s);
    bg::difference(input, cutter, expected);
    BOOST_TEST(s.copied == expected_sources);
    BOOST_TEST_EQ(actual.size(), expected.size());
    for (std::size_t i = 0; i < actual.size() && i < expected.size(); ++i)
    {
        check_range(actual[i], expected[i]);
    }
}

void test_whole_line_conversion()
{
    line const cutter{{0, -1}, {0, 1}};
    line const isolated{{4, 0}, {5, 0}};
    check_difference_conversion(isolated, cutter, {&isolated});
    using multi = bg::model::multi_linestring<line>;
    multi const disjoint{{{-5, 0}, {-4, 0}}, {{4, 0}, {5, 0}}};
    check_difference_conversion(disjoint, cutter, {&disjoint[0], &disjoint[1]});
    multi const partially_intersecting{{{-5, 0}, {-4, 0}}, {{-1, 0}, {1, 0}}, {{4, 0}, {5, 0}}};
    check_difference_conversion(partially_intersecting, cutter,
                                {&partially_intersecting[0], &partially_intersecting[2]});
}

int main()
{
    test_defaults();
    test_strategy_hooks();
    test_source_vertex_hook();
    test_graph_preparation();
    test_clipping_preparation();
    test_whole_line_conversion();
    return boost::report_errors();
}
