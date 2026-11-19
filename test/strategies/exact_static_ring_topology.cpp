// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.

#include <boost/core/lightweight_test.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/ring.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/algorithms/detail/overlay/ring_properties.hpp>
#include <boost/geometry/algorithms/detail/overlay/assign_parents.hpp>
#include <boost/geometry/algorithms/detail/overlay/add_rings.hpp>
#include <boost/geometry/strategies/cartesian/exact/static_ring_topology.hpp>
#include <boost/geometry/strategies/relate/cartesian.hpp>

#include <array>
#include <cfenv>
#include <cstdint>
#include <cstdlib>
#include <map>
#include <new>
#include <vector>

namespace bg = boost::geometry;
namespace ec = bg::strategies::detail::exact_cartesian;
namespace ov = bg::detail::overlay;
using id = bg::ring_identifier;

bool counting = false;
std::size_t allocations = 0;
void* operator new(std::size_t size)
{
    if (counting) ++allocations;
    if (void* result = std::malloc(size == 0 ? 1 : size)) return result;
    throw std::bad_alloc();
}
void operator delete(void* value) noexcept { std::free(value); }
void* operator new[](std::size_t size) { return ::operator new(size); }
void operator delete[](void* value) noexcept { ::operator delete(value); }
#if defined(__cpp_sized_deallocation)
void operator delete(void* value, std::size_t) noexcept { std::free(value); }
void operator delete[](void* value, std::size_t) noexcept { std::free(value); }
#endif

template <typename C> using point = bg::model::d2::point_xy<C>;
template <typename C> using ring = bg::model::ring<point<C>>;
template <typename C> struct strategy
{
    auto ring_topology() const { return ec::ring_topology<C>(); }
};

template <typename C>
ring<C> make_ring(std::initializer_list<std::array<int, 2>> points)
{
    ring<C> result;
    for (auto const& p : points) result.emplace_back(C(p[0]), C(p[1]));
    return result;
}

template <typename C>
void test_boundaries()
{
    ec::point_cache_scope<C> operation;
    ec::ring_topology<C> topology;
    auto square = make_ring<C>({{{0,0}},{{0,4}},{{4,4}},{{4,0}},{{0,0}}});
    ec::exact_point<C> representative;
    allocations = 0;
    counting = true;
    auto const simple = topology.analyze(square, representative);
    counting = false;
    BOOST_TEST_EQ(allocations, 0u);
    BOOST_TEST(simple.valid && simple.clockwise && simple.direction == 1);
    bg::model::polygon<point<C>> single_polygon;
    single_polygon.outer() = square;
    bg::model::box<point<C>> const single_box(point<C>(0,0), point<C>(4,4));
    allocations = 0;
    counting = true;
    auto const polygon_status = topology.analyze(single_polygon, representative);
    auto const box_status = topology.analyze(single_box, representative);
    counting = false;
    BOOST_TEST_EQ(allocations, 0u);
    BOOST_TEST(polygon_status.valid && polygon_status.clockwise);
    BOOST_TEST(box_status.valid && box_status.clockwise);
    std::reverse(square.begin(), square.end());
    auto const reverse = topology.analyze(square);
    BOOST_TEST(reverse.valid && !reverse.clockwise && reverse.direction == -1);

    auto const spike = make_ring<C>({{{0,0}},{{-4,0}},{{-2,0}},{{0,0}},{{0,4}},{{4,4}},{{4,0}},{{0,0}}});
    auto const spiked = topology.analyze(spike, representative);
    BOOST_TEST(spiked.valid && spiked.clockwise);
    BOOST_TEST_EQ(ec::compare_x(representative, ec::exact_point<C>(0,0)), 0);
    auto const bottom_spike = make_ring<C>({{{0,0}},{{3,0}},{{1,0}},{{0,0}},{{0,4}},{{4,4}},{{4,0}},{{0,0}}});
    BOOST_TEST(topology.analyze(bottom_spike).clockwise);
    auto const canceled = make_ring<C>({{{0,0}},{{0,4}},{{4,4}},{{4,0}},{{0,0}},
        {{2,0}},{{4,0}},{{4,2}},{{4,4}},{{2,4}},{{0,4}},{{0,2}},{{0,0}}});
    BOOST_TEST(!topology.analyze(canceled).valid);
    auto const line = make_ring<C>({{{0,0}},{{4,0}},{{2,0}},{{0,0}}});
    BOOST_TEST(!topology.analyze(line).valid);
    auto const touch = make_ring<C>({{{0,0}},{{0,2}},{{2,2}},{{2,0}},{{0,0}},
        {{0,-2}},{{-2,-2}},{{-2,0}},{{0,0}}});
    BOOST_TEST(topology.analyze(touch).valid);
    bg::model::ring<point<C>, true, false> open;
    for (auto const& p : make_ring<C>({{{0,0}},{{0,4}},{{4,4}},{{4,0}}})) open.push_back(p);
    BOOST_TEST(topology.analyze(open).clockwise);
    BOOST_TEST(!topology.analyze(make_ring<C>({{{0,0}},{{0,4}},{{4,4}},{{4,0}}})).valid);
    bg::model::ring<point<C>, false> ccw(square.begin(), square.end());
    BOOST_TEST_EQ(topology.analyze(ccw).direction, 1);
    bg::model::polygon<point<C>> polygon;
    polygon.outer() = square;
    polygon.inners().push_back(square);
    std::reverse(polygon.inners().back().begin(), polygon.inners().back().end());
    BOOST_TEST(!topology.analyze(polygon).valid);
    bg::model::box<point<C>> box(point<C>(0,0), point<C>(4,4));
    BOOST_TEST(topology.analyze(box).valid);

    bg::model::ring<ec::exact_point<C>> implicit;
    for (auto const& p : spike) implicit.push_back(ec::make_point<C>(p));
#ifdef BOOST_GEOMETRY_EXACT_FILTER_STATISTICS
    auto const q = ec::filtering_statistics().point_coordinate_values;
#endif
    BOOST_TEST(topology.analyze(implicit).clockwise);
#ifdef BOOST_GEOMETRY_EXACT_FILTER_STATISTICS
    BOOST_TEST_EQ(ec::filtering_statistics().point_coordinate_values, q);
#endif
}

template <typename C>
void test_parents()
{
    ec::point_cache_scope<C> operation;
    using properties = ov::ring_properties<point<C>, ec::ring_metadata<C>>;
    std::vector<ring<C>> rings;
    for (int inset : {0,2,4,6})
        rings.push_back(make_ring<C>({{{inset,inset}},{{inset,20-inset}},{{20-inset,20-inset}},
                                     {{20-inset,inset}},{{inset,inset}}}));
    std::reverse(rings[1].begin(), rings[1].end());
    std::reverse(rings[3].begin(), rings[3].end());
    std::map<id, properties> map;
    strategy<C> strategy;
    for (std::size_t i = 0; i != rings.size(); ++i) map.try_emplace(id(2,i,-1), rings[i], strategy);
    bg::model::polygon<point<C>> empty;
    {
        std::map<id, properties> single;
        single.try_emplace(id(2,0,-1), rings[0], strategy);
        allocations = 0;
        counting = true;
        ov::assign_parents<bg::overlay_union>(empty, empty, rings, single, strategy);
        counting = false;
        BOOST_TEST_EQ(allocations, 0u);
        BOOST_TEST(single.begin()->second.parent.source_index == -1);
        BOOST_TEST(single.begin()->second.children.empty());
        BOOST_TEST(!single.begin()->second.discarded);
    }
    ov::assign_parents<bg::overlay_union>(empty, empty, rings, map, strategy);
    BOOST_TEST(map.at(id(2,1,-1)).parent == id(2,0,-1));
    BOOST_TEST(map.at(id(2,3,-1)).parent == id(2,2,-1));
    BOOST_TEST(map.at(id(2,0,-1)).parent.source_index == -1);
    BOOST_TEST(map.at(id(2,2,-1)).parent.source_index == -1);

    rings.resize(2);
    rings[1] = rings[0];
    map.clear();
    for (int i = 0; i != 2; ++i) map.try_emplace(id(2,i,-1), rings[i], strategy);
    ov::assign_parents<bg::overlay_union>(empty, empty, rings, map, strategy);
    BOOST_TEST(!map.at(id(2,0,-1)).discarded);
    BOOST_TEST(map.at(id(2,1,-1)).discarded);
    std::reverse(rings[1].begin(), rings[1].end());
    map.clear();
    for (int i = 0; i != 2; ++i) map.try_emplace(id(2,i,-1), rings[i], strategy);
    ov::assign_parents<bg::overlay_union>(empty, empty, rings, map, strategy);
    BOOST_TEST(map.at(id(2,0,-1)).discarded && map.at(id(2,1,-1)).discarded);
    BOOST_TEST(map.at(id(2,0,-1)).parent.source_index == -1);
    BOOST_TEST(map.at(id(2,1,-1)).parent.source_index == -1);
}

void test_first_boundary_point_contract()
{
    using open_ring = bg::model::ring<point<double>, false, false>;
    open_ring inner;
    for (auto const& p : make_ring<double>({{{5,2}},{{4,4}},{{0,4}},{{0,0}},{{4,0}}}))
        inner.push_back(p);
    auto const outer = make_ring<double>({{{0,0}},{{0,4}},{{4,4}},{{4,0}},{{0,0}}});
    bg::strategies::relate::cartesian<> ordinary;
    // The noncanonical sample is on boundary. skip_first then omits the only
    // outside vertex of this open ring and reports "all on boundary".
    BOOST_TEST_EQ(ov::range_in_geometry(point<double>(0,0), inner, outer, ordinary), 0);
    BOOST_TEST_EQ(ov::range_in_geometry(point<double>(5,2), inner, outer, ordinary), -1);
    using properties = ov::ring_properties<point<double>, ec::ring_metadata<double>>;
    properties const metadata(inner, strategy<double>());
    BOOST_TEST(metadata.valid);
    BOOST_TEST_EQ(ec::compare_x(metadata.point, ec::exact_point<double>(5,2)), 0);
    BOOST_TEST_EQ(ec::compare_y(metadata.point, ec::exact_point<double>(5,2)), 0);
    BOOST_TEST_EQ(ec::compare_x(metadata.boundary_point, ec::exact_point<double>(0,0)), 0);
    point<double> selected;
    ec::assign_output(metadata.point, selected);
    BOOST_TEST_EQ(ov::range_in_geometry(selected, inner, outer, ordinary), -1);
}

void test_certified_sources()
{
    namespace rt = ec::ring_topology_detail;
    auto const square = make_ring<double>({{{0,0}},{{0,2}},{{0,4}},{{4,4}},{{4,0}},{{0,0}}});
    auto const diamond = make_ring<double>({{{0,2}},{{2,4}},{{4,2}},{{2,0}},{{0,2}}});
    auto outer = rt::selected_source<double, bg::ring_tag>(id(0,-1,-1), square, true);
    auto inner = rt::selected_source<double, bg::ring_tag>(id(1,-1,-1), diamond, true);
    auto const reference_outer = rt::make_boundary<double>(square);
    auto const reference_inner = rt::make_boundary<double>(diamond);
    BOOST_TEST(rt::certify_simple_cycle(outer));
    BOOST_TEST(rt::certify_simple_cycle(inner));

    auto const reference = rt::residual(reference_outer);
    auto const reduced = rt::residual(outer);
    BOOST_TEST_EQ(reduced.edges.size(), reference.edges.size());
    for (std::size_t i = 0; i < (std::min)(reduced.edges.size(), reference.edges.size()); ++i)
    {
        BOOST_TEST_EQ(reduced.edges[i].first, reference.edges[i].first);
        BOOST_TEST_EQ(reduced.edges[i].second, reference.edges[i].second);
        BOOST_TEST_EQ(reduced.edges[i].flow, reference.edges[i].flow);
    }
    // All inner vertices touch the outside ring. A midpoint is needed to
    // show strict containment; equality at the representative is not enough.
    for (auto const& p : inner.points) { BOOST_TEST_EQ(rt::locate(p, outer), 0); }
    auto const& sample = inner.points.front();
    BOOST_TEST_EQ(rt::classify(inner, outer, sample), 1);
    BOOST_TEST_EQ(rt::classify(reference_inner, reference_outer, sample), 1);
    BOOST_TEST_EQ(rt::classify(outer, outer, outer.points.front()), 0);

    // A stored hole can be certified with the opposite role orientation.
    bg::model::polygon<point<double>> polygon;
    polygon.outer() = square;
    polygon.inners().push_back(make_ring<double>({{{0,2}},{{1,3}},{{2,2}},{{1,1}},{{0,2}}}));
    std::reverse(polygon.inners()[0].begin(), polygon.inners()[0].end());
    auto const hole = rt::selected_source<double, bg::polygon_tag>(id(0,-1,0), polygon, true);
    BOOST_TEST(rt::certify_simple_cycle(hole));
    auto const wrong_role = rt::selected_source<double, bg::ring_tag>(id(0,-1,0), square, true);
    BOOST_TEST(!rt::certify_simple_cycle(wrong_role));

    // The same shape in a traversal collection has no original provenance.
    std::vector<ring<double>> const generated{square};
    auto const selected = rt::selected<double>(id(2,0,-1), square, square, generated, true);
    BOOST_TEST(!rt::certify_simple_cycle(std::get<2>(selected)));
    auto const ordinary = rt::selected_source<double, bg::ring_tag>(id(0,-1,-1), square);
    BOOST_TEST(!rt::certify_simple_cycle(ordinary));

    using properties = ov::ring_properties<point<double>, ec::ring_metadata<double>>;
    std::map<id, properties> coincident;
    auto const a = id(0,-1,-1), b = id(1,-1,-1);
    strategy<double> s;
    coincident.try_emplace(a, square, s);
    coincident.try_emplace(b, square, s);
    coincident.at(b).reversed = true;
    // Individually simple input rings can share their entire support. Keep
    // the bidirectional classification and signed multiplicity reduction.
    ov::assign_parents<bg::overlay_difference>(square, square, generated, coincident, s);
    BOOST_TEST(coincident.at(a).discarded && coincident.at(b).discarded);

    // A noncertifying original still needs cancellation along an excursion.
    auto const spur = make_ring<double>({{{0,0}},{{-4,0}},{{-2,0}},{{0,0}},
                                        {{0,4}},{{4,4}},{{4,0}},{{0,0}}});
    auto const spurred = rt::selected_source<double, bg::ring_tag>(id(0,-1,-1), spur, true);
    point<double> const erased(-1,0);
    BOOST_TEST(!rt::certify_simple_cycle(spurred));
    BOOST_TEST_EQ(rt::locate_raw(erased, spurred), 0);
    BOOST_TEST_EQ(rt::locate(erased, spurred), -1);
}

void test_generated_certificate()
{
    namespace rt = ec::ring_topology_detail;
    auto const triangle = make_ring<double>({{{0,0}},{{0,6}},{{6,0}},{{0,0}}});
    auto reverse = triangle;
    std::reverse(reverse.begin(), reverse.end());
    auto const flat = make_ring<double>({{{0,0}},{{1,0}},{{3,0}},{{0,0}}});
    // Nonlocal return along a differently subdivided edge. There is no
    // adjacent reversal to remove; each end of the edge goes around a lobe.
    auto const spur = make_ring<double>({{{-3,0}},{{3,0}},{{3,2}},{{5,2}},{{5,0}},{{3,0}},
                                        {{0,0}},{{-3,0}},{{-3,-2}},{{-5,-2}},{{-5,0}},{{-3,0}}});
    // Fully noded, no overlapping edges and no local spike. The two lobes
    // touch at a repeated vertex and have opposite signed areas.
    auto const zero = make_ring<double>({{{0,0}},{{-2,0}},{{-2,2}},{{0,2}},{{0,0}},
                                        {{0,-2}},{{2,-2}},{{2,0}},{{0,0}}});
    std::vector<ring<double>> const generated{triangle, reverse, flat, spur, zero};
    auto select = [&](int i, bool allow = true)
    { return std::get<2>(rt::selected<double>(id(2,i,-1), triangle, triangle, generated, allow)); };

    for (int i : {0, 1})
    {
        auto const candidate = select(i);
        auto const edge = point<double>(3,3);
        BOOST_TEST(rt::certify_simple_cycle(candidate));
        BOOST_TEST_EQ(rt::locate(edge, candidate), 0);
        auto const direct = rt::residual(candidate);
        auto const general = rt::residual_impl<false>(candidate);
        BOOST_TEST_EQ(direct.edges.size(), 3u);
        BOOST_TEST_EQ(direct.edges.size(), general.edges.size());
        for (std::size_t k = 0; k < (std::min)(direct.edges.size(), general.edges.size()); ++k)
        {
            BOOST_TEST_EQ(direct.edges[k].first, general.edges[k].first);
            BOOST_TEST_EQ(direct.edges[k].second, general.edges[k].second);
            BOOST_TEST_EQ(direct.edges[k].flow, general.edges[k].flow);
        }
    }
    BOOST_TEST(!rt::certify_simple_cycle(select(0, false)));
    auto const collinear = select(2);
    BOOST_TEST(!rt::certify_simple_cycle(collinear));
    BOOST_TEST(rt::residual(collinear).edges.empty());
    BOOST_TEST_EQ(rt::locate(point<double>(2,0), collinear), -1);

    auto const nonlocal = select(3);
    BOOST_TEST(!rt::certify_simple_cycle(nonlocal));
    BOOST_TEST_EQ(rt::locate_raw(point<double>(0,0), nonlocal), 0);
    BOOST_TEST_EQ(rt::locate(point<double>(0,0), nonlocal), -1);
    double twice_area = 0;
    for (std::size_t i = 1; i < zero.size(); ++i)
    {
        twice_area += bg::get<0>(zero[i-1]) * bg::get<1>(zero[i])
                    - bg::get<1>(zero[i-1]) * bg::get<0>(zero[i]);
    }
    BOOST_TEST_EQ(twice_area, 0);
    BOOST_TEST(!rt::certify_simple_cycle(select(4)));
    // A unique noncollinear extremum is not a certificate of nonzero area.
    for (std::size_t i = 0; i + 1 < zero.size(); ++i)
        if (i != 1) BOOST_TEST(rt::less(zero[1], zero[i]));
    BOOST_TEST(ec::static_orientation<double>(zero[0], zero[1], zero[2]) != 0);

    using properties = ov::ring_properties<point<double>, ec::ring_metadata<double>>;
    std::map<id, properties> coincident;
    auto const a = id(2,0,-1), b = id(2,1,-1);
    strategy<double> s;
    coincident.try_emplace(a, generated[0], s);
    coincident.try_emplace(b, generated[1], s);
    // Both triangles were individually certified. Coincident-ring
    // multiplicities still have to be reduced across the two walks.
    ov::assign_parents<bg::overlay_difference>(triangle, triangle, generated, coincident, s);
    BOOST_TEST(coincident.at(a).discarded && coincident.at(b).discarded);
}

int main()
{
    ec::point_cache_scope operation;
    test_first_boundary_point_contract();
    test_certified_sources();
    test_generated_certificate();
    int const previous = std::fegetround();
    for (int mode : {FE_TONEAREST})
    {
        if (std::fesetround(mode) != 0) continue;
        test_boundaries<std::int16_t>();
        test_boundaries<std::int64_t>();
        test_boundaries<double>();
        test_parents<std::int64_t>();
        test_parents<double>();
    }
    std::fesetround(previous);
    return boost::report_errors();
}
