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

int main()
{
    ec::point_cache_scope operation;
    test_first_boundary_point_contract();
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
