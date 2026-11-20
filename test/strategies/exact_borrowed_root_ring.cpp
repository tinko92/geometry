// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.

#include <boost/core/lightweight_test.hpp>
#include <boost/geometry/algorithms/detail/overlay/assign_parents.hpp>
#include <boost/geometry/algorithms/detail/overlay/ring_properties.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/strategies/cartesian/exact/static_ring_topology.hpp>

#include <cstdint>
#include <map>
#include <type_traits>
#include <vector>

namespace bg = boost::geometry;
namespace ec = bg::strategies::detail::exact_cartesian;
namespace rt = ec::ring_topology_detail;
namespace ov = bg::detail::overlay;

template <typename Root>
struct root_topology : ec::ring_topology<typename Root::coordinate_type>
{
    using point_type = ec::borrowed_point<Root>;
    using properties_tag = bg::strategies::detail::topological_ring_tag<point_type>;
};
template <typename Root>
struct topology_strategy
{
    static auto ring_topology() { return root_topology<Root>(); }
};

template <typename Policy>
void test_parenting()
{
    using C = std::int64_t;
    using native = bg::model::d2::point_xy<C>;
    using polygon = bg::model::polygon<native>;
    using root = ec::borrowed_root<C, Policy, polygon, polygon>;
    using point = ec::borrowed_point<root>;
    using query = ec::query_point<root>;
    using working_ring = typename ec::borrowed_geometry_types<root>::template ring<true, true>;
    using properties = ov::ring_properties<native, typename root_topology<root>::properties_tag>;
    static_assert(ec::static_predicates_detail::exact_input<point>::value, "Preserve Root geometry references");
    static_assert(!ec::static_predicates_detail::exact_input<query>::value, "Queries are not geometry samples");
    static_assert(std::is_same<typename ec::static_predicates_detail::common_root<native, point, query>::type,
        root>::value, "A Root query participates in predicate dispatch");

    polygon outer, inner;
    outer.outer() = {{0, 0}, {0, 4}, {4, 4}, {4, 0}, {0, 0}};
    inner.outer() = {{0, 0}, {4, 0}, {0, 4}, {0, 0}};
    ec::borrowed_context<root> context;
    working_ring generated;
    for (auto it = inner.outer().cbegin(); it != inner.outer().cend(); ++it)
        generated.push_back(point::from_iterator(it));
    std::vector<working_ring> const collection{generated};
    topology_strategy<root> strategy;
    std::map<bg::ring_identifier, properties> metadata;
    metadata.try_emplace(bg::ring_identifier(0, -1, -1), outer.outer(), strategy);
    metadata.try_emplace(bg::ring_identifier(2, 0, -1), collection[0], strategy);
    BOOST_TEST(metadata.at(bg::ring_identifier(0, -1, -1)).valid);
    BOOST_TEST(metadata.at(bg::ring_identifier(2, 0, -1)).valid);
    ov::assign_parents<bg::overlay_union>(outer, inner, collection, metadata, strategy);
    BOOST_TEST(metadata.at(bg::ring_identifier(2, 0, -1)).parent == bg::ring_identifier(0, -1, -1));
    BOOST_TEST(!metadata.at(bg::ring_identifier(2, 0, -1)).discarded);
    BOOST_TEST(!ec::operation_detail::require<C>().storage().auxiliary.empty());

    auto const root_inner = rt::make_boundary<C>(generated);
    auto const native_outer = rt::make_boundary<C>(outer.outer());
    static_assert(std::is_same<typename decltype(root_inner)::point_type, point>::value,
                  "Boundary construction retains Root references");
    auto const middle = rt::midpoint<C>(root_inner.points[1], root_inner.points[2]);
    static_assert(std::is_same<std::decay_t<decltype(middle)>, query>::value, "Root midpoint is a distinct query");
    BOOST_TEST_EQ(rt::locate(middle, native_outer), 1);
    BOOST_TEST_EQ((ec::static_compare_coordinate<0, C>(middle, native(2, 2))), 0);
    BOOST_TEST_EQ(ec::static_direction<C>(native(0, 0), middle, native(4, 4)), 1);
    BOOST_TEST(ec::static_equal_points<C>(middle, native(2, 2)));

    working_ring root_outer;
    for (auto it = outer.outer().cbegin(); it != outer.outer().cend(); ++it)
        root_outer.push_back(point::from_iterator(it));
    auto const native_inner = rt::make_boundary<C>(inner.outer());
    auto const borrowed_outer = rt::make_boundary<C>(root_outer);
    auto const legacy_middle = rt::midpoint<C>(native_inner.points[1], native_inner.points[2]);
    static_assert(std::is_same<std::decay_t<decltype(legacy_middle)>, ec::exact_point<C>>::value,
                  "An all-native midpoint retains the existing bounded legacy query");
    BOOST_TEST_EQ(rt::locate(legacy_middle, borrowed_outer), 1);
    BOOST_TEST_EQ(rt::classify(native_inner, borrowed_outer, native_inner.points.front()), 1);
    BOOST_TEST_EQ(rt::classify(root_inner, native_outer, root_inner.points.front()), 1);
}

int main()
{
    test_parenting<bg::strategies::relate::approximate_output>();
    test_parenting<bg::strategies::relate::nearest_output>();
    return boost::report_errors();
}
