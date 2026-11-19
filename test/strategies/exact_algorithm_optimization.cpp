// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.

#include <boost/core/lightweight_test.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/algorithms/detail/overlay/colocate_clusters.hpp>
#include <boost/geometry/algorithms/detail/overlay/get_clusters.hpp>
#include <boost/geometry/algorithms/detail/overlay/get_turn_info.hpp>
#include <boost/geometry/strategies/relate/cartesian_exact.hpp>

#include <array>
#include <cmath>
#include <iterator>
#include <map>
#include <set>
#include <type_traits>
#include <vector>

namespace bg = boost::geometry;
namespace ov = bg::detail::overlay;
namespace ex = bg::strategies::detail::exact_cartesian;
namespace services = bg::strategies::relate::services;

using exact_strategy = bg::strategies::relate::cartesian_exact<double>;
using point = ex::exact_point<double>;
using turn = ov::turn_info<point, ex::segment_position_tag>;
using clusters = std::map<bg::signed_size_type, ov::cluster_info>;

struct counters
{
    std::size_t distances = 0;
    std::size_t equalities = 0;
};

struct counted_topology : ex::topology<double>
{
    counters* counts;

    explicit counted_topology(counters* value) : counts(value) {}

    template <typename A, typename B, typename P>
    auto distance_measure(A const& a, B const& b, P const& p) const
    {
        if (counts) { ++counts->distances; }
        return ex::topology<double>::distance_measure(a, b, p);
    }

    template <typename A, typename B>
    bool cluster_equal(A const& a, B const& b) const
    {
        if (counts) { ++counts->equalities; }
        return ex::topology<double>::cluster_equal(a, b);
    }
};

// Keep the exact kernel identical while independently enabling algorithm paths.
template <bool Enabled>
struct algorithm_strategy : exact_strategy
{
    static constexpr bool exact_predicates = Enabled;
    using exact_cluster_compare = ex::point_compare<double, bg::strategy::compare::less>;

    counters* counts = nullptr;
    counted_topology topology() const { return counted_topology(counts); }
};

struct unspecified_strategy {};
static_assert(!services::exact_predicates<unspecified_strategy>::value, "No implicit side guarantee");
static_assert(!services::exact_predicates<bg::strategies::relate::cartesian<>>::value,
              "Default Cartesian side verification is unchanged");
static_assert(std::is_void<services::exact_cluster_compare<unspecified_strategy>::type>::value,
              "No implicit clustering guarantee");
static_assert(std::is_void<services::exact_cluster_compare<bg::strategies::relate::cartesian<>>::type>::value,
              "Default Cartesian clustering is unchanged");
static_assert(!services::exact_predicates<algorithm_strategy<false>>::value, "Explicit opt out");
static_assert(services::exact_predicates<algorithm_strategy<true>>::value, "Explicit opt in");

struct range
{
    using point_type = point;
    std::array<point, 3> points;
    bool last = false;

    std::size_t size() const { return last ? 2 : 3; }
    bool is_first_segment() const { return false; }
    bool is_last_segment() const { return last; }
    point const& at(std::size_t i) const { return points.at(i); }
};

void check_verified_side(range const& p, point const& q, int expected)
{
    range const other{{{q, q, q}}};
    int const side = exact_strategy::side().apply(p.at(0), p.at(1), q);
    BOOST_TEST_EQ(side, expected);
    using fun = ov::turn_info_verification_functions<ov::policy_verify_all>;
    counters old_counts, new_counts;
    algorithm_strategy<false> old_strategy;
    algorithm_strategy<true> new_strategy;
    old_strategy.counts = &old_counts;
    new_strategy.counts = &new_counts;
    BOOST_TEST_EQ(fun::verified_side(side, p, other, old_strategy, 0, 0), expected);
    BOOST_TEST_EQ(fun::verified_side(side, p, other, new_strategy, 0, 0), expected);
    BOOST_TEST_EQ(old_counts.distances, side == 0 ? 1u : 0u);
    BOOST_TEST_EQ(new_counts.distances, 0u);
}

void test_verified_side()
{
    check_verified_side(range{{{point(0, 0), point(1, 1), point(2, 2)}}}, point(2, 2), 0);
    check_verified_side(range{{{point(0, 0), point(0, 0), point(1, 1)}}}, point(2, 1), 0);
    double const tiny = std::ldexp(1.0, -200);
    range const p{{{point(0, 0), point(tiny, tiny), point(1, 1)}}};
    check_verified_side(p, point(tiny, std::nextafter(tiny, 2 * tiny)), 1);
    check_verified_side(p, point(tiny, std::nextafter(tiny, 0.0)), -1);

    point const third = ex::construct_intersection(point(0, 0), point(1, 1),
                                                   point(0, 1), point(1, -1));
    check_verified_side(range{{{point(0, 0), point(1, 1), point(2, 2)}}}, third, 0);

    using fun = ov::turn_info_verification_functions<ov::policy_verify_all>;
    range const last{{{point(0, 0), point(1, 1), point()}}, true};
    BOOST_TEST_EQ(fun::verified_side(0, last, last, algorithm_strategy<false>(), 1, 2), 0);
    BOOST_TEST_EQ(fun::verified_side(0, last, last, algorithm_strategy<true>(), 1, 2), 0);

    // A strategy without the capability must still recover a numerically zero
    // side even when the unrounded determinant is nonzero.
    range const off_line{{{point(0, 1), point(0, 1), point(0, 1)}}};
    range const horizontal{{{point(0, 0), point(1, 0), point(2, 0)}}};
    BOOST_TEST_EQ(fun::verified_side(0, horizontal, off_line, algorithm_strategy<false>(), 0, 0), 1);
}

template <bool Enabled>
std::vector<turn> get_turns(range const& p, range const& q)
{
    std::vector<turn> result;
    ov::get_turn_info<ov::assign_null_policy>::apply(p, q, turn(), algorithm_strategy<Enabled>(),
                                                   std::back_inserter(result));
    return result;
}

void compare_turns(range const& p, range const& q)
{
    auto const old_turns = get_turns<false>(p, q);
    auto const new_turns = get_turns<true>(p, q);
    BOOST_TEST_EQ(new_turns.size(), old_turns.size());
    if (new_turns.size() != old_turns.size()) { return; }
    for (std::size_t i = 0; i < old_turns.size(); ++i)
    {
        auto const& a = old_turns[i];
        auto const& b = new_turns[i];
        BOOST_TEST_EQ(a.method, b.method);
        BOOST_TEST_EQ(a.touch_only, b.touch_only);
        BOOST_TEST(ex::equal_points(a.point, b.point));
        for (std::size_t j = 0; j < 2; ++j)
        {
            BOOST_TEST_EQ(a.operations[j].operation, b.operations[j].operation);
            BOOST_TEST(a.operations[j].seg_id == b.operations[j].seg_id);
        }
    }
}

void test_turns()
{
    // Shared endpoints, spikes, zero-length segments and both exact and nearly
    // collinear continuations. Tiny determinants may underflow native double.
    for (int exponent : {0, -52, -600})
    {
        double const y = std::ldexp(1.0, exponent);
        std::array<point, 7> const points{{point(-1, 0), point(0, 0), point(1, 0),
            point(-1, y), point(1, y), point(-1, -y), point(1, -y)}};
        for (auto const& pi : points)
        for (auto const& pk : points)
        for (auto const& qi : points)
        for (auto const& qk : points)
        {
            compare_turns(range{{{pi, point(0, 0), pk}}}, range{{{qi, point(0, 0), qk}}});
        }
        // Touch interiors and unequal collinear overlap, not just shared ends.
        compare_turns(range{{{point(-2, 0), point(2, 0), point(3, y)}}},
                      range{{{point(0, y), point(0, 0), point(1, 0)}}});
        compare_turns(range{{{point(0, 0), point(2, 0), point(3, y)}}},
                      range{{{point(1, 0), point(3, 0), point(4, -y)}}});
    }
}

std::set<std::set<bg::signed_size_type>> memberships(clusters const& input)
{
    std::set<std::set<bg::signed_size_type>> result;
    for (auto const& pair : input) { result.insert(pair.second.turn_indices); }
    return result;
}

void compare_clusters(std::vector<turn> turns)
{
    auto old_turns = turns;
    clusters old_clusters, new_clusters;
    ov::get_clusters(old_turns, old_clusters, algorithm_strategy<false>());
    ov::get_clusters(turns, new_clusters, algorithm_strategy<true>());
    BOOST_TEST(memberships(old_clusters) == memberships(new_clusters));
    ov::colocate_clusters(new_clusters, turns, algorithm_strategy<true>());
    for (auto const& pair : new_clusters)
    {
        auto const& indices = pair.second.turn_indices;
        point const& canonical = turns[*indices.begin()].point;
        for (auto index : indices)
        {
            BOOST_TEST(turns[index].point.shares_node(canonical));
        }
    }
}

void test_clusters()
{
    ex::point_cache_scope operation;
    compare_clusters({});
    compare_clusters(std::vector<turn>(1));
    std::vector<turn> turns;
    for (int i = 0; i < 256; ++i)
    {
        turn item;
        item.point = point(i, 0);
        turns.push_back(item);
    }
    counters old_counts, new_counts;
    algorithm_strategy<false> old_strategy;
    algorithm_strategy<true> new_strategy;
    old_strategy.counts = &old_counts;
    new_strategy.counts = &new_counts;
    clusters old_clusters, new_clusters;
    ov::get_clusters(turns, old_clusters, old_strategy);
    ov::get_clusters(turns, new_clusters, new_strategy);
    BOOST_TEST(old_clusters.empty());
    BOOST_TEST(new_clusters.empty());
    BOOST_TEST_EQ(old_counts.equalities, 256u * 255u / 2u);
    BOOST_TEST_EQ(new_counts.equalities, 255u);

    for (int i = 0; i < 256; ++i)
    {
        turn item;
        item.point = point(i, 0); // Equal values with independent nodes.
        item.discarded = i % 3 == 0;
        turns.push_back(item);
    }
    compare_clusters(turns);
    for (auto& item : turns) { item.discarded = true; }
    compare_clusters(turns);

    // Distinct exact intersections can round to one double coordinate. They
    // must not share a cluster, while independent constructions of 1/3 must.
    turns.clear();
    ex::exact_number const third = ex::exact_number(1) / 3;
    ex::exact_number const delta = ex::exact_number(1)
        / ex::exact_number(boost::multiprecision::cpp_int(1) << 200);
    std::array<point, 4> const points{{
        ex::construct_intersection(point(0, 0), point(1, 1), point(0, 1), point(1, -1)),
        point::literal(third, third), point::literal(third + delta, third),
        point::literal(third, third + delta)}};
    for (auto const& p : points)
    {
        turn item;
        item.point = p;
        turns.push_back(item);
    }
    new_clusters.clear();
    ov::get_clusters(turns, new_clusters, algorithm_strategy<true>());
    BOOST_TEST_EQ(new_clusters.size(), 1u);
    if (!new_clusters.empty())
    {
        BOOST_TEST(new_clusters.begin()->second.turn_indices == std::set<bg::signed_size_type>({0, 1}));
    }
    compare_clusters(turns);
}

using native_point = bg::model::d2::point_xy<double>;
using polygon = bg::model::polygon<native_point>;
using exact_multi = bg::model::multi_polygon<bg::model::polygon<point>>;

template <bool Enabled>
exact_multi overlay(polygon const& a, polygon const& b, int operation)
{
    exact_multi result;
    algorithm_strategy<Enabled> strategy;
    switch (operation)
    {
        case 0: bg::union_(a, b, result, strategy); break;
        case 1: bg::intersection(a, b, result, strategy); break;
        case 2: bg::difference(a, b, result, strategy); break;
        case 3: bg::sym_difference(a, b, result, strategy); break;
    }
    return result;
}

void compare_overlays(polygon const& a, polygon const& b)
{
    for (int operation = 0; operation < 4; ++operation)
    {
        auto const old_result = overlay<false>(a, b, operation);
        auto const new_result = overlay<true>(a, b, operation);
        BOOST_TEST_EQ(new_result.size(), old_result.size());
        BOOST_TEST_EQ(bg::num_interior_rings(new_result), bg::num_interior_rings(old_result));
        BOOST_TEST(bg::area(new_result, bg::strategies::relate::cartesian<ex::exact_number>())
            == bg::area(old_result, bg::strategies::relate::cartesian<ex::exact_number>()));
        if (!new_result.empty() && !old_result.empty())
        {
            BOOST_TEST(bg::relate(new_result, old_result, bg::de9im::mask("T*F**FFF*"),
                                  algorithm_strategy<false>()));
        }
    }
}

void test_overlays()
{
    ex::point_cache_scope<double> operation_scope;
    polygon a, b;
    bg::exterior_ring(a) = {{0, 0}, {0, 2}, {2, 2}, {2, 0}, {0, 0}};
    // Shared edge, corner contact, exact identity, and a sub-ulp-scale gap.
    for (double x : {0.0, 1.0, 2.0, std::nextafter(2.0, 3.0)})
    {
        bg::exterior_ring(b) = {{x, 0}, {x, 2}, {4, 2}, {4, 0}, {x, 0}};
        compare_overlays(a, b);
    }
    compare_overlays(a, a);
    bg::exterior_ring(b) = {{2, 2}, {2, 3}, {3, 3}, {3, 2}, {2, 2}};
    compare_overlays(a, b);

    double const tiny = std::ldexp(1.0, -200);
    bg::exterior_ring(a) = {{0, 0}, {0, tiny}, {1, 0}, {0, 0}};
    bg::exterior_ring(b) = {{0, 0}, {0, tiny}, {2, 0}, {0, 0}};
    compare_overlays(a, b);

    // A genuine 2^-107 hole whose vertices round together in native output.
    bg::exterior_ring(a) = {{0, 0}, {0, 3}, {3, 3}, {3, 1.25},
        {1, 1}, {3, 0.75}, {3, 0}, {0, 0}};
    double const left = std::nextafter(1.0, 2.0);
    bg::exterior_ring(b) = {{left, 0}, {left, 3}, {4, 3}, {4, 0}, {left, 0}};
    compare_overlays(a, b);
    auto const united = overlay<true>(a, b, 0);
    BOOST_TEST_EQ(united.size(), 1u);
    BOOST_TEST_EQ(bg::num_interior_rings(united), 1u);
    ex::exact_number const hole_area = ex::exact_number(1)
        / ex::exact_number(boost::multiprecision::cpp_int(1) << 107);
    BOOST_TEST(bg::area(united, bg::strategies::relate::cartesian<ex::exact_number>())
        == ex::exact_number(12) - hole_area);
}

void test_independent_xor_parts()
{
    using p = bg::model::d2::point_xy<double>;
    using poly = bg::model::polygon<p>;
    auto read = [](char const* w) { poly v; bg::read_wkt(w, v); return v; };
    auto a = read("POLYGON((0 0,0 4,4 4,4 0,0 0))");
    auto check = [&](poly const& b, bool independent, double area) {
        exact_strategy const s;
        bg::strategies::detail::with_operation_sources(s, a, b, [&](auto const& bound) {
            BOOST_TEST_EQ(bound.no_common_boundary_edge(a,b), independent);
        });
        bg::model::multi_polygon<poly> result;
        bg::sym_difference(a,b,result,s);
        BOOST_TEST(bg::is_valid(result)); BOOST_TEST_EQ(bg::area(result), area);
    };
    check(read("POLYGON((2 -1,2 3,5 3,5 -1,2 -1))"), true, 16);
    check(read("POLYGON((4 4,4 5,5 5,5 4,4 4))"), true, 17);
    check(read("POLYGON((4 1,4 2,4 3,5 3,5 1,4 1))"), false, 18);
    check(a, false, 0);
    a = read("POLYGON((0 0,0 6,6 6,6 0,0 0),(2 2,4 2,4 4,2 4,2 2))");
    check(read("POLYGON((2 2,2 4,4 4,4 2,2 2))"), false, 36);
}

int main()
{
    test_verified_side();
    test_turns();
    test_clusters();
    test_overlays();
    test_independent_xor_parts();
    return boost::report_errors();
}
