// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Licensed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <initializer_list>
#include <map>
#include <type_traits>
#include <utility>
#include <vector>

#include <boost/core/lightweight_test.hpp>
#include <boost/geometry/algorithms/detail/overlay/get_clusters.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/strategies/relate/cartesian_exact.hpp>

namespace bg = boost::geometry;
namespace ec = bg::strategies::detail::exact_cartesian;

template <typename Point>
struct counted_point
{
    Point value;
    inline static std::size_t copies = 0;
    inline static std::size_t moves = 0;

    explicit counted_point(Point const& point) : value(point) {}
    counted_point(counted_point const& other) : value(other.value) { ++copies; }
    counted_point(counted_point&& other)
        noexcept(std::is_nothrow_move_constructible_v<Point>)
        : value(std::move(other.value)) { ++moves; }

    counted_point& operator=(counted_point const& other)
    {
        value = other.value;
        ++copies;
        return *this;
    }

    counted_point& operator=(counted_point&& other)
        noexcept(std::is_nothrow_move_assignable_v<Point>)
    {
        value = std::move(other.value);
        ++moves;
        return *this;
    }
};

template <typename Point>
struct test_turn
{
    using point_type = counted_point<Point>;
    point_type point;
    bool discarded;

    explicit test_turn(Point const& point, bool discarded = false)
        : point(point), discarded(discarded) {}
};

// Count builder copies while keeping the production ordering and equality.
template <typename Kernel>
struct counted_strategy
{
    static constexpr bool exact_predicates = true;

    struct exact_cluster_compare
    {
        template <typename Point>
        static bool apply(counted_point<Point> const& a, counted_point<Point> const& b)
        {
            using compare = typename Kernel::exact_cluster_compare;
            return compare::apply(a.value, b.value);
        }
    };

    struct topology_type
    {
        template <typename Point>
        static bool cluster_equal(counted_point<Point> const& a, counted_point<Point> const& b)
        {
            return Kernel::topology().cluster_equal(a.value, b.value);
        }

        template <typename Point>
        static int compare_y(counted_point<Point> const& a, counted_point<Point> const& b)
        {
            return Kernel::topology().compare_y(a.value, b.value);
        }

        template <typename Point>
        static bool cluster_window_exceeded(counted_point<Point> const& a,
                                            counted_point<Point> const& b)
        {
            return Kernel::topology().cluster_window_exceeded(a.value, b.value);
        }
    };

    static topology_type topology() { return {}; }
};

using index_type = bg::signed_size_type;
using expected_clusters = std::initializer_list<std::initializer_list<index_type>>;

template <typename Coordinate, typename Turns>
void check_clusters(Turns& turns, expected_clusters expected)
{
    using point_type = typename Turns::value_type::point_type;
    std::map<index_type, bg::detail::overlay::cluster_info> clusters;
    counted_strategy<bg::strategies::relate::cartesian_exact<Coordinate>> strategy;

    // Fixture construction is outside the measured builder call.
    point_type::copies = point_type::moves = 0;
    bg::detail::overlay::get_clusters(turns, clusters, strategy);
    BOOST_TEST_EQ(point_type::copies, 0u);
    BOOST_TEST_EQ(point_type::moves, 0u);

    BOOST_TEST_EQ(clusters.size(), expected.size());
    index_type id = 1;
    for (auto const& indices : expected)
    {
        auto const found = clusters.find(id++);
        BOOST_TEST(found != clusters.end());
        if (found != clusters.end())
        {
            auto const& actual = found->second.turn_indices;
            BOOST_TEST(std::equal(actual.begin(), actual.end(), indices.begin(), indices.end()));
            BOOST_TEST_EQ(found->second.open_count, 0u);
        }
    }
}

template <typename Coordinate, typename Point>
void check_containers(std::vector<test_turn<Point>>& turns, expected_clusters expected)
{
    std::deque<test_turn<Point>> deque(turns.begin(), turns.end());
    check_clusters<Coordinate>(turns, expected);
    check_clusters<Coordinate>(deque, expected);
}

template <typename Coordinate>
void test_distinct_producers()
{
    ec::point_cache_scope<Coordinate> scope;
    using point = ec::exact_point<Coordinate>;
    std::array<point, 8> const sources{{
        point(0, 0), point(1, 1), point(0, 1), point(1, -1),
        point(0, 0), point(2, 2), point(0, 2), point(1, -3)
    }};
    auto const first = ec::construct_intersection(sources[0], sources[1], sources[2], sources[3]);
    auto const second = ec::construct_intersection(sources[4], sources[5], sources[6], sources[7]);
    BOOST_TEST(!first.shares_node(second));

    // Both constructions are (1/3, 1/3); input order deliberately differs from
    // cluster order. Discarded matches and singletons must not consume IDs.
    std::vector<test_turn<point>> turns;
    turns.emplace_back(second);
    turns.emplace_back(point(2, 3));
    turns.emplace_back(point(-10, -10), true);
    turns.emplace_back(point(-2, 1));
    turns.emplace_back(first);
    turns.emplace_back(point(-2, 1));
    turns.emplace_back(second, true);
    turns.emplace_back(point(4, 0));
    turns.emplace_back(point(4, 0));
    turns.emplace_back(point(-2, 1));
    turns.emplace_back(point(-10, -10), true);
    check_containers<Coordinate>(turns, {{3, 5, 9}, {0, 4}, {7, 8}});
}

void test_distinct_exact_points_with_equal_approximations()
{
    ec::point_cache_scope<double> scope;
    using point = ec::exact_point<double>;
    double const next = std::nextafter(1.0, 2.0);
    std::array<point, 4> const sources{{point(1, 0), point(next, 1), point(0, 0.25), point(2, 0.25)}};
    auto const crossing = ec::construct_intersection(sources[0], sources[1], sources[2], sources[3]);
    point const native(1, 0.25);
    // Exact crossing x = 1 + 2^-54, but its double approximation is 1.
    BOOST_TEST_EQ(crossing.approximation(0).approximate, native.approximation(0).approximate);
    std::vector<test_turn<point>> turns;
    turns.emplace_back(crossing);
    turns.emplace_back(native);
    turns.emplace_back(point(next, 0.25));
    turns.emplace_back(native);
    turns.emplace_back(crossing);
    check_containers<double>(turns, {{1, 3}, {0, 4}});
}

void test_native_points()
{
    using point = bg::model::d2::point_xy<double>;
    double const next = std::nextafter(1.0, 2.0);
    std::vector<test_turn<point>> turns;
    turns.emplace_back(point(0.0, -0.0));
    turns.emplace_back(point(next, 2));
    turns.emplace_back(point(-0.0, 0.0));
    turns.emplace_back(point(1, 2));
    turns.emplace_back(point(0.0, 0.0));
    turns.emplace_back(point(next, 2));
    turns.emplace_back(point(1, 2));
    turns.emplace_back(point(-2, 3));
    turns.emplace_back(point(-0.0, -0.0), true);
    check_containers<double>(turns, {{0, 2, 4}, {3, 6}, {1, 5}});
    BOOST_TEST(std::signbit(bg::get<1>(turns[0].point.value)));
    BOOST_TEST(std::signbit(bg::get<0>(turns[2].point.value)));
}

void test_no_clusters()
{
    using point = bg::model::d2::point_xy<double>;
    std::vector<test_turn<point>> turns;
    check_containers<double>(turns, {});
    turns.emplace_back(point(0, 0));
    check_containers<double>(turns, {});
    turns.emplace_back(point(0, 0), true);
    turns.emplace_back(point(1, 0));
    check_containers<double>(turns, {});
    for (auto& turn : turns) turn.discarded = true;
    check_containers<double>(turns, {});
}

int main()
{
    test_distinct_producers<double>();
    test_distinct_producers<std::int_least64_t>();
    test_distinct_exact_points_with_equal_approximations();
    test_native_points();
    test_no_clusters();
    return boost::report_errors();
}
