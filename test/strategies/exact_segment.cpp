// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Licensed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt

#define BOOST_GEOMETRY_EXACT_FILTER_STATISTICS
#include <boost/core/lightweight_test.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/policies/relate/direction.hpp>
#include <boost/geometry/policies/relate/intersection_points.hpp>
#include <boost/geometry/strategies/cartesian/exact/segment.hpp>
#include <boost/geometry/strategies/intersection_result.hpp>
#include <boost/geometry/strategy/cartesian/side_by_triangle.hpp>
#include <boost/geometry/strategy/cartesian/side_robust.hpp>
#include <boost/geometry/strategies/cartesian/intersection.hpp>
#include <boost/geometry/algorithms/detail/overlay/turn_info.hpp>

#include <algorithm>
#include <array>
#include <cfenv>
#include <cmath>
#include <cstdint>
#include <limits>
#include <random>

namespace bg = boost::geometry;
namespace ec = bg::strategies::detail::exact_cartesian;

template <typename Point>
struct segment_range
{
    using point_type = Point;
    std::array<Point, 2> points;
    Point const& at(std::size_t index) const { return points.at(index); }
};

template <typename Coordinate>
void check_pair(segment_range<bg::model::d2::point_xy<Coordinate>> const& a,
                segment_range<bg::model::d2::point_xy<Coordinate>> const& b)
{
    using rational_point = bg::model::d2::point_xy<ec::exact_number>;
    auto lift = [](auto const& range)
    {
        segment_range<rational_point> result;
        for (std::size_t i = 0; i != 2; ++i)
        {
            result.points[i] = rational_point(ec::exact_value(bg::get<0>(range.at(i))),
                                              ec::exact_value(bg::get<1>(range.at(i))));
        }
        return result;
    };
    auto const ra = lift(a), rb = lift(b);
    using exact_result = ec::segment_intersection_result<ec::exact_point<Coordinate>>;
    using reference_result = bg::segment_intersection_points<rational_point>;
    auto const actual = ec::segment_intersection<Coordinate>::apply(a, b,
        bg::policies::relate::segments_intersection_points<exact_result>());
    auto const expected = bg::strategy::intersection::cartesian_segments<ec::exact_number>::apply(ra, rb,
        bg::policies::relate::segments_intersection_points<reference_result>());
    BOOST_TEST_EQ(actual.count, expected.count);
    for (std::size_t i = 0; i < (std::min)(actual.count, expected.count); ++i)
    {
        BOOST_TEST_EQ(actual.intersections[i].exact_coordinate(0), bg::get<0>(expected.intersections[i]));
        BOOST_TEST_EQ(actual.intersections[i].exact_coordinate(1), bg::get<1>(expected.intersections[i]));

    }
    auto const direction = ec::segment_intersection<Coordinate>::apply(a, b,
        bg::policies::relate::segments_direction());
    auto const reference_direction = bg::strategy::intersection::cartesian_segments<ec::exact_number>
        ::apply(ra, rb, bg::policies::relate::segments_direction());
    BOOST_TEST_EQ(direction.how, reference_direction.how);
    BOOST_TEST_EQ(direction.opposite, reference_direction.opposite);
    for (std::size_t i = 0; i != 2; ++i)
    {
        BOOST_TEST_EQ(direction.arrival[i], reference_direction.arrival[i]);
    }
}

template <typename Coordinate>
void test_small_grid()
{
    using point = bg::model::d2::point_xy<Coordinate>;
    std::array<point, 9> points;
    for (int i = 0; i != 9; ++i)
    {
        points[i] = point(Coordinate(i % 3 - 1), Coordinate(i / 3 - 1));
    }
    // Includes endpoint touches, overlaps, reversals and zero-length segments.
    for (auto const& a : points)
    for (auto const& b : points)
    for (auto const& c : points)
    for (auto const& d : points)
    {
        check_pair<Coordinate>({{{a, b}}}, {{{c, d}}});
    }
}

void test_binary_inputs()
{
    using point = bg::model::d2::point_xy<double>;
    std::mt19937_64 random(0x5345474d454e54);
    auto coordinate = [&]
    {
        return std::ldexp(double(std::int64_t(random() % 2000001) - 1000000),
                          int(random() % 361) - 180);
    };
    for (int i = 0; i != 500; ++i)
    {
        check_pair<double>({{{point(coordinate(), coordinate()), point(coordinate(), coordinate())}}},
                           {{{point(coordinate(), coordinate()), point(coordinate(), coordinate())}}});
    }
}

void test_point_positions()
{
    using point = ec::exact_point<double>;
    using position = ec::segment_position<double>;
    static_assert(std::is_empty<ec::segment_position_tag>::value, "No exact ratio representation");
    using operation = bg::detail::overlay::turn_operation<point, ec::segment_position_tag>;
    static_assert(!bg::detail::overlay::stores_fraction<ec::segment_position_tag>::value,
                  "Exact turn operations have no fraction storage");
    (void)sizeof(operation);
    for (bool reverse : {false, true})
    {
        point const a = reverse ? point(3, 0) : point(0, 0);
        point const b = reverse ? point(0, 0) : point(3, 0);
        for (int i = -1; i <= 4; ++i)
        {
            point const p(i, 0);
            ec::exact_number const t = reverse ? (ec::exact_number(3) - i) / 3 : ec::exact_number(i) / 3;
            int const expected = t < 0 ? 0 : t == 0 ? 1 : t < 1 ? 2 : t == 1 ? 3 : 4;
            BOOST_TEST_EQ(position::classify(a, b, p), expected);
            for (int j = -1; j <= 4; ++j)
                BOOST_TEST_EQ(position::apply(a, b, p, point(j, 0)), reverse ? (j > i) - (j < i) : (i > j) - (i < j));
        }
        BOOST_TEST_EQ(position::classify(a, b, point(1, 1)), 7);
    }
    BOOST_TEST_EQ(position::classify(point(0, 0), point(0, 0), point(0, 0)), 5);
    BOOST_TEST_EQ(position::classify(point(0, 0), point(0, 0), point(0, 1)), 6);
    point const third = ec::construct_intersection(point(0, 0), point(1, 1), point(0, 1), point(1, -1));
    BOOST_TEST_EQ(position::classify(point(0, 0), point(1, 1), third), 2);
    BOOST_TEST_EQ(position::apply(point(0, 0), point(1, 1), third, point(0.5, 0.5)), -1);
    BOOST_TEST_EQ(position::apply(point(0, 3), point(0, 0), point(0, 2), point(0, 1)), -1);
    using int_point = ec::exact_point<std::int64_t>;
    auto const low = (std::numeric_limits<std::int64_t>::min)();
    auto const high = (std::numeric_limits<std::int64_t>::max)();
    BOOST_TEST_EQ(ec::segment_position<std::int64_t>::classify(int_point(low, low), int_point(high, high),
                                                            int_point(0, 0)), 2);
}

void test_rational_inputs()
{
    ec::point_cache_scope operation;
    using q = ec::exact_number;
    using p = bg::model::d2::point_xy<q>;
    using result = ec::segment_intersection_result<ec::exact_point<double>>;
    segment_range<p> const a{{{p(q(1) / 7, 0), p(1, 0)}}};
    segment_range<p> const crossing{{{p(q(1) / 3, -1), p(q(1) / 3, 1)}}};
    auto const intersection = ec::segment_intersection<double>::apply(a, crossing,
        bg::policies::relate::segments_intersection_points<result>());
    BOOST_TEST_EQ(intersection.count, 1u);
    if (intersection.count == 1)
        BOOST_TEST(intersection.intersections[0].exact_coordinate(0) == q(1) / 3);
    segment_range<p> const overlap{{{p(q(1) / 3, 0), p(q(1) / 2, 0)}}};
    auto const collinear = ec::segment_intersection<double>::apply(a, overlap,
        bg::policies::relate::segments_intersection_points<result>());
    BOOST_TEST_EQ(collinear.count, 2u);
    if (collinear.count == 2)
    {
        BOOST_TEST(collinear.intersections[0].exact_coordinate(0) == q(1) / 3);
        BOOST_TEST(collinear.intersections[1].exact_coordinate(0) == q(1) / 2);
    }
}

int main()
{
    test_point_positions();
    test_rational_inputs();
    test_small_grid<std::int64_t>();
    test_small_grid<double>();
    test_binary_inputs();
    return boost::report_errors();
}
