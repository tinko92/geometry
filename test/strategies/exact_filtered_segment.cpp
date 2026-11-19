// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.

#include <boost/core/lightweight_test.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/strategy/cartesian/side_by_triangle.hpp>
#include <boost/geometry/strategy/cartesian/side_robust.hpp>
#include <boost/geometry/strategies/cartesian/intersection.hpp>
#include <boost/geometry/strategies/cartesian/exact/segment.hpp>
#include <array>
#include <cfenv>
#include <cstdint>

namespace bg = boost::geometry;
namespace ec = bg::strategies::detail::exact_cartesian;
using point = ec::exact_point<double>;
using rational_point = bg::model::d2::point_xy<ec::exact_number>;

template <typename P>
struct range
{
    using point_type = P;
    std::array<P, 2> points;
    P const& at(std::size_t i) const { return points[i]; }
};

void check(range<point> const& a, range<point> const& b)
{
    auto lift = [](range<point> const& input)
    {
        range<rational_point> output;
        for (unsigned i = 0; i < 2; ++i)
            output.points[i] = rational_point(input.at(i).exact_coordinate(0), input.at(i).exact_coordinate(1));
        return output;
    };
    using result = ec::segment_intersection_result<point>;
    using reference = bg::segment_intersection_points<rational_point>;
    auto const ar = lift(a), br = lift(b);
    auto const actual = ec::filtered_collinear_intersection<double>::apply(a, b,
        bg::policies::relate::segments_intersection_policy<result>());
    auto const expected = bg::strategy::intersection::cartesian_segments<ec::exact_number>::apply(ar, br,
        bg::policies::relate::segments_intersection_policy<reference>());
    BOOST_TEST_EQ(actual.intersection_points.count, expected.intersection_points.count);
    for (std::size_t i = 0; i < (std::min)(actual.intersection_points.count, expected.intersection_points.count); ++i)
    {
        BOOST_TEST(actual.intersection_points.intersections[i].exact_coordinate(0) == bg::get<0>(expected.intersection_points.intersections[i]));
        BOOST_TEST(actual.intersection_points.intersections[i].exact_coordinate(1) == bg::get<1>(expected.intersection_points.intersections[i]));
    }
    BOOST_TEST_EQ(actual.direction.how, expected.direction.how);
    BOOST_TEST_EQ(actual.direction.opposite, expected.direction.opposite);
    BOOST_TEST_EQ(actual.direction.arrival[0], expected.direction.arrival[0]);
    BOOST_TEST_EQ(actual.direction.arrival[1], expected.direction.arrival[1]);
}

void test_axis(bool vertical)
{
    auto p = [vertical](double x) { return vertical ? point(0, x) : point(x, 0); };
    point const third = vertical
        ? ec::construct_intersection(point(0, 0), point(0, 1), point(-1, 0), point(2, 1))
        : ec::construct_intersection(point(0, 0), point(1, 0), point(0, -1), point(1, 2));
    std::array<point, 6> const points{{p(-1), p(0), third, p(0.5), p(1), p(2)}};
    for (auto const& a : points) for (auto const& b : points)
    for (auto const& c : points) for (auto const& d : points)
        check({{{a, b}}}, {{{c, d}}});
}

int main()
{
    int const original = std::fegetround();
    for (int mode : {FE_TONEAREST})
    {
        if (std::fesetround(mode) != 0) continue;
        test_axis(false);
        test_axis(true);
    }
    BOOST_TEST_EQ(std::fesetround(original), 0);
    auto retained = []
    {
        point const third = ec::construct_intersection(point(0, 0), point(1, 0), point(0, -1), point(1, 2));
        using result = ec::segment_intersection_result<point>;
        return ec::filtered_collinear_intersection<double>::apply(range<point>{{{point(0, 0), point(1, 0)}}},
            range<point>{{{third, point(0.5, 0)}}}, bg::policies::relate::segments_intersection_points<result>());
    }();
    BOOST_TEST_EQ(retained.count, 2u);
    BOOST_TEST(retained.intersections[0].exact_coordinate(0) == ec::exact_number(1) / 3);
    BOOST_TEST(retained.intersections[1].exact_coordinate(0) == ec::exact_number(1) / 2);
    return boost::report_errors();
}
