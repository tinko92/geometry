// Boost.Geometry (aka GGL, Generic Geometry Library)
// Unit Test

// Copyright (c) 2018 Barend Gehrels, Amsterdam, the Netherlands.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)


#if defined(_MSC_VER)
// We deliberately mix float/double's here so turn off warning
#pragma warning( disable : 4244 )
#endif // defined(_MSC_VER)

#include <geometry_test_common.hpp>

#include <boost/geometry/algorithms/assign.hpp>

#include <boost/geometry/strategies/cartesian/intersection.hpp>
#include <boost/geometry/strategies/cartesian/general_intersection.hpp>
#include <boost/geometry/strategies/geographic/intersection.hpp>
#include <boost/geometry/strategies/intersection_result.hpp>

#include <boost/geometry/policies/relate/intersection_points.hpp>
#include <boost/geometry/policies/relate/direction.hpp>
#include <boost/geometry/policies/relate/tupled.hpp>

#include <boost/geometry/algorithms/intersection.hpp>


#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/segment.hpp>
#include <boost/geometry/geometries/adapted/boost_tuple.hpp>

BOOST_GEOMETRY_REGISTER_BOOST_TUPLE_CS(cs::cartesian);


template <typename R>
inline double get_fraction(R const& ratio)
{
    return ratio.numerator() / double(ratio.denominator());
}


template <typename Point>
static void test_segment_intersection(std::string const& caseid,
                double x1, double y1, double x2, double y2,
                double x3, double y3, double x4, double y4,
                char expected_how,
                double expected_x1 = -99, double expected_y1 = -99)
{
//    std::cout << caseid << std::endl;

    typedef typename bg::coordinate_type<Point>::type coordinate_type;
    typedef bg::model::referring_segment<const Point> segment_type;

    bool const is_cart = boost::is_same
        <
            typename bg::coordinate_system<Point>::type,
            bg::cs::cartesian
        >::value;
    bool const is_geo = boost::is_same
        <
            typename bg::coordinate_system<Point>::type,
            bg::cs::geographic<bg::degree>
        >::value;


    Point p1, p2, p3, p4;
    bg::assign_values(p1, x1, y1);
    bg::assign_values(p2, x2, y2);
    bg::assign_values(p3, x3, y3);
    bg::assign_values(p4, x4, y4);

    segment_type s12(p1,p2);
    segment_type s34(p3,p4);

    std::size_t expected_count = 0;

    if (expected_x1 != -99 && expected_y1 != -99)
    {
        expected_count++;
    }

    // Using intersection_insert

    std::vector<Point> out;
    bg::detail::intersection::intersection_insert<Point>(s12, s34, std::back_inserter(out));

    // Using strategy
    typedef bg::detail::no_rescale_policy rescale_policy_type;
    rescale_policy_type rescale_policy;
    typedef typename bg::segment_ratio_type<Point, rescale_policy_type>::type ratio_type;
    typedef bg::segment_intersection_points
        <
            Point,
            ratio_type
        > result_type;

    typedef bg::policies::relate::segments_intersection_points
        <
            result_type
        > points_policy_type;

    result_type ip_cartesian, ip_geographic, ip_general_form;

    if (is_cart)
    {
        ip_cartesian = bg::strategy::intersection::cartesian_segments<>
                ::apply(s12, s34, points_policy_type(), rescale_policy, p1, p2, p3, p4);
    }

#if 1
    // Geographic
    if (is_geo)
    {
        bg::strategy::intersection::geographic_segments
                <bg::strategy::andoyer, 1> geographic_strategy;

        ip_geographic = geographic_strategy.apply(s12, s34,
                points_policy_type(), rescale_policy, p1, p2, p3, p4);
    }

    // Using general form
    typedef bg::strategy::intersection::cartesian_general_segments<> cgs;
    ip_general_form = cgs::apply(s12, s34, points_policy_type(),
                            rescale_policy, p1, p2, p3, p4);

#endif
//    BOOST_CHECK_EQUAL(boost::size(out), expected_count);

    if (expected_x1 != -99 && expected_y1 != -99)
    {
        if (boost::size(out) >= 1)
        {
            BOOST_CHECK_CLOSE(bg::get<0>(out[0]), coordinate_type(expected_x1), 0.001);
            BOOST_CHECK_CLOSE(bg::get<1>(out[0]), coordinate_type(expected_y1), 0.001);
        }

        if (is_cart)
        {
            BOOST_CHECK_EQUAL(ip_cartesian.count, expected_count);

            BOOST_CHECK_CLOSE(bg::get<0>(ip_cartesian.intersections[0]), expected_x1, 0.001);
            BOOST_CHECK_CLOSE(bg::get<1>(ip_cartesian.intersections[0]), expected_y1, 0.001);

#if 1
            BOOST_CHECK_CLOSE(bg::get<0>(ip_general_form.intersections[0]), expected_x1, 0.001);
            BOOST_CHECK_CLOSE(bg::get<1>(ip_general_form.intersections[0]), expected_y1, 0.001);
#endif
        }

#if 1
        if (is_geo)
        {
            BOOST_CHECK_CLOSE(bg::get<0>(ip_geographic.intersections[0]), expected_x1, 0.001);
            BOOST_CHECK_CLOSE(bg::get<1>(ip_geographic.intersections[0]), expected_y1, 0.001);
        }
#endif
    }
}



int test_main(int, char* [])
{
    // Test the available segment intersection strategies

    typedef bg::model::point<double, 2, bg::cs::cartesian> cartesian_point;
    typedef bg::model::point<double, 2, bg::cs::spherical_equatorial<bg::degree> > spherical_point;
    typedef bg::model::point<double, 2, bg::cs::geographic<bg::degree> > geographic_point;

//    test_segment_intersection<cartesian_point>("cartesian",
//                                               0, 2, 2, 0,
//                                               0, 0, 2, 2,
//                                               'i',
//                                               1, 1);
//    test_segment_intersection<cartesian_point>("cartesian",
//                                               0, 1, 2, 5,
//                                               0, 3, 4, 5,
//                                               'i',
//                                               1.33333333, 3.6666667);

    test_segment_intersection<cartesian_point>("cartesian",
                                               3, 4, 1, 2,
                                               4, 3, 2, 1,
                                               'd',
                                               -99.0, -99.0);

//    test_segment_intersection<spherical_point>("spherical",
//                                               10, -1, 20, 1,
//                                               12.5, -1, 12.5, 1,
//                                               'i',
//                                               12.5 -0.50051443471392);

//    test_segment_intersection<geographic_point>("geographic",
//                                               10, -1, 20, 1,
//                                               12.5, -1, 12.5, 1,
//                                               'i',
//                                               12.5 -0.50051443471392);

    return 0;
}
