// Boost.Geometry (aka GGL, Generic Geometry Library)
// Unit Test

// Copyright (c) 2007-2012 Barend Gehrels, Amsterdam, the Netherlands.
// Copyright (c) 2008-2012 Bruno Lalande, Paris, France.
// Copyright (c) 2009-2012 Mateusz Loskot, London, UK.

// Parts of Boost.Geometry are redesigned from Geodan's Geographic Library
// (geolib/GGL), copyright (c) 1995-2010 Geodan, Amsterdam, the Netherlands.

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

template <typename Collection>
void check_intersection(int caseid,
                Collection const& points,
                int expected_x1 = -99, int expected_y1 = -99,
                int expected_x2 = -99, int expected_y2 = -99)
{
    double const eps = 1.0e-7;
    bool const first_matches =
            std::fabs(bg::get<0>(points[0]) - expected_x1) < eps
            && std::fabs(bg::get<1>(points[0]) - expected_y1) < eps;

    int const index1 = first_matches ? 0 : 1;
    int const index2 = 1 - index1;

    BOOST_CHECK_CLOSE(bg::get<0>(points[index1]), expected_x1, 0.001);
    BOOST_CHECK_CLOSE(bg::get<1>(points[index1]), expected_y1, 0.001);
    BOOST_CHECK_CLOSE(bg::get<0>(points[index2]), expected_x2, 0.001);
    BOOST_CHECK_CLOSE(bg::get<1>(points[index2]), expected_y2, 0.001);
}


template <typename P, typename T>
static void test_segment_intersection(int caseid,
                T x1, T y1, T x2, T y2,
                T x3, T y3, T x4, T y4,
                char expected_how,
                T expected_x1 = -99, T expected_y1 = -99,
                T expected_x2 = -99, T expected_y2 = -99)
{
//    std::cout << caseid << std::endl;

    typedef typename bg::coordinate_type<P>::type coordinate_type;
    typedef bg::model::referring_segment<const P> segment_type;

    P p1, p2, p3, p4;
    bg::assign_values(p1, x1, y1);
    bg::assign_values(p2, x2, y2);
    bg::assign_values(p3, x3, y3);
    bg::assign_values(p4, x4, y4);

    segment_type s12(p1,p2);
    segment_type s34(p3,p4);

    std::size_t expected_count = 0;

    if (expected_x1 > -98 && expected_y1 > -98)
    {
        expected_count++;
    }
    if (expected_x2 > -98 && expected_y2 > -98)
    {
        expected_count++;
    }

    // Using intersection_insert

    std::vector<P> out;
    bg::detail::intersection::intersection_insert<P>(s12, s34, std::back_inserter(out));

    // Using strategy
    typedef bg::detail::no_rescale_policy rescale_policy_type;
    rescale_policy_type rescale_policy;
    typedef typename bg::segment_ratio_type<P, rescale_policy_type>::type ratio_type;
    typedef bg::segment_intersection_points
    <
        P,
        ratio_type
    > result_type;
    typedef bg::policies::relate::segments_intersection_points
        <
            result_type
        > points_policy_type;

    result_type ir_cartesian
        = bg::strategy::intersection::cartesian_segments<>
            ::apply(s12, s34, points_policy_type(), rescale_policy, p1, p2, p3, p4);

    bg::policies::relate::direction_type dir
        = bg::strategy::intersection::cartesian_segments<>
            ::apply(s12, s34, bg::policies::relate::segments_direction(),
                    rescale_policy, p1, p2, p3, p4);

    result_type ir_general;

    // Using general form

    {
        typedef bg::strategy::intersection::cartesian_general_segments<> cgs;

#if 0
        cgs::general_form<int> pu = cgs::construct_line<int>(x1, y1, x2, y2);
        cgs::general_form<int> qu = cgs::construct_line<int>(x3, y3, x4, y4);

        cgs::general_form<double> p = cgs::normalize_line(pu);
        cgs::general_form<double> q = cgs::normalize_line(qu);

        P ip;
        bool intersecting = cgs::get_intersection(ip, p, q);
        if (intersecting)
        {
            double f1, f2;
            bool const on1 = cgs::on_segment(ip, s12, f1, p);
            bool const on2 = cgs::on_segment(ip, s34, f2, q);

            // Verify if sides are as expected
            int const side1 = cgs::side(x1, y1, q);
            int const side2 = cgs::side(x2, y2, q);
            int const side3 = cgs::side(x3, y3, p);
            int const side4 = cgs::side(x4, y4, p);

            int const sd1 = cgs::side_strategy_type::apply(p3, p4, p1);
            int const sd2 = cgs::side_strategy_type::apply(p3, p4, p2);
            int const sd3 = cgs::side_strategy_type::apply(p1, p2, p3);
            int const sd4 = cgs::side_strategy_type::apply(p1, p2, p4);

            BOOST_CHECK_EQUAL(side1, sd1);
            BOOST_CHECK_EQUAL(side2, sd2);
            BOOST_CHECK_EQUAL(side3, sd3);
            BOOST_CHECK_EQUAL(side4, sd4);

            if (expected_count == 1)
            {
                BOOST_CHECK(on1);
                BOOST_CHECK(on2);
                BOOST_CHECK_CLOSE(bg::get<0>(ip), expected_x1, 0.001);
                BOOST_CHECK_CLOSE(bg::get<1>(ip), expected_y1, 0.001);

                BOOST_CHECK_CLOSE(f1, get_fraction(ir_cartesian.fractions[0].robust_ra), 0.001);
                BOOST_CHECK_CLOSE(f2, get_fraction(ir_cartesian.fractions[0].robust_rb), 0.001);
            }
            else
            {
                BOOST_CHECK_EQUAL(0, expected_count);
                BOOST_CHECK(on1 != on2);
            }
        }
        else
        {

        }

        if (expected_count == 2)
        {
            // Segments are parallel or collinear
            BOOST_CHECK(!intersecting);
        }
#endif

        ir_general = cgs::apply(s12, s34, points_policy_type(),
                       rescale_policy, p1, p2, p3, p4);
    }


    BOOST_CHECK_EQUAL(boost::size(out), expected_count);
    BOOST_CHECK_EQUAL(ir_cartesian.count, expected_count);
    BOOST_CHECK_MESSAGE(dir.how == expected_how,
            caseid
            << " how: detected: " << dir.how
            << " expected: "  << expected_how);

    if (expected_count == 2)
    {
        // Check if either first or second matches
        check_intersection(caseid, out, expected_x1, expected_y1, expected_x2, expected_y2);
        check_intersection(caseid, ir_cartesian.intersections, expected_x1, expected_y1, expected_x2, expected_y2);
        check_intersection(caseid, ir_general.intersections, expected_x1, expected_y1, expected_x2, expected_y2);
    }
    else if (expected_count == 1 && !out.empty())
    {
        BOOST_CHECK_CLOSE(bg::get<0>(out[0]), coordinate_type(expected_x1), 0.001);
        BOOST_CHECK_CLOSE(bg::get<1>(out[0]), coordinate_type(expected_y1), 0.001);

        BOOST_CHECK_CLOSE(bg::get<0>(ir_cartesian.intersections[0]), expected_x1, 0.001);
        BOOST_CHECK_CLOSE(bg::get<1>(ir_cartesian.intersections[0]), expected_y1, 0.001);

        BOOST_CHECK_CLOSE(bg::get<0>(ir_general.intersections[0]), expected_x1, 0.001);
        BOOST_CHECK_CLOSE(bg::get<1>(ir_general.intersections[0]), expected_y1, 0.001);
    }
}


template <typename P>
void test_all()
{
goto wrong;
    test_segment_intersection<P>( 1, 0,2, 2,0, 0,0, 2,2, 'i', 1, 1);
    test_segment_intersection<P>( 2, 2,2, 3,1, 0,0, 2,2, 'a', 2, 2);
    test_segment_intersection<P>( 3, 3,1, 2,2, 0,0, 2,2, 't', 2, 2);
    test_segment_intersection<P>( 4, 0,2, 1,1, 0,0, 2,2, 'm', 1, 1);

    test_segment_intersection<P>( 5, 1,1, 0,2, 0,0, 2,2, 's', 1, 1);
    test_segment_intersection<P>( 6, 0,2, 2,0, 0,0, 1,1, 'm', 1, 1);
    test_segment_intersection<P>( 7, 2,0, 0,2, 0,0, 1,1, 'm', 1, 1);
    test_segment_intersection<P>( 8, 2,3, 3,2, 0,0, 2,2, 'd');

    test_segment_intersection<P>( 9, 0,0, 2,2, 0,0, 2,2, 'e', 0, 0, 2, 2);
    test_segment_intersection<P>(10, 2,2, 0,0, 0,0, 2,2, 'e', 2, 2, 0, 0);
    test_segment_intersection<P>(11, 1,1, 3,3, 0,0, 2,2, 'c', 1, 1, 2, 2);
    test_segment_intersection<P>(12, 3,3, 1,1, 0,0, 2,2, 'c', 1, 1, 2, 2);

    test_segment_intersection<P>(13, 0,2, 2,2, 2,1, 2,3, 'm', 2, 2);
    test_segment_intersection<P>(14, 2,2, 2,4, 2,0, 2,2, 'a', 2, 2);
    test_segment_intersection<P>(15, 2,2, 2,4, 2,0, 2,1, 'd');
    test_segment_intersection<P>(16, 2,4, 2,2, 2,0, 2,1, 'd');

    test_segment_intersection<P>(17, 2,1, 2,3, 2,2, 2,4, 'c', 2, 2, 2, 3);
    test_segment_intersection<P>(18, 2,3, 2,1, 2,2, 2,4, 'c', 2, 3, 2, 2);
    test_segment_intersection<P>(19, 0,2, 2,2, 4,2, 2,2, 't', 2, 2);
    test_segment_intersection<P>(20, 0,2, 2,2, 2,2, 4,2, 'a', 2, 2);

    test_segment_intersection<P>(21, 1,2, 3,2, 2,1, 2,3, 'i', 2, 2);
    test_segment_intersection<P>(22, 2,4, 2,1, 2,1, 2,3, 'c', 2, 1, 2, 3);
    test_segment_intersection<P>(23, 2,4, 2,1, 2,3, 2,1, 'c', 2, 3, 2, 1);
    test_segment_intersection<P>(24, 1,1, 3,3, 0,0, 3,3, 'c', 1, 1, 3, 3);

    test_segment_intersection<P>(25, 2,0, 2,4, 2,1, 2,3, 'c', 2, 1, 2, 3);
    test_segment_intersection<P>(26, 2,0, 2,4, 2,3, 2,1, 'c', 2, 3, 2, 1);
    test_segment_intersection<P>(27, 0,0, 4,4, 1,1, 3,3, 'c', 1, 1, 3, 3);
    test_segment_intersection<P>(28, 0,0, 4,4, 3,3, 1,1, 'c', 3, 3, 1, 1);

    test_segment_intersection<P>(29, 1,1, 3,3, 0,0, 4,4, 'c', 1, 1, 3, 3);
    test_segment_intersection<P>(30, 0,0, 2,2, 2,2, 3,1, 'a', 2, 2);
    test_segment_intersection<P>(31, 0,0, 2,2, 2,2, 1,3, 'a', 2, 2);
    test_segment_intersection<P>(32, 0,0, 2,2, 1,1, 2,0, 's', 1, 1);

    test_segment_intersection<P>(33, 0,0, 2,2, 1,1, 0,2, 's', 1, 1);
    test_segment_intersection<P>(34, 2,2, 1,3, 0,0, 2,2, 'a', 2, 2);
    test_segment_intersection<P>(35, 2,2, 3,1, 0,0, 2,2, 'a', 2, 2);
    test_segment_intersection<P>(36, 0,0, 2,2, 0,2, 1,1, 'm', 1, 1);

    test_segment_intersection<P>(37, 0,0, 2,2, 2,0, 1,1, 'm', 1, 1);
    test_segment_intersection<P>(38, 1,1, 0,2, 0,0, 2,2, 's', 1, 1);
    test_segment_intersection<P>(39, 1,1, 2,0, 0,0, 2,2, 's', 1, 1);
    test_segment_intersection<P>(40, 2,0, 1,1, 0,0, 2,2, 'm', 1, 1);

    test_segment_intersection<P>(41, 1,2, 0,2, 2,2, 0,2, 'c', 1, 2, 0, 2);
    test_segment_intersection<P>(42, 2,1, 1,1, 2,2, 0,2, 'd');
    test_segment_intersection<P>(43, 4,1, 3,1, 2,2, 0,2, 'd');
    test_segment_intersection<P>(44, 4,2, 3,2, 2,2, 0,2, 'd');

    test_segment_intersection<P>(45, 2,0, 0,2, 0,0, 2,2, 'i', 1, 1);

    // In figure: times 2
    test_segment_intersection<P>(46, 8,2, 4,6, 0,0, 8, 8, 'i', 5, 5);

    // Specific cases for general form
//    test_segment_intersection<P>(47, 0.0, 1.0, 2.0, 5.0, 0.0, 3.0, 4.0, 5.0, 'i', 1.33333333, 3.6666667);
//    test_segment_intersection<P>(48, 3, 4, 1, 2, 4, 3, 2, 1, 'd');
//    test_segment_intersection<P>(49, 4.1, 2.5,5.3, 2.5, 1.5, 2.5,4.5, 2.5, 'c', 4.1, 2.5, 4.5, 2.5);

    test_segment_intersection<P>(50, 7,1, 8,2, 7,0, 7,1, 'a', 7, 1);
    test_segment_intersection<P>(51, -16.121616363525390625, -18.300218582153320312,
                                 -14.870415687561035156, -6.9987998008728027344,
                                 -22.508068084716796875, -27.9248046875,
                                 -14.870415687561035156, -6.9987998008728027344, 't',
                                 -14.870415687561035156, -6.9987998008728027344);
    test_segment_intersection<P>(52, 0.0, 3.3554432e-07,
                                 3.3554432e-07, 3.3554431e-07,
                                 2.5165824e-07, 2.5165824e-07,
                                 2.5165824e-07, 4.194304e-07, 't',
                                 2.5165824e-07, 3.3554432e-07);

    // They actually do not really intersect - but side value is still 0
    test_segment_intersection<P>(53,
                                 0.5489109414010371335, 0.57748351100509265343,
                                 0.40996112820544472477, 0.4644351568071597991,
                                 0.40996112820544511335, 0.46443515680716007665,
                                 0.39842498650182062159, 0.45263359648085582654,
                                 'a',
                                 0.40996112820544527988, 0.46443515680716024319);

    // Lines are nearly collinear (is_zero=1), but denominator!=0. All calculation
    // go wrong if considered as intersecting. But, seen normalized, they
    // are NOT considered as collinear. However, the two end-points (IP) are
    // on the other segment. This should be either considered as intersecting,
    // or as collinear (2 points)
    test_segment_intersection<P>(54,
                                 0.40996112820544472477, 0.4644351568071597991,  // =IP
                                 0.4294011278595494252, 0.48432242367292388519,
                                 0.4183516736935783964, 0.47301874838333629603,
                                 0.40996112820544511335, 0.46443515680716007665, // =IP
                                 't',
                                 0.40996112820544511335, 0.46443515680716007665);


wrong:
    // Lines are nearly collinear (is_zero=1), but denominator!=0. All calculation
    // go wrong if considered as intersecting. But, seen normalized, they
    // are NOT considered as collinear. However, the two end-points (IP) are
    // on the other segment. This should be either considered as intersecting,
    // or as collinear (2 points)
    test_segment_intersection<P>(55,
                                 50.97488895714744217, -30.277285722290763204,  // =IP
                                 37.1402130126953125, 1.3446992635726928711,
                                 57.524304765518266436 ,-37.807195733984784169,
                                 41.988733475572281861, -19.94583874943721824, // =IP
                                 'i',
                                 52.268058064579619781, -31.764051457399769873);

}

int test_main(int, char* [])
{
#if !defined(BOOST_GEOMETRY_TEST_ONLY_ONE_TYPE)
    test_all<boost::tuple<double, double> >();
    test_all<bg::model::point<float, 2, bg::cs::cartesian> >();
#endif
    test_all<bg::model::point<double, 2, bg::cs::cartesian> >();

    return 0;
}
