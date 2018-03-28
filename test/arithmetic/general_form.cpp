// Boost.Geometry (aka GGL, Generic Geometry Library)
// Unit Test

// Copyright (c) 2018 Barend Gehrels, Amsterdam, the Netherlands.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include <geometry_test_common.hpp>

#include <boost/geometry/arithmetic/general_form.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/io/wkt/wkt.hpp>


template <typename T, typename C>
void verify_point_on_line(bg::arithmetic::general_form<T> const& f,
                          C const& x, C const& y)
{
    T const epsilon = 1.0e-5;
    BOOST_CHECK_SMALL(f.a * x + f.b * y + f.c, epsilon);
}

template <typename T>
void test_construct_line()
{
    // Horizontal through origin
    bg::arithmetic::general_form<T> p = bg::arithmetic::construct_line<T>(0, 0, 10, 0);
    verify_point_on_line(p, 0, 0);
    verify_point_on_line(p, 10, 0);
    bg::arithmetic::general_form<T> n = bg::arithmetic::normalize_line<T>(p);
    verify_point_on_line(n, 0, 0);
    verify_point_on_line(n, 10, 0);

    // Horizontal line above origin
    p = bg::arithmetic::construct_line<T>(0, 5, 10, 5);
    verify_point_on_line(p, 0, 5);
    verify_point_on_line(p, 10, 5);

    // Vertical through origin
    p = bg::arithmetic::construct_line<T>(0, 0, 0, 10);
    verify_point_on_line(p, 0, 0);
    verify_point_on_line(p, 0, 10);

    // Vertical line left from origin
    p = bg::arithmetic::construct_line<T>(5, 0, 5, 10);
    verify_point_on_line(p, 5, 0);
    verify_point_on_line(p, 5, 10);

    // Diagonal through origin
    p = bg::arithmetic::construct_line<T>(0, 0, 8, 10);
    verify_point_on_line(p, 0, 0);
    verify_point_on_line(p, 8, 10);

    // Diagonal not through origin
    p = bg::arithmetic::construct_line<T>(5, 2, -8, 10);
    verify_point_on_line(p, 5, 2);
    verify_point_on_line(p, -8, 10);
}

template <typename T>
void test_normalize_line()
{
    // Diagonal not through origin
    bg::arithmetic::general_form<T> p = bg::arithmetic::construct_line<T>(5, 2, -8, 10);
    BOOST_CHECK(! p.normalized);

    bg::arithmetic::general_form<T> n = bg::arithmetic::normalize_line<T>(p);
    BOOST_CHECK(n.normalized);
    verify_point_on_line(n, 5, 2);
    verify_point_on_line(n, -8, 10);

    BOOST_CHECK(p.a != n.a);
    BOOST_CHECK(p.b != n.b);
    BOOST_CHECK(p.c != n.c);
}

template <typename T>
void test_distance_measure()
{
    T const epsilon = 1.0e-5;

    // Horizontal line going right
    bg::arithmetic::general_form<T> p = bg::arithmetic::construct_line<T>(0, 0, 10, 0);

    // Point above (= on left side)
    T d = bg::arithmetic::signed_comparable_distance(p, 5, 5);
    BOOST_CHECK_CLOSE(d, 25.0, epsilon);

    // Point below (= on right side)
    d = bg::arithmetic::signed_comparable_distance(p, 5, -5);
    BOOST_CHECK_CLOSE(d, -25.0, epsilon);

    // Diagonal not through origin, from right (down) to left (up)
    p = bg::arithmetic::construct_line<T>(5, 2, -7, 10);
    d = bg::arithmetic::signed_comparable_distance(p, 5, 2);
    BOOST_CHECK_SMALL(d, epsilon);
    d = bg::arithmetic::signed_comparable_distance(p, -7, 10);
    BOOST_CHECK_SMALL(d, epsilon);

    // vector is (-12, 8), move (-3,2) on the line from (5,2)
    d = bg::arithmetic::signed_comparable_distance(p, 2, 4);
    BOOST_CHECK_SMALL(d, epsilon);

    // Go perpendicular (2,3) from (2,4) up, so right of the line (negative)
    d = bg::arithmetic::signed_comparable_distance(p, 4, 7);
    BOOST_CHECK_CLOSE(d, -(4 + 9), epsilon);

    // Go perpendicular (2,3) from (2,4) down, so left of the line (positive)
    d = bg::arithmetic::signed_comparable_distance(p, 0, 1);
    BOOST_CHECK_CLOSE(d, 4 + 9, epsilon);
}


template <typename T>
void test_get_intersection()
{
    // Diagonal lines (first is same as in distance measure,
    // second is perpendicular and used there for distance measures)
    bg::arithmetic::general_form<T> p = bg::arithmetic::construct_line<T>(5, 2, -7, 10);
    bg::arithmetic::general_form<T> q = bg::arithmetic::construct_line<T>(4, 7, 0, 1);

    typedef bg::model::point<T, 2, bg::cs::cartesian> point_type;
    point_type ip;
    bool doubt = false;
    BOOST_CHECK(bg::arithmetic::get_intersection(ip, doubt, p, q));
    BOOST_CHECK(!doubt);

    T const epsilon = 1.0e-5;
    BOOST_CHECK_CLOSE(bg::get<0>(ip), 2, epsilon);
    BOOST_CHECK_CLOSE(bg::get<1>(ip), 4, epsilon);

    verify_point_on_line(p, bg::get<0>(ip), bg::get<1>(ip));
    verify_point_on_line(q, bg::get<0>(ip), bg::get<1>(ip));
}

void increase_counts(std::size_t& count,
                     std::size_t& count_crossing,
                     std::size_t& count_doubt,
                     bool crossing,
                     bool doubt)
{
    count++;
    if (crossing)
    {
        count_crossing++;
    }
    if (doubt)
    {
        count_doubt++;
    }
}

template <typename T>
void test_random(std::size_t count,
                 T const& xmin, T const& ymin,
                 T const& xmax, T const& ymax,
                 T const& epsilon)
{
    srand(12345);
    typedef bg::model::point<T, 2, bg::cs::cartesian> point_type;

    bg::arithmetic::general_form<T> p
            = bg::arithmetic::construct_line<T>(xmin, ymin, xmax, ymax);

    T const range = xmax - xmin;

    std::size_t count_all = 0;
    std::size_t count_doubt = 0;
    std::size_t count_crossing = 0;

    for (std::size_t i = 0; i < count; i++)
    {
        T const multiplier1 = rand() % 4;
        T const multiplier2 = rand() % 4;

        // Generate between -1.0 and 1.0
        double const value1 = -1.0 + (rand() % 100) / 50.0;
        double const value2 = -1.0 + (rand() % 100) / 50.0;

        T const x = range * multiplier1 * value1;
        T const y = range * multiplier2 * value2;

        bg::arithmetic::general_form<T> q
                = bg::arithmetic::construct_line<T>(xmin, ymin, x, y);

        point_type ip;
        bg::set<0>(ip, -99);
        bg::set<1>(ip, -99);

        bool doubt = false;
        bool crossing = bg::arithmetic::get_intersection(ip, doubt, p, q);

        increase_counts(count_all, count_crossing, count_doubt, crossing, doubt);

        if (crossing && !doubt)
        {
            BOOST_CHECK_CLOSE(bg::get<0>(ip), xmin, epsilon);
            BOOST_CHECK_CLOSE(bg::get<1>(ip), ymin, epsilon);
        }
    }
}

template <typename T>
void test_nearly_collinear(T const& threshold,
                           T const& xmin, T const& ymin,
                           T const& xmax, T const& ymax,
                           T const& epsilon)
{
    typedef bg::model::point<T, 2, bg::cs::cartesian> point_type;

    bg::arithmetic::general_form<T> p
            = bg::arithmetic::construct_line<T>(xmin, ymin, xmax, ymax);

    std::size_t count_all = 0;
    std::size_t count_doubt = 0;
    std::size_t count_crossing = 0;

    for (T i = 1.0; i > threshold; i /= 10.0)
    {
        bg::arithmetic::general_form<T> q
                = bg::arithmetic::construct_line<T>(xmin, ymin, xmax, ymax - i);

        point_type ip;
        bg::set<0>(ip, -99);
        bg::set<1>(ip, -99);

        bool doubt = false;
        bool crossing = bg::arithmetic::get_intersection(ip, doubt, p, q);

        if (crossing && !doubt)
        {
            BOOST_CHECK_CLOSE(bg::get<0>(ip), xmin, epsilon);
            BOOST_CHECK_CLOSE(bg::get<1>(ip), ymin, epsilon);
        }

        increase_counts(count_all, count_crossing, count_doubt, crossing, doubt);
    }
}

template <typename T>
void test_component()
{
    bg::arithmetic::general_form<T> p = bg::arithmetic::construct_line<T>(0, 0, 9, 1);

    BOOST_CHECK(bg::arithmetic::more_horizontal(p));
    BOOST_CHECK(bg::arithmetic::has_horizontal_component(p));
    BOOST_CHECK(bg::arithmetic::has_vertical_component(p));

    p = bg::arithmetic::construct_line<T>(0, 0, 1, 8);
    BOOST_CHECK(! bg::arithmetic::more_horizontal(p));
    BOOST_CHECK(bg::arithmetic::has_horizontal_component(p));
    BOOST_CHECK(bg::arithmetic::has_vertical_component(p));

    p = bg::arithmetic::construct_line<T>(0, 0, 9, 0);
    BOOST_CHECK(bg::arithmetic::more_horizontal(p));
    BOOST_CHECK(bg::arithmetic::has_horizontal_component(p));
    BOOST_CHECK(! bg::arithmetic::has_vertical_component(p));

    p = bg::arithmetic::construct_line<T>(0, 0, 0, 8);
    BOOST_CHECK(! bg::arithmetic::more_horizontal(p));
    BOOST_CHECK(! bg::arithmetic::has_horizontal_component(p));
    BOOST_CHECK(bg::arithmetic::has_vertical_component(p));
}

template <typename T>
void test_same_direction()
{
    bg::arithmetic::general_form<T> p, q;

    p = bg::arithmetic::construct_line<T>(2, 1, 12, 11);
    q = bg::arithmetic::construct_line<T>(12, 11, 2, 1);
    BOOST_CHECK(! bg::arithmetic::similar_direction(p, q));

    p = bg::arithmetic::construct_line<T>(0, 0, 10, 0);
    q = bg::arithmetic::construct_line<T>(10, 0, 0, 0);
    BOOST_CHECK(! bg::arithmetic::similar_direction(p, q));

    p = bg::arithmetic::construct_line<T>(0, 0, 0, 10);
    q = bg::arithmetic::construct_line<T>(0, 10, 0, 0);
    BOOST_CHECK(! bg::arithmetic::similar_direction(p, q));

    p = bg::arithmetic::construct_line<T>(0, 0, 10, 10);
    q = bg::arithmetic::construct_line<T>(0, 0, 10, 10);
    BOOST_CHECK(bg::arithmetic::similar_direction(p, q));

    p = bg::arithmetic::construct_line<T>(0, 0, 10, 0);
    q = bg::arithmetic::construct_line<T>(0, 0, 10, 0);
    BOOST_CHECK(bg::arithmetic::similar_direction(p, q));

    p = bg::arithmetic::construct_line<T>(0, 0, 0, 10);
    q = bg::arithmetic::construct_line<T>(0, 0, 0, 10);
    BOOST_CHECK(bg::arithmetic::similar_direction(p, q));

    // (Nearly) perpendicular lines:
    p = bg::arithmetic::construct_line<T>(0, 0, 10, 10);
    q = bg::arithmetic::construct_line<T>(0, 0, -10, 10);
    BOOST_CHECK(! bg::arithmetic::similar_direction(p, q));

    // 45 deg
    p = bg::arithmetic::construct_line<T>(0, 0, 10, 10);
    q = bg::arithmetic::construct_line<T>(0, 0, 0, 10);
    BOOST_CHECK(bg::arithmetic::similar_direction(p, q));

    // a bit more than 45 deg
    p = bg::arithmetic::construct_line<T>(0, 0, 10, 10);
    q = bg::arithmetic::construct_line<T>(0, 0, -1, 10);
    BOOST_CHECK(! bg::arithmetic::similar_direction(p, q));

    // 135 deg
    p = bg::arithmetic::construct_line<T>(0, 0, 10, 10);
    q = bg::arithmetic::construct_line<T>(0, 0, -10, 0);
    BOOST_CHECK(! bg::arithmetic::similar_direction(p, q));
}


template <typename T>
void test_all(T threshold, T const& epsilon)
{
    test_construct_line<T>();
    test_normalize_line<T>();
    test_distance_measure<T>();
    test_get_intersection<T>();
    test_component<T>();
    test_same_direction<T>();

    // Random test, currently turned off (but still working)
    if (false)
    {
        test_random(2000, T(100.0), T(100.0), T(400.0), T(400.0), epsilon);
    }

    test_nearly_collinear(threshold, T(0.01), T(0.01), T(0.02), T(0.01), epsilon);
    test_nearly_collinear(threshold, T(0.1), T(0.1), T(0.2), T(0.1), epsilon);
    test_nearly_collinear(threshold, T(1.0), T(1.0), T(2.0), T(1.0), epsilon);
    test_nearly_collinear(threshold, T(10.0), T(10.0), T(20.0), T(10.0), epsilon);
    test_nearly_collinear(threshold, T(100.0), T(100.0), T(200.0), T(100.0), epsilon);
    test_nearly_collinear(threshold, T(1000.0), T(1000.0), T(2000.0), T(1000.0), epsilon);
    test_nearly_collinear(threshold, T(10000.0), T(10000.0), T(20000.0), T(10000.0), epsilon);
    test_nearly_collinear(threshold, T(100000.0), T(100000.0), T(200000.0), T(100000.0), epsilon);
    test_nearly_collinear(threshold, T(1000000.0), T(1000000.0), T(2000000.0), T(1000000.0), epsilon);
}

int test_main(int, char* [])
{
    test_all<long double>(1.0e-16, 1.0e-7);
    test_all<double>(1.0e-10, 1.0e-7);
    test_all<float>(1.0e-6, 1.0e-3);
    return 0;
}
