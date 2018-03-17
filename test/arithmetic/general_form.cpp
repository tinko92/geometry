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
//    std::cout << "=== RANDOM " << string_from_type<T>::name()
//              << " " << xmin << " " << ymin
//              << " " << xmax << " " << ymax
//              << std::endl;
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
        double value1 = -1.0 + (rand() % 100) / 50.0; // between -1.0 and 1.0
        double value2 = -1.0 + (rand() % 100) / 50.0;
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

//        if (!crossing || doubt || i < 10)
//        {
//            std::cout << std::boolalpha << string_from_type<T>::name()
//                      << " " << i
//                      << " " << x
//                      << " " << y
//                      << " " << crossing
//                      << " " << doubt
//                      << std::setprecision(12)
//                      << " " << bg::wkt(ip)
//                      << std::endl;
//        }
    }
    std::cout << string_from_type<T>::name()
              << " " << count_all
              << " " << count_crossing
              << " " << count_doubt
              << std::endl;
}

template <typename T>
void test_nearly_collinear(T const& threshold,
                           T const& xmin, T const& ymin,
                           T const& xmax, T const& ymax,
                           T const& epsilon)
{
//    std::cout << "=== " << string_from_type<T>::name()
//              << " " << xmin << " " << ymin
//              << " " << xmax << " " << ymax
//              << std::endl;
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

        std::cout << std::boolalpha << string_from_type<T>::name()
                  << " " << i
                  << " " << crossing
                  << " " << doubt
                  << std::setprecision(12)
                  << " " << bg::wkt(ip)
                  << std::endl;
    }
    std::cout << string_from_type<T>::name()
              << " " << count_all
              << " " << count_crossing
              << " " << count_doubt
              << std::endl;
}

template <typename T>
void test_component()
{
    bg::arithmetic::general_form<T> p = bg::arithmetic::construct_line<T>(0, 0, 9, 1);

    std::cout << p.a << " " << p.b << " " << p.c << " " << p.a / p.b << std::endl;

    BOOST_CHECK(p.more_horizontal());
    BOOST_CHECK(p.has_horizontal_component());
    BOOST_CHECK(p.has_vertical_component());

    p = bg::arithmetic::construct_line<T>(0, 0, 1, 8);
    BOOST_CHECK(! p.more_horizontal());
    BOOST_CHECK(p.has_horizontal_component());
    BOOST_CHECK(p.has_vertical_component());

    std::cout << p.a << " " << p.b << " " << p.c << " " << p.b / p.a << std::endl;


    p = bg::arithmetic::construct_line<T>(0, 0, 9, 0);
    BOOST_CHECK(p.more_horizontal());
    BOOST_CHECK(p.has_horizontal_component());
    BOOST_CHECK(! p.has_vertical_component());

    p = bg::arithmetic::construct_line<T>(0, 0, 0, 8);
    BOOST_CHECK(! p.more_horizontal());
    BOOST_CHECK(! p.has_horizontal_component());
    BOOST_CHECK(p.has_vertical_component());
}

template <typename T>
void test_same_direction()
{
    bg::arithmetic::general_form<T> p, q;

    p = bg::arithmetic::construct_line<T>(2, 1, 12, 11);
    q = bg::arithmetic::construct_line<T>(12, 11, 2, 1);
    BOOST_CHECK(! p.similar_direction(q));

    p = bg::arithmetic::construct_line<T>(0, 0, 10, 0);
    q = bg::arithmetic::construct_line<T>(10, 0, 0, 0);
    BOOST_CHECK(! p.similar_direction(q));

    p = bg::arithmetic::construct_line<T>(0, 0, 0, 10);
    q = bg::arithmetic::construct_line<T>(0, 10, 0, 0);
    BOOST_CHECK(! p.similar_direction(q));

    p = bg::arithmetic::construct_line<T>(0, 0, 10, 10);
    q = bg::arithmetic::construct_line<T>(0, 0, 10, 10);
    BOOST_CHECK(p.similar_direction(q));

    p = bg::arithmetic::construct_line<T>(0, 0, 10, 0);
    q = bg::arithmetic::construct_line<T>(0, 0, 10, 0);
    BOOST_CHECK(p.similar_direction(q));

    p = bg::arithmetic::construct_line<T>(0, 0, 0, 10);
    q = bg::arithmetic::construct_line<T>(0, 0, 0, 10);
    BOOST_CHECK(p.similar_direction(q));

    // (Nearly) perpendicular lines:
    p = bg::arithmetic::construct_line<T>(0, 0, 10, 10);
    q = bg::arithmetic::construct_line<T>(0, 0, -10, 10);
    BOOST_CHECK(! p.similar_direction(q));

    // 45 deg
    p = bg::arithmetic::construct_line<T>(0, 0, 10, 10);
    q = bg::arithmetic::construct_line<T>(0, 0, 0, 10);
    BOOST_CHECK(p.similar_direction(q));

    // a bit more than 45 deg
    p = bg::arithmetic::construct_line<T>(0, 0, 10, 10);
    q = bg::arithmetic::construct_line<T>(0, 0, -1, 10);
    BOOST_CHECK(! p.similar_direction(q));

    // 135 deg
    p = bg::arithmetic::construct_line<T>(0, 0, 10, 10);
    q = bg::arithmetic::construct_line<T>(0, 0, -10, 0);
    BOOST_CHECK(! p.similar_direction(q));
}


template <typename T>
void test_all(T threshold, T const& epsilon)
{
    test_component<T>();
    test_same_direction<T>();

    test_random(2000, T(100.0), T(100.0), T(400.0), T(400.0), epsilon);
//    return;

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
