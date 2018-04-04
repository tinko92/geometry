// Boost.Geometry (aka GGL, Generic Geometry Library)
// Unit Test

// Copyright (c) 2018 Barend Gehrels, Amsterdam, the Netherlands.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include "debug_sort_by_side_svg.hpp"


#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>

#include <boost/type_traits/is_same.hpp>

#include <geometry_test_common.hpp>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/geometries.hpp>


static std::string case_a[2] =
    {
    "MULTIPOLYGON(((0 0,0 4,2 4,2 3,4 3,4 0,0 0)))",
    "MULTIPOLYGON(((2 7,4 7,4 3,2 3,2 7)))"
    };



template <typename Geometry, bg::overlay_type OverlayType>
void test_overlay(std::string const& caseid,
        Geometry const& g1, Geometry const& g2,
        double expected_area)
{

    typedef typename boost::range_value<Geometry>::type geometry_out;
    typedef bg::detail::overlay::overlay
        <
            Geometry, Geometry,
            bg::detail::overlay::do_reverse<bg::point_order<Geometry>::value>::value,
            OverlayType == bg::overlay_difference
            ? ! bg::detail::overlay::do_reverse<bg::point_order<Geometry>::value>::value
            : bg::detail::overlay::do_reverse<bg::point_order<Geometry>::value>::value,
            bg::detail::overlay::do_reverse<bg::point_order<Geometry>::value>::value,
            geometry_out,
            OverlayType
        > overlay;

    typedef typename bg::strategy::intersection::services::default_strategy
        <
            typename bg::cs_tag<Geometry>::type
        >::type strategy_type;

    strategy_type strategy;

    typedef typename bg::rescale_overlay_policy_type
    <
        Geometry,
        Geometry
    >::type rescale_policy_type;

    rescale_policy_type robust_policy
        = bg::get_rescale_policy<rescale_policy_type>(g1, g2);

    Geometry result;
    bg::detail::overlay::overlay_null_visitor visitor;
    overlay::apply(g1, g2, robust_policy, std::back_inserter(result),
                   strategy, visitor);

    const double detected_area = bg::area(result);
    if (std::fabs(detected_area - expected_area) > 0.1)
    {
        std::cout << caseid << " area=" << std::setprecision(6) << detected_area
                  << " " << std::setprecision(12) << bg::wkt(g2) << std::endl;
    }
}

template <typename Ring>
void update(Ring& ring, double x, double y, std::size_t index)
{
    bg::set<0>(ring[index], bg::get<0>(ring[index]) + x);
    bg::set<1>(ring[index], bg::get<1>(ring[index]) + y);
}

template <typename T, bool Clockwise>
void test_all()
{
    typedef bg::model::point<T, 2, bg::cs::cartesian> point_type;
    typedef bg::model::polygon<point_type, Clockwise> polygon;
    typedef bg::model::multi_polygon<polygon> multi_polygon;

    multi_polygon poly1;
    bg::read_wkt(case_a[0], poly1);

    multi_polygon poly2;
    bg::read_wkt(case_a[1], poly2);

    int i = 0;
    double const step = 1.0e-6;
    for (double x = -1.0e-5; x <= +1.0e-5; x += step, ++i)
    {
        int j = 0;
        for (double y = -1.0e-5; y <= +1.0e-5; y += step, ++j)
        {
            multi_polygon poly_adapted = poly2;

            update(bg::exterior_ring(poly_adapted.front()), x, y, 2);

            std::ostringstream out;
            out << "case_a_union_" << i << "_" << j;
            test_overlay<multi_polygon, bg::overlay_union>(out.str(), poly1, poly_adapted, 22.0);
        }
    }
}

int test_main(int, char* [])
{
    print_configuration();
    test_all<double, true>();
    return 0;
 }
