// Boost.Geometry (aka GGL, Generic Geometry Library)
// Unit Test

// Copyright (c) 2018 Barend Gehrels, Amsterdam, the Netherlands.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

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

static std::string case_b[2] =
    {
    "MULTIPOLYGON(((0 0,0 4,2 4,2 3,4 3,4 0,0 0)))",
    "MULTIPOLYGON(((-1 -1,-1 8,8 8,8 -1,-1 -1),(2 7,2 3,4 3,4 7,2 7)))"
    };



template <typename Geometry, bg::overlay_type OverlayType>
bool test_overlay(std::string const& caseid,
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
    if (std::fabs(detected_area - expected_area) > 0.01)
    {
        std::cout << std::endl
                  << "ERROR: " << caseid << " area=" << std::setprecision(6) << detected_area
                  << " " << std::setprecision(12) << bg::wkt(g2) << std::endl;
        return false;
    }
    return true;
}

template <typename Ring>
void update(Ring& ring, double x, double y, std::size_t index)
{
    bg::set<0>(ring[index], bg::get<0>(ring[index]) + x);
    bg::set<1>(ring[index], bg::get<1>(ring[index]) + y);
    if (index == 0)
    {
        ring.back() = ring.front();
    }
}

template <typename T, bool Clockwise>
void test_all(bool with_holes)
{
    typedef bg::model::point<T, 2, bg::cs::cartesian> point_type;
    typedef bg::model::polygon<point_type, Clockwise> polygon;
    typedef bg::model::multi_polygon<polygon> multi_polygon;

    const std::string& first = with_holes ? case_b[0] : case_a[0];
    const std::string& second = with_holes ? case_b[1] : case_a[1];
    const double expectation = with_holes ? 73 : 22;

    multi_polygon poly1;
    bg::read_wkt(first, poly1);

    multi_polygon poly2;
    bg::read_wkt(second, poly2);

    int error_count = 0;
    int n = 0;
    for (double factor = 1.0; factor > 1.0e-10; factor /= 10.0)
    {
        int i = 0;
        double const bound = 1.0e-5 * factor;
        double const step = 1.0e-6 * factor;
        for (double x = -bound; x <= bound; x += step, ++i)
        {
            int j = 0;
            for (double y = -bound; y <= bound; y += step, ++j, ++n)
            {
                for (int k = 0; k < 4; k++, ++n)
                {
                    multi_polygon poly_adapted = poly2;

                    if (with_holes)
                        update(bg::interior_rings(poly_adapted.front()).front(), x, y, k);
                    else
                        update(bg::exterior_ring(poly_adapted.front()), x, y, k);


                    std::ostringstream out;
                    out << "case_a_union_" << i << "_" << j << "_" << k;
                    if (!test_overlay<multi_polygon, bg::overlay_union>(out.str(), poly1, poly_adapted, expectation))
                    {
                        error_count++;
                    }

                    if (i == 1 && j == 1)
                    {
                        std::cout << std::endl
                                  << "SAMPLE: " << out.str()
                                  << " " << std::setprecision(18)
                                  << bg::wkt(poly_adapted) << std::endl;
                    }
                }
            }
        }
    }
    std::cout
            << std::endl
            << (with_holes ? "with interior rings" : "only exterior rings")
            << " #cases: " << n << " #errors: " << error_count << std::endl;
    BOOST_CHECK_EQUAL(error_count, 0);
}

int test_main(int, char* [])
{
    typedef double coordinate_type;
    print_configuration();
    test_all<coordinate_type, true>(false);
    test_all<coordinate_type, true>(true);
    return 0;
 }
