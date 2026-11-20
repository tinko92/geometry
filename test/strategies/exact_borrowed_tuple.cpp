// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.

#include <boost/core/lightweight_test.hpp>
#include <boost/geometry/geometries/register/point.hpp>
#include <boost/geometry/strategies/detail/segment_as_subrange.hpp>

struct native_point
{
    double x = 0, y = 0;
    int metadata = 17;
};
BOOST_GEOMETRY_REGISTER_POINT_2D(native_point, double, boost::geometry::cs::cartesian, x, y)

namespace bg = boost::geometry;

void test_native_fallback()
{
    native_point const a{1, 2, 99}, b{3, 4, 88};
    bg::model::pointing_segment<native_point const> segment(a, b);
    bg::strategies::detail::segment_as_subrange<decltype(segment)> range(segment);
    BOOST_TEST_EQ(range.at(0).x, 1);
    BOOST_TEST_EQ(range.at(1).y, 4);
    BOOST_TEST_EQ(range.at(0).metadata, 17);
    BOOST_TEST_EQ(range.at(1).metadata, 17);
}

#if __cplusplus >= 201703L
#include <boost/geometry.hpp>
#include <boost/geometry/strategies/relate/cartesian_exact.hpp>
#include "exact_additional_issue_cases.hpp"
#include <cstring>
#include <tuple>

namespace ec = bg::strategies::detail::exact_cartesian;

void test_record_endpoints()
{
    using p = bg::model::d2::point_xy<double>;
    using segment = bg::model::segment<p>;
    using root = ec::borrowed_root<double, bg::strategies::relate::approximate_output, segment, segment>;
    using point = ec::borrowed_point<root>;
    ec::borrowed_record<root> first, second;
    point const a(first), b(second);
    bg::model::pointing_segment<point const> segment_view(a, b);
    bg::strategies::detail::segment_as_subrange<decltype(segment_view)> range(segment_view);
    // No coordinate access is valid here: these records have no producers.
    BOOST_TEST(&range.at(0).record() == &first);
    BOOST_TEST(&range.at(1).record() == &second);
}

template <typename Policy>
void test_issue_1235()
{
    using p = bg::model::d2::point_xy<double>;
    using rational_point = bg::model::d2::point_xy<ec::exact_number>;
    using rational_polygons = bg::model::multi_polygon<bg::model::polygon<rational_point>>;
    bg::model::polygon<p> a;
    bg::model::multi_polygon<bg::model::polygon<p>> b;
    bool found = false;
    for (auto const& item : exact_additional_issue_cases::pairs())
    {
        if (std::strcmp(item.data.name, "issue_1235_tuple") != 0) continue;
        bg::read_wkt(item.data.first, a);
        bg::read_wkt(item.data.second, b);
        found = true;
        break;
    }
    BOOST_TEST(found);
    bg::strategies::relate::cartesian_exact<double, Policy> strategy;
    bg::correct(a, strategy);
    bg::correct(b, strategy);
    std::tuple<bg::model::multi_point<rational_point>,
               bg::model::multi_linestring<bg::model::linestring<rational_point>>,
               rational_polygons> result;
    bg::intersection(a, b, result, strategy);
    BOOST_TEST(std::get<0>(result).empty());
    BOOST_TEST(std::get<1>(result).empty());
    BOOST_TEST_EQ(std::get<2>(result).size(), 2u);
    bg::strategies::relate::cartesian<ec::exact_number> validation;
    BOOST_TEST(bg::is_valid(std::get<2>(result), validation));
    rational_polygons areal;
    bg::intersection(a, b, areal, strategy);
    BOOST_TEST(bg::relate(areal, std::get<2>(result), bg::de9im::mask("T*F**FFF*"), validation));
}
#endif

int main()
{
    test_native_fallback();
#if __cplusplus >= 201703L
    test_record_endpoints();
    test_issue_1235<bg::strategies::relate::approximate_output>();
    test_issue_1235<bg::strategies::relate::nearest_output>();
#endif
    return boost::report_errors();
}
