// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.

#include <boost/core/lightweight_test.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/strategies/relate/cartesian_exact.hpp>
#include <type_traits>

namespace bg = boost::geometry;
namespace ec = bg::strategies::detail::exact_cartesian;
using native = bg::model::d2::point_xy<double>;
using box = bg::model::box<native>;
using boxes = bg::model::geometry_collection<boost::variant<box>>;
using exact = bg::model::d2::point_xy<ec::exact_number>;
using points = bg::model::multi_point<exact>;
using lines = bg::model::multi_linestring<bg::model::linestring<exact>>;
using polygons = bg::model::multi_polygon<bg::model::polygon<exact>>;
using output = bg::model::geometry_collection<boost::variant<points, lines, polygons>>;

ec::exact_number area(output const& collection)
{
    ec::exact_number result = 0;
    for (auto const& member : collection)
        bg::traits::visit<output::value_type>::apply([&](auto const& geometry)
        {
            using G = std::decay_t<decltype(geometry)>;
            if constexpr (bg::util::is_multi_polygon<G>::value)
                result += bg::area(geometry, bg::strategies::relate::cartesian<ec::exact_number>());
        }, member);
    return result;
}

void test_original_boxes()
{
    boxes const source{box(native(0, 0), native(4, 4))};
    boxes const disjoint{box(native(8, 0), native(9, 1))};
    boxes const overlapping{box(native(2, -1), native(5, 3))};
    boxes const inner{box(native(1, 1), native(3, 3))};
    bg::strategies::relate::cartesian_exact<double> strategy;
    output result;
    bg::union_(source, disjoint, result, strategy);
    BOOST_TEST_EQ(area(result), ec::exact_number(17));
    result.clear();
    bg::union_(source, overlapping, result, strategy);
    BOOST_TEST_EQ(area(result), ec::exact_number(22));
    result.clear();
    bg::difference(source, inner, result, strategy);
    BOOST_TEST_EQ(area(result), ec::exact_number(12));
    result.clear();
    bg::difference(source, disjoint, result, strategy);
    BOOST_TEST_EQ(area(result), ec::exact_number(16));
}

#if defined(BOOST_GEOMETRY_TEST_GC_BOX_INTERSECTION) || defined(BOOST_GEOMETRY_TEST_GC_BOX_BOX_INTERSECTION)
void test_box_intersections()
{
    boxes const source{box(native(0, 0), native(4, 4))};
    boxes const overlapping{box(native(2, -1), native(5, 3))};
    bg::strategies::relate::cartesian_exact<double> strategy;
    output result;
    bg::intersection(source, overlapping, result, strategy);
    BOOST_TEST_EQ(area(result), ec::exact_number(6));
}
#endif

#if defined(BOOST_GEOMETRY_TEST_GC_BOX_INTERSECTION) || defined(BOOST_GEOMETRY_TEST_GC_BOX_LINE_INTERSECTION)
void test_box_line_intersection()
{
    boxes const source{box(native(0, 0), native(4, 4))};
    bg::strategies::relate::cartesian_exact<double> strategy;
    output result;
    using line = bg::model::linestring<native>;
    bg::model::geometry_collection<boost::variant<line>> const crossing{line{{-1, 0.5}, {5, 0.5}}};
    bg::intersection(crossing, source, result, strategy);
    BOOST_TEST_EQ(result.size(), 1u);
    if (result.empty()) return;
    auto const* clipped = boost::get<lines>(&result.front());
    BOOST_TEST(clipped != nullptr);
    if (!clipped || clipped->empty()) return;
    BOOST_TEST_EQ(bg::get<0>(clipped->front().front()), ec::exact_number(0));
    BOOST_TEST_EQ(bg::get<0>(clipped->front().back()), ec::exact_number(4));
    BOOST_TEST_EQ(bg::get<1>(clipped->front().front()), ec::exact_number(1) / 2);
    BOOST_TEST_EQ(bg::get<1>(clipped->front().back()), ec::exact_number(1) / 2);
}
#endif

int main()
{
    test_original_boxes();
#if defined(BOOST_GEOMETRY_TEST_GC_BOX_INTERSECTION) || defined(BOOST_GEOMETRY_TEST_GC_BOX_BOX_INTERSECTION)
    test_box_intersections();
#endif
#if defined(BOOST_GEOMETRY_TEST_GC_BOX_INTERSECTION) || defined(BOOST_GEOMETRY_TEST_GC_BOX_LINE_INTERSECTION)
    test_box_line_intersection();
#endif
    return boost::report_errors();
}
