// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Licensed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt

#include <boost/geometry.hpp>
#include <boost/geometry/strategies/relate/cartesian_exact.hpp>
#include <boost/core/lightweight_test.hpp>

#include <cmath>
#include <type_traits>

namespace bg = boost::geometry;
namespace ec = bg::strategies::detail::exact_cartesian;

template <typename Coordinate>
struct geometry_types
{
    using point = bg::model::d2::point_xy<Coordinate>;
    using polygon = bg::model::polygon<point>;
    using polygons = bg::model::multi_polygon<polygon>;
    using collection = bg::model::geometry_collection<boost::variant<polygons>>;
};

template <typename Collection>
std::size_t polygon_count(Collection const& collection)
{
    std::size_t count = 0;
    for (auto const& variant : collection)
    {
        bg::traits::visit<typename Collection::value_type>::apply([&](auto const& polygons)
        {
            count += polygons.size();
        }, variant);
    }
    return count;
}

void test_final_conversion()
{
    using types = geometry_types<double>;
    types::polygons triangle, subtractor, enclosing;
    bg::read_wkt("MULTIPOLYGON(((0 0,0 0.25,0.25 0,0 0)))", triangle);
    bg::read_wkt("MULTIPOLYGON(((-1 -1,-1 0.125,1 0.125,1 -1,-1 -1)))", subtractor);
    bg::read_wkt("MULTIPOLYGON(((-1 -1,-1 1,1 1,1 -1,-1 -1)))", enclosing);
    types::collection const a{triangle}, b{subtractor}, c{enclosing};
    bg::strategies::relate::cartesian_exact<double> strategy;
    bg::strategies::relate::cartesian<ec::exact_number> validation;
    using output = geometry_types<int>::collection;
    static_assert(!std::is_same<output, decltype(strategy)::overlay_output_type<output>>::value,
                  "Collection intermediates must retain exact coordinates");
    output rounded;
    bg::difference(a, b, rounded, strategy);
    BOOST_TEST_EQ(polygon_count(rounded), 1u);
    geometry_types<ec::exact_number>::collection exact;
    bg::difference(a, b, exact, strategy);
    BOOST_TEST_EQ(polygon_count(exact), 1u);
    bg::traits::visit<decltype(exact)::value_type>::apply([&](auto const& result)
    {
        BOOST_TEST_EQ(bg::area(result, validation), ec::exact_number(1) / 128);
    }, exact.front());
    // GC intersection requires all three multi-geometry output alternatives.
    using point = geometry_types<int>::point;
    using full_output = bg::model::geometry_collection<boost::variant
        <bg::model::multi_point<point>, bg::model::multi_linestring<bg::model::linestring<point>>,
         geometry_types<int>::polygons>>;
    full_output intersection;
    bg::intersection(a, c, intersection, strategy);
    BOOST_TEST_EQ(polygon_count(intersection), 1u);
    // Existing members are retained, not reconverted or included in computation.
    bg::difference(a, b, intersection, strategy);
    BOOST_TEST_EQ(polygon_count(intersection), 2u);
}

void test_relate_working_geometry()
{
    using types = geometry_types<double>;
    types::polygon notched, bridge;
    bg::exterior_ring(notched) = {{0, 0}, {0, 3}, {3, 3}, {3, 1.25},
        {1, 1}, {3, 0.75}, {3, 0}, {0, 0}};
    double const left = std::nextafter(1.0, 2.0);
    bg::exterior_ring(bridge) = {{left, 0}, {left, 3}, {4, 3}, {4, 0}, {left, 0}};
    types::collection const input{types::polygons{notched}, types::polygons{bridge}};
    using point = geometry_types<ec::exact_number>::point;
    using query_collection = bg::model::geometry_collection<boost::variant<point>>;
    ec::exact_number const x = ec::exact_number(1) + ec::exact_number(1)
        / (boost::multiprecision::cpp_int(1) << 53);
    query_collection const query{point(x, 1)};
    bg::strategies::relate::cartesian_exact<double> strategy;
    // The query lies in the true triangular hole, not on its rounded boundary.
    BOOST_TEST_EQ(bg::relation(input, query, strategy).str(), "FF2FF10F2");
    BOOST_TEST_EQ(bg::relation(query, input, strategy).str(), "FF0FFF212");
}

void test_gc_copy_support_lifetime()
{
    using types = geometry_types<double>;
    using exact_types = geometry_types<ec::exact_number>;
    types::polygons source, first_cut, second_cut, disjoint;
    bg::read_wkt("MULTIPOLYGON(((0 0,0 4,4 4,4 0,0 0)))", source);
    bg::read_wkt("MULTIPOLYGON(((-1 -1,-1 1,5 1,5 -1,-1 -1)))", first_cut);
    bg::read_wkt("MULTIPOLYGON(((-1 3,-1 5,5 5,5 3,-1 3)))", second_cut);
    bg::read_wkt("MULTIPOLYGON(((8 0,8 1,9 1,9 0,8 0)))", disjoint);
    types::collection const original{source}, cuts{first_cut, second_cut}, other{disjoint};
    bg::strategies::relate::cartesian_exact<double> strategy;
    bg::strategies::relate::cartesian<ec::exact_number> validation;
    // Owning rational output remains usable after the public operation's Root
    // context dies. Repeated cuts exercise retained original support labels.
    exact_types::collection result;
    bg::difference(original, cuts, result, strategy);
    BOOST_TEST_EQ(polygon_count(result), 1u);
    ec::exact_number area = 0;
    for (auto const& member : result)
        bg::traits::visit<exact_types::collection::value_type>::apply(
            [&](auto const& polygons) { area += bg::area(polygons, validation); }, member);
    BOOST_TEST_EQ(area, ec::exact_number(8));

    using p = exact_types::point;
    using full_output = bg::model::geometry_collection<boost::variant
        <bg::model::multi_point<p>, bg::model::multi_linestring<bg::model::linestring<p>>,
         exact_types::polygons>>;
    full_output combined;
    bg::union_(original, other, combined, strategy);
    BOOST_TEST_EQ(combined.size(), 2u);
    area = 0;
    for (auto const& member : combined)
        bg::traits::visit<full_output::value_type>::apply([&](auto const& geometry)
        {
            using geometry_type = std::decay_t<decltype(geometry)>;
            if constexpr (bg::util::is_multi_polygon<geometry_type>::value)
                area += bg::area(geometry, validation);
            else BOOST_TEST(false);
        }, member);
    BOOST_TEST_EQ(area, ec::exact_number(17));
}

int main()
{
    test_final_conversion();
    test_relate_working_geometry();
    test_gc_copy_support_lifetime();
    return boost::report_errors();
}
