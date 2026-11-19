// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.

#include <boost/core/lightweight_test.hpp>
#include <boost/geometry.hpp>
#include <algorithm>
#include <sstream>
#include <string>
#include <vector>

namespace bg = boost::geometry;
using point = bg::model::d2::point_xy<double>;
using line = bg::model::linestring<point>;
using polygon = bg::model::polygon<point>;
using points = bg::model::multi_point<point>;
using lines = bg::model::multi_linestring<line>;
using polygons = bg::model::multi_polygon<polygon>;
using input = bg::model::geometry_collection<boost::variant<point, line, polygon, points, lines, polygons>>;
using output = bg::model::geometry_collection<boost::variant<points, lines, polygons>>;

struct tracking_strategy : bg::strategies::relate::cartesian<>
{
    mutable std::vector<void const*> converted;
    template <typename Source, typename Destination>
    void convert_geometry(Source const& source, Destination& destination) const
    {
        converted.push_back(&source);
        bg::convert(source, destination);
    }
};

std::vector<std::string> normalized(output const& collection)
{
    std::vector<std::string> result;
    for (auto const& member : collection)
        bg::traits::visit<output::value_type>::apply([&](auto const& geometry)
        {
            std::ostringstream text;
            text << bg::wkt(geometry);
            result.push_back(text.str());
        }, member);
    std::sort(result.begin(), result.end());
    return result;
}

void check_sources(std::vector<void const*> actual, std::vector<void const*> expected)
{
    std::sort(actual.begin(), actual.end(), std::less<void const*>());
    std::sort(expected.begin(), expected.end(), std::less<void const*>());
    BOOST_TEST(actual == expected);
}

void test_original_members()
{
    input source;
    bg::read_wkt("GEOMETRYCOLLECTION(POINT(1 1),LINESTRING(10 0,11 0),"
        "POLYGON((20 0,20 1,21 1,21 0,20 0)),MULTIPOINT((30 0),(32 0)),"
        "MULTILINESTRING((40 0,41 0),(43 0,44 0)),"
        "MULTIPOLYGON(((50 0,50 1,51 1,51 0,50 0)),((53 0,53 1,54 1,54 0,53 0))))", source);
    input const empty;
    std::vector<void const*> const expected{
        &boost::get<point>(source[0]), &boost::get<line>(source[1]), &boost::get<polygon>(source[2]),
        &boost::get<points>(source[3]), &boost::get<lines>(source[4]), &boost::get<polygons>(source[5])};
    for (bool difference : {false, true})
    {
        output actual, reference;
        tracking_strategy const strategy;
        if (difference)
        {
            bg::difference(source, empty, actual, strategy);
            bg::difference(source, empty, reference);
        }
        else
        {
            bg::union_(source, empty, actual, strategy);
            bg::union_(source, empty, reference);
        }
        check_sources(strategy.converted, expected);
        BOOST_TEST(normalized(actual) == normalized(reference));
    }
}

void test_default_native_box()
{
    using box = bg::model::box<point>;
    using boxes = bg::model::geometry_collection<boost::variant<box>>;
    boxes const source{box(point(0, 0), point(2, 3))}, empty;
    for (bool difference : {false, true})
    {
        output actual, reference;
        tracking_strategy const strategy;
        if (difference)
        {
            bg::difference(source, empty, actual, strategy);
            bg::difference(source, empty, reference);
        }
        else
        {
            bg::union_(source, empty, actual, strategy);
            bg::union_(source, empty, reference);
        }
        check_sources(strategy.converted, {&boost::get<box>(source.front())});
        BOOST_TEST(normalized(actual) == normalized(reference));
        BOOST_TEST_EQ(actual.size(), 1u);
        if (!actual.empty())
            BOOST_TEST_EQ(bg::area(boost::get<polygons>(actual.front())), 6.0);
    }
}

int main()
{
    test_original_members();
    test_default_native_box();
    return boost::report_errors();
}
