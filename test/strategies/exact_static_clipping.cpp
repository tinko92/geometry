// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.

#define BOOST_GEOMETRY_EXACT_FILTER_STATISTICS
#include <boost/core/lightweight_test.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/strategies/relate/cartesian_exact.hpp>
#include <array>
#include <cstdint>
#include <random>

namespace bg = boost::geometry;
namespace ec = bg::strategies::detail::exact_cartesian;
using rational = ec::exact_number;
using rational_point = bg::model::d2::point_xy<rational>;
using result_type = bg::model::multi_linestring<bg::model::linestring<rational_point>>;

template <typename Coordinate>
void compare_clip(std::array<Coordinate, 8> const& values)
{
    using point = bg::model::d2::point_xy<Coordinate>;
    bg::model::linestring<point> const line{{values[0], values[1]}, {values[2], values[3]}};
    bg::model::box<point> const box({values[4], values[5]}, {values[6], values[7]});
    result_type expected, actual;
    bg::intersection(line, box, expected, bg::strategies::relate::cartesian<rational>());
    bg::intersection(line, box, actual, bg::strategies::relate::cartesian_exact<Coordinate>());
    BOOST_TEST_EQ(actual.size(), expected.size());
    if (actual.size() != expected.size()) return;
    for (std::size_t i = 0; i != actual.size(); ++i)
    {
        BOOST_TEST_EQ(actual[i].size(), expected[i].size());
        if (actual[i].size() != expected[i].size()) continue;
        for (std::size_t j = 0; j != actual[i].size(); ++j)
        {
            BOOST_TEST(bg::get<0>(actual[i][j]) == bg::get<0>(expected[i][j]));
            BOOST_TEST(bg::get<1>(actual[i][j]) == bg::get<1>(expected[i][j]));
        }
    }
}

template <typename Coordinate, typename OutputPolicy>
void test_isolated_contacts()
{
    ec::point_cache_scope<Coordinate> operation;
    using point = bg::model::d2::point_xy<Coordinate>;
    using line = bg::model::linestring<point>;
    bg::model::box<point> const box({0, 0}, {1, 2});
    bg::strategies::relate::cartesian_exact<Coordinate, OutputPolicy> strategy;
    for (line const& input : {line{{1, 1}, {2, 1}}, line{{2, 1}, {1, 1}},
            line{{-1, 1}, {0, 2}, {-1, 3}}})
    {
        bg::model::multi_linestring<line> output;
        bg::intersection(input, box, output, strategy);
        BOOST_TEST(output.empty());
        bg::intersection(input, box, output);
        BOOST_TEST(output.empty());
    }
    {
        line const input{{4096, 4000}, {4100, 4010}};
        bg::model::box<point> const large({0, 0}, {4096, 4096});
        result_type output;
        bg::intersection(input, large, output, strategy);
        BOOST_TEST(output.empty());
    }
    // Discard an isolated contact, then retain the later positive-length piece.
    line const input{{1, 1}, {2, 1}, {2, 0}, {0, 0}};
    result_type output;
    bg::intersection(input, box, output, strategy);
    BOOST_TEST_EQ(output.size(), 1u);
    BOOST_TEST_EQ(output.front().size(), 2u);
    BOOST_TEST_EQ(bg::get<0>(output.front().front()), rational(1));
    BOOST_TEST_EQ(bg::get<0>(output.front().back()), rational(0));
    if constexpr (std::is_same<Coordinate, double>::value)
    {
        line const short_line{{0.5, 1}, {std::nextafter(0.5, 1.0), 1}};
        using single_line = bg::model::linestring<bg::model::d2::point_xy<float>>;
        bg::model::multi_linestring<single_line> rounded;
        bg::intersection(short_line, box, rounded, strategy);
        BOOST_TEST_EQ(rounded.size(), 1u);
        BOOST_TEST_EQ(rounded.front().size(), 2u);
        BOOST_TEST_EQ(bg::get<0>(rounded.front().front()), bg::get<0>(rounded.front().back()));
    }
}

template <typename Coordinate>
void test_clipping()
{
    for (auto values : {std::array<Coordinate, 8>{{-3, 1, 4, 2, 0, 0, 2, 2}},
            {{4, 2, -3, 1, 0, 0, 2, 2}}, {{-2, 0, 0, 0, 0, 0, 2, 2}},
            {{-2, 2, 2, -2, 0, 0, 2, 2}}, {{0, -3, 0, 4, 0, 0, 2, 2}},
            {{1, 1, 1, 1, 0, 0, 2, 2}}, {{-1, -1, -1, -1, 0, 0, 2, 2}}})
        compare_clip(values);
    std::mt19937 random(7293);
    for (int i = 0; i != 2000; ++i)
    {
        std::array<Coordinate, 8> values{{0, 0, 0, 0, -10, -10, 10, 10}};
        for (int j = 0; j != 4; ++j) values[j] = Coordinate(int(random() % 81) - 40);
        compare_clip(values);
    }

    ec::point_cache_scope scope;
    using exact = ec::exact_point<Coordinate>;
    exact first(-3, 1), second(4, 2);
    bg::model::referring_segment<exact> segment(first, second);
    bg::model::box<exact> const bounds(exact(0, 0), exact(2, 2));
    bool a = false, b = false;
    ec::filtering_statistics() = {};
    BOOST_TEST(ec::clipping<Coordinate>().clip_segment(bounds, segment, a, b));
    BOOST_TEST(a && b);
    BOOST_TEST_EQ(ec::filtering_statistics().native_rationals, 0u);
    BOOST_TEST_EQ(ec::filtering_statistics().point_coordinate_values, 0u);
    BOOST_TEST_EQ(ec::filtering_statistics().homogeneous_values, 0u);
}

int main()
{
    test_isolated_contacts<double, bg::strategies::relate::approximate_output>();
    test_isolated_contacts<double, bg::strategies::relate::nearest_output>();
    test_isolated_contacts<std::int64_t, bg::strategies::relate::approximate_output>();
    test_isolated_contacts<std::int64_t, bg::strategies::relate::nearest_output>();
    test_clipping<double>();
    test_clipping<std::int64_t>();
    return boost::report_errors();
}
