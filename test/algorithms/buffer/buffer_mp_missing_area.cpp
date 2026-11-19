// Boost.Geometry
// Distributed under the Boost Software License, Version 1.0.

#include <cmath>
#include <type_traits>
#include <utility>
#include <boost/core/lightweight_test.hpp>
#include <boost/geometry/strategies/buffer/cartesian_mp.hpp>
#include <boost/geometry.hpp>

namespace bg = boost::geometry;
using mp = boost::multiprecision::cpp_bin_float_100;

template <typename... Args>
void buffer(std::true_type, Args&&... args)
{
    bg::buffer(std::forward<Args>(args)..., bg::strategies::buffer::cartesian_mp{});
}

template <typename... Args>
void buffer(std::false_type, Args&&... args)
{
    bg::buffer(std::forward<Args>(args)...);
}

template <typename C, typename Join>
void check(char const* text, Join const& join, std::size_t components, double area, double tolerance)
{
    using point = bg::model::d2::point_xy<C>;
    using polygon = bg::model::polygon<point>;
    using multi = bg::model::multi_polygon<polygon>;
    static_assert(std::is_same<typename bg::strategies::buffer::services::default_strategy<multi>::type,
                              bg::strategies::buffer::cartesian<>>::value,
                  "Including the MP header must not change the default strategy");
    multi input, output;
    bg::read_wkt(text, input);
    bg::correct(input);
    BOOST_TEST(bg::is_valid(input));
    buffer(std::is_same<C, mp>{}, input, output, bg::strategy::buffer::distance_symmetric<C>(1),
        bg::strategy::buffer::side_straight{}, join, bg::strategy::buffer::end_flat{},
        bg::strategy::buffer::point_circle(88));
    BOOST_TEST_EQ(output.size(), components);
    BOOST_TEST(bg::is_valid(output));
    BOOST_TEST(std::abs(static_cast<double>(bg::area(output)) - area) < tolerance);
    for (auto const& p : input)
        for (auto const& v : p.outer())
            BOOST_TEST(bg::covered_by(v, output));
}

int main()
{
    char const* const crossing
        = "MULTIPOLYGON(((6 5,6 6,7 6,6 5)),((7 4,6 3,6 4,7 4)),((7 4,7 5,8 5,7 4)))";
    char const* const cardinal
        = "MULTIPOLYGON(((0 4,0 5,1 4,0 4)),((0 7,0 8,1 8,0 7)),"
          "((3 2,2 2,2 3,3 2)),((2 5,1 5,1 6,2 5)))";
    double const miter_area = 11.5 + 5 * std::sqrt(2.0);
    check<mp>(crossing, bg::strategy::buffer::join_miter{}, 1, miter_area, 1e-10);
    check<double>(crossing, bg::strategy::buffer::join_miter{}, 1, miter_area, 1e-10);
    check<mp>(cardinal, bg::strategy::buffer::join_round(32), 1, 24.44698493496707, 1e-10);
    // Round fan point counts can differ by one with ordinary binary64
    // atan2/ceil. Keep the area and coverage checks on that unchanged path.
    check<double>(cardinal, bg::strategy::buffer::join_round(32), 1, 24.449, 0.02);
    return boost::report_errors();
}
