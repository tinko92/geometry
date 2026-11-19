// Boost.Geometry
// Distributed under the Boost Software License, Version 1.0.

#include <cmath>
#include <string>
#include <type_traits>
#include <utility>
#include <boost/core/lightweight_test.hpp>
#include <boost/geometry/strategies/buffer/cartesian_mp.hpp>
#include <boost/geometry.hpp>
#include "buffer_cases.hpp"

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
void check(std::string const& text, Join const& join, std::size_t components, double area,
           double tolerance, int holes = -1)
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
    if (holes >= 0)
    {
        std::size_t count = 0;
        for (auto const& p : output)
            count += bg::num_interior_rings(p);
        BOOST_TEST_EQ(count, std::size_t(holes));
    }
    for (auto const& p : input)
        for (auto const& v : p.outer())
            BOOST_TEST(bg::covered_by(v, output));
}

void check_touch_certificate()
{
    using point = bg::model::d2::point_xy<mp>;
    using turn_type = bg::detail::overlay::turn_info<point, bg::segment_ratio<mp>>;
    using strategy = bg::strategies::buffer::cartesian_mp;
    mp const step = ldexp(mp(1), -332);
    point const b(1, 1), c(1, 0), d(1, 3), a(2, 1);
    point const misplaced(1, 1 + step);
    turn_type t;
    auto& first = t.operations[0].fraction;
    auto& second = t.operations[1].fraction;
    first.assign(1, 1);
    second.assign(1, 2);
    t.point = misplaced;
    strategy::recover_touch_point(t, a, b, c, d);
    BOOST_TEST_EQ(bg::get<1>(t.point), mp(1));
    BOOST_TEST_EQ(second.numerator(), mp(1));
    BOOST_TEST_EQ(second.denominator(), mp(3));

    // Swapped operands and a reversed containing segment.
    first.assign(1, 2);
    second.assign(1, 1);
    t.point = misplaced;
    strategy::recover_touch_point(t, d, c, a, b);
    BOOST_TEST_EQ(bg::get<1>(t.point), mp(1));
    BOOST_TEST_EQ(first.numerator(), mp(2));
    BOOST_TEST_EQ(first.denominator(), mp(3));

    first.assign(1, 1);
    second.assign(1, 2);
    auto declines = [&](point const& p, point const& q, point const& r, point const& s)
    {
        t.point = misplaced;
        strategy::recover_touch_point(t, p, q, r, s);
        BOOST_TEST_EQ(bg::get<1>(t.point), bg::get<1>(misplaced));
        BOOST_TEST_EQ(second.numerator(), mp(1));
        BOOST_TEST_EQ(second.denominator(), mp(2));
    };
    declines(a, point(1 + step, 1), c, d); // Near, but not on the support.
    declines(point(1, 2), b, c, d); // Collinear overlap is not a singleton.
    declines(a, b, point(1, 2), d); // On the support, outside the segment.
    declines(a, b, point(1, ldexp(mp(1), -512)), d); // Projection cannot be stored exactly.
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

    char const* const rt_w6_subset
        = "MULTIPOLYGON(((7 0,7 1,8 0,7 0)),((2 8,3 9,3 8,2 8)),"
          "((4 3,4 4,5 4,4 3)),((6 5,6 6,7 6,6 5)),((1 5,0 4,0 5,1 5)),"
          "((7 4,6 3,6 4,7 4)),((7 4,7 5,8 5,7 4)))";
    check_touch_certificate();
    check<mp>(rt_w6_subset, bg::strategy::buffer::join_miter{}, 3, 53.79163056034262, 1e-10, 0);
    check<mp>(rt_w6, bg::strategy::buffer::join_miter{}, 1, 79.79898987322333, 1e-10, 1);
    check<double>(rt_w6_subset, bg::strategy::buffer::join_miter{}, 3, 53.79163056034262, 1e-10, 0);
    check<double>(rt_w6, bg::strategy::buffer::join_miter{}, 1, 79.79898987322333, 1e-10, 1);
    return boost::report_errors();
}
