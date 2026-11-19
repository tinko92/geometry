// Boost.Geometry
// Unit Test
// Copyright (c) 2026
// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include <boost/geometry/strategies/cartesian/exact/round_buffer.hpp>

#include <boost/core/lightweight_test.hpp>

namespace ec = boost::geometry::strategies::detail::exact_cartesian;
namespace rb = ec::round_buffer;

using ec::exact_number;
using rb::circle_direction;
using rb::offset_side;
using rb::rational_xy;
using rb::straight_offset;

exact_number q(int n, int d)
{
    return exact_number(n) / d;
}

void test_irrational_offset()
{
    rational_xy const a{q(1, 3), q(-2, 7)}, b{a.x + q(1, 2), a.y + q(1, 2)};
    straight_offset const left(a, b, q(2, 3), offset_side::left);
    straight_offset const right(a, b, q(2, 3), offset_side::right);
    BOOST_TEST_EQ(left.normal().squared_length(), q(1, 2));
    BOOST_TEST_EQ(left.normal().numerator().x, q(-1, 2));
    BOOST_TEST_EQ(left.normal().numerator().y, q(1, 2));

    // Each left displacement has magnitude sqrt(2)/3, between 1/3 and 1/2.
    auto const start = left.start();
    auto const end = left.end();
    BOOST_TEST_EQ(start.compare_x(a.x - q(1, 2)), 1);
    BOOST_TEST_EQ(start.compare_x(a.x - q(1, 3)), -1);
    BOOST_TEST_EQ(start.compare_y(a.y + q(1, 3)), 1);
    BOOST_TEST_EQ(start.compare_y(a.y + q(1, 2)), -1);
    BOOST_TEST_EQ(end.compare_x(b.x - q(1, 2)), 1);
    BOOST_TEST_EQ(end.compare_y(b.y + q(1, 2)), -1);
    BOOST_TEST_EQ(right.start().compare_x(a.x + q(1, 3)), 1);
    BOOST_TEST_EQ(right.start().compare_y(a.y - q(1, 3)), -1);

    rational_xy const inside{a.x - q(1, 3), a.y + q(1, 3)};
    rational_xy const outside{a.x - q(1, 2), a.y + q(1, 2)};
    BOOST_TEST_EQ(left.side_of_line(inside), -1);
    BOOST_TEST_EQ(left.side_of_line(outside), 1);
    BOOST_TEST_EQ(right.side_of_line(outside), -1);
    BOOST_TEST_EQ(right.side_of_line({a.x + q(1, 2), a.y - q(1, 2)}), 1);
    BOOST_TEST_EQ(straight_offset(b, a, q(2, 3), offset_side::right).side_of_line(outside), 1);
}

void test_rational_offset()
{
    // A rational (3,4,5) direction, and a translation too large for a double
    // to distinguish the two endpoints.
    exact_number const large(boost::multiprecision::cpp_int(1) << 200);
    rational_xy const a{large + q(1, 3), q(1, 5)}, b{a.x + q(3, 7), a.y + q(4, 7)};
    straight_offset const left(a, b, q(5, 6), offset_side::left);
    straight_offset const right(a, b, q(5, 6), offset_side::right);
    BOOST_TEST_EQ(left.normal().squared_length(), q(25, 49));

    rational_xy const start{a.x - q(2, 3), a.y + q(1, 2)};
    rational_xy const end{b.x - q(2, 3), b.y + q(1, 2)};
    BOOST_TEST_EQ(left.start().compare_x(start.x), 0);
    BOOST_TEST_EQ(left.start().compare_y(start.y), 0);
    BOOST_TEST_EQ(left.end().compare_x(end.x), 0);
    BOOST_TEST_EQ(left.end().compare_y(end.y), 0);
    BOOST_TEST_EQ(left.side_of_line(start), 0);
    BOOST_TEST_EQ(left.side_of_line(end), 0);
    BOOST_TEST_EQ(left.side_of_line({end.x + q(3, 7), end.y + q(4, 7)}), 0);
    BOOST_TEST_EQ(right.side_of_line({b.x + q(2, 3), b.y - q(1, 2)}), 0);

    exact_number const tiny = 1 / large;
    BOOST_TEST_EQ(left.side_of_line({start.x - tiny, start.y}), 1);
    BOOST_TEST_EQ(left.side_of_line({start.x + tiny, start.y}), -1);
    BOOST_TEST_EQ(left.start().compare_x(start.x - tiny), 1);
    BOOST_TEST_EQ(left.start().compare_x(start.x + tiny), -1);

    straight_offset const zero(a, b, 0, offset_side::left);
    BOOST_TEST_EQ(zero.start().compare_x(a.x), 0);
    BOOST_TEST_EQ(zero.side_of_line(b), 0);
    straight_offset const axis({0, 0}, {2, 0}, 1, offset_side::left);
    BOOST_TEST_EQ(axis.end().compare_x(2), 0);
    BOOST_TEST_EQ(axis.end().compare_y(1), 0);
}

void test_rational_samples()
{
    rational_xy const axes[] = {{1, 0}, {0, 1}, {-1, 0}, {0, -1}};
    for (unsigned i = 0; i < 4; ++i)
    {
        auto const start = rb::rational_arc_sample(i, 0);
        auto const end = rb::rational_arc_sample(i, 1);
        BOOST_TEST_EQ(start.x, axes[i].x);
        BOOST_TEST_EQ(start.y, axes[i].y);
        BOOST_TEST_EQ(end.x, axes[(i + 1) % 4].x);
        BOOST_TEST_EQ(end.y, axes[(i + 1) % 4].y);

        auto previous = start;
        for (int j = 1; j <= 7; ++j)
        {
            auto const sample = rb::rational_arc_sample(i, q(j, 7));
            BOOST_TEST_EQ(sample.x * sample.x + sample.y * sample.y, 1);
            BOOST_TEST(previous.x * sample.y - previous.y * sample.x > 0);
            previous = sample;
        }
    }
    auto const sample = rb::rational_arc_sample(0, q(1, 2));
    BOOST_TEST_EQ(sample.x, q(3, 5));
    BOOST_TEST_EQ(sample.y, q(4, 5));
}

void test_chord_certificates()
{
    using rb::minor_arc_radial_error_within;
    circle_direction const x({1, 0}), opposite({-2, 0}), y({0, 1});
    BOOST_TEST(minor_arc_radial_error_within(x, circle_direction({4, 0}), 1, 0));
    BOOST_TEST_NOT(minor_arc_radial_error_within(x, y, 1, 0));
    BOOST_TEST(minor_arc_radial_error_within(x, y, 1, q(3, 10)));
    BOOST_TEST_NOT(minor_arc_radial_error_within(x, y, 1, q(29, 100)));
    BOOST_TEST(minor_arc_radial_error_within(x, opposite, 1, 1));
    BOOST_TEST_NOT(minor_arc_radial_error_within(x, opposite, 1, q(999, 1000)));
    BOOST_TEST(minor_arc_radial_error_within(x, opposite, 1, 3));
    BOOST_TEST(minor_arc_radial_error_within(x, opposite, 0, 0));

    // Opposite y coordinates give an independent midpoint oracle. Use both
    // signs of the dot product and a rational radius that is not dyadic.
    exact_number const r = q(5, 7);
    exact_number const tiny = exact_number(1) / (boost::multiprecision::cpp_int(1) << 200);
    for (int n : {2, 3})
    {
        auto const a = rb::rational_arc_sample(0, q(1, n));
        auto const b = rb::rational_arc_sample(3, q(n - 1, n + 1));
        BOOST_TEST_EQ(a.x, b.x);
        BOOST_TEST_EQ(a.y, -b.y);
        exact_number const error = r * (1 - a.x);
        circle_direction const u(a), v(b);
        BOOST_TEST(minor_arc_radial_error_within(u, v, r, error));
        BOOST_TEST(minor_arc_radial_error_within(v, u, r, error + tiny));
        BOOST_TEST_NOT(minor_arc_radial_error_within(u, v, r, error - tiny));
    }

    // Mixed chord from the irrational endpoint of an offset to a rational
    // circle sample at the same source vertex, and a two-radical chord.
    straight_offset const diagonal({0, 0}, {1, 1}, 1, offset_side::right);
    circle_direction const axis(rb::rational_arc_sample(0, 0));
    BOOST_TEST(minor_arc_radial_error_within(diagonal.normal(), axis, 1, q(1, 10)));
    BOOST_TEST_NOT(minor_arc_radial_error_within(diagonal.normal(), axis, 1, q(1, 20)));
    circle_direction const other({2, 1}); // Dot = 1/sqrt(10).
    BOOST_TEST(minor_arc_radial_error_within(diagonal.normal(), other, 1, q(1, 5)));
    BOOST_TEST_NOT(minor_arc_radial_error_within(diagonal.normal(), other, 1, q(9, 50)));
    // Both directions have irrational coordinates but the midpoint radius is 4/5.
    circle_direction const u({7, 1}), v({1, 7});
    BOOST_TEST(minor_arc_radial_error_within(u, v, 1, q(1, 5)));
    BOOST_TEST_NOT(minor_arc_radial_error_within(u, v, 1, q(1, 5) - tiny));
}

void test_invalid_inputs()
{
    BOOST_TEST_THROWS(straight_offset({0, 0}, {0, 0}, 1, offset_side::left), std::domain_error);
    BOOST_TEST_THROWS(straight_offset({0, 0}, {0, 1}, -1, offset_side::left), std::domain_error);
    BOOST_TEST_THROWS(circle_direction({0, 0}), std::domain_error);
    BOOST_TEST_THROWS(rb::rational_arc_sample(4, 0), std::domain_error);
    BOOST_TEST_THROWS(rb::rational_arc_sample(0, q(-1, 7)), std::domain_error);
    BOOST_TEST_THROWS(rb::rational_arc_sample(0, q(8, 7)), std::domain_error);
    circle_direction const x({1, 0});
    BOOST_TEST_THROWS(rb::minor_arc_radial_error_within(x, x, -1, 0), std::domain_error);
    BOOST_TEST_THROWS(rb::minor_arc_radial_error_within(x, x, 0, -1), std::domain_error);
}

int main()
{
    test_irrational_offset();
    test_rational_offset();
    test_rational_samples();
    test_chord_certificates();
    test_invalid_inputs();
    return boost::report_errors();
}
