// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.
// Directed point ordering replaces numerical exact segment parameters.

#include <boost/core/lightweight_test.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/linestring.hpp>
#include <boost/geometry/strategies/cartesian/exact/segment_position.hpp>
#include <boost/geometry/algorithms/detail/overlay/turn_info.hpp>
#include <boost/geometry/strategies/relate/intersection_types.hpp>
#include <boost/geometry/algorithms/detail/turns/compare_turns.hpp>
#include <array>
#include <cfenv>
#include <cstdint>
#include <limits>
#include <type_traits>
#include <utility>

namespace bg = boost::geometry;
namespace ec = bg::strategies::detail::exact_cartesian;
using rational = ec::exact_number;

template <typename...> using void_t = void;
template <typename T, typename = void> struct has_fraction : std::false_type {};
template <typename T> struct has_fraction<T, void_t<decltype(std::declval<T>().fraction)>> : std::true_type {};
template <typename T, typename = void> struct has_fractions : std::false_type {};
template <typename T> struct has_fractions<T, void_t<decltype(std::declval<T>().fractions)>> : std::true_type {};

using point = ec::exact_point<double>;
struct exact_strategy_types { using segment_position_type = ec::segment_position_tag; };
struct position_strategy : exact_strategy_types
{
    static auto segment_position() { return ec::segment_position<double>(); }
};
struct legacy_strategy_types { template <typename> using segment_ratio_type = bg::segment_ratio<int>; };
static_assert(std::is_same<bg::strategies::relate::services::segment_ratio_type<exact_strategy_types, point>::type,
                          ec::segment_position_tag>::value, "New exact strategy hook selects no fraction");
static_assert(std::is_same<bg::strategies::relate::services::segment_ratio_type<legacy_strategy_types, point>::type,
                          bg::segment_ratio<int>>::value, "Legacy strategy ratio customization is unchanged");
static_assert(std::is_empty<ec::segment_position_tag>::value, "Only a type-level exact position tag");
static_assert(!has_fraction<bg::detail::overlay::turn_operation<point, ec::segment_position_tag>>::value,
              "No placeholder fraction in exact operations");
static_assert(!has_fractions<bg::segment_intersection_points<point, ec::segment_position_tag>>::value,
              "Exact intersection results contain points and count only");
static_assert(has_fraction<bg::detail::overlay::turn_operation<bg::model::d2::point_xy<double>,
              bg::segment_ratio<double>>>::value, "Legacy fractions remain available");

template <typename C>
void native_positions()
{
    using p = bg::model::d2::point_xy<C>;
    std::array<p, 9> points;
    for (int i = 0; i < 9; ++i) points[i] = p(C(i % 3 - 1), C(i / 3 - 1));
    for (auto const& a : points) for (auto const& b : points) for (auto const& v : points)
    {
        rational const ax = ec::exact_value(bg::get<0>(a)), ay = ec::exact_value(bg::get<1>(a));
        rational const bx = ec::exact_value(bg::get<0>(b)), by = ec::exact_value(bg::get<1>(b));
        rational const x = ec::exact_value(bg::get<0>(v)), y = ec::exact_value(bg::get<1>(v));
        int expected;
        if (ax == bx && ay == by) expected = x == ax && y == ay ? 5 : 6;
        else if ((bx - ax) * (y - ay) != (by - ay) * (x - ax)) expected = 7;
        else
        {
            rational const t = ax != bx ? rational((x - ax) / (bx - ax)) : rational((y - ay) / (by - ay));
            expected = t < 0 ? 0 : t == 0 ? 1 : t < 1 ? 2 : t == 1 ? 3 : 4;
            BOOST_TEST_EQ(ec::segment_position<C>::apply(a, b, v, a), t.sign());
            BOOST_TEST_EQ(ec::segment_position<C>::apply(a, b, a, v), -t.sign());
        }
        BOOST_TEST_EQ(ec::segment_position<C>::classify(a, b, v), expected);
    }
}

int main()
{
    using native_point = bg::model::d2::point_xy<double>;
    using legacy_operation = bg::detail::overlay::turn_operation<native_point, bg::segment_ratio<double>>;
    legacy_operation aggregate{bg::detail::overlay::operation_union,
        bg::segment_identifier(0, -1, -1, 3), bg::segment_ratio<double>(1, 2)};
    BOOST_TEST_EQ(aggregate.operation, bg::detail::overlay::operation_union);
    BOOST_TEST_EQ(aggregate.seg_id.segment_index, 3);
    BOOST_TEST(aggregate.fraction == bg::segment_ratio<double>(1, 2));
    using line = bg::model::linestring<native_point>;
    using turn = bg::detail::overlay::turn_info<point, ec::segment_position_tag>;
    turn left, right;
    left.point = point(0, 1); right.point = point(3, 1);
    left.operations[0].seg_id = right.operations[0].seg_id = bg::segment_identifier(0, -1, -1, 0);
    line const forward{{-1, 1}, {4, 1}}, backward{{4, 1}, {-1, 1}};
    auto const order = bg::detail::turns::make_less_seg_fraction_other_op(position_strategy());
    // A one-segment open line must not wrap its last endpoint to its first.
    BOOST_TEST(order.on_geometry(forward)(left, right));
    BOOST_TEST(order.on_geometry(backward)(right, left));
    int const original = std::fegetround();
    for (int mode : {FE_TONEAREST})
    {
        if (std::fesetround(mode) != 0) continue;
        native_positions<double>();
        native_positions<std::int64_t>();
    }
    BOOST_TEST_EQ(std::fesetround(original), 0);
    point const third = ec::construct_intersection(point(0, 0), point(1, 1), point(0, 1), point(1, -1));
    auto retained = third;
    BOOST_TEST_EQ(ec::segment_position<double>::classify(point(0, 0), point(1, 1), retained), 2);
    BOOST_TEST_EQ(ec::segment_position<double>::apply(point(0, 0), point(1, 1), retained, point(0.5, 0.5)), -1);
    bg::segment_ratio<int> const legacy(1, 3);
    BOOST_TEST(legacy.numerator() == 1 && legacy.denominator() == 3 && legacy.in_segment());
    return boost::report_errors();
}
