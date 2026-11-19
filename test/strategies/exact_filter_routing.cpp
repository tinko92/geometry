// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.
#define BOOST_GEOMETRY_EXACT_FILTER_STATISTICS
#include <boost/core/lightweight_test.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/strategies/relate/cartesian_exact.hpp>
#include <array>
#include <cfenv>
#include <iostream>

namespace bg = boost::geometry;
namespace ec = bg::strategies::detail::exact_cartesian;
using point = bg::model::d2::point_xy<double>;
using polygon = bg::model::polygon<point>;
using polygons = bg::model::multi_polygon<polygon>;

template <typename Coordinate, typename Points>
void test_predicates(Points const& points)
{
    for (auto const& a : points)
    for (auto const& b : points)
    for (auto const& p : points)
    {
        ec::exact_number const ax(bg::get<0>(a)), ay(bg::get<1>(a));
        ec::exact_number const bx(bg::get<0>(b)), by(bg::get<1>(b));
        ec::exact_number const px(bg::get<0>(p)), py(bg::get<1>(p));
        ec::exact_number const dot = (bx - ax) * (px - bx) + (by - ay) * (py - by);
        ec::exact_number const determinant = (bx - ax) * (py - ay) - (by - ay) * (px - ax);
        ec::exact_number const squared = (bx - ax) * (bx - ax) + (by - ay) * (by - ay);
        BOOST_TEST_EQ(ec::static_direction<Coordinate>(a, b, p), dot.sign());
        auto const measure = ec::static_distance_measure<Coordinate>(a, b, p);
        BOOST_TEST_EQ(measure.side, determinant.sign());
        BOOST_TEST_EQ(measure.is_zero(), determinant == 0);
        auto const distance = ec::static_squared_distance<Coordinate>(a, b);
        BOOST_TEST_EQ(ec::compare_distance(distance, 0), squared.sign());
        auto oracle_value = [](auto const& pair)
        {
            return ec::exact_number(boost::multiprecision::cpp_int(pair.N))
                 / ec::exact_number(boost::multiprecision::cpp_int(pair.D));
        };
        BOOST_TEST((distance.native() ? oracle_value(distance.fixed_pair())
                                     : oracle_value(distance.precision_pair())) == squared);
    }
}

void test_predicate_oracles()
{
    std::array<point, 7> const native{{{0, 0}, {1, 0}, {2, 0}, {-1, 1}, {1, 1},
                                      {0.1, 0.3}, {0.3, 0.9}}};
    using integer_point = bg::model::d2::point_xy<std::int64_t>;
    auto const high = (std::numeric_limits<std::int64_t>::max)();
    auto const low = (std::numeric_limits<std::int64_t>::min)();
    std::array<integer_point, 5> const integral{{{low, low}, {high, high},
                                               {high, high - 1}, {0, 0}, {1, -1}}};
    std::array<ec::exact_point<double>, 4> const implicit{{ec::make_point(point(0, 0)),
        ec::make_point(point(1, 1)),
        ec::construct_intersection(ec::make_point(point(0, 0)), ec::make_point(point(1, 1)),
            ec::make_point(point(0, 1)), ec::make_point(point(1, -1))),
        ec::construct_intersection(ec::make_point(point(0, 0)), ec::make_point(point(1, 1)),
            ec::make_point(point(0, 2)), ec::make_point(point(1, -1)))}};
    int const rounding = std::fegetround();
    for (int mode : {FE_TONEAREST})
    {
        std::fesetround(mode);
        test_predicates<double>(native);
        test_predicates<std::int64_t>(integral);
        test_predicates<double>(implicit);
    }
    std::fesetround(rounding);
}

void statistics(char const* operation)
{
    auto const& s = ec::filtering_statistics();
    std::cout << operation << ": native=" << s.native_exact_values << " rational=" << s.native_rationals
              << " lazy=" << s.lazy_exact_values << " point=" << s.point_coordinate_values
              << " homogeneous=" << s.homogeneous_values << " ratio=" << s.ratio_values
              << " construction=" << s.construction_filters << '\n';
}

int main()
{
    test_predicate_oracles();
    bg::strategies::relate::cartesian_exact<double> strategy;
    polygon a, b, separate;
    bg::read_wkt("POLYGON((0 0,0 10,10 0,0 0))", a);
    bg::read_wkt("POLYGON((2 -1,-1 7,11 7,2 -1))", b);
    bg::read_wkt("POLYGON((20 20,20 30,30 20,20 20))", separate);
    polygons out;
    ec::filtering_statistics() = {};
    bg::union_(a, separate, out, strategy);
    statistics("disjoint union");
    BOOST_TEST_EQ(out.size(), 2u);
    BOOST_TEST_EQ(ec::filtering_statistics().native_rationals, 0u);
    BOOST_TEST_EQ(ec::filtering_statistics().point_coordinate_values, 0u);
    BOOST_TEST_EQ(ec::filtering_statistics().homogeneous_values, 0u);
    BOOST_TEST_EQ(ec::filtering_statistics().lazy_exact_values, 0u);

    for (int operation = 0; operation < 4; ++operation)
    {
        out.clear();
        ec::filtering_statistics() = {};
        if (operation == 0) { bg::union_(a, b, out, strategy); }
        if (operation == 1) { bg::intersection(a, b, out, strategy); }
        if (operation == 2) { bg::difference(a, b, out, strategy); }
        if (operation == 3) { bg::sym_difference(a, b, out, strategy); }
        statistics(operation == 0 ? "union" : operation == 1 ? "intersection"
            : operation == 2 ? "difference" : "symmetric difference");
        BOOST_TEST(!out.empty());
        BOOST_TEST_EQ(ec::filtering_statistics().native_rationals, 0u);
        BOOST_TEST_EQ(ec::filtering_statistics().point_coordinate_values, 0u);
        BOOST_TEST_EQ(ec::filtering_statistics().homogeneous_values, 0u);
        BOOST_TEST_EQ(ec::filtering_statistics().lazy_exact_values, 0u);
    }
    ec::filtering_statistics() = {};
    BOOST_TEST_EQ(bg::relation(a, b, strategy).str(), "212101212");
    statistics("relation");
    BOOST_TEST_EQ(ec::filtering_statistics().native_rationals, 0u);
    BOOST_TEST_EQ(ec::filtering_statistics().point_coordinate_values, 0u);
    BOOST_TEST_EQ(ec::filtering_statistics().homogeneous_values, 0u);
    BOOST_TEST_EQ(ec::filtering_statistics().lazy_exact_values, 0u);
    return boost::report_errors();
}
