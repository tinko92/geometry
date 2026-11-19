// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.

#include <boost/geometry.hpp>
#include <boost/geometry/strategies/relate/cartesian_exact.hpp>
#include <boost/core/lightweight_test.hpp>
#include <cstdint>

namespace bg = boost::geometry;

template <typename C, typename Policy>
void check()
{
    using point = bg::model::d2::point_xy<C>;
    using line = bg::model::linestring<point>;
    using polygon = bg::model::polygon<point>;
    using collection = bg::model::geometry_collection<boost::variant<point, line, polygon>>;
    polygon square;
    bg::read_wkt("POLYGON((0 0,0 2,2 2,2 0,0 0))", square);
    collection const first{square, line{{3, 0}, {3, 2}}};
    collection const second{line{{0, 1}, {4, 1}}};
    bg::strategies::relate::cartesian_exact<C, Policy> strategy;
    // Exercise both strategy-preserving L/A splits. This is a routing control,
    // not an independent GC oracle: the native and rational baselines omit
    // B/I at (2,1) here (1F2F01102 instead of 1F2001102).
    BOOST_TEST_EQ(bg::relation(first, second, strategy).str(), bg::relation(first, second).str());
    BOOST_TEST_EQ(bg::relation(second, first, strategy).str(), bg::relation(second, first).str());
}

int main()
{
    check<double, bg::strategies::relate::approximate_output>();
    check<double, bg::strategies::relate::eager_nearest_output>();
    check<double, bg::strategies::relate::deferred_nearest_output>();
    check<std::int64_t, bg::strategies::relate::approximate_output>();
    return boost::report_errors();
}
