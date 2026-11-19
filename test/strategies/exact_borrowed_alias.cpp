// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.

#include <sstream>

#include <boost/core/lightweight_test.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/strategies/relate/cartesian_exact.hpp>

namespace bg = boost::geometry;
using polygon = bg::model::polygon<bg::model::d2::point_xy<double>>;
using multi = bg::model::multi_polygon<polygon>;

template <typename G>
std::string signature(G const& g)
{
    std::ostringstream out;
    out << std::hexfloat << bg::wkt(g);
    return out.str();
}

int main()
{
    polygon triangle, crossing, disjoint;
    bg::read_wkt("POLYGON((0 0,0 3,9 0,0 0))", triangle);
    bg::read_wkt("POLYGON((1 -1,1 3,4 3,4 -1,1 -1))", crossing);
    bg::read_wkt("POLYGON((10 10,10 12,12 12,12 10,10 10))", disjoint);
    bg::strategies::relate::cartesian_exact<double> const strategy;

    multi aliased{triangle, disjoint};
    aliased.shrink_to_fit();
    multi separate = aliased;
    bg::union_(aliased, crossing, separate, strategy);
    BOOST_TEST_EQ(separate.size(), 4u); // Two appended results, not a single conversion.
    bg::union_(aliased, crossing, aliased, strategy);
    BOOST_TEST_EQ(signature(aliased), signature(separate));

    // Object identity alone is insufficient: the output can own an operand.
    multi parent{triangle};
    parent.shrink_to_fit();
    multi rhs{crossing, disjoint};
    separate = parent;
    bg::union_(parent.front(), rhs, separate, strategy);
    BOOST_TEST_EQ(separate.size(), 3u);
    bg::union_(parent.front(), rhs, parent, strategy);
    BOOST_TEST_EQ(signature(parent), signature(separate));

    return boost::report_errors();
}
