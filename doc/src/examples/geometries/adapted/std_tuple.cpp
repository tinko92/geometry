// Boost.Geometry (aka GGL, Generic Geometry Library)
// QuickBook Example

// Copyright (c) 2011-2012 Barend Gehrels, Amsterdam, the Netherlands.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//[std_tuple
//` Shows how to use C++ tuples as points in Boost.Geometry

#include <iostream>
#include <tuple>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/adapted/std_tuple.hpp>
#include <boost/geometry/geometries/polygon.hpp>

BOOST_GEOMETRY_REGISTER_STD_TUPLE_CS(cs::cartesian)

int main()
{
    boost::geometry::model::polygon<std::tuple<double, double>> poly;
    poly.outer().push_back(std::make_tuple(1.0, 2.0));
    poly.outer().push_back(std::make_tuple(6.0, 4.0));
    poly.outer().push_back(std::make_tuple(5.0, 1.0));
    poly.outer().push_back(std::make_tuple(1.0, 2.0));

    std::cout << "Area: " << boost::geometry::area(poly) << std::endl;
    std::cout << "Contains (1.5, 2.5): "
        << std::boolalpha
        << boost::geometry::within(std::make_tuple(1.5, 2.5), poly)
        << std::endl;

    return 0;
}

//]

//[std_tuple_output
/*`
Output:
[pre
Area: 6.5
Contains (1.5, 2.5): false
]
*/
//]
