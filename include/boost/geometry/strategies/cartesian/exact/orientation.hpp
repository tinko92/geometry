// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Licensed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_GEOMETRY_STRATEGIES_CARTESIAN_EXACT_ORIENTATION_HPP
#define BOOST_GEOMETRY_STRATEGIES_CARTESIAN_EXACT_ORIENTATION_HPP

#include <array>
#include <cmath>
#include <cstddef>
#include <boost/geometry/strategies/cartesian/exact/arithmetic.hpp>

namespace boost { namespace geometry { namespace strategies { namespace detail
{
namespace exact_cartesian { namespace orientation_detail
{

// Error-free sums/differences and Grow-Expansion, following Shewchuk:
// https://www.cs.cmu.edu/~quake/robust.html .
// Requires binary64 evaluation without reassociation or implicit contraction.
inline std::array<double, 2> two_sum(double a, double b)
{
    using namespace arithmetic_detail;
    double const sum = rounded_add(a, b);
    double const virtual_b = rounded_sub(sum, a);
    double const virtual_a = rounded_sub(sum, virtual_b);
    return {{rounded_add(rounded_sub(a, virtual_a), rounded_sub(b, virtual_b)), sum}};
}

inline std::array<double, 2> two_diff(double a, double b)
{
    using namespace arithmetic_detail;
    double const difference = rounded_sub(a, b);
    double const virtual_b = rounded_sub(a, difference);
    double const virtual_a = rounded_add(difference, virtual_b);
    return {{rounded_add(rounded_sub(a, virtual_a), rounded_sub(virtual_b, b)), difference}};
}

}} // namespace exact_cartesian::orientation_detail
}}}} // namespace boost::geometry::strategies::detail

#endif
