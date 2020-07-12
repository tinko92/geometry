// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2020 Tinko Bartels, Berlin, Germany.

// Contributed and/or modified by Tinko Bartels,
//   as part of Google Summer of Code 2020 program.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_SIMPLE_ORIENT2D_FILTERS_HPP
#define BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_SIMPLE_ORIENT2D_FILTERS_HPP

#include <limits>

#include <boost/mp11/list.hpp>
#include <boost/mp11/algorithm.hpp>
#include <boost/mp11/set.hpp>

#include <boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/expression_tree.hpp>
#include <boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/expressions.hpp>
#include <boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/coefficient_list.hpp>
#include <boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/error_bound.hpp>

#include <boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/semi_static_filter.hpp>
#include <boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/almost_static_filter.hpp>
#include <boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/static_filter.hpp>

namespace boost { namespace geometry
{

namespace detail { namespace generic_robust_predicates
{

// The following error expressions are based on "Simple Floating-Point Filters 
// for the Two-Dimensional Orientation Problem" by Katsuhisa Ozaki, Florian
// Bunger, Takeshi Ogita, Shin'ichi Oishi, Siegfried M. Rump.
// See: https://doi.org/10.1007/s10543-015-0574-9
//
// They are applicable to the 2D orientation problem.

template <typename Real>
struct u_n : public static_constant_interface<Real>
{
    static constexpr Real value = std::numeric_limits<Real>::min();
};

//The following could be defined generally and statically if sqrt and floor
//were constexpr functions.
template <std::size_t MantissaLength>
struct small_phi {};

template <>
struct small_phi<64>
{
    static constexpr std::size_t value = 4294967294;
};

template <>
struct small_phi<53>
{
    static constexpr std::size_t value = 94906264;
};

template <>
struct small_phi<22>
{
    static constexpr std::size_t value = 4096;
};

template <typename Real>
struct theta : public static_constant_interface<Real>
{
private:
    static constexpr std::size_t sphi = small_phi
        <
            std::numeric_limits<Real>::digits
        >::value;
    static constexpr Real u = std::numeric_limits<Real>::epsilon() / 2;
public:
    static constexpr Real value = 3 * u - (sphi - 22) * u * u;
};

template <typename Real>
struct simple_orient2d_semi_static_error_expression_impl
{
private:
    using left  = orient2d::left;
    using right = orient2d::right;
    using magnitude = sum < abs < sum < left, right > >, u_n<Real> >;
public:
    using type = product < theta<Real> , magnitude >;
};

template
<
    typename Real
>
using simple_orient2d_semi_static = semi_static_filter
        <
            orient2d,
            Real,
            typename simple_orient2d_semi_static_error_expression_impl<Real>
                ::type
        >;

}} // namespace detail::generic_robust_predicates

}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_SIMPLE_ORIENT2D_FILTERS_HPP
