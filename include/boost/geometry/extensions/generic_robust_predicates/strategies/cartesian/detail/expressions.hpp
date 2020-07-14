// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2020 Tinko Bartels, Berlin, Germany.

// Contributed and/or modified by Tinko Bartels,
//   as part of Google Summer of Code 2020 program.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_EXPRESSIONS_HPP
#define BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_EXPRESSIONS_HPP

#include <boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/expression_tree.hpp>

namespace boost { namespace geometry
{

namespace detail { namespace generic_robust_predicates
{

using orient2d = difference
        <
            product
                <
                    difference < _1, _5 >,
                    difference < _4, _6 >
                >,
            product
                <
                    difference < _2, _6 >,
                    difference < _3, _5 >
                >
        >;

struct incircle_helper
{
    using adx = difference<_1, _7>;
    using ady = difference<_2, _8>;
    using bdx = difference<_3, _7>;
    using bdy = difference<_4, _8>;
    using cdx = difference<_5, _7>;
    using cdy = difference<_6, _8>;
    using abdet = difference<product<adx, bdy>, product<bdx, ady>>;
    using bcdet = difference<product<bdx, cdy>, product<cdx, bdy>>;
    using cadet = difference<product<cdx, ady>, product<adx, cdy>>;
    using alift = sum<product<adx, adx>, product<ady, ady>>;
    using blift = sum<product<bdx, bdx>, product<bdy, bdy>>;
    using clift = sum<product<cdx, cdx>, product<cdy, cdy>>;
    using expression = sum
        <
            product
                <
                    alift,
                    bcdet
                >,
            sum
                <
                    product
                        <
                            blift,
                            cadet
                        >,
                    product
                        <
                            clift,
                            abdet
                        >
                >
        >;
};

using incircle = incircle_helper::expression;


}} // namespace detail::generic_robust_predicates

}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_EXPRESSIONS_HPP
