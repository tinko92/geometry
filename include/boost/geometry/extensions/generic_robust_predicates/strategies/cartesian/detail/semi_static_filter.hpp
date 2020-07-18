// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2020 Tinko Bartels, Berlin, Germany.

// Contributed and/or modified by Tinko Bartels,
//   as part of Google Summer of Code 2020 program.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_SEMI_STATIC_FILTER_HPP
#define BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_SEMI_STATIC_FILTER_HPP

#include <array>

#include <boost/mp11/list.hpp>
#include <boost/mp11/algorithm.hpp>
#include <boost/mp11/set.hpp>

#include <boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/expression_tree.hpp>
#include <boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/approximate.hpp>

namespace boost { namespace geometry
{

namespace detail { namespace generic_robust_predicates
{

template <typename Expression, typename Real, typename ErrorExpression>
struct semi_static_filter
{
private:
    using stack = typename boost::mp11::mp_unique<post_order<Expression>>;
    using evals = typename boost::mp11::mp_remove_if<stack, is_leaf>;
    using error_eval_stack = boost::mp11::mp_unique
        <
            post_order<ErrorExpression>
        >;
    using error_eval_stack_remainder = boost::mp11::mp_set_difference
        <
            error_eval_stack,
            evals
        >;
    using all_evals = boost::mp11::mp_append
        <
            evals,
            error_eval_stack_remainder
        >;
public:
    template <typename ...Reals>
    static inline int apply(const Reals&... args)
    {
        std::array<Real, sizeof...(Reals)> input {{ static_cast<Real>(args)... }};
        std::array<Real, boost::mp11::mp_size<all_evals>::value> results;
        approximate_interim<all_evals, all_evals, Real>(results, input);
        const Real error_bound =
            get_approx<all_evals, ErrorExpression, Real>(results, input);
        const Real det =
            get_approx<all_evals, Expression, Real>(results, input);
        if (det > error_bound)
        {
            return 1;
        }
        else if (det < -error_bound)
        {
            return -1;
        }
        else if (error_bound == 0 && det == 0)
        {
            return 0;
        }
        else
        {
            return sign_uncertain;
        }
    }
};

}} // namespace detail::generic_robust_predicates

}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_SEMI_STATIC_FILTER_HPP
