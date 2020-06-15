// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2020 Tinko Bartels, Berlin, Germany.

// Contributed and/or modified by Tinko Bartels,
//   as part of Google Summer of Code 2020 program.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_STAGE_A_HPP
#define BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_STAGE_A_HPP

#include <array>

#include <boost/mp11/list.hpp>
#include <boost/mp11/algorithm.hpp>
#include <boost/mp11/map.hpp>
#include <boost/mp11/set.hpp>

#include <boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/expression_tree.hpp>
#include <boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/approximate.hpp>
#include <boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/coefficient_list.hpp>
#include <boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/error_bound.hpp>

namespace boost { namespace geometry
{

namespace detail { namespace generic_robust_predicates
{

template<typename Expression, typename Real, typename ...Reals>
inline int stage_a(const Reals&... args) {
    using root = Expression;
    using stack = typename boost::mp11::mp_unique<post_order<Expression>>;
    using evals = typename boost::mp11::mp_remove_if<stack, is_leaf>;
    using interim_evals = typename boost::mp11::mp_remove<boost::mp11::mp_remove_if<stack, is_leaf>, root>;
    using interim_errors = evals_error<interim_evals>;
    using final_children = add_children<
        boost::mp11::mp_second<boost::mp11::mp_map_find<interim_errors, typename root::left>>,
        boost::mp11::mp_second<boost::mp11::mp_map_find<interim_errors, typename root::right>>
    >;
    using final_children_ltp = error_map_list_to_product<final_children>;
    using final_children_ltp_abs = abs_all<final_children_ltp>;
    using final_children_sum_fold = error_map_sum_up<final_children_ltp_abs>;
    using final_coeff = coeff_round<div_by_1_m_eps<mult_by_1_p_eps<boost::mp11::mp_second<final_children_sum_fold>>>>;
    using error_expression = boost::mp11::mp_front<final_children_sum_fold>;
    using error_eval_stack = boost::mp11::mp_unique<post_order<error_expression>>;
    using error_eval_stack_remainder = boost::mp11::mp_set_difference<error_eval_stack, evals>;
    using all_evals = boost::mp11::mp_append<evals, error_eval_stack_remainder>;

    std::array<Real, boost::mp11::mp_size<all_evals>::value> results;
    approximate_interim<all_evals, all_evals, Real>(results, args...);

    const Real stage_a_bound = 
          eval_eps_polynomial<Real, final_coeff>::value 
        * get_approx<all_evals, error_expression, Real>(results, args...);
    const Real det = get_approx<all_evals, root, Real>(results, args...);
    if(det > stage_a_bound) return 1;
    if(det < -stage_a_bound) return -1;
    if(stage_a_bound == 0) return 0;
    return sign_uncertain;
}

}} // namespace detail::generic_robust_predicates

}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_STAGE_A_HPP
