// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2020 Tinko Bartels, Berlin, Germany.

// Contributed and/or modified by Tinko Bartels,
//   as part of Google Summer of Code 2020 program.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_APPROXIMATE_HPP
#define BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_APPROXIMATE_HPP

#include <cstddef>
#include <cmath>
#include <array>

#include <boost/mp11/integral.hpp>
#include <boost/mp11/list.hpp>
#include <boost/mp11/algorithm.hpp>

#include <boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/expression_tree.hpp>

namespace boost { namespace geometry
{

namespace detail { namespace generic_robust_predicates
{


template<std::size_t N, typename Real, typename ...Reals>
inline Real get_nth_real(const Real& arg, const Reals&... args) {
    if constexpr(N == 1) return arg;
    else return get_nth_real<N - 1, Real>(args...);
}

template<typename All, typename Node, typename Real, typename Arr, typename ...Reals>
inline Real get_approx(Arr& interim_results, const Reals&...args) {
    using node = Node;
    if constexpr(is_leaf<node>::value) {
        return get_nth_real<node::argn, Real>(args...);
    } else {
        return interim_results[boost::mp11::mp_find<All, Node>::value];
    }
}

template<typename All, typename Remaining, typename Real, typename Arr, typename ...Reals>
inline void approximate_interim(Arr& interim_results, const Reals&... args) {
    if constexpr( !boost::mp11::mp_empty<Remaining>::value ) {
        using node = boost::mp11::mp_front<Remaining>;
        if constexpr(node::operator_type == operator_types::product)
            interim_results[boost::mp11::mp_find<All, node>::value] =
                  get_approx<All, typename node::left, Real>(interim_results, args...)
                * get_approx<All, typename node::right, Real>(interim_results, args...);
        else if constexpr(node::operator_type == operator_types::sum)
            interim_results[boost::mp11::mp_find<All, node>::value] =
                  get_approx<All, typename node::left, Real>(interim_results, args...)
                + get_approx<All, typename node::right, Real>(interim_results, args...);
        else if constexpr(node::operator_type == operator_types::difference)
            interim_results[boost::mp11::mp_find<All, node>::value] =
                  get_approx<All, typename node::left, Real>(interim_results, args...)
                - get_approx<All, typename node::right, Real>(interim_results, args...);
        else if constexpr(node::operator_type == operator_types::abs) {
            interim_results[boost::mp11::mp_find<All, node>::value] =
                std::abs(get_approx<All, typename node::child, Real>(interim_results, args...));
        }
        approximate_interim<All, boost::mp11::mp_pop_front<Remaining>, Real>(interim_results, args...);
    }
}

template<typename Expression, typename Real, typename ...Reals>
inline int approximate_sign(const Reals&... args) {
    using root = Expression;
    using stack = typename boost::mp11::mp_unique<post_order<Expression>>;
    using interim_evals = typename boost::mp11::mp_remove<boost::mp11::mp_remove_if<stack, is_leaf>, root>;
    std::array<Real, boost::mp11::mp_size<interim_evals>::value> interim_results;
    approximate_interim<interim_evals, interim_evals, Real>(interim_results, args...);
    if constexpr(root::operator_type == operator_types::product)
        return   get_approx<interim_evals, typename root::left, Real>(interim_results, args...)
                    * get_approx<interim_evals, typename root::right, Real>(interim_results, args...) > 0 ? 1 :
                       ( get_approx<interim_evals, typename root::left, Real>(interim_results, args...)
                               * get_approx<interim_evals, typename root::right, Real>(interim_results, args...) < 0 ? -1 : 0 );
    else if constexpr(root::operator_type == operator_types::sum)
        return   get_approx<interim_evals, typename root::left, Real>(interim_results, args...)
                       > -get_approx<interim_evals, typename root::right, Real>(interim_results, args...) ? 1 : 
                                ( get_approx<interim_evals, typename root::left, Real>(interim_results, args...)
                                < -get_approx<interim_evals, typename root::right, Real>(interim_results, args...) ? -1 : 0 );
    else if constexpr(root::operator_type == operator_types::difference)
                return   get_approx<interim_evals, typename root::left, Real>(interim_results, args...)
                       > get_approx<interim_evals, typename root::right, Real>(interim_results, args...) ? 1 :    
                                ( get_approx<interim_evals, typename root::left, Real>(interim_results, args...)
                                < get_approx<interim_evals, typename root::right, Real>(interim_results, args...) ? -1 : 0 );
}

template<typename Expression, typename Real, typename ...Reals>
inline double approximate_value(const Reals&... args) {
    using root = Expression;
    using stack = typename boost::mp11::mp_unique<post_order<Expression>>;
    using evals = typename boost::mp11::mp_remove_if<stack, is_leaf>;
    std::array<Real, boost::mp11::mp_size<evals>::value> results;
        approximate_interim<evals, evals, Real>(results, args...);
    return results.back();
}

}} // namespace detail::generic_robust_predicates

}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_APPROXIMATE_HPP
