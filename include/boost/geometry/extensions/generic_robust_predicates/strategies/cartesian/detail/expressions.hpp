// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2020 Tinko Bartels, Berlin, Germany.

// Contributed and/or modified by Tinko Bartels,
//   as part of Google Summer of Code 2020 program.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_EXPRESSIONS_HPP
#define BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_EXPRESSIONS_HPP

#include <cstddef>
#include <array>

#include <boost/mp11/integral.hpp>
#include <boost/mp11/list.hpp>
#include <boost/mp11/algorithm.hpp>

namespace boost { namespace geometry
{

namespace detail { namespace generic_robust_predicates
{

enum class operator_types { sum, difference, product, no_op };

template<typename Left, typename Right>
struct internal_node
{
    using left  = Left;
    using right = Right;
    static constexpr bool is_leaf = false;
};

template<typename Left, typename Right>
struct sum : public internal_node<Left, Right>
{
    static constexpr bool sign_exact = Left::is_leaf && Right::is_leaf;
    static constexpr operator_types operator_type = operator_types::sum;
};

template<typename Left, typename Right>
struct difference : public internal_node<Left, Right>
{
    static constexpr bool sign_exact = Left::is_leaf && Right::is_leaf;
    static constexpr operator_types operator_type = operator_types::difference;
};

template<typename Left, typename Right>
struct product : public internal_node<Left, Right>
{
    static constexpr bool sign_exact = Left::sign_exact && Right::sign_exact;
    static constexpr operator_types operator_type = operator_types::product;
};

template<std::size_t Argn>
struct leaf
{
    static constexpr std::size_t argn = Argn;
    static constexpr bool is_leaf = true;
    static constexpr bool sign_exact = true;
    static constexpr operator_types operator_type = operator_types::no_op;
};

template<typename In, typename Out, bool at_leaf = In::is_leaf>
struct post_order_impl;

template<typename In, typename Out>
struct post_order_impl<In, Out, true>
{
    using type = boost::mp11::mp_push_back<Out, In>;
};

template<typename In, typename Out>
struct post_order_impl<In, Out, false>
{
    using leftl  = typename post_order_impl<typename In::left, boost::mp11::mp_list<>>::type;
    using rightl = typename post_order_impl<typename In::right, boost::mp11::mp_list<>>::type;
    using merged = boost::mp11::mp_append<Out, leftl, rightl>;
    using type   = boost::mp11::mp_push_back<merged, In>;
};

template<typename In>
using post_order = typename post_order_impl<In, boost::mp11::mp_list<>>::type;

template<typename Node>
using is_leaf = boost::mp11::mp_bool<Node::is_leaf>;

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

#endif // BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_EXPRESSIONS_HPP
