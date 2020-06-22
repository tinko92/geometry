// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2020 Tinko Bartels, Berlin, Germany.

// Contributed and/or modified by Tinko Bartels,
//   as part of Google Summer of Code 2020 program.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_EXPRESSION_TREE_HPP
#define BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_EXPRESSION_TREE_HPP

#include <cstddef>

#include <boost/mp11/integral.hpp>
#include <boost/mp11/list.hpp>

namespace boost { namespace geometry
{

namespace detail { namespace generic_robust_predicates
{

enum class operator_types {
    sum, difference, product, abs, no_op,
    two_sum_tail, two_difference_tail, two_product_tail, //check
    fast_two_sum_tail, fast_two_difference_tail, //check
    grow_expansion, grow_expansion_ze, //check
    grow_diff_expansion, grow_diff_expansion_ze, //check
    expansion_sum, expansion_sum_ze, //check
    expansion_diff, expansion_diff_ze, //check
    fast_expansion_sum, fast_expansion_sum_ze, //check
    fast_expansion_diff, fast_expansion_diff_ze, //check
    scale_expansion, scale_expansion_ze //check
};
enum class operator_arities { unary, binary };

constexpr int sign_uncertain = -2;

struct sum_error_type {};
struct product_error_type {};
struct no_error_type {};

template<typename... Children>
struct internal_node
{
    static constexpr bool is_leaf = false;
};

template<typename Left, typename Right>
struct internal_binary_node : internal_node<Left, Right>
{
    using left  = Left;
    using right = Right;
    static constexpr operator_arities operator_arity = operator_arities::binary;
};

template<typename Child>
struct internal_unary_node : internal_node<Child>
{
    using child = Child;
    static constexpr operator_arities operator_arity = operator_arities::unary;
};

template<typename Left, typename Right>
struct sum : public internal_binary_node<Left, Right>
{
    static constexpr bool sign_exact = Left::is_leaf && Right::is_leaf;
    static constexpr operator_types operator_type = operator_types::sum;
    using error_type = sum_error_type;
};

template<typename Left, typename Right>
struct difference : public internal_binary_node<Left, Right>
{
    static constexpr bool sign_exact = Left::is_leaf && Right::is_leaf;
    static constexpr operator_types operator_type = operator_types::difference;
    using error_type = sum_error_type;
};

template<typename Left, typename Right>
struct product : public internal_binary_node<Left, Right>
{
    static constexpr bool sign_exact = Left::sign_exact && Right::sign_exact;
    static constexpr operator_types operator_type = operator_types::product;
    using error_type = product_error_type;
};

template<typename Child>
struct abs : public internal_unary_node<Child>
{
    using error_type = no_error_type;
    static constexpr operator_types operator_type = operator_types::abs;
    static constexpr bool sign_exact = Child::sign_exact;
};

template<std::size_t Argn>
struct leaf
{
    static constexpr std::size_t argn = Argn;
    static constexpr bool is_leaf = true;
    static constexpr bool sign_exact = true;
    static constexpr operator_types operator_type = operator_types::no_op;

    static constexpr operator_arities operator_arity = operator_arities::unary;
};

template<
    typename In,
    typename Out,
    bool at_leaf = In::is_leaf,
    bool is_binary = In::operator_arity == operator_arities::binary
>
struct post_order_impl;

template<typename In, typename Out>
struct post_order_impl<In, Out, true, false>
{
    using type = boost::mp11::mp_push_back<Out, In>;
};

template<typename In, typename Out>
struct post_order_impl<In, Out, false, true>
{
    using leftl  = typename post_order_impl<typename In::left, boost::mp11::mp_list<>>::type;
    using rightl = typename post_order_impl<typename In::right, boost::mp11::mp_list<>>::type;
    using merged = boost::mp11::mp_append<Out, leftl, rightl>;
    using type   = boost::mp11::mp_push_back<merged, In>;
};

template<typename In, typename Out>
struct post_order_impl<In, Out, false, false>
{
    using childl  = typename post_order_impl<typename In::child, boost::mp11::mp_list<>>::type;
    using merged = boost::mp11::mp_append<Out, childl>;
    using type   = boost::mp11::mp_push_back<merged, In>;
};

template<typename In>
using post_order = typename post_order_impl<In, boost::mp11::mp_list<>>::type;

template<typename Node>
using is_leaf = boost::mp11::mp_bool<Node::is_leaf>;



using  _1 = leaf<1>;
using  _2 = leaf<2>;
using  _3 = leaf<3>;
using  _4 = leaf<4>;
using  _5 = leaf<5>;
using  _6 = leaf<6>;
using  _7 = leaf<7>;
using  _8 = leaf<8>;
using  _9 = leaf<9>;
using _10 = leaf<10>;
using _11 = leaf<11>;
using _12 = leaf<12>;

}} // namespace detail::generic_robust_predicates

}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_EXPRESSION_TREE_HPP
