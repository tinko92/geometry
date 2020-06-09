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
#include <cmath>
#include <limits>

#include <boost/mp11/integral.hpp>
#include <boost/mp11/list.hpp>
#include <boost/mp11/algorithm.hpp>
#include <boost/mp11/map.hpp>
#include <boost/mp11/set.hpp>

namespace boost { namespace geometry
{

namespace detail { namespace generic_robust_predicates
{

enum class operator_types { sum, difference, product, abs, no_op };
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

template<typename T> using inc = boost::mp11::mp_int<T::value + 1>;

template<typename L> using app_zero_b = boost::mp11::mp_push_front<L, boost::mp11::mp_int<0>>;
template<typename L> using app_zero_f = boost::mp11::mp_push_back<L, boost::mp11::mp_int<0>>;
template<typename L> using mult_by_1_p_eps = 
    boost::mp11::mp_transform<
        boost::mp11::mp_plus,
        app_zero_b<L>,
        app_zero_f<L>
    >;

template<typename L, typename N, typename done = boost::mp11::mp_bool<N::value == 0>>
struct mult_by_1_p_eps_pow_impl
{
private:
    using next = mult_by_1_p_eps<L>;
public:
    using type = typename mult_by_1_p_eps_pow_impl<next, boost::mp11::mp_int<N::value - 1>>::type;
};

template<typename L, typename N>
struct mult_by_1_p_eps_pow_impl<L, N, boost::mp11::mp_true>
{  
public: 
    using type = L;
};

template<typename L, typename N> using mult_by_1_p_eps_pow = typename mult_by_1_p_eps_pow_impl<L, N>::type;

template<typename L> using div_by_1_m_eps_helper = boost::mp11::mp_partial_sum<L, boost::mp11::mp_int<0>, boost::mp11::mp_plus>;
template<typename L> using div_by_1_m_eps = boost::mp11::mp_push_back
    <
        boost::mp11::mp_pop_back<div_by_1_m_eps_helper<L>>,
        inc<boost::mp11::mp_back<div_by_1_m_eps_helper<L>>>
    >;

template<typename LOV, typename Void = boost::mp11::mp_same<LOV, void>>
struct empty_or_void { using type = boost::mp11::mp_empty<LOV>; };

template<typename LOV>
struct empty_or_void<LOV, boost::mp11::mp_true> { using type = boost::mp11::mp_true; };


template
<
    typename L1,
    typename L2,
    typename L,
    typename L1_Empty = typename empty_or_void<L1>::type,
    typename L2_Empty = typename empty_or_void<L2>::type
>
struct coeff_merge_impl {};

template
<
    typename L1,
    typename L2,
    typename L
>
struct coeff_merge_impl<L1, L2, L, boost::mp11::mp_false, boost::mp11::mp_false>
{
    using type = typename coeff_merge_impl<
        boost::mp11::mp_pop_front<L1>,
        boost::mp11::mp_pop_front<L2>,
        boost::mp11::mp_push_back<
            L,
            boost::mp11::mp_plus<
                boost::mp11::mp_front<L1>,
                boost::mp11::mp_front<L2>
            >
        >
    >::type;
};

template
< 
    typename L1,
    typename L2,
    typename L
>   
struct coeff_merge_impl<L1, L2, L, boost::mp11::mp_true, boost::mp11::mp_true>
{
    using type = L;
};

template
< 
    typename L1,
    typename L2,
    typename L
>   
struct coeff_merge_impl<L1, L2, L, boost::mp11::mp_true, boost::mp11::mp_false>
{
    using type = boost::mp11::mp_append<L, L2>;
};

template
<
    typename L1,
    typename L2,
    typename L
>   
struct coeff_merge_impl<L1, L2, L, boost::mp11::mp_false, boost::mp11::mp_true>
{
    using type = boost::mp11::mp_append<L, L1>;
};

template<typename L1, typename L2> using coeff_merge =
    typename coeff_merge_impl<L1, L2, boost::mp11::mp_list<>>::type;

template<typename N>
struct log_2_floor_impl {
    using type = typename boost::mp11::mp_plus<typename log_2_floor_impl<boost::mp11::mp_int<N::value / 2>>::type, boost::mp11::mp_int<1>>;
};

template<>
struct log_2_floor_impl<boost::mp11::mp_int<1>> {
    using type = boost::mp11::mp_int<0>;
};

template<typename N> using log_2_floor = typename log_2_floor_impl<N>::type;

template<typename N>
struct log_2_ceil_impl {
private:
    using floor = log_2_floor<N>;
public:
    using type = boost::mp11::mp_int< (1 << floor::value) == N::value ? floor::value : floor::value + 1 >;
};

template<typename N> using log_2_ceil = typename log_2_ceil_impl<N>::type;

template
<
    typename L1,
    typename L2,
    typename L,
    typename L1_Empty = typename empty_or_void<L1>::type,
    typename L2_Empty = typename empty_or_void<L2>::type
>
struct coeff_max_impl {};

template
<
    typename L1,
    typename L2,
    typename L
>
struct coeff_max_impl<L1, L2, L, boost::mp11::mp_false, boost::mp11::mp_false>
{
    using type = typename coeff_max_impl<
        boost::mp11::mp_if<
            boost::mp11::mp_less<boost::mp11::mp_front<L1>, boost::mp11::mp_front<L2>>,
            boost::mp11::mp_list<>,
            boost::mp11::mp_pop_front<L1>
        >,
        boost::mp11::mp_if<
            boost::mp11::mp_less<boost::mp11::mp_front<L2>, boost::mp11::mp_front<L1>>,
            boost::mp11::mp_list<>,
            boost::mp11::mp_pop_front<L2>
        >,
        boost::mp11::mp_push_back<L, boost::mp11::mp_max<boost::mp11::mp_front<L1>, boost::mp11::mp_front<L2>>>>::type;
};

template
< 
    typename L1,
    typename L2,
    typename L
>   
struct coeff_max_impl<L1, L2, L, boost::mp11::mp_true, boost::mp11::mp_true>
{
    using type = L;
};

template
< 
    typename L1,
    typename L2,
    typename L
>   
struct coeff_max_impl<L1, L2, L, boost::mp11::mp_true, boost::mp11::mp_false>
{
    using type = boost::mp11::mp_append<L, L2>;
};

template
<
    typename L1,
    typename L2,
    typename L
>   
struct coeff_max_impl<L1, L2, L, boost::mp11::mp_false, boost::mp11::mp_true>
{
    using type = boost::mp11::mp_append<L, L1>;
};

template<typename L1, typename L2> using coeff_max = typename coeff_max_impl<L1, L2, boost::mp11::mp_list<>>::type;

template<typename S> using is_not_zero = boost::mp11::mp_bool<S::value != 0>;

template
<
    typename L, 
    std::size_t Tail_Size =
        boost::mp11::mp_size<L>::value - boost::mp11::mp_find_if<L, is_not_zero>::value
>
struct coeff_round_impl
{
private:
    using first_nz = boost::mp11::mp_find_if<L, is_not_zero>;
    using tail = boost::mp11::mp_erase_c<L, 0, first_nz::value + 2>;
    using zero_tail = boost::mp11::mp_same<boost::mp11::mp_find_if<tail, is_not_zero>, boost::mp11::mp_size<tail>>;
    using head = boost::mp11::mp_erase_c<L, first_nz::value + 2, boost::mp11::mp_size<L>::value>;
    using minor = boost::mp11::mp_if<zero_tail, boost::mp11::mp_back<head>, inc<boost::mp11::mp_back<head>>>;
    using major = boost::mp11::mp_at<head, first_nz>;
    using major_rounded = boost::mp11::mp_int<1 << log_2_ceil<major>::value>;
    using minor_rounded = boost::mp11::mp_int<
        (minor::value / major_rounded::value) * major_rounded::value < minor::value ?
         (minor::value / major_rounded::value + 1) * major_rounded::value
       : (minor::value / major_rounded::value) * major_rounded::value >;
public:
    using type = boost::mp11::mp_push_back< boost::mp11::mp_pop_back<head>, minor_rounded >;
};

template<typename L>
struct coeff_round_impl<L, 0> { using type = L; };

template<typename L>
struct coeff_round_impl<L, 1> { using type = L; };

template<typename L> using coeff_round = typename coeff_round_impl<L>::type;

template<typename M, typename K, typename Contained = boost::mp11::mp_map_contains<M, K>>
struct mp_map_at_second_or_void
{
    using type = void;
};

template<typename M, typename K>
struct mp_map_at_second_or_void<M, K, boost::mp11::mp_true>
{
    using type = boost::mp11::mp_second<boost::mp11::mp_map_find<M, K>>;
};

template<typename M1, typename M2>
struct add_fold_operator
{
public:
    template<typename M, typename K> using fn = boost::mp11::mp_map_insert<
        M,
        boost::mp11::mp_list<K, 
            coeff_merge<
                typename mp_map_at_second_or_void<M1, K>::type,
                typename mp_map_at_second_or_void<M2, K>::type
            >
        >
    >;
};

template<typename M1, typename M2> using add_children = boost::mp11::mp_fold<
    boost::mp11::mp_set_union<boost::mp11::mp_map_keys<M1>, boost::mp11::mp_map_keys<M2>>,
    boost::mp11::mp_list<>,
    add_fold_operator<M1, M2>::template fn>;

template<typename L> using indexed = boost::mp11::mp_transform<
    boost::mp11::mp_list,
    boost::mp11::mp_iota<boost::mp11::mp_size<L>>,
    L>;

template<typename L> using strip_index = boost::mp11::mp_transform<boost::mp11::mp_second, L>;
template<typename IV1, typename IV2> using indexed_value_product =
    boost::mp11::mp_list<
	boost::mp11::mp_plus< boost::mp11::mp_first<IV1>, boost::mp11::mp_first<IV2>>,
	boost::mp11::mp_int< boost::mp11::mp_second<IV1>::value * boost::mp11::mp_second<IV2>::value >
    >;

template<typename IV1, typename IV2> using indexed_value_sum =
    boost::mp11::mp_list<
        boost::mp11::mp_first<IV1>,
        boost::mp11::mp_plus< boost::mp11::mp_second<IV1>, boost::mp11::mp_second<IV2> >
    >;

template<typename V1> struct add_nested
{
    template<typename K, typename V2> using fn = boost::mp11::mp_plus<V1, V2>;
};

template<typename M, typename IV> using list_product_fold =
    boost::mp11::mp_map_update_q<M, IV, add_nested<boost::mp11::mp_second<IV>>>;

template<typename IV1, typename IV2> using index_less =
    boost::mp11::mp_less<boost::mp11::mp_first<IV1>, boost::mp11::mp_first<IV2>>;

template<typename L1, typename L2> using list_product =
    strip_index<
        boost::mp11::mp_sort<
            boost::mp11::mp_fold<
                boost::mp11::mp_product<
                    indexed_value_product,
                    indexed<L1>, indexed<L2>
                >,
                boost::mp11::mp_list<>,
                list_product_fold
            >,
            index_less
        >
    >;

template
<
    typename Exp,
    typename LErr,
    typename RErr,
    typename Children_Empty = 
        boost::mp11::mp_and<
            typename empty_or_void<LErr>::type,
            typename empty_or_void<RErr>::type
        >
>
struct sum_err_impl
{
private:
    using children = add_children<LErr, RErr>;
public:
    using type = boost::mp11::mp_map_insert<
        children,
        boost::mp11::mp_list<Exp, boost::mp11::mp_list<boost::mp11::mp_int<1>>>
    >;
};

template
<
    typename Exp,
    typename LErr,
    typename RErr
>
struct sum_err_impl<Exp, LErr, RErr, boost::mp11::mp_true>
{
    using type = boost::mp11::mp_list<
        boost::mp11::mp_list<Exp, boost::mp11::mp_list<boost::mp11::mp_int<1>>>
    >;
};

template<typename Exp, typename LErr, typename RErr> using sum_err = typename sum_err_impl<Exp, LErr, RErr>::type;

template<typename L> using pad_second = boost::mp11::mp_list<
    boost::mp11::mp_front<L>,
    boost::mp11::mp_push_front<boost::mp11::mp_second<L>, boost::mp11::mp_int<0>>>;

template<typename L> using pop_front_second = boost::mp11::mp_list<
    boost::mp11::mp_front<L>,
    boost::mp11::mp_pop_front<boost::mp11::mp_second<L>>>;

template<typename K, typename V> using increment_first_of_second =
    boost::mp11::mp_transform_front<V, inc>;

template<typename KV1, typename KV2> using prod_entry_merge =
    boost::mp11::mp_list<
        boost::mp11::mp_list<
            boost::mp11::mp_first<KV1>,
            boost::mp11::mp_first<KV2>
        >,
        list_product<boost::mp11::mp_second<KV1>, boost::mp11::mp_second<KV2>>>;


template
<
    typename Exp,
    typename LErr,
    typename RErr
>
struct prod_children_impl
{
private:
    using left = typename Exp::left;
    using right = typename Exp::right;
    using padded_lerr = boost::mp11::mp_map_update<
        boost::mp11::mp_transform<pad_second, LErr>,
        boost::mp11::mp_list<left, boost::mp11::mp_list<boost::mp11::mp_int<1>>>,
        increment_first_of_second
    >;
    using padded_rerr = boost::mp11::mp_map_update<
        boost::mp11::mp_transform<pad_second, RErr>, 
        boost::mp11::mp_list<right, boost::mp11::mp_list<boost::mp11::mp_int<1>>>, 
        increment_first_of_second
    >;
    using prod = boost::mp11::mp_product<prod_entry_merge, padded_lerr, padded_rerr>;
    using stripped_prod = boost::mp11::mp_transform<pop_front_second, prod>;
public:
    using type = stripped_prod;
};

template<typename Exp, typename LErr, typename RErr> using prod_children =
    typename prod_children_impl<Exp, LErr, RErr>::type;

template
<
    typename Exp,
    typename LErr,
    typename RErr
>
struct product_err_impl
{
private:
    using children = prod_children<Exp, LErr, RErr>;
public:
    using type = boost::mp11::mp_map_insert<
        children,
        boost::mp11::mp_list<Exp, boost::mp11::mp_list<boost::mp11::mp_int<1>>>
    >;
};

template<typename Exp, typename LErr, typename RErr> using product_err = typename product_err_impl<Exp, LErr, RErr>::type;

template
<
    typename Map,
    typename Key,
    typename Contains = boost::mp11::mp_map_contains<Map, Key>
>
struct val_or_empty_list
{
    using type = boost::mp11::mp_second<boost::mp11::mp_map_find<Map, Key>>;
};

template
<
    typename Map,
    typename Key
>
struct val_or_empty_list<Map, Key, boost::mp11::mp_false>
{
    using type = boost::mp11::mp_list<>;
};

template
<
    typename Errors,
    typename Exp
>
struct error_fold_impl
{
private:
    using lerr = typename val_or_empty_list<Errors, typename Exp::left>::type;
    using rerr = typename val_or_empty_list<Errors, typename Exp::right>::type;
    using err = boost::mp11::mp_if<
        boost::mp11::mp_same<typename Exp::error_type, sum_error_type>,
        sum_err<Exp, lerr, rerr>,
        product_err<Exp, lerr, rerr>
    >;
public:
    using type = boost::mp11::mp_map_insert<Errors, boost::mp11::mp_list<Exp, err>>;
};


template<typename Errors, typename Exp> using error_fold = typename error_fold_impl<Errors, Exp>::type;

template<typename Evals> using evals_error = boost::mp11::mp_fold<Evals, boost::mp11::mp_list<>, error_fold>;

template<typename T> using is_mp_list = boost::mp11::mp_similar<boost::mp11::mp_list<>, T>;

template<typename KV>
struct list_to_product_impl
{
private:
    using key = boost::mp11::mp_front<KV>;
    using value = boost::mp11::mp_second<KV>;
    using multiplications = boost::mp11::mp_int<boost::mp11::mp_size<key>::value - 1>;
    using nvalue = mult_by_1_p_eps_pow<value, multiplications>;
    using nkey = boost::mp11::mp_fold<
        boost::mp11::mp_pop_front<key>,
        boost::mp11::mp_front<key>,
        product
    >;
public:
    using type = boost::mp11::mp_list<nkey, nvalue>;
};

template<typename KV> using list_to_product = typename list_to_product_impl<KV>::type;

template<typename KV, typename M>
struct error_map_insert_impl
{
private:
    using key = boost::mp11::mp_front<KV>;
    using value = boost::mp11::mp_second<KV>;
    using other_value = typename mp_map_at_second_or_void<M, key>::type;
    using merged_value = coeff_merge<value, other_value>;
    using nkv = boost::mp11::mp_list<key, merged_value>;
public:
    using type = boost::mp11::mp_map_replace<M, nkv>;
};

template<typename KV, typename M> using error_map_insert = typename error_map_insert_impl<KV, M>::type;
template<typename M, typename KV, typename KeyMPList = is_mp_list<boost::mp11::mp_front<KV>>>
struct error_map_list_to_product_fold_impl
{
    using type = error_map_insert<list_to_product<KV>, M>;
};

template<typename M, typename KV>
struct error_map_list_to_product_fold_impl<M, KV, boost::mp11::mp_false>
{
    using type = error_map_insert<KV, M>;
};

template<typename M, typename KV> using error_map_list_to_product_fold =
    typename error_map_list_to_product_fold_impl<M, KV>::type;

template<typename M>
struct error_map_list_to_product_impl
{
    using type = boost::mp11::mp_fold<
        M,
        boost::mp11::mp_list<>,
        error_map_list_to_product_fold
    >;
};

template<typename M> using error_map_list_to_product = typename error_map_list_to_product_impl<M>::type;

template<
    typename KV1,
    typename KV2
>
struct abs_sum_error_term_impl
{
private:
    using key1 = boost::mp11::mp_front<KV1>;
    using key2 = boost::mp11::mp_front<KV2>;
    using nkey = sum<abs<key1>, abs<key2>>;
    using val1 = boost::mp11::mp_second<KV2>;
    using val2 = boost::mp11::mp_second<KV2>;
    using mval = coeff_max<val1, val2>;
    using nval = mult_by_1_p_eps<mval>;
public:
    using type = boost::mp11::mp_list<nkey, nval>;
};

template<typename KV1, typename KV2> using abs_sum_error_term = typename abs_sum_error_term_impl<KV1, KV2>::type;

//TODO: The following could be probably optimized for potentially produce better error bounds in some cases 
//      if the error map is treated as a minheap by ordering of epsilon-polynomial.
template<typename M> using error_map_sum_up =
    boost::mp11::mp_fold<
        boost::mp11::mp_pop_front<M>,
        boost::mp11::mp_first<M>,
        abs_sum_error_term
    >;

template<
    typename Real,
    typename Exp
>
struct eps_pow
{
    static constexpr Real value = 
          std::numeric_limits<Real>::epsilon()/2.0
        * eps_pow<Real, boost::mp11::mp_size_t<Exp::value - 1>>::value;
};

template<typename Real> struct eps_pow<Real, boost::mp11::mp_size_t<0>>
{
    static constexpr Real value = 1.0;
};

template<
    typename Real,
    typename L,
    typename S = boost::mp11::mp_size<L>
>
struct eval_eps_polynomial
{
private:
    using last = boost::mp11::mp_back<L>;
    using s2last = boost::mp11::mp_back<boost::mp11::mp_pop_back<L>>;
public:
    static constexpr Real value =
          s2last::value * eps_pow<Real, boost::mp11::mp_size_t<S::value - 1>>::value
        + last::value * eps_pow<Real, S>::value;
};

template<
    typename Real,
    typename L
>
struct eval_eps_polynomial<Real, L, boost::mp11::mp_size_t<1>>
{
    static constexpr Real value = boost::mp11::mp_front<L>::value * std::numeric_limits<Real>::epsilon()/2.0;
};

template<
    typename Real,
    typename L                          
>
struct eval_eps_polynomial<Real, L, boost::mp11::mp_size_t<0>>
{
    static constexpr Real value = 0;
};

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
    using final_children_sum_fold = error_map_sum_up<final_children_ltp>;
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

#endif // BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_EXPRESSIONS_HPP
