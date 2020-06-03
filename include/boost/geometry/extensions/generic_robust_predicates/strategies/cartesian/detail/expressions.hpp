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
#include <boost/mp11/map.hpp>
#include <boost/mp11/set.hpp>

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

template<typename T> using inc = boost::mp11::mp_int<T::value + 1>;

template<typename L> using app_zero_b = boost::mp11::mp_push_front<L, boost::mp11::mp_int<0>>;
template<typename L> using app_zero_f = boost::mp11::mp_push_back<L, boost::mp11::mp_int<0>>;
template<typename L> using mult_by_1_p_eps = boost::mp11::mp_transform<boost::mp11::mp_plus, app_zero_b<L>, app_zero_f<L>>;

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
        boost::mp11::mp_push_back<L, boost::mp11::mp_plus<boost::mp11::mp_front<L1>, boost::mp11::mp_front<L2>>>>::type;
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

template<typename L1, typename L2> using coeff_merge = typename coeff_merge_impl<L1, L2, boost::mp11::mp_list<>>::type;

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
        boost::mp11::mp_pop_front<L1>,
        boost::mp11::mp_pop_front<L2>,
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
        boost::mp11::mp_list<K, coeff_merge<
            typename mp_map_at_second_or_void<M1, K>::type,
	    typename mp_map_at_second_or_void<M2, K>::type>>>;
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


//TODO: magnitude_expressions

}} // namespace detail::generic_robust_predicates

}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_EXPRESSIONS_HPP
