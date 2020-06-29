// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2020 Tinko Bartels, Berlin, Germany.

// Contributed and/or modified by Tinko Bartels,
//   as part of Google Summer of Code 2020 program.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_STAGE_D_HPP
#define BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_STAGE_D_HPP

#include <array>
#include <iterator>

#include <boost/mp11/list.hpp>
#include <boost/mp11/algorithm.hpp>

#include <boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/approximate.hpp>
#include <boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/expression_tree.hpp>
#include <boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/expansion_arithmetic.hpp>

namespace boost { namespace geometry
{

namespace detail { namespace generic_robust_predicates
{

template<
    typename Expression,
    operator_types Op = Expression::operator_type
>
struct expansion_size_impl
{};



template<
    typename Expression
>
struct expansion_size_impl<Expression, operator_types::no_op>
{
    static constexpr std::size_t value = 1;
};

template<
    typename Expression
>
struct expansion_size_impl<Expression, operator_types::sum>
{
private:
    static constexpr std::size_t left_size = expansion_size_impl<typename Expression::left>::value;
    static constexpr std::size_t right_size = expansion_size_impl<typename Expression::right>::value;
public:
    static constexpr std::size_t value = left_size + right_size;
};

template<
    typename Expression
>
struct expansion_size_impl<Expression, operator_types::difference>
{
private:
    static constexpr std::size_t left_size = expansion_size_impl<typename Expression::left>::value;
    static constexpr std::size_t right_size = expansion_size_impl<typename Expression::right>::value;
public:
    static constexpr std::size_t value = left_size + right_size;
};

template<
    typename Expression
>
struct expansion_size_impl<Expression, operator_types::product>
{
private:
    static constexpr std::size_t left_size = expansion_size_impl<typename Expression::left>::value;
    static constexpr std::size_t right_size = expansion_size_impl<typename Expression::right>::value;
public:
    static constexpr std::size_t value = 2 * left_size * right_size;
};

template<typename Expression> using expansion_size = boost::mp11::mp_size_t<expansion_size_impl<Expression>::value>;

template<
    operator_types Op,
    std::size_t LeftLength,
    std::size_t RightLength,
    bool Inplace,
    typename Iter
>
struct perform_op_impl {};

template<
    std::size_t LeftLength,
    std::size_t RightLength,
    bool Inplace,
    typename Iter
>
struct perform_op_impl<operator_types::sum, LeftLength, RightLength, Inplace, Iter>
{
    template<typename ...Args>
    static Iter apply(Args...args) {
        return expansion_plus<LeftLength, RightLength, Inplace>(args...);
    }
};

template<
    std::size_t LeftLength,
    std::size_t RightLength,
    bool Inplace,
    typename Iter
>
struct perform_op_impl<operator_types::difference, LeftLength, RightLength, Inplace, Iter>
{
    template<typename ...Args>
    static Iter apply(Args...args) {
        return expansion_minus<LeftLength, RightLength, Inplace>(args...);
    }
};

template<
    std::size_t LeftLength,
    std::size_t RightLength,
    bool Inplace,
    typename Iter
>
struct perform_op_impl<operator_types::product, LeftLength, RightLength, Inplace, Iter>
{
    template<typename ...Args>
    static Iter apply(Args...args) {
        return expansion_times<LeftLength, RightLength, Inplace>(args...);
    }
};


template<
    typename Evals,
    typename Eval,
    typename Sizes,
    typename AccumulatedSizes,
    typename Iter,
    typename Real,
    operator_types Op = Eval::operator_type,
    bool LeftLeaf = Eval::left::is_leaf,
    bool RightLeaf = Eval::right::is_leaf
>
struct eval_expansion_impl
{};

template<
    typename Evals,
    typename Eval,
    typename Sizes,
    typename AccumulatedSizes,
    typename Iter,
    typename Real,
    operator_types Op
>
struct eval_expansion_impl<Evals, Eval, Sizes, AccumulatedSizes, Iter, Real, Op, true, true>
{
private:
    using left = Eval::left;
    using right = Eval::right;
    using eval_index = boost::mp11::mp_find<Evals, Eval>;
    static constexpr std::size_t size = boost::mp11::mp_at<Sizes, eval_index>::value;
    static constexpr std::size_t start = boost::mp11::mp_at<AccumulatedSizes, eval_index>::value;
public:
    template<typename ...Reals>
    static Iter apply(Iter begin, Iter end, const Reals&... args) {
        Real left_val = get_nth_real<left::argn, Real>(args...);
        Real right_val = get_nth_real<right::argn, Real>(args...);
        return perform_op_impl<Op, 1, 1, false, Iter>::apply(left_val, right_val, begin + start, begin + start + size);
    }
};

template<   
    typename Evals,
    typename Eval,
    typename Sizes,
    typename AccumulatedSizes,
    typename Iter,
    typename Real,
    operator_types Op
>
struct eval_expansion_impl<Evals, Eval, Sizes, AccumulatedSizes, Iter, Real, Op, true, false>
{
private:
    using left = Eval::left;
    using right = Eval::right;
    using eval_index = boost::mp11::mp_find<Evals, Eval>;
    static constexpr std::size_t size = boost::mp11::mp_at<Sizes, eval_index>::value;
    static constexpr std::size_t start = boost::mp11::mp_at<AccumulatedSizes, eval_index>::value;
    using right_eval_index = boost::mp11::mp_find<Evals, right>;
    static constexpr std::size_t right_size = boost::mp11::mp_at<Sizes, right_eval_index>::value;
    static constexpr std::size_t right_start = boost::mp11::mp_at<AccumulatedSizes, right_eval_index>::value;
public:
    template<typename ...Reals>
    static Iter apply(Iter begin, Iter end, const Reals&... args) {
        Real left_val = get_nth_real<left::argn, Real>(args...);
        return perform_op_impl<Op, 1, right_size, false, Iter>::apply(
            left_val,
            begin + right_start,
            begin + right_start + right_size,
            begin + start,
            begin + start + size);
    }
};

template<
    typename Evals,
    typename Eval,
    typename Sizes,
    typename AccumulatedSizes,
    typename Iter,
    typename Real,
    operator_types Op
>
struct eval_expansion_impl<Evals, Eval, Sizes, AccumulatedSizes, Iter, Real, Op, false, true>
{
private:
    using left = Eval::left;
    using right = Eval::right;
    using eval_index = boost::mp11::mp_find<Evals, Eval>;
    static constexpr std::size_t size = boost::mp11::mp_at<Sizes, eval_index>::value;
    static constexpr std::size_t start = boost::mp11::mp_at<AccumulatedSizes, eval_index>::value;
    using left_eval_index = boost::mp11::mp_find<Evals, left>;
    static constexpr std::size_t left_size = boost::mp11::mp_at<Sizes, left_eval_index>::value;
    static constexpr std::size_t left_start = boost::mp11::mp_at<AccumulatedSizes, left_eval_index>::value;
public:
    template<typename ...Reals>
    static Iter apply(Iter begin, Iter end, const Reals&... args) {
        Real right_val = get_nth_real<right::argn, Real>(args...);
        return perform_op_impl<Op, left_size, 1, false, Iter>::apply(
            begin + left_start,
            begin + left_start + left_size,
            right_val,
            begin + start,
            begin + start + size);
    }
};

template<
    typename Evals,
    typename Eval,
    typename Sizes,
    typename AccumulatedSizes,
    typename Iter,
    typename Real,
    operator_types Op
>
struct eval_expansion_impl<Evals, Eval, Sizes, AccumulatedSizes, Iter, Real, Op, false, false>
{
private:
    using left = Eval::left;
    using right = Eval::right;
    using eval_index = boost::mp11::mp_find<Evals, Eval>;
    static constexpr std::size_t size = boost::mp11::mp_at<Sizes, eval_index>::value;
    static constexpr std::size_t start = boost::mp11::mp_at<AccumulatedSizes, eval_index>::value;
    using left_eval_index = boost::mp11::mp_find<Evals, left>;
    static constexpr std::size_t left_size = boost::mp11::mp_at<Sizes, left_eval_index>::value;
    static constexpr std::size_t left_start = boost::mp11::mp_at<AccumulatedSizes, left_eval_index>::value;
    using right_eval_index = boost::mp11::mp_find<Evals, right>;
    static constexpr std::size_t right_size = boost::mp11::mp_at<Sizes, right_eval_index>::value;
    static constexpr std::size_t right_start = boost::mp11::mp_at<AccumulatedSizes, right_eval_index>::value;
public:
    template<typename ...Reals>
    static Iter apply(Iter begin, Iter end, const Reals&... args) {
        return perform_op_impl<Op, left_size, right_size, false, Iter>::apply(
            begin + left_start,
            begin + left_start + left_size,
            begin + right_start,
            begin + right_start + right_size,
            begin + start,
            begin + start + size);
    }
};

template<
    typename Evals,
    typename RemainingEvals,
    typename Sizes,
    typename AccumulatedSizes,
    typename Iter,
    typename Real,
    std::size_t = boost::mp11::mp_size<RemainingEvals>::value
>
struct eval_expansions_impl
{
    template<typename ...Reals>
    static Iter apply(Iter begin, Iter end, const Reals&... args) {
        eval_expansion_impl
            <
                Evals,
                boost::mp11::mp_front<RemainingEvals>,
                Sizes,
                AccumulatedSizes,
                Iter,
                Real
            >::apply(begin, end, args...);
        return eval_expansions_impl
            <
                Evals,
                boost::mp11::mp_pop_front<RemainingEvals>,
                Sizes,
                AccumulatedSizes,
                Iter,
                Real
            >::apply(begin, end, args...);
    }
};

template<
    typename Evals,
    typename RemainingEvals,
    typename Sizes,
    typename AccumulatedSizes,
    typename Iter,
    typename Real
>
struct eval_expansions_impl<Evals, RemainingEvals, Sizes, AccumulatedSizes, Iter, Real, 1>
{
    template<typename ...Reals>
    static Iter apply(Iter begin, Iter end, const Reals&... args) {
        return eval_expansion_impl
            <
                Evals,
                boost::mp11::mp_front<RemainingEvals>,
                Sizes,
                AccumulatedSizes,
                Iter,
                Real
            >::apply(begin, end, args...);
    }
};

template<typename Expression, typename Real, typename ...Reals>
inline int stage_d(const Reals&... args) {
    using root = Expression;
    using stack = typename boost::mp11::mp_unique<post_order<Expression>>;
    using evals = typename boost::mp11::mp_remove_if<stack, is_leaf>;
    using sizes = boost::mp11::mp_transform<expansion_size, evals>;
    using accumulated_sizes = boost::mp11::mp_push_front<
        boost::mp11::mp_partial_sum<sizes, boost::mp11::mp_size_t<0>, boost::mp11::mp_plus>,
        boost::mp11::mp_size_t<0>
    >;

    using result_array = std::array<Real, boost::mp11::mp_back<accumulated_sizes>::value>;
    result_array results;
    
    auto final_exp_end = eval_expansions_impl<evals, evals, sizes, accumulated_sizes, decltype(results.begin()), Real>
        ::apply(results.begin(), results.end(), args...);

    
    //TODO: If the last value is not computed with zero elimination, we need to either search
    //      for the last non-zero entry (multiplate branches) or add them all up (multiple 
    //      float ops). Evaluate which one is faster (if the case of no zero elimination
    //      even matters at all).
    constexpr std::size_t final_exp_size = boost::mp11::mp_back<sizes>::value;
    auto is_zero = [](Real d) { return d == Real(0); };
    auto most_significant =
        std::find_if_not(
            results.crbegin(),
            results.crbegin() + final_exp_size, is_zero);
    if( most_significant == results.crbegin() + final_exp_size) return 0;
    else if( *most_significant > 0 ) return 1;
    else return -1;
}

}} // namespace detail::generic_robust_predicates

}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_STAGE_D_HPP
