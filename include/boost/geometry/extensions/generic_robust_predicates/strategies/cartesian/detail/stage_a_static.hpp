// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2020 Tinko Bartels, Berlin, Germany.

// Contributed and/or modified by Tinko Bartels,
//   as part of Google Summer of Code 2020 program.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_STAGE_A_STATIC_HPP
#define BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_STAGE_A_STATIC_HPP

#include <array>

#include <boost/mp11/list.hpp>
#include <boost/mp11/algorithm.hpp>
#include <boost/mp11/map.hpp>
#include <boost/mp11/set.hpp>

#include <boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/stage_a.hpp>

namespace boost { namespace geometry
{

namespace detail { namespace generic_robust_predicates
{


template
<
    typename All,
    typename Remaining,
    typename Real,
    typename Arr,
    operator_types Op,
    typename ...Reals
>
struct maximize_abs_impl {};

template
<
    typename All,
    typename Remaining,
    typename Real,
    typename Arr,
    bool Empty,
    typename ...Reals
>
struct maximize_abs_remainder_impl
{
    static inline void apply(Arr& interim_results, const Reals&... args)
    {
        using node = boost::mp11::mp_front<Remaining>;
        maximize_abs_impl
            <
                All,
                Remaining,
                Real,
                Arr,
                node::operator_type,
                Reals...
            >::apply(interim_results, args...);
    }
};

template
<
    typename All,
    typename Remaining,
    typename Real,
    typename Arr,
    typename ...Reals
>
struct maximize_abs_remainder_impl
    <
        All,
        Remaining,
        Real,
        Arr,
        true,
        Reals...
    >
{
    static inline void apply(Arr& interim_results, const Reals&... args) {}
};

template
<
    typename All,
    typename Remaining,
    typename Real,
    typename Arr,
    typename ...Reals
>
inline void maximize_abs_remainder(Arr& interim_results, const Reals&... args)
{
    maximize_abs_remainder_impl
        <
            All,
            Remaining,
            Real,
            Arr,
            boost::mp11::mp_empty<Remaining>::value,
            Reals...
        >::apply(interim_results, args...);
}

template
<
    typename All,
    typename Node,
    typename Real,
    typename Arr,
    bool is_leaf,
    bool max,
    typename ...Reals
>
struct get_min_max_impl
{
    static inline Real apply(Arr& interim_results, const Reals&... args)
    {
        return interim_results[boost::mp11::mp_find<All, Node>::value];
    }
};

template
<
    typename All,
    typename Node,
    typename Real,
    typename Arr,
    typename ...Reals
>
struct get_min_max_impl<All, Node, Real, Arr, true, true, Reals...>
{
    static inline Real apply(Arr& interim_results, const Reals&... args)
    {
        return get_nth_real<Node::argn, Real>(args...);
    }
};

template
<
    typename All,
    typename Node,
    typename Real,
    typename Arr,
    typename ...Reals
>
struct get_min_max_impl<All, Node, Real, Arr, true, false, Reals...>
{
    static inline Real apply(Arr& interim_results, const Reals&... args)
    {
        return get_nth_real<sizeof...(Reals) / 2 + Node::argn, Real>(args...);
    }
};

template
<
    typename All,
    typename Node,
    typename Real,
    typename Arr,
    bool max,
    typename ...Reals
>
inline Real get_min_max(Arr& interim_results, const Reals&...args)
{
    return get_min_max_impl
        <
            All,
            Node,
            Real,
            Arr,
            is_leaf<Node>::value,
            max,
            Reals...
        >::apply(interim_results, args...);
}

template
<
    typename All,
    typename Remaining,
    typename Real,
    typename Arr,
    typename ...Reals
>
struct maximize_abs_impl
    <
        All,
        Remaining,
        Real,
        Arr,
        operator_types::product,
        Reals...
    >
{
    static inline void apply(Arr& interim_results, const Reals&... args)
    {
        using node = boost::mp11::mp_front<Remaining>;
        const Real l = std::max(
                std::abs(get_min_max
                    <
                        All,
                        typename node::left,
                        Real,
                        Arr,
                        false
                    >(interim_results, args...)),
                std::abs(get_min_max
                    <
                        All,
                        typename node::left,
                        Real,
                        Arr,
                        true
                    >(interim_results, args...)));
        const Real r = std::max(
                std::abs(get_min_max
                    <
                        All,
                        typename node::right,
                        Real,
                        Arr,
                        false
                    >(interim_results, args...)),
                std::abs(get_min_max
                    <
                        All,
                        typename node::right,
                        Real,
                        Arr,
                        true
                    >(interim_results, args...)));
        interim_results[boost::mp11::mp_find<All, node>::value] = l * r;
        maximize_abs_remainder
            <
                All,
                boost::mp11::mp_pop_front<Remaining>,
                Real
            >(interim_results, args...);
    }
};

template
<
    typename All,
    typename Remaining,
    typename Real,
    typename Arr,
    typename ...Reals
>
struct maximize_abs_impl
    <
        All,
        Remaining,
        Real,
        Arr,
        operator_types::sum,
        Reals...
    >
{
    static inline void apply(Arr& interim_results, const Reals&... args)
    {
        using node = boost::mp11::mp_front<Remaining>;
        Real maxsum = get_min_max
                    <
                        All,
                        typename node::left,
                        Real,
                        Arr,
                        true
                    >(interim_results, args...)
                    + get_min_max
                    <
                        All,
                        typename node::right,
                        Real,
                        Arr,
                        true
                    >(interim_results, args...);
        Real minsum = get_min_max
                    <
                        All,
                        typename node::left,
                        Real,
                        Arr,
                        false
                    >(interim_results, args...)
                    + get_min_max
                    <
                        All,
                        typename node::right,
                        Real,
                        Arr,
                        false
                    >(interim_results, args...);
        interim_results[boost::mp11::mp_find<All, node>::value] =
            std::max(std::abs(maxsum), std::abs(minsum));
        maximize_abs_remainder
            <
                All,
                boost::mp11::mp_pop_front<Remaining>,
                Real
            >(interim_results, args...);
    }
};

template
<
    typename All,
    typename Remaining,
    typename Real,
    typename Arr,
    typename ...Reals
>
struct maximize_abs_impl
    <
        All,
        Remaining,
        Real,
        Arr,
        operator_types::difference,
        Reals...
    >
{
    static inline void apply(Arr& interim_results, const Reals&... args)
    {
        using node = boost::mp11::mp_front<Remaining>;
        Real maxmin = get_min_max
                    <
                        All,
                        typename node::left,
                        Real,
                        Arr,
                        true
                    >(interim_results, args...)
                    - get_min_max
                    <
                        All,
                        typename node::right,
                        Real,
                        Arr,
                        false
                    >(interim_results, args...);
        Real minmax = get_min_max
                    <
                        All,
                        typename node::left,
                        Real,
                        Arr,
                        false
                    >(interim_results, args...)
                    + get_min_max
                    <
                        All,
                        typename node::right,
                        Real,
                        Arr,
                        true
                    >(interim_results, args...);
        interim_results[boost::mp11::mp_find<All, node>::value] =
            std::max(std::abs(maxmin), std::abs(minmax));

        maximize_abs_remainder
            <
                All,
                boost::mp11::mp_pop_front<Remaining>,
                Real
            >(interim_results, args...);
    }
};

template
<
    typename All,
    typename Remaining,
    typename Real,
    typename Arr,
    typename ...Reals
>
struct maximize_abs_impl
    <
        All,
        Remaining,
        Real,
        Arr,
        operator_types::abs,
        Reals...
    >
{
    static inline void apply(Arr& interim_results, const Reals&... args)
    {
        using node = boost::mp11::mp_front<Remaining>;
        interim_results[boost::mp11::mp_find<All, node>::value] =
            std::abs(get_approx
                <
                    All,
                    typename node::child,
                    Real
                >(interim_results, args...));
        maximize_abs_remainder
            <
                All,
                boost::mp11::mp_pop_front<Remaining>,
                Real
            >(interim_results, args...);
    }
};

template
<
    typename All,
    typename Remaining,
    typename Real,
    typename Arr,
    typename ...Reals
>
struct maximize_abs_impl
    <
        All,
        Remaining,
        Real,
        Arr,
        operator_types::no_op,
        Reals...
    >
{
    static inline void apply(Arr& interim_results, const Reals&... args)
    {
        maximize_abs_remainder
            <
                All,
                boost::mp11::mp_pop_front<Remaining>,
                Real
            >(interim_results, args...);
    }
};

template
<
    typename All,
    typename Remaining,
    typename Real,
    typename Arr,
    typename ...Reals
>
inline void maximize_abs(Arr& interim_results, const Reals&... args)
{
    maximize_abs_remainder
        <
            All,
            Remaining,
            Real
        >(interim_results, args...);
}

template
<
    typename Real,
    typename ErrorExpression,
    typename ErrorEvalStack,
    typename ArgumentN,
    typename MaxLeaf
>
struct error_bound_impl {};

template
<
    typename Real,
    typename ErrorExpression,
    typename ErrorEvalStack,
    typename MaxLeaf
>
struct error_bound_impl
    <
        Real,
        ErrorExpression,
        ErrorEvalStack,
        boost::mp11::mp_size_t<1>,
        MaxLeaf
    >
{
    template<typename Real1>
    static constexpr Real apply(const Real1& arg)
    {
        return arg;
    }
};

template
<
    typename Real,
    typename ErrorExpression,
    typename ErrorEvalStack,
    typename MaxLeaf
>
struct error_bound_impl
    <
        Real,
        ErrorExpression,
        ErrorEvalStack,
        boost::mp11::mp_size_t<MaxLeaf::value * 2>,
        MaxLeaf
    >
{
    template<typename ...Reals>
    static constexpr Real apply(const Reals&... args)
    {
        std::array<Real, boost::mp11::mp_size<ErrorEvalStack>::value> evals;
        maximize_abs
            <
                ErrorEvalStack,
                ErrorEvalStack,
                Real
            >(evals, args...);
        return
            get_approx<ErrorEvalStack, ErrorExpression, Real>(evals, args...);
    }
};

template <typename Expression, typename Real>
class stage_a_static
{
private:
    using base_predicate = stage_a<Expression, Real>;
    using evals = base_predicate::evals;
    using error_expression = base_predicate::error_expression;
    using error_eval_stack = base_predicate::error_eval_stack;
    using final_coeff = base_predicate::final_coeff;

public:
    const Real error_bound;
    template <typename ...Reals>
    inline stage_a_static(const Reals&... args)
        : error_bound(
              error_bound_impl
                <
                    Real,
                    error_expression,
                    error_eval_stack,
                    boost::mp11::mp_size_t<sizeof...(Reals)>,
                    max_leaf<root>
                >::apply(args...)
            * eval_eps_polynomial<Real, final_coeff>::value) {}

    template <typename ...Reals>
    inline int apply(const Reals&... args)
    {
        if (error_bound == 0)
        {
            return 0;
        }
        std::array<Real, boost::mp11::mp_size<evals>::value> results;
        approximate_interim<evals, evals, Real>(results, args...);
        const Real det = get_approx<evals, Expression, Real>(results, args...);
        if (det > error_bound)
        {
            return 1;
        }
        else if (det < -error_bound)
        {
            return -1;
        }
        else
        {
            return sign_uncertain;
        }
    }

    template <typename ...Reals>
    inline int operator()(const Reals&... args)
    {
        return apply(args...);
    }
};

}} // namespace detail::generic_robust_predicates

}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_STAGE_A_STATIC_HPP
