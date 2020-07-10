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

template
<
    typename Node,
    std::size_t N,
    typename Real,
    typename ...Reals
>
struct get_nth_real_impl
{
    static inline Real apply(const Real& arg, const Reals&... args)
    {
        return get_nth_real_impl<Node, N - 1, Reals...>::apply(args...);
    }
};

template
<
    typename Node,
    typename Real,
    typename ...Reals
>
struct get_nth_real_impl<Node, 1, Real, Reals...>
{
    static inline Real apply(const Real& arg, const Reals&... args)
    {
        return arg;
    }
};

template
<
    typename Node,
    typename Real,
    typename ...Reals
>
struct get_nth_real_impl<Node, 0, Real, Reals...>
{
    static inline Real apply(const Real& arg, const Reals&... args)
    {
        return Node::value;
    }
};

template
<
    typename Node,
    std::size_t N,
    typename Real,
    typename ...Reals
>
inline Real get_nth_real(const Real& arg, const Reals&... args)
{
    return get_nth_real_impl<Node, N, Real, Reals...>::apply(arg, args...);
}

template
<
    typename All,
    typename Node,
    typename Real,
    typename Arr,
    bool is_leaf,
    typename ...Reals
>
struct get_approx_impl
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
struct get_approx_impl<All, Node, Real, Arr, true, Reals...>
{
    static inline Real apply(Arr& interim_results, const Reals&... args)
    {
        return get_nth_real<Node, Node::argn, Real>(args...);
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
inline Real get_approx(Arr& interim_results, const Reals&...args)
{
    return get_approx_impl
        <
            All,
            Node,
            Real,
            Arr,
            is_leaf<Node>::value,
            Reals...
        >::apply(interim_results, args...);
}

template
<
    typename All,
    typename Remaining,
    typename Real,
    typename Arr,
    operator_types Op,
    typename ...Reals
>
struct approximate_interim_impl {};

template
<
    typename All,
    typename Remaining,
    typename Real,
    typename Arr,
    bool Empty,
    typename ...Reals
>
struct approximate_remainder_impl
{
    static inline void apply(Arr& interim_results, const Reals&... args)
    {
        using node = boost::mp11::mp_front<Remaining>;
        approximate_interim_impl
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
struct approximate_remainder_impl
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
inline void approximate_remainder(Arr& interim_results, const Reals&... args)
{
    approximate_remainder_impl
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
    typename Remaining,
    typename Real,
    typename Arr,
    typename ...Reals
>
struct approximate_interim_impl
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
        interim_results[boost::mp11::mp_find<All, node>::value] =
                  get_approx<All, typename node::left, Real>(interim_results,
                                                             args...)
                * get_approx<All, typename node::right, Real>(interim_results,
                                                              args...);
        approximate_remainder
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
struct approximate_interim_impl
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
        interim_results[boost::mp11::mp_find<All, node>::value] =
                  get_approx<All, typename node::left, Real>(interim_results,
                                                             args...)
                + get_approx<All, typename node::right, Real>(interim_results,
                                                              args...);
        approximate_remainder
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
struct approximate_interim_impl
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
        interim_results[boost::mp11::mp_find<All, node>::value] =
                  get_approx<All, typename node::left, Real>(interim_results,
                                                             args...)
                - get_approx<All, typename node::right, Real>(interim_results,
                                                              args...);
        approximate_remainder
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
struct approximate_interim_impl
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
        approximate_remainder
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
struct approximate_interim_impl
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
        approximate_remainder
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
inline void approximate_interim(Arr& interim_results, const Reals&... args)
{
    approximate_remainder
        <
            All,
            Remaining,
            Real
        >(interim_results, args...);
}

template<typename Expression, typename Real, typename ...Reals>
inline double approximate_value(const Reals&... args)
{
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
