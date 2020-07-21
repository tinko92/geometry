// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2020 Tinko Bartels, Berlin, Germany.

// Contributed and/or modified by Tinko Bartels,
//   as part of Google Summer of Code 2020 program.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_SIGNS_ONLY_FILTER_HPP
#define BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_SIGNS_ONLY_FILTER_HPP

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

template <typename Expression> using is_sign_exact =
    boost::mp11::mp_bool<Expression::sign_exact>;

template
<
    typename Exacts,
    typename All,
    typename Remaining,
    typename SignArr,
    typename ApproxArr,
    operator_types Op,
    bool LeftExact,
    bool RightExact,
    typename InputArr
>
struct deduce_sign_binary_impl;

template
<
    typename Exacts,
    typename All,
    typename Remaining,
    typename SignArr,
    typename ApproxArr,
    typename InputArr
>
struct deduce_sign_binary_impl
    <
        Exacts,
        All,
        Remaining,
        SignArr,
        ApproxArr,
        operator_types::product,
        true,
        true,
        InputArr
    >
{
    static inline void apply(SignArr& signs,
                             const ApproxArr& approx,
                             const InputArr& input)
    {
        using node = boost::mp11::mp_front<Remaining>;
        int& out = signs[boost::mp11::mp_find<All, node>::value];
        const auto& l = get_approx
                <
                    Exacts,
                    typename node::left,
                    typename ApproxArr::value_type
                >(approx, input);
        const auto& r = get_approx
                <
                    Exacts,
                    typename node::right,
                    typename ApproxArr::value_type
                >(approx, input);
        out =   ( l > 0 ? 1 : ( l < 0 ? -1 : 0 ) )
              * ( r > 0 ? 1 : ( r < 0 ? -1 : 0 ) );
    }
};

template
<
    typename Exacts,
    typename All,
    typename Remaining,
    typename SignArr,
    typename ApproxArr,
    typename InputArr
>
struct deduce_sign_binary_impl
    <
        Exacts,
        All,
        Remaining,
        SignArr,
        ApproxArr,
        operator_types::product,
        true,
        false,
        InputArr
    >
{
    static inline void apply(SignArr& signs,
                             const ApproxArr& approx,
                             const InputArr& input)
    {
        using node = boost::mp11::mp_front<Remaining>;
        int& out = signs[boost::mp11::mp_find<All, node>::value];
        const auto& l = get_approx
                <
                    Exacts,
                    typename node::left,
                    typename ApproxArr::value_type
                >(approx, input);
        const int& sr =
            signs[boost::mp11::mp_find<All, typename node::right>::value];
        if ( sr == sign_uncertain && l != 0 )
        {
            out = sign_uncertain;
        }
        else
        {
            out = ( l > 0 ? 1 : ( l < 0 ? -1 : 0 ) ) * sr;
        }
    }
};

template
<
    typename Exacts,
    typename All,
    typename Remaining,
    typename SignArr,
    typename ApproxArr,
    typename InputArr
>
struct deduce_sign_binary_impl
    <
        Exacts,
        All,
        Remaining,
        SignArr,
        ApproxArr,
        operator_types::product,
        false,
        true,
        InputArr
    >
{
    static inline void apply(SignArr& signs,
                             const ApproxArr& approx,
                             const InputArr& input)
    {
        using node = boost::mp11::mp_front<Remaining>;
        int& out = signs[boost::mp11::mp_find<All, node>::value];
        const auto& r = get_approx
                <
                    Exacts,
                    typename node::right,
                    typename ApproxArr::value_type
                >(approx, input);
        const int& sl =
            signs[boost::mp11::mp_find<All, typename node::left>::value];
        if ( sl == sign_uncertain && r != 0 )
        {
            out = sign_uncertain;
        }
        else
        {
            out = sl * ( r > 0 ? 1 : ( r < 0 ? -1 : 0 ) );
        }
    }
};

template
<
    typename Exacts,
    typename All,
    typename Remaining,
    typename SignArr,
    typename ApproxArr,
    typename InputArr
>
struct deduce_sign_binary_impl
    <
        Exacts,
        All,
        Remaining,
        SignArr,
        ApproxArr,
        operator_types::product,
        false,
        false,
        InputArr
    >
{
    static inline void apply(SignArr& signs,
                             const ApproxArr& approx,
                             const InputArr&)
    {
        using node = boost::mp11::mp_front<Remaining>;
        int& out = signs[boost::mp11::mp_find<All, node>::value];
        const int& sl =
            signs[boost::mp11::mp_find<All, typename node::left>::value];
        const int& sr =
            signs[boost::mp11::mp_find<All, typename node::right>::value];
        if ( sl == 0 || sr == 0 )
        {
            out = 0;
        }
        else if ( sl == sign_uncertain || sr == sign_uncertain )
        {
            out = sign_uncertain;
        }
        else
        {
            out = sl * sr;
        }
    }
};

template
<
    typename Exacts,
    typename All,
    typename Remaining,
    typename SignArr,
    typename ApproxArr,
    typename InputArr
>
struct deduce_sign_binary_impl
    <
        Exacts,
        All,
        Remaining,
        SignArr,
        ApproxArr,
        operator_types::sum,
        true,
        true,
        InputArr
    >
{
    static inline void apply(SignArr& signs,
                             const ApproxArr& approx,
                             const InputArr& input)
    {
        using node = boost::mp11::mp_front<Remaining>;
        int& out = signs[boost::mp11::mp_find<All, node>::value];
        const auto& l = get_approx
                <
                    Exacts,
                    typename node::left,
                    typename ApproxArr::value_type
                >(approx, input);
        const auto& r = get_approx
                <
                    Exacts,
                    typename node::right,
                    typename ApproxArr::value_type
                >(approx, input);
        if ( (l > 0 && r >= 0) || (l >= 0 && r > 0) )
        {
            out = 1;
        }
        else if ( (l < 0 && r <= 0) || (l <= 0 && r < 0) )
        {
            out = -1;
        }
        else if ( l == 0 && r == 0 )
        {
            out = 0;
        }
        else
        {
            out = sign_uncertain;
        }
    }
};

template
<
    typename Exacts,
    typename All,
    typename Remaining,
    typename SignArr,
    typename ApproxArr,
    typename InputArr
>
struct deduce_sign_binary_impl
    <
        Exacts,
        All,
        Remaining,
        SignArr,
        ApproxArr,
        operator_types::sum,
        true,
        false,
        InputArr
    >
{
    static inline void apply(SignArr& signs,
                             const ApproxArr& approx,
                             const InputArr& input)
    {
        using node = boost::mp11::mp_front<Remaining>;
        int& out = signs[boost::mp11::mp_find<All, node>::value];
        const auto& l = get_approx
                <
                    Exacts,
                    typename node::left,
                    typename ApproxArr::value_type
                >(approx, input);
        const int& sr =
            signs[boost::mp11::mp_find<All, typename node::right>::value];
        if ( sr == sign_uncertain )
        {
            out = sign_uncertain;
        }
        else if ( sr == 1 && l >= 0 )
        {
            out = 1;
        }
        else if ( sr == -1 && l <= 0 )
        {
            out = -1;
        }
        else if ( sr == 0 && l == 0 )
        {
            out = 0;
        }
        else
        {
            out = sign_uncertain;
        }
    }
};

template
<
    typename Exacts,
    typename All,
    typename Remaining,
    typename SignArr,
    typename ApproxArr,
    typename InputArr
>
struct deduce_sign_binary_impl
    <
        Exacts,
        All,
        Remaining,
        SignArr,
        ApproxArr,
        operator_types::sum,
        false,
        true,
        InputArr
    >
{
    static inline void apply(SignArr& signs,
                             const ApproxArr& approx,
                             const InputArr& input)
    {
        using node = boost::mp11::mp_front<Remaining>;
        int& out = signs[boost::mp11::mp_find<All, node>::value];
        const auto& r = get_approx
                <
                    Exacts,
                    typename node::right,
                    typename ApproxArr::value_type
                >(approx, input);
        const int& sl =
            signs[boost::mp11::mp_find<All, typename node::left>::value];
        if ( sl == sign_uncertain )
        {
            out = sign_uncertain;
        }
        else if ( sl == 1 && r >= 0 )
        {
            out = 1;
        }
        else if ( sl == -1 && r <= 0 )
        {
            out = -1;
        }
        else if ( sl == 0 && r == 0 )
        {
            out = 0;
        }
        else
        {
            out = sign_uncertain;
        }
    }
};

template
<
    typename Exacts,
    typename All,
    typename Remaining,
    typename SignArr,
    typename ApproxArr,
    typename InputArr
>
struct deduce_sign_binary_impl
    <
        Exacts,
        All,
        Remaining,
        SignArr,
        ApproxArr,
        operator_types::sum,
        false,
        false,
        InputArr
    >
{
    static inline void apply(SignArr& signs,
                             const ApproxArr& approx,
                             const InputArr&)
    {
        using node = boost::mp11::mp_front<Remaining>;
        int& out = signs[boost::mp11::mp_find<All, node>::value];
        const int& sl =
            signs[boost::mp11::mp_find<All, typename node::left>::value];
        const int& sr =
            signs[boost::mp11::mp_find<All, typename node::right>::value];
        if ( sl == 0 && sr == 0)
        {
            out = 0;
        }
        else if ( sl == sign_uncertain || sr == sign_uncertain || sl == -sr )
        {
            out = sign_uncertain;
        }
        else if ( sl == sr )
        {
            out = sl;
        }
        else
        {
            out = sl + sr;
        }
    }
};

template
<
    typename Exacts,
    typename All,
    typename Remaining,
    typename SignArr,
    typename ApproxArr,
    typename InputArr
>
struct deduce_sign_binary_impl
    <
        Exacts,
        All,
        Remaining,
        SignArr,
        ApproxArr,
        operator_types::difference,
        true,
        true,
        InputArr
    >
{
    static inline void apply(SignArr& signs,
                             const ApproxArr& approx,
                             const InputArr& input)
    {
        using node = boost::mp11::mp_front<Remaining>;
        int& out = signs[boost::mp11::mp_find<All, node>::value];
        const auto& l = get_approx
                <
                    Exacts,
                    typename node::left,
                    typename ApproxArr::value_type
                >(approx, input);
        const auto& r = get_approx
                <
                    Exacts,
                    typename node::right,
                    typename ApproxArr::value_type
                >(approx, input);
        if ( (l > 0 && r <= 0) || (l >= 0 && r < 0) )
        {
            out = 1;
        }
        else if ( (l < 0 && r >= 0) || (l <= 0 && r > 0) )
        {
            out = -1;
        }
        else if ( l == 0 && r == 0 )
        {
            out = 0;
        }
        else
        {
            out = sign_uncertain;
        }
    }
};

template
<
    typename Exacts,
    typename All,
    typename Remaining,
    typename SignArr,
    typename ApproxArr,
    typename InputArr
>
struct deduce_sign_binary_impl
    <
        Exacts,
        All,
        Remaining,
        SignArr,
        ApproxArr,
        operator_types::difference,
        true,
        false,
        InputArr
    >
{
    static inline void apply(SignArr& signs,
                             const ApproxArr& approx,
                             const InputArr& input)
    {
        using node = boost::mp11::mp_front<Remaining>;
        int& out = signs[boost::mp11::mp_find<All, node>::value];
        const auto& l = get_approx
                <
                    Exacts,
                    typename node::left,
                    typename ApproxArr::value_type
                >(approx, input);
        const int& sr =
            signs[boost::mp11::mp_find<All, typename node::right>::value];
        if ( sr == sign_uncertain )
        {
            out = sign_uncertain;
        }
        else if ( sr == -1 && l >= 0 )
        {
            out = 1;
        }
        else if ( sr == 1 && l <= 0 )
        {
            out = -1;
        }
        else if ( sr == 0 && l == 0 )
        {
            out = 0;
        }
        else
        {
            out = sign_uncertain;
        }
    }
};

template
<
    typename Exacts,
    typename All,
    typename Remaining,
    typename SignArr,
    typename ApproxArr,
    typename InputArr
>
struct deduce_sign_binary_impl
    <
        Exacts,
        All,
        Remaining,
        SignArr,
        ApproxArr,
        operator_types::difference,
        false,
        true,
        InputArr
    >
{
    static inline void apply(SignArr& signs,
                             const ApproxArr& approx,
                             const InputArr& input)
    {
        using node = boost::mp11::mp_front<Remaining>;
        int& out = signs[boost::mp11::mp_find<All, node>::value];
        const auto& r = get_approx
                <
                    Exacts,
                    typename node::right,
                    typename ApproxArr::value_type
                >(approx, input);
        const int& sl =
            signs[boost::mp11::mp_find<All, typename node::left>::value];
        if ( sl == sign_uncertain )
        {
            out = sign_uncertain;
        }
        else if ( sl == 1 && r <= 0 )
        {
            out = 1;
        }
        else if ( sl == -1 && r >= 0 )
        {
            out = -1;
        }
        else if ( sl == 0 && r == 0 )
        {
            out = 0;
        }
        else
        {
            out = sign_uncertain;
        }
    }
};

template
<
    typename Exacts,
    typename All,
    typename Remaining,
    typename SignArr,
    typename ApproxArr,
    typename InputArr
>
struct deduce_sign_binary_impl
    <
        Exacts,
        All,
        Remaining,
        SignArr,
        ApproxArr,
        operator_types::difference,
        false,
        false,
        InputArr
    >
{
    static inline void apply(SignArr& signs,
                             const ApproxArr& approx,
                             const InputArr&)
    {
        using node = boost::mp11::mp_front<Remaining>;
        int& out = signs[boost::mp11::mp_find<All, node>::value];
        const int& sl =
            signs[boost::mp11::mp_find<All, typename node::left>::value];
        const int& sr =
            signs[boost::mp11::mp_find<All, typename node::right>::value];
        if ( sl == 0 && sr == 0 )
        {
            out = 0;
        }
        else if ( sl == sr || sl == sign_uncertain || sr == sign_uncertain )
        {
            out = sign_uncertain;
        }
        else if ( sl == 0 || sr == 0 )
        {
            out = sl - sr;
        }
        else
        {
            out = sl;
        }
    }
};

template
<
    typename Exacts,
    typename All,
    typename Remaining,
    typename SignArr,
    typename ApproxArr,
    operator_arities Arity,
    typename InputArr
>
struct deduce_sign_arity_helper_impl {};

template
<
    typename Exacts,
    typename All,
    typename Remaining,
    typename SignArr,
    typename ApproxArr,
    typename InputArr
>
struct deduce_sign_arity_helper_impl
    <
        Exacts,
        All,
        Remaining,
        SignArr,
        ApproxArr,
        operator_arities::binary,
        InputArr
    >
{
    static inline void apply(SignArr& s, const ApproxArr& a, const InputArr& i)
    {
        using node = boost::mp11::mp_front<Remaining>;
        deduce_sign_binary_impl
            <
                Exacts,
                All,
                Remaining,
                SignArr,
                ApproxArr,
                node::operator_type,
                node::left::sign_exact,
                node::right::sign_exact,
                InputArr
            >::apply(s, a, i);
    }
};

template
<
    typename Exacts,
    typename All,
    typename Remaining,
    typename SignArr,
    typename ApproxArr,
    std::size_t RemainingSize,
    typename InputArr
>
struct deduce_signs_remainder_impl
{
    static inline void apply(SignArr& signs,
                             const ApproxArr& approx,
                             const InputArr& input)
    {
        using node = boost::mp11::mp_front<Remaining>;
        deduce_sign_arity_helper_impl
            <
                Exacts,
                All,
                Remaining,
                SignArr,
                ApproxArr,
                node::operator_arity,
                InputArr
            >::apply(signs, approx, input);
        deduce_signs_remainder_impl
            <
                Exacts,
                All,
                boost::mp11::mp_pop_front<Remaining>,
                SignArr,
                ApproxArr,
                boost::mp11::mp_size<Remaining>::value - 1,
                InputArr
            >::apply(signs, approx, input);
    }
};

template
<
    typename Exacts,
    typename All,
    typename Remaining,
    typename SignArr,
    typename ApproxArr,
    typename InputArr
>
struct deduce_signs_remainder_impl
    <
        Exacts,
        All,
        Remaining,
        SignArr,
        ApproxArr,
        1,
        InputArr
    >
{
    static inline void apply(SignArr& signs,
                             const ApproxArr& approx,
                             const InputArr& input)
    {
        using node = boost::mp11::mp_front<Remaining>;
        deduce_sign_arity_helper_impl
            <
                Exacts,
                All,
                Remaining,
                SignArr,
                ApproxArr,
                node::operator_arity,
                InputArr
            >::apply(signs, approx, input);
    }
};

template
<
    typename Exacts,
    typename All,
    typename Remaining,
    typename SignArr,
    typename ApproxArr,
    typename InputArr
>
inline void deduce_signs_remainder(SignArr& signs,
                                  const ApproxArr& approx,
                                  const InputArr& input)
{
    deduce_signs_remainder_impl
        <
            Exacts,
            All,
            Remaining,
            SignArr,
            ApproxArr,
            boost::mp11::mp_size<Remaining>::value,
            InputArr
        >::apply(signs, approx, input);
}

template
<
    typename Exacts,
    typename All,
    typename Remaining,
    typename SignArr,
    typename ApproxArr,
    typename InputArr
>
inline void deduce_signs(SignArr& signs,
                         const ApproxArr& approx,
                         const InputArr& input)
{
    deduce_signs_remainder
        <
            Exacts,
            All,
            Remaining
        >(signs, approx, input);
}

template <typename Expression, typename Real>
struct signs_only_filter
{
private:
    using non_exact_signs_po =
        typename boost::mp11::mp_unique<post_order<Expression, is_sign_exact>>;
    using non_exact_signs =
        typename boost::mp11::mp_remove_if<non_exact_signs_po, is_sign_exact>;
    using stack = typename boost::mp11::mp_unique<post_order<Expression>>;
    using evals = typename boost::mp11::mp_remove_if<stack, is_leaf>;
    using evals_sign_exact =
        typename boost::mp11::mp_copy_if<evals, is_sign_exact>;
public:
    static constexpr bool stateful = false;
    static constexpr bool updates = false;
    using computations = evals_sign_exact;
    template <typename ...Reals>
    static inline int apply(const Reals&... args)
    {
        std::array<Real, sizeof...(args)> input {{ static_cast<Real>(args)... }};
        std::array<Real, boost::mp11::mp_size<evals_sign_exact>::value>
            results_sign_exact;
        approximate_interim
            <
                evals_sign_exact,
                evals_sign_exact,
                Real
            >(results_sign_exact, input);
        std::array<int, boost::mp11::mp_size<non_exact_signs>::value>
            remainder_signs;
        deduce_signs
            <
                evals_sign_exact,
                non_exact_signs,
                non_exact_signs
            >(remainder_signs, results_sign_exact, input);
        return remainder_signs[
            boost::mp11::mp_find<non_exact_signs, Expression>::value];
    }
};

}} // namespace detail::generic_robust_predicates

}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_SIGNS_ONLY_FILTER_HPP
