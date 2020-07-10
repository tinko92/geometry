// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2020 Tinko Bartels, Berlin, Germany.

// Contributed and/or modified by Tinko Bartels,
//   as part of Google Summer of Code 2020 program.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_ALMOST_STATIC_FILTER_HPP
#define BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_ALMOST_STATIC_FILTER_HPP

#include <array>
#include <algorithm>
#include <limits>

#include <boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/expression_tree.hpp>

namespace boost { namespace geometry
{

namespace detail { namespace generic_robust_predicates
{

template
<
    std::size_t N,
    typename ExtremaArray,
    typename Real,
    typename ...Reals
>
inline void update_extremum_check(bool& check,
                                  ExtremaArray& extrema,
                                  const Reals&... args)
{
    Real arg = get_nth_real<void, N + 1, Real>(args...);
    if (arg > extrema[N])
    {
        check = true;
        extrema[N] = arg;
    }
    else if (arg < extrema[N + extrema.size() / 2])
    {
        check = true;
        extrema[N + extrema.size() / 2] = arg;
    }
}

template
<
    std::size_t N,
    typename Real,
    typename ExtremaArray,
    typename ...Reals
>
inline void update_extremum(ExtremaArray& extrema, const Reals&... args)
{
    Real arg = get_nth_real<void, N + 1, Real>(args...);
    extrema[N] = std::max(extrema[N], arg);
    extrema[N + extrema.size() / 2] =
        std::min(extrema[N + extrema.size() / 2], arg);
}

template <std::size_t I, typename Real>
struct update_extrema_check_impl
{
    template <typename ExtremaArray, Real, typename ...Reals>
    static inline void apply(bool& changed,
                             ExtremaArray& extrema,
                             const Reals&... args)
    {
        update_extremum_check<I, Real>(changed, extrema, args...);
        update_extrema_check_impl<I - 1, Real>
            ::apply(changed, extrema, args...);
    }
};

template <typename Real>
struct update_extrema_check_impl<0, Real>
{
    template <typename ExtremaArray, typename ...Reals>
    static inline void apply(bool& changed,
                             ExtremaArray& extrema,
                             const Reals&... args)
    {
        update_extremum_check<0, Real>(changed, extrema, args...);
    }
};

template <std::size_t I, typename Real>
struct update_extrema_impl
{
    template <typename ExtremaArray, typename ...Reals>
    static inline void apply(ExtremaArray& extrema, const Reals&... args)
    {
        update_extremum<I, Real>(extrema, args...);
        update_extrema_impl<I - 1, Real>::apply(extrema, args...);
    }
};

template <typename Real>
struct update_extrema_impl<0, Real>
{
    template <typename ExtremaArray, typename ...Reals>
    static inline void apply(ExtremaArray& extrema, const Reals&... args)
    {
            update_extremum<0, Real>(extrema, args...);
    }
};

template <typename Filter, std::size_t N, std::size_t End>
struct make_filter_impl
{
    template <typename ExtremaArray, typename ...Reals>
    static Filter apply(const ExtremaArray& extrema, const Reals&... args)
    {
        return make_filter_impl<Filter, N + 1, End>
            ::apply(extrema, args..., extrema[N]);
    }
};

template <typename Filter, std::size_t End>
struct make_filter_impl<Filter, End, End>
{
    template <typename ExtremaArray, typename ...Reals>
    static Filter apply(const ExtremaArray& extrema, const Reals&... args)
    {
        Filter f(args...);
        return f;
    }
};

template
<
    typename Expression,
    typename Real,
    typename StaticFilter
>
class almost_static_filter
{
private:
    using expression_max_argn = max_argn<Expression>;
    using extrema_array = std::array<Real, 2 * expression_max_argn::value>;
    extrema_array m_extrema;
    StaticFilter m_filter;

public:
    const StaticFilter& filter() const { return m_filter; }
    inline almost_static_filter()
    {
        std::fill(m_extrema.begin(),
                  m_extrema.begin() + m_extrema.size() / 2,
                  -std::numeric_limits<Real>::infinity());
        std::fill(m_extrema.begin() + m_extrema.size() / 2,
                  m_extrema.end(),
                  std::numeric_limits<Real>::infinity());
    }
    template <typename ...Reals>
    int apply(const Reals&... args) const
    {
        return m_filter.apply(args...);
    }

    template <typename ...Reals>
    inline void update_extrema(const Reals&... args)
    {
        update_extrema_impl
            <
                expression_max_argn::value - 1,
                Real
            >::apply(m_extrema, args...);
    }

    template <typename ...Reals>
    inline bool update_extrema_check(const Reals&... args)
    {
        bool changed = false;
        update_extrema_check_impl<expression_max_argn::value - 1, Real>
            ::apply(changed, m_extrema, args...);
        return changed;
    }

    inline void update_filter()
    {
        m_filter = make_filter_impl
                <
                    StaticFilter,
                    0,
                    2 * expression_max_argn::value
                >::apply(m_extrema);
    }
};

}} // namespace detail::generic_robust_predicates

}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_ALMOST_STATIC_FILTER_HPP
