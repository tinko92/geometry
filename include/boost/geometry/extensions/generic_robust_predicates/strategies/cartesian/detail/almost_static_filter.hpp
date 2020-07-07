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

#include <boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/expression_tree.hpp>

namespace boost { namespace geometry
{

namespace detail { namespace generic_robust_predicates
{

template
<
    typename Expression,
    typename Real,
    template <typename, typename> StaticFilter
>
class almost_static_filter
{
private:
    using static_filter = StaticFilter<Expression, Real>;
    using expression_max_leaf = max_leaf<Expression>;
    using index_set = boost::mp11::mp_iota<expression_max_leaf>;
    std::array<Real, 2 * expression_max_leaf::value> extrema;
    static_filter filter;

    template <std::size_t N, typename ...Reals>
    inline void update_extremum(const Reals&... args)
    {
        Real arg = get_nth_real<N, Real>(args);
        extrema[N] = std::max(extrema[N], arg);
        extrema[N + expression_max_leaf::value] =
            std::min(extrema[N + expression_max_leaf::value], arg);
    }

    template <std::size_t N, typename ...Reals>
    inline void update_extremum_check(bool& check, const Reals&... args)
    {
        Real arg = get_nth_real<N, Real>(args);
        if (arg > extrema[N])
        {
            check = true;
            extrema[N] = arg;
        }
        else if (arg < extrema[N + expression_max_leaf::value])
        {
            check = true;
            extrema[N + expression_max_leaf::value] = arg;
        }
    }
public:
    inline almost_static_filter() {}
    template <typename ...Reals>
    int operator()(const Reals&... args)
    {
        return filter(args...);
    }

    template <typename ...Reals>
    inline void update_extrema(const Reals&... args)
    {
        update_extremum<index_set>(args...)...;
    }

    template <typename ...Reals>
    inline bool update_extrema_check(const Reals&... args)
    {
        bool changed = false;
        update_extremum_check<index_set>(changed, args...)...;
        return changed;
    }

    inline void update_filter()
    {
        filter = static_filter(std::get<index_set>(extrema)...);
    }
};

}} // namespace detail::generic_robust_predicates

}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_ALMOST_STATIC_FILTER_HPP
