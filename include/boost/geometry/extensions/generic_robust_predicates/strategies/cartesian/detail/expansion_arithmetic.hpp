// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2020 Tinko Bartels, Berlin, Germany.

// Contributed and/or modified by Tinko Bartels,
//   as part of Google Summer of Code 2020 program.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_EXPANSION_ARITHMETIC_HPP
#define BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_EXPANSION_ARITHMETIC_HPP

#include <cassert>
#include <cmath>
#include <array>
#include <algorithm>
#include <iterator>
#include <type_traits>

#include <boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/expression_tree.hpp>

namespace boost { namespace geometry
{

namespace detail { namespace generic_robust_predicates
{

namespace debug_expansion
{

struct abs_comp
{
    template<typename Real>
    constexpr bool operator()(Real a, Real b) const {
        return std::abs(a) < std::abs(b) || b == 0;
    }
};

unsigned long long r2pow2(unsigned long long num)
{
    if (!num)
        return 0;
    unsigned long long ret = 1;
    while (num >>= 1)
        ret <<= 1;
    return ret;
};

template<bool Negate, typename Real>
struct negate_impl {
    constexpr Real apply(Real a) {
        return a;
    }
};

template<typename Real>
struct negate_impl<true, Real> {
    constexpr Real apply(Real a) {
        return -a;
    }
};

template<bool Negate, typename Real>
constexpr Real negate(Real a) {
    return negate_impl<Negate, Real>::apply(a);
}

template<typename Real>
inline bool nonoverlapping(Real a, Real b) {
    int a_exp, b_exp, min_exp, max_exp;
    Real a_mant = std::frexp(a, &a_exp);
    Real b_mant = std::frexp(b, &b_exp);
    unsigned long long scale = 1ULL<<63;
    unsigned long long a_mantll = a_mant * scale;
    unsigned long long b_mantll = b_mant * scale;
    if(a_mantll == 0 || b_mantll == 0) return true;
    unsigned long long min_mantll, max_mantll;
    if(a_exp < b_exp) {
        min_exp = a_exp;
        max_exp = b_exp;
        min_mantll = a_mantll;
        max_mantll = b_mantll;
    } else {
        min_exp = b_exp;
        max_exp = a_exp;
        min_mantll = b_mantll;
        max_mantll = a_mantll;
    }
    int scale_down = max_exp - min_exp;
    unsigned long long min_mantll_sc = min_mantll >> scale_down;
    auto min_mantll_sc_rd = r2pow2(min_mantll_sc);
    return (max_mantll % (2*min_mantll_sc_rd)) == 0;
}

template<typename Real>
inline bool nonadjacent(Real a, Real b) {
    return nonoverlapping(a, b) && nonoverlapping(a, 2*b) && nonoverlapping(2*a, b);
}

template<typename Iter>
inline bool expansion_nonoverlapping(Iter begin, Iter end) {
    if(!std::is_sorted(begin, end, abs_comp{})) return false;
    while(std::next(begin) != end) {
        if( !nonoverlapping(*begin, *std::next(begin)) ) return false;
        ++begin;
    }
    return true;
}

template<typename Iter>
inline bool expansion_nonadjacent(Iter begin, Iter end) {
    if(!std::is_sorted(begin, end, abs_comp{})) return false;
    while(std::next(begin) != end) {
        if( !nonadjacent(*begin, *std::next(begin)) ) return false;
        ++begin;
    }
    return true;
}

template<typename Iter>
inline bool expansion_strongly_nonoverlapping(Iter begin, Iter end) {
    if(!std::is_sorted(begin, end, abs_comp{})) return false;
    int null;
    while(std::next(begin) != end) {
        if( !nonoverlapping(*begin, *std::next(begin)) ) return false;
        if( nonadjacent(*begin, *std::next(begin)) ) {}
        else {
            if( std::frexp(*begin, &null) != 0.5 || std::frexp(*std::next(begin), &null) != 0.5 ) return false;
            if( std::next(std::next(begin)) != end && !nonadjacent( *std::next(begin), *std::next(std::next(begin)) ) ) return false;
        }
        ++begin;
    }
    return true;
}

} // namespace debug_expansion

template<typename Real>
constexpr Real two_sum_tail(Real a, Real b, Real x) {
    Real b_virtual = x - a;
    Real a_virtual = x - b_virtual;
    Real b_rounded = b - b_virtual;
    Real a_rounded = a - a_virtual;
    Real y = a_rounded + b_rounded;
    assert(debug_expansion::nonadjacent(x, y));
    return y;
}

template<typename Real>
constexpr Real fast_two_sum_tail(Real a, Real b, Real x) {
    assert(std::abs(a) > std::abs(b) || a == 0);
    Real b_virtual = x - a;
    Real y = b - b_virtual;
    assert(debug_expansion::nonadjacent(x, y));
    return y;
}

template<typename Real>
constexpr Real two_difference_tail(Real a, Real b, Real x) {
    Real b_virtual = a - x;
    Real a_virtual = x + b_virtual;
    Real b_rounded = b_virtual - b;
    Real a_rounded = a - a_virtual;
    Real y = a_rounded + b_rounded;
    assert(debug_expansion::nonadjacent(x, y));
    return y;
}

template<typename Real>
constexpr Real fast_two_difference_tail(Real a, Real b, Real x) {
    assert(std::abs(a) > std::abs(b) || a == 0);
    Real b_virtual = a - x;
    Real y = b_virtual - b;
    assert(debug_expansion::nonadjacent(x, y));
    return y;
}

template<typename Real>
constexpr Real two_product_tail(Real a, Real b, Real x) {
    Real y = std::fma(a, b, -x);
    assert(debug_expansion::nonadjacent(x, y));
    return y;
}

template<typename InIter, typename OutIter, typename Real, bool NegateE = false, bool NegateB = false>
inline OutIter grow_expansion(InIter e_begin, InIter e_end, Real b, OutIter h_begin, OutIter h_end) {
    assert(std::distance(e_begin, e_end) + 1 == std::distance(h_begin, h_end));
    static_assert(std::is_same<typename std::iterator_traits<InIter>::value_type, Real>::value);
    static_assert(std::is_same<typename std::iterator_traits<OutIter>::value_type, Real>::value);
    assert(debug_expansion::expansion_nonoverlapping(e_begin, e_end));  
    Real Q = negate<NegateB>(b);
    auto h_it = h_begin;
    for(auto e_it = e_begin; e_it != e_end; ++e_it) {
        Real Q_new = negate<NegateE>(*e_it) + Q;
        *h_it = two_sum_tail(negate<NegateE>(*e_it), Q, Q_new);
        Q = Q_new;
        ++h_it;
    }
    *h_it = Q;
    assert(debug_expansion::expansion_nonoverlapping(h_begin, h_end));
    assert( !debug_expansion::expansion_nonadjacent(e_begin, e_end) || debug_expansion::expansion_nonadjacent(h_begin, h_end));
    return h_it;
}

template<typename InIter, typename OutIter, typename Real>
inline OutIter grow_expansion_difference(InIter e_begin, InIter e_end, Real b, OutIter h_begin, OutIter h_end) {
    return grow_expansion_sum<InIter, OutIter, Real, false, true>(e_begin, e_end, b, h_begin, h_end);
}

template<typename InIter, typename OutIter, typename Real, bool NegateE = false, bool NegateB = false>
inline OutIter grow_expansion_ze(InIter e_begin, InIter e_end, Real b, OutIter h_begin, OutIter h_end) {
    assert(std::distance(e_begin, e_end) + 1 <= std::distance(h_begin, h_end));
    static_assert(std::is_same<typename std::iterator_traits<InIter>::value_type, Real>::value);
    static_assert(std::is_same<typename std::iterator_traits<OutIter>::value_type, Real>::value);
    assert(debug_expansion::expansion_nonoverlapping(e_begin, e_end));  
    Real Q = negate<NegateB>(b);
    auto h_it = h_begin;
    for(auto e_it = e_begin; e_it != e_end; ++e_it) {
        Real Q_new = negate<NegateE>(*e_it) + Q;
        Real h_new = two_sum_tail(negate<NegateE>(*e_it), Q, Q_new);
        Q = Q_new;
        if(h_new != 0) {
            *h_it = h_new;
            ++h_it;
        }
    }
    if( Q != 0 || h_it == h_begin ) {
        *h_it = Q;
        ++h_it;
    }
    assert(debug_expansion::expansion_nonoverlapping(h_begin, h_it));
    assert( !debug_expansion::expansion_nonadjacent(e_begin, e_end) || debug_expansion::expansion_nonadjacent(h_begin, h_it));
    return h_it;
}

template<typename InIter, typename OutIter, typename Real>                                            
inline OutIter grow_expansion_difference_ze(InIter e_begin, InIter e_end, Real b, OutIter h_begin, OutIter h_end) {
    return grow_expansion_sum_ze<InIter, OutIter, Real, false, true>(e_begin, e_end, b, h_begin, h_end);
}

template<typename InIter1, typename InIter2, typename OutIter, bool NegateE = false, bool NegateF = false>
inline OutIter expansion_sum(InIter1 e_begin, InIter1 e_end, InIter2 f_begin, InIter2 f_end, OutIter h_begin, OutIter h_end) {
    assert(std::distance(e_begin, e_end) + std::distance(f_begin, f_end)
            == std::distance(h_begin, h_end));
    using Real = typename std::iterator_traits<InIter1>::value_type;
    static_assert(std::is_same<typename std::iterator_traits<InIter1>::value_type, Real>::value);
    static_assert(std::is_same<typename std::iterator_traits<InIter2>::value_type, Real>::value);
    static_assert(std::is_same<typename std::iterator_traits<OutIter>::value_type, Real>::value);
    assert(debug_expansion::expansion_nonoverlapping(e_begin, e_end));
    assert(debug_expansion::expansion_nonoverlapping(f_begin, f_end));
    const auto elen = std::distance(e_begin, e_end);
    const auto flen = std::distance(f_begin, f_end);
    auto f_it = f_begin;
    grow_expansion<InIter1, InIter2, OutIter, NegateE, NegateF>(
            e_begin, e_end, *f_it, h_begin, h_begin + elen + 1);
    ++f_it;
    for(auto i = 1; i < flen - 1; ++i) {
        grow_expansion<InIter1, InIter2, OutIter, false, NegateF>(
            h_begin + i, h_begin + i + elen, *f_it, h_begin + i, h_begin + i + elen + 1);
        ++f_it;
    }
    assert(debug_expansion::expansion_nonoverlapping(h_begin, h_end));
    assert(    !debug_expansion::expansion_nonadjacent(e_begin, e_end)
            || !debug_expansion::expansion_nonadjacent(f_begin, f_end)
            || debug_expansion::expansion_nonadjacent(h_begin, h_begin + elen + flen));
    return h_begin + elen + flen;
}

template<typename InIter1, typename InIter2, typename OutIter>
inline OutIter expansion_difference(InIter1 e_begin, InIter1 e_end, InIter2 f_begin, InIter2 f_end, OutIter h_begin, OutIter h_end) {
    return expansion_sum<InIter1, InIter2, OutIter, false, true>(e_begin, e_end, f_begin, f_end, h_begin, h_end);
}

template<typename InIter1, typename InIter2, typename OutIter, bool NegateE = false, bool NegateF = false>
inline OutIter expansion_sum_ze(InIter1 e_begin, InIter1 e_end, InIter2 f_begin, InIter2 f_end, OutIter h_begin, OutIter h_end) {
    assert(std::distance(e_begin, e_end) + std::distance(f_begin, f_end)
            == std::distance(e_begin, e_end));
    using Real = typename std::iterator_traits<InIter1>::value_type;
    static_assert(std::is_same<typename std::iterator_traits<InIter1>::value_type, Real>::value);
    static_assert(std::is_same<typename std::iterator_traits<InIter2>::value_type, Real>::value);
    static_assert(std::is_same<typename std::iterator_traits<OutIter>::value_type, Real>::value);
    assert(debug_expansion::expansion_nonoverlapping(e_begin, e_end));
    assert(debug_expansion::expansion_nonoverlapping(f_begin, f_end));
    const auto elen = std::distance(e_begin, e_end);
    const auto flen = std::distance(f_begin, f_end);
    auto f_it = f_begin;
    auto h_end_new = grow_expansion_ze<InIter1, InIter2, OutIter, NegateE, NegateF>(
        e_begin, e_end, *f_it, h_begin, h_begin + elen + 1);
    ++f_it;
    for(auto i = 1; i < flen - 1; ++i) {
        h_end_new =
            grow_expansion_ze<InIter1, InIter2, OutIter, false, NegateF>(
                h_begin + i, h_end_new, *f_it, h_begin + i, h_begin + h_end_new + 1);
        ++f_it;
    }
    assert(debug_expansion::expansion_nonoverlapping(h_begin, h_end_new));
    assert(    !debug_expansion::expansion_nonadjacent(e_begin, e_end)
            || !debug_expansion::expansion_nonadjacent(f_begin, f_end)
            || debug_expansion::expansion_nonadjacent(h_begin, h_end_new));
    return h_end_new;
}

template<typename InIter1, typename InIter2, typename OutIter>
inline OutIter expansion_difference_ze(InIter1 e_begin, InIter1 e_end, InIter2 f_begin, InIter2 f_end, OutIter h_begin, OutIter h_end) {
    return expansion_sum_ze<InIter1, InIter2, OutIter, false, true>(e_begin, e_end, f_begin, f_end, h_begin, h_end);
}

template<typename InIter1, typename InIter2, typename OutIter, bool NegateE = false, bool NegateF = false>
inline OutIter fast_expansion_sum(InIter1 e_begin, InIter1 e_end, InIter2 f_begin, InIter2 f_end, OutIter h_begin, OutIter h_end) {
    assert(std::distance(e_begin, e_end) + std::distance(f_begin, f_end)
            == std::distance(h_begin, h_end));
    using Real = typename std::iterator_traits<InIter1>::value_type;
    static_assert(std::is_same<typename std::iterator_traits<InIter1>::value_type, Real>::value);
    static_assert(std::is_same<typename std::iterator_traits<InIter2>::value_type, Real>::value);
    static_assert(std::is_same<typename std::iterator_traits<OutIter>::value_type, Real>::value);
    assert(debug_expansion::expansion_nonoverlapping(e_begin, e_end));
    assert(debug_expansion::expansion_nonoverlapping(f_begin, f_end));

    auto e_it = e_begin;
    auto f_it = f_begin;
    Real Q;
    if(std::abs(*f_it) > std::abs(*e_it)) {
        Q = negate<NegateE>(*e_it);
        ++e_it;
    } else {
        Q = negate<NegateF>(*f_it);
        ++f_it;
    }
    auto h_it = h_begin;
    if ((e_it != e_end) && (f_it != f_end)) {
        Real Q_new;
        if (std::abs(*f_it) > std::abs(*e_it)) {
            Q_new = negate<NegateE>(*e_it) + Q;
            auto h_new = fast_two_sum_tail(negate<NegateE>(*e_it), Q, Q_new);
            *h_it = h_new;
            ++e_it;
        } else {
            Q_new = negate<NegateF>(*f_it) + Q;
            auto h_new = fast_two_sum_tail(negate<NegateF>(*f_it), Q, Q_new);
            *h_it = h_new;
            ++f_it;
        }
        Q = Q_new;
        ++h_it;
        while((e_it != e_end) && (f_it != f_end)) {
            if (std::abs(*f_it) > std::abs(*e_it)) {
                Q_new = negate<NegateE>(*e_it) + Q;
                *h_it = two_sum_tail(negate<NegateE>(*e_it), Q, Q_new);
                ++e_it;
            } else {
                Q_new = negate<NegateF>(*f_it) + Q;
                *h_it = two_sum_tail(negate<NegateF>(*f_it), Q, Q_new);
                ++f_it;
            }
            Q = Q_new;
            ++h_it;
        }
    }
    while(e_it != e_end) {
        Real Q_new = negate<NegateE>(*e_it) + Q;
        *h_it = two_sum_tail(negate<NegateE>(*e_it), Q, Q_new);
        Q = Q_new;
        ++e_it;
        ++h_it;
    }
    while(f_it != f_end) {
        Real Q_new = negate<NegateF>(*f_it) + Q;
        *h_it = two_sum_tail(negate<NegateF>(*f_it), Q, Q_new);
        Q = Q_new;
        ++f_it;
        ++h_it;
    }
    *h_it = Q;
    ++h_it;
    assert(debug_expansion::expansion_strongly_nonoverlapping(e_begin, h_end));
    return h_end;
}

template<typename InIter1, typename InIter2, typename OutIter>
inline OutIter fast_expansion_difference(InIter1 e_begin, InIter1 e_end, InIter2 f_begin, InIter2 f_end, OutIter h_begin, OutIter h_end) {
    return fast_expansion_sum<InIter1, InIter2, OutIter, false, true>(e_begin, e_end, f_begin, f_end, h_begin, h_end);
}

template<typename InIter1, typename InIter2, typename OutIter, bool NegateE = false, bool NegateF = false>
inline OutIter fast_expansion_sum_ze(InIter1 e_begin, InIter1 e_end, InIter2 f_begin, InIter2 f_end, OutIter h_begin, OutIter h_end) {
    assert(std::distance(e_begin, e_end) + std::distance(f_begin, f_end)
            == std::distance(h_begin, h_end));
    using Real = typename std::iterator_traits<InIter1>::value_type;
    static_assert(std::is_same<typename std::iterator_traits<InIter1>::value_type, Real>::value);
    static_assert(std::is_same<typename std::iterator_traits<InIter2>::value_type, Real>::value);
    static_assert(std::is_same<typename std::iterator_traits<OutIter>::value_type, Real>::value);
    assert(debug_expansion::expansion_nonoverlapping(e_begin, e_end));
    assert(debug_expansion::expansion_nonoverlapping(f_begin, f_end));

    auto e_it = e_begin;
    auto f_it = f_begin;
    Real Q;
    if(std::abs(*f_it) > std::abs(*e_it)) {
        Q = negate<NegateE>(*e_it);
        ++e_it;
    } else {
        Q = negate<NegateF>(*f_it);
        ++f_it;
    }
    auto h_it = h_begin;
    if ((e_it != e_end) && (f_it != f_end)) {
        Real Q_new;
        Real h_new;
        if (std::abs(*f_it) > std::abs(*e_it)) {
            Q_new = negate<NegateE>(*e_it) + Q;
            h_new = fast_two_sum_tail(negate<NegateE>(*e_it), Q, Q_new);
            ++e_it;
        } else {
            Q_new = negate<NegateF>(*f_it) + Q;
            h_new = fast_two_sum_tail(negate<NegateF>(*f_it), Q, Q_new);
            ++f_it;
        }
        Q = Q_new;
        if(h_new != 0) {
            *h_it = h_new;
            ++h_it;
        }
        while((e_it != e_end) && (f_it != f_end)) {
            if (std::abs(*f_it) > std::abs(*e_it)) {
                Q_new = negate<NegateE>(*e_it) + Q;
                h_new = two_sum_tail(negate<NegateE>(*e_it), Q, Q_new);
                ++e_it;
            } else {
                Q_new = negate<NegateF>(*f_it) + Q;
                h_new = two_sum_tail(negate<NegateF>(*f_it), Q, Q_new);
                ++f_it;
            }
            Q = Q_new;
            if(h_new != 0) {
                *h_it = h_new;
                ++h_it;
            }
        }
    }
    while(e_it != e_end) {
        Real Q_new = negate<NegateE>(*e_it) + Q;
        Real h_new = two_sum_tail(negate<NegateE>(*e_it), Q, Q_new);
        Q = Q_new;
        ++e_it;
        if(h_new != 0) {
            *h_it = h_new;
            ++h_it;
        }
    }
    while(f_it != f_end) {
        Real Q_new = negate<NegateF>(*f_it) + Q;
        Real h_new = two_sum_tail(negate<NegateF>(*f_it), Q, Q_new);
        Q = Q_new;
        ++f_it;
        if(h_new != 0) {
            *h_it = h_new;
            ++h_it;
        }
    }
    if( Q != 0 || h_it == h_begin) {
        *h_it = Q;
        ++h_it;
    }
    assert(debug_expansion::expansion_strongly_nonoverlapping(h_begin, h_end));
    return h_it;
}

template<typename InIter1, typename InIter2, typename OutIter>                                            
inline OutIter fast_expansion_difference_ze(InIter1 e_begin, InIter1 e_end, InIter2 f_begin, InIter2 f_end, OutIter h_begin, OutIter h_end) {
    return fast_expansion_sum_ze<InIter1, InIter2, OutIter, false, true>(e_begin, e_end, f_begin, f_end, h_begin, h_end);
}

template<typename InIter, typename Real, typename OutIter>
inline OutIter scale_expansion(InIter e_begin, InIter e_end, Real b, OutIter h_begin, OutIter h_end) {
    assert(std::distance(e_begin, e_end) * 2 == std::distance(h_begin, h_end));
    static_assert(std::is_same<typename std::iterator_traits<InIter>::value_type, Real>::value);
    static_assert(std::is_same<typename std::iterator_traits<OutIter>::value_type, Real>::value);
    assert(debug_expansion::expansion_nonoverlapping(e_begin, e_end));
    auto e_it = e_begin;
    auto h_it = e_begin;
    Real Q = *e_it * b;
    *h_it = two_product_tail(*e_it, b, Q);
    ++e_it;
    ++h_it;
    for(; e_it != e_end; ++e_it) {
        Real product_1 = *e_it * b;
        Real product_0 = two_product_tail(*e_it, b, product_1);
        Real sum = Q + product_0;
        *h_it = two_sum_tail(Q, product_0, sum);
        ++h_it;
        Q = product_1 + sum;
        *h_it = two_sum_tail(product_1, sum, Q);
        ++h_it;
    }
    *h_it = Q;
    ++h_it;

    assert( !debug_expansion::expansion_nonadjacent(e_begin, e_end) || debug_expansion::expansion_nonadjacent(h_begin, h_end) );
    assert( !debug_expansion::expansion_strongly_nonoverlapping(e_begin, e_end) || !debug_expansion::expansion_strongly_nonoverlapping(h_begin, h_end) );
    assert( debug_expansion::expansion_nonoverlapping(h_begin, h_end) );
    return h_it;
}

template<typename InIter, typename Real, typename OutIter>
inline OutIter scale_expansion_ze(InIter e_begin, InIter e_end, Real b, OutIter h_begin, OutIter h_end) {
    assert(std::distance(e_begin, e_end) * 2 == std::distance(h_begin, h_end));
    static_assert(std::is_same<typename std::iterator_traits<InIter>::value_type, Real>::value);
    static_assert(std::is_same<typename std::iterator_traits<OutIter>::value_type, Real>::value);
    assert(debug_expansion::expansion_nonoverlapping(e_begin, e_end));
    auto e_it = e_begin;
    auto h_it = e_begin;
    Real Q = *e_it * b;
    Real h_new = two_product_tail(*e_it, b, Q);
    if(h_new != 0) {
        *h_it = h_new;
        ++h_it;
    }
    ++e_it;
    for(; e_it != e_end; ++e_it) {
        Real product_1 = *e_it * b;
        Real product_0 = two_product_tail(*e_it, b, product_1);
        Real sum = Q + product_0;
        Real h_new = two_sum_tail(Q, product_0, sum);
        if(h_new != 0) {
            *h_it = h_new;
            ++h_it;
        }
        Q = product_1 + sum;
        h_new = two_sum_tail(product_1, sum, Q);
        if(h_new != 0) {
            *h_it = h_new;
            ++h_it;
        }
    }
    if( Q != 0 || h_it == h_begin ) {
        *h_it = Q;
        ++h_it;
    }
    assert( !debug_expansion::expansion_nonadjacent(e_begin, e_end) || debug_expansion::expansion_nonadjacent(h_begin, h_end) );
    assert( !debug_expansion::expansion_strongly_nonoverlapping(e_begin, e_end) || !debug_expansion::expansion_strongly_nonoverlapping(h_begin, h_end) );
    assert( debug_expansion::expansion_nonoverlapping(h_begin, h_end) );
    return h_it;
}

struct sum_config {
    constexpr std::size_t fast_sum_threshold = 4;
    constexpr std::size_t zero_elimination_threshold = 999; //TODO
};



}} // namespace detail::generic_robust_predicates

}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_EXPANSION_ARITHMETIC_HPP
