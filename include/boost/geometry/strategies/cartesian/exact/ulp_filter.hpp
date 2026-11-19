// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.

#ifndef BOOST_GEOMETRY_STRATEGIES_CARTESIAN_EXACT_ULP_FILTER_HPP
#define BOOST_GEOMETRY_STRATEGIES_CARTESIAN_EXACT_ULP_FILTER_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

#include <boost/geometry/strategies/cartesian/exact/arithmetic.hpp>

namespace boost { namespace geometry { namespace strategies { namespace detail
{ namespace exact_cartesian { namespace ulp_filter
{

// First-stage binary64 predicates. N = exact native coordinate, I = certified
// implicit approximation a with prev(a) <= exact <= next(a), independently in
// each dimension. No point ownership, stored errors, interval evaluation,
// precision allocation, or fallback is hidden here. False leaves the result
// unchanged and requests the caller's exact/indirect predicate.
//
// Contract: finite inputs, IEEE binary64, round-to-nearest/ties-even, gradual
// underflow (no FTZ/DAZ), no reassociation/fast-math. Value arithmetic
// overflow/underflow is out of scope.
// No release-mode floating-environment checks. Rounded operations below also
// prevent accidental contraction. Magnitude arithmetic may become subnormal
// at zero coordinates; the floors below explicitly cover THAT case.
//
// PROOF (u=2^-53, eta=2^-1074, m=2^-1022, so eta=2u*m).
//
// 1. H(a)=max(|a|,m). Adjacent spacing is <=2u*|a| for a normal a, including
//    either side of a binade. At subnormal a and zero, spacing is eta=2u*m.
//    Thus every implicit coordinate has |delta a| <=2u*H(a). Native error is
//    zero. At +/-max-double one neighbour is infinite: decline certification.
//    This is a bound on both sides, not a claim that a=0 implies exact=0.
//
// 2. For d=fl(b-a), |(exact_b-exact_a)-d| <=u/(1-u)*|d| plus endpoint errors,
//    hence <=2u*S, S=|d|+[I(a)]H(a)+[I(b)]H(b). Cancellation is covered by the
//    ENDPOINT magnitudes, not by eps*|d| alone. S uses <=2 positive additions.
//
// 3. For coordinate comparison, D=exact_a-exact_b, d=fl(a-b), put
//    M=|d|+[I(a)]H(a)+[I(b)]H(b). Then |D-d| <=2u*M.
//    Computed Mhat >=(1-u)^2 M and Mhat>=m. B=fl(4u*Mhat) >=3u*Mhat: scaling
//    by the power of two is exact normally; if subnormal, the loss <=eta/2
//    is <=one quarter of 4u*Mhat. Since 3(1-u)^2>2, B bounds the total error.
//
// 4. Orientation: uvec=fl(b-a), vvec=fl(c-a), p=fl(ux*vy), q=fl(uy*vx),
//    d=fl(p-q). Define the four S's as in (2), and
//       M=Sux*Svy+Suy*Svx.
//    Perturbing vector coordinates costs <=(4u+4u^2)M, including the quadratic
//    error terms. Rounding the products/subtraction costs <=(2u+u^2)M.
//    Therefore |exact_orientation-d| <=(6u+5u^2)M.
//
//    Each computed S >=(1-u)^2 S. Each positive magnitude product may lose
//    eta/2 to gradual underflow. After the final sum,
//       Mhat0 >=(1-u)^6 M-(1-u)eta.
//    Set Mhat=max(Mhat0,m). Then
//       (1-u)^6 M <=Mhat+eta <=(1+2u)Mhat.
//    B=fl(8u*Mhat) >=7u*Mhat: its unrounded value is >=4eta, so rounding can
//    lose at most one eighth, even subnormally. The strict binary64 inequality
//       7(1-u)^6 > (6+5u)(1+2u)
//    proves that B bounds all errors, including rounded magnitude evaluation.
//    Artificial tiny magnitude products at zero cannot erase the certificate.
//
// 5. Only d>B or d<-B certifies a numerical sign. d=0 NEVER certifies zero for
//    uncertain inputs. A coincident native pair in NNI does prove exact zero.
//    Coordinate equality NI/II can certify inequality, not equality, without
//    additional exact information from the caller.
//
// NNI/NII/III specialize the magnitude expression at compile time; no runtime
// operand-kind dispatch. Normalize other permutations by cyclic permutation
// (same sign) and exchange (opposite sign). NNN is deliberately absent: use
// side_robust<double, fp_equals_policy, 3> for native floating orientation.

using coordinates = std::array<double, 2>;
static_assert(std::numeric_limits<double>::is_iec559
    && std::numeric_limits<double>::radix == 2
    && std::numeric_limits<double>::digits == 53
    && std::numeric_limits<double>::min_exponent == -1021
    && std::numeric_limits<double>::max_exponent == 1024,
    "ULP filters require IEEE binary64");

namespace detail
{
namespace arithmetic = arithmetic_detail;
constexpr double unit_roundoff = std::numeric_limits<double>::epsilon() / 2;
constexpr double minimum_normal = (std::numeric_limits<double>::min)();

inline double implicit_magnitude(double value)
{
    double const absolute = std::abs(value);
    return absolute == (std::numeric_limits<double>::max)()
        ? (std::numeric_limits<double>::infinity)()
        : (std::max)(absolute, minimum_normal);
}

inline bool sign(double value, double bound, int& result)
{
    if (!std::isfinite(value) || !std::isfinite(bound)) return false;
    if (value > bound) { result = 1; return true; }
    if (value < -bound) { result = -1; return true; }
    return false;
}

template <bool FirstImplicit, bool SecondImplicit>
inline double difference_magnitude(double difference, double first, double second)
{
    double result = std::abs(difference);
    if constexpr (FirstImplicit)
        result = arithmetic::rounded_add(result, implicit_magnitude(first));
    if constexpr (SecondImplicit)
        result = arithmetic::rounded_add(result, implicit_magnitude(second));
    return result;
}

template <bool FirstImplicit, bool SecondImplicit>
inline bool compare(double a, double b, int& result)
{
    double const difference = arithmetic::rounded_sub(a, b);
    double const magnitude = difference_magnitude<FirstImplicit, SecondImplicit>(difference, a, b);
    double const bound = arithmetic::rounded_mul(4 * unit_roundoff,
        (std::max)(magnitude, minimum_normal));
    return sign(difference, bound, result);
}

template <bool AImplicit, bool BImplicit>
inline bool orientation(coordinates const& a, coordinates const& b, coordinates const& c, int& result)
{
    // C is implicit in every canonical combination provided here.
    double const ux = arithmetic::rounded_sub(b[0], a[0]);
    double const uy = arithmetic::rounded_sub(b[1], a[1]);
    double const vx = arithmetic::rounded_sub(c[0], a[0]);
    double const vy = arithmetic::rounded_sub(c[1], a[1]);
    double const p = arithmetic::rounded_mul(ux, vy);
    double const q = arithmetic::rounded_mul(uy, vx);
    double const value = arithmetic::rounded_sub(p, q);
    double const sux = difference_magnitude<AImplicit, BImplicit>(ux, a[0], b[0]);
    double const suy = difference_magnitude<AImplicit, BImplicit>(uy, a[1], b[1]);
    double const svx = difference_magnitude<AImplicit, true>(vx, a[0], c[0]);
    double const svy = difference_magnitude<AImplicit, true>(vy, a[1], c[1]);
    double const magnitude = arithmetic::rounded_add(arithmetic::rounded_mul(sux, svy),
        arithmetic::rounded_mul(suy, svx));
    double const bound = arithmetic::rounded_mul(8 * unit_roundoff,
        (std::max)(magnitude, minimum_normal));
    return sign(value, bound, result);
}
} // namespace detail

inline bool try_compare_nn(double a, double b, int& result)
{
    result = int(a > b) - int(a < b);
    return true;
}

inline bool try_compare_ni(double native, double implicit, int& result)
{ return detail::compare<false, true>(native, implicit, result); }

inline bool try_compare_ii(double first, double second, int& result)
{ return detail::compare<true, true>(first, second, result); }

inline bool try_equal_nn(double a, double b, bool& result)
{
    result = a == b;
    return true;
}

inline bool try_equal_ni(double native, double implicit, bool& result)
{
    int comparison;
    if (!try_compare_ni(native, implicit, comparison)) return false;
    result = false;
    return true;
}

inline bool try_equal_ii(double first, double second, bool& result)
{
    int comparison;
    if (!try_compare_ii(first, second, comparison)) return false;
    result = false;
    return true;
}

inline bool try_orientation_nni(coordinates const& a, coordinates const& b,
                                coordinates const& implicit, int& result)
{
    if (a[0] == b[0] && a[1] == b[1]) { result = 0; return true; }
    return detail::orientation<false, false>(a, b, implicit, result);
}

inline bool try_orientation_nii(coordinates const& native, coordinates const& b,
                                coordinates const& c, int& result)
{ return detail::orientation<false, true>(native, b, c, result); }

inline bool try_orientation_iii(coordinates const& a, coordinates const& b,
                                coordinates const& c, int& result)
{ return detail::orientation<true, true>(a, b, c, result); }

}}}}}} // namespace boost::geometry::strategies::detail::exact_cartesian::ulp_filter

#endif
