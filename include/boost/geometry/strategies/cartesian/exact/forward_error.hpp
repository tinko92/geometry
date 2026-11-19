// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Licensed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_GEOMETRY_STRATEGIES_CARTESIAN_EXACT_FORWARD_ERROR_HPP
#define BOOST_GEOMETRY_STRATEGIES_CARTESIAN_EXACT_FORWARD_ERROR_HPP

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <type_traits>

#include <boost/geometry/strategies/cartesian/exact/arithmetic.hpp>

namespace boost { namespace geometry { namespace strategies { namespace detail
{
namespace exact_cartesian { namespace forward_error
{

// Absolute forward errors: |exact_coordinate - value[d]| <= error[d].
// Errors include any preceding native-integer conversion. No interval
// arithmetic, expression ownership, allocation, or exact fallback is hidden here.
struct point_approximation
{
    std::array<double, 2> value{};
    std::array<double, 2> error{};

    bool valid() const
    {
        return std::isfinite(value[0]) && std::isfinite(value[1])
            && std::isfinite(error[0]) && std::isfinite(error[1])
            && error[0] >= 0 && error[1] >= 0;
    }
};

struct scalar_approximation
{
    double value = 0;
    double error = 0;

    bool valid() const
    { return std::isfinite(value) && std::isfinite(error) && error >= 0; }
};

namespace detail
{
namespace arithmetic = arithmetic_detail;
constexpr double u = std::numeric_limits<double>::epsilon() / 2;
constexpr double rho = 2 * u;

inline bool nearest()
{
    return std::numeric_limits<double>::is_iec559
        && std::numeric_limits<double>::radix == 2
        && std::numeric_limits<double>::digits == 53;
}

// FP overflow/underflow (including in residual/error evaluation) is outside
// this first-stage contract. All formulas use fl(z)=z(1+delta), |delta|<=u.
// Posterior single-operation error is <=u/(1-u)*|fl(z)| < rho*|fl(z)|.
// Each fixed nonnegative error polynomial below has a path of <=16 rounded
// operations. Its computed value is >=(1-u)^16 times the real polynomial.
// The final multiplication is also rounded: (1+64u)(1-u)^17 > 1 for u=2^-53.
// Hence this STATIC inflation covers positive-bound evaluation roundoff too.
inline double inflate(double value)
{ return arithmetic::rounded_mul(1 + 64 * u, value); }

inline bool sign(scalar_approximation const& a, int& result)
{
    if (!a.valid()) return false;
    if (a.value > a.error) { result = 1; return true; }
    if (a.value < -a.error) { result = -1; return true; }
    if (a.value == 0 && a.error == 0) { result = 0; return true; }
    return false;
}

template <typename T>
struct native_coordinate : std::integral_constant<bool,
    std::is_same<T, double>::value || std::is_same<T, float>::value
    || (std::is_integral<T>::value && std::numeric_limits<T>::digits <= 64)> {};

template <typename T>
inline bool native(T value, double& approximate, double& error, std::true_type)
{
    double const converted = static_cast<double>(value);
    if (!std::isfinite(converted)) return false;
    bool const negative = std::is_signed<T>::value && value < 0;
    std::uint64_t const magnitude = negative ? std::uint64_t(0) - std::uint64_t(value)
                                             : std::uint64_t(value);
    double const absolute = std::abs(converted);
    // Guard the backconversion, including uint64_max -> double(2^64).
    bool const exact = absolute < std::ldexp(1.0, 64)
        && static_cast<std::uint64_t>(absolute) == magnitude;
    approximate = converted;
    // An integer-to-double conversion selects an adjacent binary64 value.
    // Its full adjacent-ULP error is <=2u*|converted|. This power-of-two
    // scaling is exact under the no-under/overflow contract, even if the
    // implementation selects either neighbour instead of nearest-even.
    error = exact ? 0 : arithmetic::rounded_mul(rho, absolute);
    return std::isfinite(error);
}

template <typename T>
inline bool native(T value, double& approximate, double& error, std::false_type)
{
    if (std::numeric_limits<T>::radix != 2 || std::numeric_limits<T>::digits > 53) return false;
    double const converted = static_cast<double>(value);
    if (!std::isfinite(converted)) return false;
    approximate = converted;
    error = 0;
    return true;
}

template <typename T>
inline typename std::enable_if<native_coordinate<T>::value, bool>::type
extract(T value, double& approximate, double& error)
{ return native(value, approximate, error, std::is_integral<T>()); }

template <typename T>
inline typename std::enable_if<!native_coordinate<T>::value, bool>::type
extract(T const&, double&, double&)
{ return false; }

inline bool exact_coordinates(point_approximation const& point)
{ return point.error[0] == 0 && point.error[1] == 0; }

inline std::array<double, 2> difference_values(point_approximation const& a, point_approximation const& b)
{
    return {{arithmetic::rounded_sub(a.value[0], b.value[0]),
             arithmetic::rounded_sub(a.value[1], b.value[1])}};
}

// Exact-converted endpoints give each product at most three roundings (two
// differences and one product), including repeated differences in a square.
// Its reverse relative error is <=(1-u)^-3-1 <4u. The final sum/difference
// adds <=u*M, M=|p|+|q|. Computed M >=(1-u)*M, so 8*(1-u)>5 proves that
// ONE final 8u coefficient covers the error and bound evaluation. Scaling by
// this power of two is exact under the first-stage no-under/overflow contract.
inline scalar_approximation native_bilinear(std::array<double, 2> const& a,
                                            std::array<double, 2> const& b, bool determinant)
{
    double const p = arithmetic::rounded_mul(a[0], b[determinant ? 1 : 0]);
    double const q = arithmetic::rounded_mul(a[1], b[determinant ? 0 : 1]);
    return {determinant ? arithmetic::rounded_sub(p, q) : arithmetic::rounded_add(p, q),
        arithmetic::rounded_mul(8 * u, arithmetic::rounded_add(std::abs(p), std::abs(q)))};
}

inline point_approximation difference(point_approximation const& a, point_approximation const& b)
{
    point_approximation result;
    for (std::size_t d = 0; d < 2; ++d)
    {
        result.value[d] = arithmetic::rounded_sub(a.value[d], b.value[d]);
        // |(A-B)-fl(a-b)| <= Ea+Eb+rho*|fl(a-b)|.
        result.error[d] = inflate(arithmetic::rounded_add(
            arithmetic::rounded_add(a.error[d], b.error[d]),
            arithmetic::rounded_mul(rho, std::abs(result.value[d]))));
    }
    return result;
}

// Determinant and dot product share this nonnegative polynomial, with either
// (ax,by),(ay,bx) or (ax,bx),(ay,by) as the two products. Input perturbation of
// x*y is bounded by |x|Ey+|y|Ex+Ex*Ey (the second-order term is not dropped).
// For computed products p,q and their computed sum/difference s, arithmetic
// error is <=rho*(|p|+|q|+|s|), by the posterior single-operation bound above.
inline double bilinear_error(double ax, double ay, double bx, double by,
                              double eax, double eay, double ebx, double eby,
                              double p, double q, double s)
{
    double error = arithmetic::rounded_mul(rho, arithmetic::rounded_add(
        arithmetic::rounded_add(std::abs(p), std::abs(q)), std::abs(s)));
    error = arithmetic::rounded_add(error, arithmetic::rounded_mul(std::abs(ax), ebx));
    error = arithmetic::rounded_add(error, arithmetic::rounded_mul(std::abs(bx), eax));
    error = arithmetic::rounded_add(error, arithmetic::rounded_mul(eax, ebx));
    error = arithmetic::rounded_add(error, arithmetic::rounded_mul(std::abs(ay), eby));
    error = arithmetic::rounded_add(error, arithmetic::rounded_mul(std::abs(by), eay));
    error = arithmetic::rounded_add(error, arithmetic::rounded_mul(eay, eby));
    return inflate(error);
}

inline scalar_approximation determinant(point_approximation const& a, point_approximation const& b)
{
    if (exact_coordinates(a) && exact_coordinates(b)) { return native_bilinear(a.value, b.value, true); }
    double const p = arithmetic::rounded_mul(a.value[0], b.value[1]);
    double const q = arithmetic::rounded_mul(a.value[1], b.value[0]);
    double const value = arithmetic::rounded_sub(p, q);
    return {value, bilinear_error(a.value[0], a.value[1], b.value[1], b.value[0],
        a.error[0], a.error[1], b.error[1], b.error[0], p, q, value)};
}

inline scalar_approximation dot(point_approximation const& a, point_approximation const& b)
{
    if (exact_coordinates(a) && exact_coordinates(b)) { return native_bilinear(a.value, b.value, false); }
    double const p = arithmetic::rounded_mul(a.value[0], b.value[0]);
    double const q = arithmetic::rounded_mul(a.value[1], b.value[1]);
    double const value = arithmetic::rounded_add(p, q);
    return {value, bilinear_error(a.value[0], a.value[1], b.value[0], b.value[1],
        a.error[0], a.error[1], b.error[0], b.error[1], p, q, value)};
}

inline bool native_intersection(std::array<double, 8> const& p, bool exact, point_approximation* output,
                                scalar_approximation* parameter = nullptr, bool second = false)
{
    double const ux = arithmetic::rounded_sub(p[2], p[0]);
    double const uy = arithmetic::rounded_sub(p[3], p[1]);
    double const vx = arithmetic::rounded_sub(p[6], p[4]);
    double const vy = arithmetic::rounded_sub(p[7], p[5]);
    double const wx = arithmetic::rounded_sub(p[4], p[0]);
    double const wy = arithmetic::rounded_sub(p[5], p[1]);
    double const dp = arithmetic::rounded_mul(ux, vy), dq = arithmetic::rounded_mul(uy, vx);
    double const tp = arithmetic::rounded_mul(wx, second ? uy : vy);
    double const tq = arithmetic::rounded_mul(wy, second ? ux : vx);
    double const denominator = arithmetic::rounded_sub(dp, dq);
    double const numerator = arithmetic::rounded_sub(tp, tq);
    double sx = std::abs(ux), sy = std::abs(uy);
    double md = arithmetic::rounded_add(std::abs(dp), std::abs(dq));
    double mt = arithmetic::rounded_add(std::abs(tp), std::abs(tq));
    if (!exact)
    {
        // Rounded-integer endpoint errors are <=2u*|p[i]|. A computed
        // difference then has error <=3u*(|endpoint1|+|endpoint2|).
        sx = arithmetic::rounded_add(std::abs(p[2]), std::abs(p[0]));
        sy = arithmetic::rounded_add(std::abs(p[3]), std::abs(p[1]));
        double const svx = arithmetic::rounded_add(std::abs(p[6]), std::abs(p[4]));
        double const svy = arithmetic::rounded_add(std::abs(p[7]), std::abs(p[5]));
        double const swx = arithmetic::rounded_add(std::abs(p[4]), std::abs(p[0]));
        double const swy = arithmetic::rounded_add(std::abs(p[5]), std::abs(p[1]));
        md = arithmetic::rounded_add(arithmetic::rounded_mul(sx, svy), arithmetic::rounded_mul(sy, svx));
        mt = arithmetic::rounded_add(arithmetic::rounded_mul(swx, second ? sy : svy),
            arithmetic::rounded_mul(swy, second ? sx : svx));
    }
    // Exact-converted native data: every product has three rounded operations
    // (two differences, one product). For M=|product1|+|product2|, determinant
    // error <=[gamma3+u*(1+u)/(1-u)]*M <5u*M. The computed M may understate by
    // (1-u), so K=8 bounds both D and T, with power-of-two error scaling exact.
    // Rounded integers: writing S=|endpoint1|+|endpoint2|, a product's error is
    // <=(7u+17u^2+u^3)*S1*S2. The final subtraction makes the determinant error
    // <9u*(S1*S2+S3*S4). Its computed magnitude has paths of <=4 operations;
    // K=16 is safe since 16*(1-u)^4>9. No generic per-operation error objects.
    double const k = exact ? 8 : 16;
    double const denominator_error = arithmetic::rounded_mul(k * u, md);
    if (!std::isfinite(denominator) || !std::isfinite(numerator)
        || !std::isfinite(md) || !std::isfinite(mt)
        || !(std::abs(denominator) > denominator_error)) return false;
    double const lower = arithmetic::rounded_mul(1 - 2 * u,
        arithmetic::rounded_sub(std::abs(denominator), denominator_error));
    if (!(lower > 0) || !std::isfinite(lower)) return false;
    double const t = arithmetic::rounded_div(numerator, denominator);
    if (!std::isfinite(t)) return false;
    double const absolute_t = std::abs(t);
    double const quotient_budget = arithmetic::rounded_add(arithmetic::rounded_add(mt,
        arithmetic::rounded_mul(absolute_t, md)), std::abs(numerator));
    if (parameter)
    {
        // |exact_t-t| <=u*(K*MT+K*|t|MD+|T|)/L. As K>=1, the
        // computed (MT+|t|MD+|T|)/L needs only the final factor 2K*u:
        // its paths have <=4 operations, and 2*(1-u)^4>1. Both segment
        // parameters use the same D; N=cross(c-a,d-c) or cross(c-a,b-a).
        scalar_approximation const result{t, arithmetic::rounded_mul(2 * k * u,
            arithmetic::rounded_div(quotient_budget, lower))};
        if (!result.valid()) return false;
        *parameter = result;
        return true;
    }
    std::array<double, 2> const scale{{sx, sy}}, direction{{ux, uy}};
    point_approximation result;
    for (std::size_t d = 0; d < 2; ++d)
    {
        double const coordinate = std::fma(t, direction[d], p[d]);
        result.value[d] = coordinate;
        double base = arithmetic::rounded_add(std::abs(coordinate), arithmetic::rounded_mul(absolute_t, scale[d]));
        if (!exact) base = arithmetic::rounded_add(base, std::abs(p[d]));
        double const budget = arithmetic::rounded_add(base,
            arithmetic::rounded_mul(arithmetic::rounded_div(scale[d], lower), quotient_budget));
        // With L<=|exact D|, |exact t-t| <=u*(K*MT+K*|t|MD+|T|)/L.
        // Endpoint error is <=0 (exact) or 2u*|a| (rounded). Difference error
        // is <=2u*S or 4u*S, respectively, and |exact U|<=(1+6u)*S.
        // One FMA contributes <=2u*|coordinate|. Therefore the real coordinate
        // error is bounded by K*u*(1+6u)*B, with
        // B=|coordinate|+|t|S+[rounded:|a|]+(S/L)*(MT+|t|MD+|T|).
        // B's evaluation paths have <=8 operations. ONE final coefficient
        // 2K*u suffices: 2*(1-u)^8 >1+6u at u=2^-53. Thus C=16 (exact-
        // converted) or C=32 (rounded integer), including bound-evaluation
        // rounding. This is not an empirically selected epsilon.
        result.error[d] = arithmetic::rounded_mul(2 * k * u, budget);
    }
    if (!result.valid()) return false;
    *output = result;
    return true;
}

inline bool segment_parameter(std::array<point_approximation, 4> const& p, bool second,
                               scalar_approximation& output)
{
    for (auto const& point : p) if (!point.valid()) return false;
    auto const u = difference(p[1], p[0]);
    auto const v = difference(p[3], p[2]);
    auto const w = difference(p[2], p[0]);
    auto const numerator = determinant(w, second ? u : v);
    auto const denominator = determinant(u, v);
    if (!u.valid() || !v.valid() || !w.valid() || !numerator.valid() || !denominator.valid()
        || !(std::abs(denominator.value) > denominator.error)) return false;
    // Same detN/detD residual proof as the construction below. Lower is
    // certified despite its two rounded operations: (1-2u)*(1+u)^2<1.
    double const lower = arithmetic::rounded_mul(1 - 2 * detail::u,
        arithmetic::rounded_sub(std::abs(denominator.value), denominator.error));
    if (!(lower > 0) || !std::isfinite(lower)) return false;
    double const t = arithmetic::rounded_div(numerator.value, denominator.value);
    double const error = inflate(arithmetic::rounded_div(arithmetic::rounded_add(
        arithmetic::rounded_add(numerator.error,
            arithmetic::rounded_mul(std::abs(t), denominator.error)),
        arithmetic::rounded_mul(detail::u, std::abs(numerator.value))), lower));
    scalar_approximation const result{t, error};
    if (!result.valid()) return false;
    output = result;
    return true;
}

inline bool intersection(std::array<point_approximation, 4> const& p, point_approximation& output)
{
    for (auto const& point : p) if (!point.valid()) return false;
    auto const u = difference(p[1], p[0]);
    auto const v = difference(p[3], p[2]);
    auto const q = difference(p[2], p[0]);
    auto const numerator = determinant(q, v);
    auto const denominator = determinant(u, v);
    if (!u.valid() || !v.valid() || !q.valid() || !numerator.valid() || !denominator.valid()
        || !(std::abs(denominator.value) > denominator.error)) return false;

    // L <= |D|-ED <= |exact denominator|. The subtraction and multiplication
    // may each round up, but (1-2u)*(1+u)^2 < 1, so no nextafter is required.
    double const lower = arithmetic::rounded_mul(1 - 2 * detail::u,
        arithmetic::rounded_sub(std::abs(denominator.value), denominator.error));
    if (!(lower > 0) || !std::isfinite(lower)) return false;
    double const t = arithmetic::rounded_div(numerator.value, denominator.value);
    // If t=fl(N/D), |N-tD| <= u*|N|. Thus
    // |exact_t-t| <= (EN+|t|ED+u|N|)/L, including native conversion error.
    double const t_error = inflate(arithmetic::rounded_div(arithmetic::rounded_add(
        arithmetic::rounded_add(numerator.error,
            arithmetic::rounded_mul(std::abs(t), denominator.error)),
        arithmetic::rounded_mul(detail::u, std::abs(numerator.value))), lower));
    if (!std::isfinite(t) || !std::isfinite(t_error)) return false;

    point_approximation result;
    for (std::size_t d = 0; d < 2; ++d)
    {
        double const coordinate = std::fma(t, u.value[d], p[0].value[d]);
        result.value[d] = coordinate;
        // One fused rounding, not separately rounded t*u+a:
        // EP <= EA+|t|EU+(|u|+EU)Et+rho*|fl(fma(t,u,a))|.
        result.error[d] = inflate(arithmetic::rounded_add(arithmetic::rounded_add(
            arithmetic::rounded_add(p[0].error[d], arithmetic::rounded_mul(std::abs(t), u.error[d])),
            arithmetic::rounded_mul(arithmetic::rounded_add(std::abs(u.value[d]), u.error[d]), t_error)),
            arithmetic::rounded_mul(rho, std::abs(result.value[d]))));
    }
    if (!result.valid()) return false;
    output = result;
    return true;
}

struct approximation_reader
{
    bool operator()(point_approximation const& input, point_approximation& output) const
    { output = input; return input.valid(); }
};
} // namespace detail

// Every public API leaves output unchanged on false. A successful approximation
// is an enclosure, not a rounded exact answer; a sign API certifies only signs.
template <typename Coordinate>
inline bool try_native_coordinate(Coordinate value, double& approximate, double& error)
{
    if (!detail::nearest()) return false;
    double a, e;
    if (!detail::extract(value, a, e)) return false;
    approximate = a; error = e;
    return true;
}

template <typename X, typename Y>
inline bool try_native_point(X x, Y y, point_approximation& output)
{
    if (!detail::nearest()) return false;
    point_approximation result;
    if (!detail::extract(x, result.value[0], result.error[0])
        || !detail::extract(y, result.value[1], result.error[1])) return false;
    output = result;
    return true;
}

inline bool try_intersection(std::array<point_approximation, 4> const& endpoints,
                              point_approximation& output)
{ return detail::nearest() && detail::intersection(endpoints, output); }

template <typename Coordinate>
inline bool native_intersection(std::array<Coordinate, 8> const& endpoints, point_approximation& output)
{
    if (!detail::nearest()) return false;
    std::array<double, 8> points;
    bool exact = true;
    for (std::size_t i = 0; i < 8; ++i)
    {
        double error;
        if (!detail::extract(endpoints[i], points[i], error)) return false;
        exact = exact && error == 0;
    }
    return detail::native_intersection(points, exact, &output);
}

template <typename Coordinate>
inline bool try_intersection(std::array<Coordinate, 8> const& endpoints, point_approximation& output)
{ return native_intersection(endpoints, output); }

// Parameters satisfy a+t*(b-a) == c+s*(d-c). second selects s instead of t;
// endpoint order is {a,b,c,d}, or {ax,ay,bx,by,cx,cy,dx,dy} for native inputs.
inline bool try_segment_parameter(std::array<point_approximation, 4> const& endpoints, bool second,
                                   scalar_approximation& output)
{ return detail::nearest() && detail::segment_parameter(endpoints, second, output); }

template <typename Coordinate>
inline bool try_segment_parameter(std::array<Coordinate, 8> const& endpoints, bool second,
                                   scalar_approximation& output)
{
    if (!detail::nearest()) return false;
    std::array<double, 8> points;
    bool exact = true;
    for (std::size_t i = 0; i < 8; ++i)
    {
        double error;
        if (!detail::extract(endpoints[i], points[i], error)) return false;
        exact = exact && error == 0;
    }
    return detail::native_intersection(points, exact, nullptr, &output, second);
}

// Already-rounded integers must supply conversion errors here; passing only a
// double array instead declares the binary64 endpoint values themselves exact.
inline bool try_intersection(std::array<double, 8> const& values,
                              std::array<double, 8> const& errors, point_approximation& output)
{
    std::array<point_approximation, 4> points;
    for (std::size_t i = 0; i < 4; ++i)
        for (std::size_t d = 0; d < 2; ++d)
        { points[i].value[d] = values[2 * i + d]; points[i].error[d] = errors[2 * i + d]; }
    return try_intersection(points, output);
}

inline bool try_determinant(point_approximation const& a, point_approximation const& b,
                             scalar_approximation& output)
{
    if (!detail::nearest() || !a.valid() || !b.valid()) return false;
    auto const result = detail::determinant(a, b);
    if (!result.valid()) return false;
    output = result;
    return true;
}

inline bool try_determinant_sign(point_approximation const& a, point_approximation const& b, int& output)
{
    if (!detail::nearest() || !a.valid() || !b.valid()) return false;
    return detail::sign(detail::determinant(a, b), output);
}

inline bool try_dot(point_approximation const& a, point_approximation const& b, scalar_approximation& output)
{
    if (!detail::nearest() || !a.valid() || !b.valid()) return false;
    auto const result = detail::dot(a, b);
    if (!result.valid()) return false;
    output = result;
    return true;
}

inline bool try_signed_determinant(point_approximation const& a, point_approximation const& b,
                                    point_approximation const& p, scalar_approximation& output)
{
    if (!detail::nearest() || !a.valid() || !b.valid() || !p.valid()) return false;
    auto const result = detail::exact_coordinates(a) && detail::exact_coordinates(b) && detail::exact_coordinates(p)
        ? detail::native_bilinear(detail::difference_values(b, a), detail::difference_values(p, a), true)
        : detail::determinant(detail::difference(b, a), detail::difference(p, a));
    if (!result.valid()) return false;
    output = result;
    return true;
}

inline bool try_orientation(point_approximation const& a, point_approximation const& b,
                             point_approximation const& c, int& output)
{
    scalar_approximation determinant;
    return try_signed_determinant(a, b, c, determinant) && detail::sign(determinant, output);
}

// Sign of (b-a).(p-b), as used by the Cartesian direction code.
inline bool try_direction(point_approximation const& a, point_approximation const& b,
                           point_approximation const& p, int& output)
{
    if (!detail::nearest() || !a.valid() || !b.valid() || !p.valid()) return false;
    if (detail::exact_coordinates(a) && detail::exact_coordinates(b) && detail::exact_coordinates(p))
    {
        return detail::sign(detail::native_bilinear(detail::difference_values(b, a),
            detail::difference_values(p, b), false), output);
    }
    return detail::sign(detail::dot(detail::difference(b, a), detail::difference(p, b)), output);
}

inline bool try_squared_distance(point_approximation const& a, point_approximation const& b,
                                  scalar_approximation& output)
{
    if (!detail::nearest() || !a.valid() || !b.valid()) return false;
    if (detail::exact_coordinates(a) && detail::exact_coordinates(b))
    {
        auto const d = detail::difference_values(b, a);
        auto const result = detail::native_bilinear(d, d, false);
        if (!result.valid()) return false;
        output = result;
        return true;
    }
    auto const d = detail::difference(b, a);
    auto const result = detail::dot(d, d);
    if (!result.valid()) return false;
    output = result;
    return true;
}

// Sign of |b-a|^2-|d-c|^2. The three-point overload has a common origin.
inline bool try_squared_distance_compare(point_approximation const& a, point_approximation const& b,
                                          point_approximation const& c, point_approximation const& d,
                                          int& output)
{
    if (!detail::nearest() || !a.valid() || !b.valid() || !c.valid() || !d.valid()) return false;
    if (detail::exact_coordinates(a) && detail::exact_coordinates(b)
        && detail::exact_coordinates(c) && detail::exact_coordinates(d))
    {
        auto const first = detail::difference_values(b, a), second = detail::difference_values(d, c);
        double const x = detail::arithmetic::rounded_add(detail::arithmetic::rounded_mul(first[0], first[0]),
                                                        detail::arithmetic::rounded_mul(first[1], first[1]));
        double const y = detail::arithmetic::rounded_add(detail::arithmetic::rounded_mul(second[0], second[0]),
                                                        detail::arithmetic::rounded_mul(second[1], second[1]));
        // Product errors <4u*M, distance-sum errors <=u*M and final
        // subtraction error <=u*(1+u)*M total <7u*M. The reused positive
        // sum x+y has paths of two adds: 16*(1-u)^2>7 covers its rounding.
        return detail::sign({detail::arithmetic::rounded_sub(x, y),
            detail::arithmetic::rounded_mul(16 * detail::u, detail::arithmetic::rounded_add(x, y))}, output);
    }
    auto const first = detail::difference(b, a), second = detail::difference(d, c);
    auto const x = detail::dot(first, first), y = detail::dot(second, second);
    double const value = detail::arithmetic::rounded_sub(x.value, y.value);
    double const error = detail::inflate(detail::arithmetic::rounded_add(
        detail::arithmetic::rounded_add(x.error, y.error),
        detail::arithmetic::rounded_mul(detail::rho, std::abs(value))));
    return detail::sign({value, error}, output);
}

inline bool try_squared_distance_compare(point_approximation const& origin, point_approximation const& a,
                                          point_approximation const& b, int& output)
{ return try_squared_distance_compare(origin, a, origin, b, output); }

inline bool try_compare_scalars(scalar_approximation const& a, scalar_approximation const& b, int& output)
{
    if (!detail::nearest() || !a.valid() || !b.valid()) return false;
    if (a.error == 0 && b.error == 0)
    {
        output = (a.value > b.value) - (a.value < b.value);
        return true;
    }
    double const value = detail::arithmetic::rounded_sub(a.value, b.value);
    // RN is monotone and symmetric. RN(|a-b|) > RN(ea+eb) implies
    // |a-b| > ea+eb, even when both rounded quantities round down.
    double const error = detail::arithmetic::rounded_add(a.error, b.error);
    if (!(std::abs(value) > error)) return false;
    output = value > 0 ? 1 : -1;
    return true;
}

inline bool try_compare_scalars(double a, double ea, double b, double eb, int& output)
{ return try_compare_scalars(scalar_approximation{a, ea}, scalar_approximation{b, eb}, output); }

inline bool try_compare_coordinate(point_approximation const& a, std::size_t dimension,
                                    scalar_approximation const& b, int& output)
{
    if (dimension > 1 || !a.valid()) return false;
    return try_compare_scalars({a.value[dimension], a.error[dimension]}, b, output);
}

inline bool try_compare_coordinates(point_approximation const& a, point_approximation const& b,
                                     std::size_t dimension, int& output)
{
    if (dimension > 1 || !a.valid() || !b.valid()) return false;
    return try_compare_scalars({a.value[dimension], a.error[dimension]},
        {b.value[dimension], b.error[dimension]}, output);
}

// Product of the two directed box extents; no implicit correction of an
// invalid box. Input errors and the product's second-order term are retained.
inline bool try_box_measure(point_approximation const& minimum, point_approximation const& maximum,
                             scalar_approximation& output)
{
    if (!detail::nearest() || !minimum.valid() || !maximum.valid()) return false;
    if (detail::exact_coordinates(minimum) && detail::exact_coordinates(maximum))
    {
        auto const extent = detail::difference_values(maximum, minimum);
        double const value = detail::arithmetic::rounded_mul(extent[0], extent[1]);
        // The same three-rounding product proof needs <4u*|value|; 8u
        // leaves spare margin without constructing intermediate error objects.
        scalar_approximation const result{value, detail::arithmetic::rounded_mul(8 * detail::u, std::abs(value))};
        if (!result.valid()) return false;
        output = result;
        return true;
    }
    auto const extent = detail::difference(maximum, minimum);
    double const value = detail::arithmetic::rounded_mul(extent.value[0], extent.value[1]);
    double const error = detail::inflate(detail::arithmetic::rounded_add(detail::arithmetic::rounded_add(
        detail::arithmetic::rounded_add(
            detail::arithmetic::rounded_mul(std::abs(extent.value[0]), extent.error[1]),
            detail::arithmetic::rounded_mul(std::abs(extent.value[1]), extent.error[0])),
        detail::arithmetic::rounded_mul(extent.error[0], extent.error[1])),
        detail::arithmetic::rounded_mul(detail::rho, std::abs(value))));
    scalar_approximation const result{value, error};
    if (!result.valid()) return false;
    output = result;
    return true;
}

// Signed measure is sum((ax+bx)*(ay-by))/2, positive for clockwise traversal.
// Reader(element, point_approximation&) permits a native ring to be streamed
// without allocating a point array. close=false uses only the supplied edges,
// preserving declared-closed/view semantics even for a missing closing vertex.
template <typename Iterator, typename Reader>
inline bool try_ring_measure(Iterator first, Iterator last, Reader const& read, scalar_approximation& output,
                              bool close = true)
{
    using namespace detail;
    if (!nearest()) return false;
    if (first == last) { output = {0, 0}; return true; }
    point_approximation start, previous;
    if (!read(*first, start) || !start.valid()) return false;
    previous = start;
    double sum = 0, magnitude = 0, input_error = 0;
    std::size_t edges = 0;
    auto edge = [&](point_approximation const& a, point_approximation const& b)
    {
        if (edges == (std::size_t(1) << 20)) return false;
        ++edges;
        double const x = arithmetic::rounded_add(a.value[0], b.value[0]);
        double const y = arithmetic::rounded_sub(a.value[1], b.value[1]);
        double const ex = inflate(arithmetic::rounded_add(arithmetic::rounded_add(a.error[0], b.error[0]),
            arithmetic::rounded_mul(rho, std::abs(x))));
        double const ey = inflate(arithmetic::rounded_add(arithmetic::rounded_add(a.error[1], b.error[1]),
            arithmetic::rounded_mul(rho, std::abs(y))));
        double const term = arithmetic::rounded_mul(x, y);
        double const error = inflate(arithmetic::rounded_add(arithmetic::rounded_add(
            arithmetic::rounded_add(arithmetic::rounded_mul(std::abs(x), ey), arithmetic::rounded_mul(std::abs(y), ex)),
            arithmetic::rounded_mul(ex, ey)), arithmetic::rounded_mul(rho, std::abs(term))));
        sum = arithmetic::rounded_add(sum, term);
        magnitude = arithmetic::rounded_add(magnitude, std::abs(term));
        input_error = arithmetic::rounded_add(input_error, error);
        return std::isfinite(sum) && std::isfinite(magnitude) && std::isfinite(input_error);
    };
    for (++first; first != last; ++first)
    {
        point_approximation next;
        if (!read(*first, next) || !next.valid() || !edge(previous, next)) return false;
        previous = next;
    }
    if (close && !edge(previous, start)) return false;
    // For n<=2^20, gamma_n/(1-gamma_n) <4nu and
    // 1/(1-gamma_n) <1+4nu. These bound BOTH positive accumulation errors
    // and the signed-sum error; 4nu and 1+4nu are exactly representable.
    double const factor = std::ldexp(static_cast<double>(edges), -51);
    double const error = inflate(arithmetic::rounded_add(
        arithmetic::rounded_mul(1 + factor, input_error), arithmetic::rounded_mul(factor, magnitude)));
    scalar_approximation const result{arithmetic::rounded_mul(0.5, sum), arithmetic::rounded_mul(0.5, error)};
    if (!result.valid()) return false;
    output = result;
    return true;
}

template <typename Iterator>
inline bool try_ring_measure(Iterator first, Iterator last, scalar_approximation& output, bool close = true)
{ return try_ring_measure(first, last, detail::approximation_reader(), output, close); }

} // namespace forward_error
} // namespace exact_cartesian
}}}} // namespace boost::geometry::strategies::detail

#endif
