// Boost.Geometry

// Copyright (c) 2018 Barend Gehrels, Amsterdam, the Netherlands.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_ARITHMETIC_GENERAL_FORM_HPP
#define BOOST_GEOMETRY_ARITHMETIC_GENERAL_FORM_HPP

#include <boost/geometry/util/math.hpp>
#include <boost/geometry/core/access.hpp>
#include <boost/geometry/util/select_most_precise.hpp>

namespace boost { namespace geometry
{

namespace arithmetic
{

// Structure containing thresholds, per type, for denominator.
// Determined with corresponding unit test.
// It should not be replaced by machine epsilon or math::equals
template <typename GeneralType>
struct general_threshold {};

template <>
struct general_threshold<long double>
{
   static long double get() { return 1.0e-10; }
};

template <>
struct general_threshold<double>
{
   static double get() { return 1.0e-7; }
};

template <>
struct general_threshold<float>
{
   static float get() { return 1.0e-2; }
};


template <typename ValueType, typename Policy>
static inline
bool is_zero(ValueType const& value, Policy const& policy)
{
    return math::detail::equals_by_policy(value, ValueType(), policy);
}

//--------------------------------------------------------------------------
// Structure containing the General Form of a line, a*x + b*y + c == 0
// Might be conceptized later. Therefore operations are implemented outside
// the structure itself.
template <typename GeneralType = double>
struct general_form
{

    general_form()
        : a(GeneralType())
        , b(GeneralType())
        , c(GeneralType())
        , normalized(false)
    {}

    // Horizontal: a == 0, for example y-3=0, y==3
    // Vertical: b == 0, for example x-2=0, x==2
    // Through origin: c == 0
    GeneralType a;
    GeneralType b;
    GeneralType c;
    bool normalized;

    GeneralType magnitude() const
    {
        // TODO: store this
        return (std::max)(std::abs(a), std::abs(b));
    }
};

template <typename GeneralType, typename Coordinate>
inline
general_form<GeneralType> construct_line(Coordinate const& x1,
    Coordinate const& y1, Coordinate const& x2, Coordinate const& y2)
{
    general_form<GeneralType> result;
    result.a = y1 - y2;
    result.b = x2 - x1;
    result.c = -result.a * x1 - result.b * y1;
    return result;
}

template <typename T, typename Point>
inline general_form<T> construct_line(Point const& a, Point const& b)
{
    return construct_line<T>(geometry::get<0>(a),
                             geometry::get<1>(a),
                             geometry::get<0>(b),
                             geometry::get<1>(b));
}

template <typename T, typename Segment>
inline general_form<T> construct_line(Segment const& segment)
{
    return construct_line<T>(geometry::get<0, 0>(segment),
                             geometry::get<0, 1>(segment),
                             geometry::get<1, 0>(segment),
                             geometry::get<1, 1>(segment));
}

//! Normalize the line. For robustness reasons it is often better to
//! NOT use normalization. It uses sqrt and therefore the intersection
//! point, if calculated with normalization, might go,
//! for example, from 7 to 6.99...997
template <typename FloatingPointType, typename InputType>
inline
general_form<FloatingPointType> normalize_line(general_form<InputType> const& p)
{
    general_form<FloatingPointType> result;
    result.a = p.a;
    result.b = p.b;
    result.c = p.c;

    FloatingPointType const norm
            = std::sqrt(result.a * result.a + result.b * result.b);

    // Compare with 0 (even very small values like 1.0e-12 are supported)
    if (norm != 0)
    {
        result.a /= norm;
        result.b /= norm;
        result.c /= norm;
        result.normalized = true;
    }
    return result;
}

template <typename GeneralType>
inline bool more_horizontal(general_form<GeneralType> const& p)
{
    // If a=0, then line is 'by+c=0' so y=-c/b  which is horizontal
    // If a=0.1 and b=0.9, then the line is quite horizontal
    return std::abs(p.a) < std::abs(p.b);
}

template <typename GeneralType>
inline bool has_horizontal_component(general_form<GeneralType> const& p)
{
    // TODO: thresholds should be determined
    return std::abs(p.a) < 1e-6 || std::abs(p.b / p.a) >= 0.1;
}

template <typename GeneralType>
inline bool has_vertical_component(general_form<GeneralType> const& p)
{
    // TODO: thresholds should be determined
    return std::abs(p.b) < 1e-6 || std::abs(p.a / p.b) >= 0.1;
}


// Returns a comparable and sortable distance measure:
// - if left then positive,
// - if collinear than zero - but NOT compared with any epsilon
// - if right then negative
// SQRT need not to be taken, therefore numerator is squared, but sign
// is preserved
template <typename GeneralType, typename CoordinateType>
inline
typename select_most_precise<GeneralType, CoordinateType>::type
signed_comparable_distance(general_form<GeneralType> const& p,
    CoordinateType const& x, CoordinateType const& y)
{
    // Distance from point to line in general form is given as:
    // (a * x + b * y + c) / sqrt(a * a + b * b);
    // In most use cases comparisons are enough, saving the sqrt
    // (performance plus making it a bit more precise)
    typedef typename select_most_precise<GeneralType, CoordinateType>::type
            ct;

    const ct num = p.a * x + p.b * y + p.c;
    if (num == 0)
    {
        return 0;
    }

    // TODO: should be precalculated
    const GeneralType denom = p.a * p.a + p.b * p.b;
    BOOST_ASSERT(denom != 0);
    return (num > 0 ? 1.0 : -1.0) * num * num / denom;
}

// Returns true for two lines which are supposed to be (close to) collinear
// (which is not checked) and have a similar direction
// (in practice up to 45 degrees, TO BE VERIFIED)
// true: -----------------> p -----------------> q
// false: -----------------> p <----------------- q
template <typename GeneralType>
inline
bool similar_direction(const general_form<GeneralType>& p,
                       const general_form<GeneralType>& q)
{
    return p.a * q.a >= 0 && p.b * q.b >= 0;
}

// Calculates intersection point of two infinite lines.
// Returns true if the lines intersect.
// Returns false if lines are considered as collinear.
// Set the doubt flag for nearly collinear lines.
template <typename Point, typename GeneralType>
inline
bool get_intersection(Point& ip,
                      bool& doubt,
                      general_form<GeneralType> const& p,
                      general_form<GeneralType> const& q)
{
    GeneralType const magnitude = (std::max)(p.magnitude(), q.magnitude());
    GeneralType const threshold = magnitude * magnitude
            * general_threshold<GeneralType>::get();
    GeneralType const denominator = p.b * q.a - p.a * q.b;
    GeneralType const abs_den = std::abs(denominator);

    doubt = false;
    bool result = true;
    if (abs_den < threshold)
    {
        doubt = denominator != 0; // If it is very small, set the doubt flag
        result = false;
    }

    // TODO: threshold is to be determined, and actually it should also
    // set a flag that intersection point is set.
    if (result || abs_den > threshold / 1000.0)
    {
        // Calculate y and x (even in the case of some doubt)
        geometry::set<1>(ip, (p.a * q.c - p.c * q.a) / denominator);
        geometry::set<0>(ip, (p.c * q.b - p.b * q.c) / denominator);
    }

    return result;
}

// TODO: verify if policy should be used here
template <typename Policy>
inline
bool lines_collinear(general_form<double> const& a,
                     general_form<double> const& b,
                     Policy const& policy)
{
    if (a.normalized && b.normalized)
    {
        bool const same_sign = more_horizontal(a)
                ? a.b * b.b > 0
                : a.a * b.a > 0;

        // c is the interception on x or y axis of normalized line
        // The normalized lign is still directed, if they have the same
        // direction (same_sign), check for intercept. If they are opposite,
        // then reverse one intercept
        return same_sign ? is_zero(a.c - b.c, policy)
                         : is_zero(a.c + b.c, policy)
                         ;

    }

    // Not (yet) implemented
    return false;
}

} // namespace arithmetic


}} // namespace boost::geometry


#endif // BOOST_GEOMETRY_ARITHMETIC_GENERAL_FORM_HPP
