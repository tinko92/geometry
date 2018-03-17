// Boost.Geometry

// Copyright (c) 2018 Barend Gehrels, Amsterdam, the Netherlands.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_ARITHMETIC_GENERAL_FORM_HPP
#define BOOST_GEOMETRY_ARITHMETIC_GENERAL_FORM_HPP

#include <algorithm>
#include <cmath>
#include <cstdlib>

#include <boost/array.hpp>
#include <boost/math/special_functions/round.hpp>

#include <boost/geometry/core/exception.hpp>

#include <boost/geometry/arithmetic/determinant.hpp>

#include <boost/geometry/util/math.hpp>
#include <boost/geometry/util/promote_integral.hpp>
#include <boost/geometry/util/select_calculation_type.hpp>

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

    bool more_horizontal() const
    {
        // If a=0, then line is 'by+c=0' so y=-c/b  which is horizontal
        // If a=0.1 and b=0.9, then the line is quite horizontal
        return std::abs(a) < std::abs(b);
    }

    bool has_horizontal_component() const
    {
        return std::abs(a) < 1e-6 || std::abs(b / a) >= 0.1;
    }

    bool has_vertical_component() const
    {
        return std::abs(b) < 1e-6 || std::abs(a / b) >= 0.1;
    }

    GeneralType magnitude() const
    {
        // TODO: store this
        return (std::max)(std::abs(a), std::abs(b));
    }

    // Returns true for two lines which are supposed to be (close to) collinear
    // (which is not checked) and have a similar direction
    // (in practice up to 45 degrees, TO BE VERIFIED)
    // true: -----------------> p -----------------> q
    // false: -----------------> p <----------------- q
    bool similar_direction(const general_form<GeneralType>& other) const
    {
        return a * other.a >= 0 && b * other.b >= 0;
    }

    // Returns a comparable and sortable distance measure:
    // - if left then positive,
    // - if collinear than zero - but NOT compared with any epsilon
    // - if right then negative
    // SQRT need not to be taken, therefore numerator is squared, but sign
    // is presever
    template <typename CoordinateType>
    typename select_most_precise<GeneralType, CoordinateType>::type
    distance_measure(CoordinateType const& x, CoordinateType const& y) const
    {
        // Distance from point to line in general form is given as:
        // (a * x + b * y + c) / sqrt(a * a + b * b);
        // In most use cases comparisons are enough, saving the sqrt
        // (performance plus making it a bit more precise)
        typedef typename select_most_precise<GeneralType, CoordinateType>::type
                ct;

        const ct num = a * x + b * y + c;
        if (num == 0)
        {
            return 0;
        }

        const GeneralType denom = a * a + b * b; // can be precalculated
        BOOST_ASSERT(denom != 0);
        return (num > 0 ? 1.0 : -1.0) * num * num / denom;
    }

    //! Normalize the line. For robustness reasons it is often better to
    //! NOT use normalization. It uses sqrt and therefore the intersection
    //! point, if calculated with normalization, might go,
    //! for example, from 7 to 6.99...997
    template <typename T>
    general_form<T> normalize() const
    {
        general_form<T> result;
        result.a = a;
        result.b = b;
        result.c = c;

        double const norm = std::sqrt(result.a * result.a + result.b * result.b);
        if (norm != 0) // TODO
        {
            result.a /= norm;
            result.b /= norm;
            result.c /= norm;
            result.normalized = true;
        }
        return result;
    }

};

template <typename GeneralType, typename CoordinateType>
inline
general_form<GeneralType> construct_line(CoordinateType const& x1, CoordinateType const& y1,
                               CoordinateType const& x2, CoordinateType const& y2)
{
    general_form<GeneralType> result;
    result.a = y1 - y2;
    result.b = x2 - x1;
    result.c = -result.a * x1 - result.b * y1;
    return result;
}

template <typename T, typename Point>
inline
general_form<T> construct_line(Point const& a, Point const& b)
{
    return construct_line<T>(geometry::get<0>(a),
                             geometry::get<1>(a),
                             geometry::get<0>(b),
                             geometry::get<1>(b));
}

template <typename T, typename Segment>
inline
general_form<T> construct_line(Segment const& segment)
{
    return construct_line<T>(geometry::get<0, 0>(segment),
                             geometry::get<0, 1>(segment),
                             geometry::get<1, 0>(segment),
                             geometry::get<1, 1>(segment));
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
    GeneralType const threshold = magnitude * magnitude * general_threshold<GeneralType>::get();
    GeneralType const denominator = p.b * q.a - p.a * q.b;
    GeneralType const abs_den = std::abs(denominator);

    doubt = false;
    bool result = true;
    if (abs_den < threshold)
    {
        doubt = denominator != 0; // If it is very small, set the doubt flag
        result = false;
    }

//    std::cout << threshold << " " << abs_den << " -> ";

    if (result || abs_den > threshold / 1000.0)
    {
        // Calculate y and x (even in the case of some doubt)
        geometry::set<1>(ip, (p.a * q.c - p.c * q.a) / denominator);
        geometry::set<0>(ip, (p.c * q.b - p.b * q.c) / denominator);
    }

    return result;
}

template <typename Policy>
static inline
bool lines_collinear(general_form<double> const& a,
                     general_form<double> const& b,
                     Policy const& policy)
{
    if (a.normalized && b.normalized)
    {
        bool const same_sign = a.more_horizontal()
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
