// Boost.Geometry

// Copyright (c) 2018 Barend Gehrels, Amsterdam, the Netherlands.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_STRATEGIES_CARTESIAN_GENERAL_INTERSECTION_HPP
#define BOOST_GEOMETRY_STRATEGIES_CARTESIAN_GENERAL_INTERSECTION_HPP

#include <algorithm>
#include <cmath>
#include <cstdlib>

#include <boost/array.hpp>
#include <boost/math/special_functions/round.hpp>

#include <boost/geometry/core/exception.hpp>

#include <boost/geometry/geometries/concepts/point_concept.hpp>
#include <boost/geometry/geometries/concepts/segment_concept.hpp>

#include <boost/geometry/arithmetic/determinant.hpp>
#include <boost/geometry/arithmetic/general_form.hpp>
#include <boost/geometry/algorithms/detail/assign_values.hpp>
#include <boost/geometry/algorithms/detail/assign_indexed_point.hpp>
#include <boost/geometry/algorithms/detail/equals/point_point.hpp>
#include <boost/geometry/algorithms/detail/recalculate.hpp>

#include <boost/geometry/util/math.hpp>
#include <boost/geometry/util/promote_integral.hpp>
#include <boost/geometry/util/select_calculation_type.hpp>

#include <boost/geometry/strategies/cartesian/area.hpp>
#include <boost/geometry/strategies/cartesian/distance_pythagoras.hpp>
#include <boost/geometry/strategies/cartesian/envelope_segment.hpp>
#include <boost/geometry/strategies/cartesian/point_in_poly_winding.hpp>
#include <boost/geometry/strategies/cartesian/side_by_triangle.hpp>
#include <boost/geometry/strategies/covered_by.hpp>
#include <boost/geometry/strategies/intersection.hpp>
#include <boost/geometry/strategies/intersection_result.hpp>
#include <boost/geometry/strategies/side.hpp>
#include <boost/geometry/strategies/side_info.hpp>
#include <boost/geometry/strategies/within.hpp>

#include <boost/geometry/policies/robustness/robust_point_type.hpp>
#include <boost/geometry/policies/robustness/segment_ratio_type.hpp>


namespace boost { namespace geometry
{

namespace strategy { namespace side
{

template <typename CalculationType = void>
struct side_by_generic_form
{

    template <typename P1, typename P2, typename P>
    static inline int apply(P1 const& p1, P2 const& p2, P const& p)
    {
        arithmetic::general_form<double> form = arithmetic::construct_line<double>(p1, p2);
        double dist = arithmetic::signed_comparable_distance(form, get<0>(p), get<1>(p));

        // TODO remove policy
        typedef math::detail::equals_factor_policy<double> equal_policy_type;
        equal_policy_type policy(geometry::get<0>(p1), geometry::get<1>(p1),
                       geometry::get<0>(p2), geometry::get<1>(p2));
        return math::detail::equals_by_policy(dist, 0, policy) ? 0
            : dist > 0 ? 1
            : -1;

    }

    template <typename P1, typename P2, typename P>
    static inline double signed_comparable_distance(P1 const& p1, P2 const& p2, P const& p)
    {
        arithmetic::general_form<double> const form = arithmetic::construct_line<double>(p1, p2);
        return arithmetic::signed_comparable_distance(form, get<0>(p), get<1>(p));
    }

};


}}

namespace strategy { namespace intersection
{


// Structure containing thresholds, per type, for denominator.
// Determined with corresponding unit test.
// It should not be replaced by machine epsilon or math::equals
template <typename Type>
struct general_distance_threshold {};

template <>
struct general_distance_threshold<long double>
{
   static long double get() { return 1.0e-10; }
};

template <>
struct general_distance_threshold<double>
{
   static double get() { return 1.0e-10; }
};

template <>
struct general_distance_threshold<float>
{
   static float get() { return 1.0e-6; }
};


// This intersection strategy is based
// on the general form of a line: ax + by + c
// Advantages above the old version (Kramer's rule) are:
// 1) no sides necessary (which makes it way more robust)
// 2) the code is simpler
// 3) it adds intersection of lines in terms of general form to the library
// 4) less code, should be faster
template
<
    typename CalculationType = void
>
struct cartesian_general_segments
{
    typedef side::side_by_generic_form<CalculationType> side_strategy_type;

    static inline side_strategy_type get_side_strategy()
    {
        return side_strategy_type();
    }

    template <typename Geometry1, typename Geometry2>
    struct point_in_geometry_strategy
    {
        typedef strategy::within::cartesian_winding
            <
                typename point_type<Geometry1>::type,
                typename point_type<Geometry2>::type,
                CalculationType
            > type;
    };

    template <typename Geometry1, typename Geometry2>
    static inline typename point_in_geometry_strategy<Geometry1, Geometry2>::type
        get_point_in_geometry_strategy()
    {
        typedef typename point_in_geometry_strategy
            <
                Geometry1, Geometry2
            >::type strategy_type;
        return strategy_type();
    }

    template <typename Geometry>
    struct area_strategy
    {
        typedef area::cartesian
            <
                CalculationType
            > type;
    };

    template <typename Geometry>
    static inline typename area_strategy<Geometry>::type get_area_strategy()
    {
        typedef typename area_strategy<Geometry>::type strategy_type;
        return strategy_type();
    }

    template <typename Geometry>
    struct distance_strategy
    {
        typedef distance::pythagoras
            <
                CalculationType
            > type;
    };

    template <typename Geometry>
    static inline typename distance_strategy<Geometry>::type get_distance_strategy()
    {
        typedef typename distance_strategy<Geometry>::type strategy_type;
        return strategy_type();
    }

    typedef envelope::cartesian_segment<CalculationType>
        envelope_strategy_type;

    static inline envelope_strategy_type get_envelope_strategy()
    {
        return envelope_strategy_type();
    }

    //--------------------------------------------------------------------------

    template <std::size_t Index, typename Point, typename Segment>
    static inline typename geometry::coordinate_type<Point>::type
            comparable_distance(Point const& p, Segment const& s)
    {
        typedef typename geometry::coordinate_type<Point>::type ctype;
        ctype const dx = get<0>(p) - get<Index, 0>(s);
        ctype const dy = get<1>(p) - get<Index, 1>(s);
        return dx * dx + dy * dy;
    }

    template
    <
        std::size_t Index1,
        std::size_t Index2,
        typename Segment1,
        typename Segment2
    >
    static inline
    typename geometry::coordinate_type<Segment1>::type
    comparable_distance_segseg(Segment1 const& a, Segment2 const& b)
    {
        typedef typename geometry::coordinate_type<Segment1>::type ct;
        ct const dx = get<Index1, 0>(a) - get<Index2, 0>(b);
        ct const dy = get<Index1, 1>(a) - get<Index2, 1>(b);
        return dx * dx + dy * dy;
    }


    struct side_data
    {
        double value;
        double abs_value;
        bool on_end;
        int side;  // -1 (right), 0 (collinear) or 1 (left)

        side_data()
            : value(0)
            , abs_value(0)
            , on_end(false)
            , side(0)
        {}
    };


    template
    <
        typename Data,
        typename Segment1,
        typename Segment2,
        typename T
     >
    static inline
    void initialize_signed_comparable_distances(Data& data,
               Segment1 const& a, Segment2 const& b,
               arithmetic::general_form<T> const& gfa,
               arithmetic::general_form<T> const& gfb)
    {
        data[0].value = arithmetic::signed_comparable_distance(gfb, geometry::get<0, 0>(a), geometry::get<0, 1>(a)); // a[0] w.r.t. b
        data[1].value = arithmetic::signed_comparable_distance(gfb, geometry::get<1, 0>(a), geometry::get<1, 1>(a)); // a[1] w.r.t. b
        data[2].value = arithmetic::signed_comparable_distance(gfa, geometry::get<0, 0>(b), geometry::get<0, 1>(b)); // b[0] w.r.t. a
        data[3].value = arithmetic::signed_comparable_distance(gfa, geometry::get<1, 0>(b), geometry::get<1, 1>(b)); // b[1] w.r.t. a

        for (int i = 0; i < 4; i++)
        {
            data[i].abs_value = std::fabs(data[i].value);
        }
    }

    template
    <
        typename Data
    >
    static inline
    void initialize_on_end(Data& data, double const& fa, double const& fb)
    {
        data[0].on_end = fa == 0;
        data[1].on_end = fa == 1;
        data[2].on_end = fb == 0;
        data[3].on_end = fb == 1;
    }

    template
    <
        std::size_t Index1,
        std::size_t Index2,
        typename Data,
        typename Segment1,
        typename Segment2
     >
    static inline
    void try_fix_two_common(Data& data, double& fa, double& fb,
                    Segment1 const& a, Segment2 const& b)
    {
        typedef typename geometry::coordinate_type<Segment1>::type ct;
        double const d = comparable_distance_segseg<Index1, Index2>(a, b);
        if (d < general_distance_threshold<ct>::get())
        {
            fa = Index1;
            fb = Index2;
            data[Index1].on_end = true;
            data[2 + Index2].on_end = true;
        }
    }

    template
    <
        typename Data,
        typename Segment1,
        typename Segment2
     >
    static inline
    void fix_two_common_endpoints(Data& data, double& fa, double& fb,
               Segment1 const& a, Segment2 const& b, int fa_i, int fb_i)
    {
        int on_end_counter = 0;
        for (int i = 0; i < 4; i++)
        {
            if (data[i].on_end)
            {
                on_end_counter++;
            }
        }
        if (on_end_counter >= 2)
        {
            return;
        }

        if (fa_i == 0 && fb_i == 0)
        {
            try_fix_two_common<0, 0>(data, fa, fb, a, b);
        }
        else if (fa_i == 0 && fb_i == 1)
        {
            try_fix_two_common<0, 1>(data, fa, fb, a, b);
        }
        else if (fa_i == 1 && fb_i == 0)
        {
            try_fix_two_common<1, 0>(data, fa, fb, a, b);
        }
        else if (fa_i == 1 && fb_i == 1)
        {
            try_fix_two_common<1, 1>(data, fa, fb, a, b);
        }
    }

    template
    <
        typename Data,
        typename Segment1,
        typename Segment2,
        typename T

     >
    static inline
    void fix_one_common_endpoint(Data& data, double& fa, double& fb,
                Segment1 const& a, Segment2 const& b,
                arithmetic::general_form<T> const& gf_a,
                arithmetic::general_form<T> const& gf_b,
                int fa_i, int fb_i)
    {
        typedef typename geometry::coordinate_type<Segment1>::type ct;
        int index = -1;
        int on_end_counter = 0;
        for (int i = 0; i < 4; i++)
        {
            if (data[i].on_end)
            {
                on_end_counter++;
            }
            else if (index == -1
                     && data[i].abs_value < general_distance_threshold<ct>::get())
            {
                index = i;
            }
        }

        if (on_end_counter >= 2 || index == -1)
        {
            return;
        }
        if (index == 0 && fa_i == 0)
        {
            data[index].on_end = true;
            fa = fa_i;
        }
        if (index == 1 && fa_i == 1)
        {
            data[index].on_end = true;
            fa = fa_i;
        }
        if (index == 2 && fb_i == 0)
        {
            data[index].on_end = true;
            fb = fb_i;
        }
        if (index == 3 && fb_i == 1)
        {
            data[index].on_end = true;
            fb = fb_i;
        }
    }

    static inline
    bool suspicious(int& fa_i, int& fb_i, double const& fa, double const& fb)
    {
        fa_i = boost::math::iround(fa);
        fb_i = boost::math::iround(fb);

        bool const perfect = (fa == 0 || fa == 1) && (fb == 0 || fb == 1);
        if (perfect)
        {
            return false;
        }

        bool const possible = (fa_i == 0 || fa_i == 1) && (fb_i == 0 || fb_i == 1);
        if (! possible)
        {
            return false;
        }

        double const threshold = 1.0e-3; // TODO
        if (fa_i == 0 && std::fabs(fa) < threshold)       { return true; }
        if (fa_i == 1 && std::fabs(fa - 1.0) < threshold) { return true; }
        if (fb_i == 0 && std::fabs(fb) < threshold)       { return true; }
        if (fb_i == 1 && std::fabs(fb - 1.0) < threshold) { return true; }

        return false;
    }

    //--------------------------------------------------------------------------

    template <std::size_t Dimension, typename Point, typename Segment>
    static inline
    bool on_segment(Point const& point, Segment const& segment,
                    double& fraction, bool& doubt)
    {
        doubt = false;
        typedef typename geometry::coordinate_type<Point>::type coordinate_type;

        coordinate_type const c1 = geometry::get<0, Dimension>(segment);
        coordinate_type const c2 = geometry::get<1, Dimension>(segment);
        coordinate_type const c = geometry::get<Dimension>(point);

        if (c1 == c)
        {
            fraction = 0.0;
            return true;
        }
        else if (c2 == c)
        {
            fraction = 1.0;
            return true;
        }

        bool const increasing = c1 < c2;
        bool result = (  increasing && c > c1 && c < c2)
                   || (! increasing && c > c2 && c < c1);

        coordinate_type const length = increasing ? c2 - c1 : c1 - c2;

        if (length == 0)
        {
            // Should not occur - segment is degenerate, point is on segment
            // but not equal. This is not possible. But avoid division by zero.
            fraction = 0;
            return false;
        }

        if (increasing)
        {
            coordinate_type const pos = c - c1;
            fraction = pos / static_cast<double>(length);
        }
        else
        {
            coordinate_type const pos = c - c2;
            fraction = 1.0 - pos / static_cast<double>(length);
        }

        bool near_zero = false;
        bool near_one = false;
        if (! result)
        {
            near_zero = std::fabs(fraction) < 0.01;
            near_one = ! near_zero && std::fabs(fraction - 1.0) < 0.01;
        }

        doubt = ! result && (near_zero || near_one);

        if (doubt)
        {
            double const d
                = near_zero ? comparable_distance<0>(point, segment)
                : near_one ? comparable_distance<1>(point, segment)
                : 1.0;

            if (d < general_distance_threshold<coordinate_type>::get())
            {
                result = true;
                doubt = false;
                fraction = near_zero ? 0.0 : near_one ? 1.0 : -1;
            }
        }

        return result;
    }

    template <typename P, typename S, typename T>
    static inline
    bool on_segment(P const& p, S const& s,
                    arithmetic::general_form<T> const& f,
                    double& fraction, bool& doubt)
    {
        // Depending on if line is more vertical than horizontal (a>b), compare
        // y or x direction
        bool const horizontalish = arithmetic::more_horizontal(f);
        bool result = horizontalish
                ? on_segment<0>(p, s, fraction, doubt)
                : on_segment<1>(p, s, fraction, doubt);

        if (! result && doubt)
        {
            // The intersection point is NOT on the segment, but it is close.
            // This might be a numerical error. First try if it can be
            // established in the other direction.
            if (horizontalish && arithmetic::has_vertical_component(f))
            {
                if (on_segment<1>(p, s, fraction, doubt))
                {
                    fraction = static_cast<double>(boost::math::iround(fraction));
                    return true;
                }
            }
            if (! horizontalish && arithmetic::has_horizontal_component(f))
            {
                if (on_segment<0>(p, s, fraction, doubt))
                {
                    fraction = static_cast<double>(boost::math::iround(fraction));
                    return true;
                }
            }
        }

        return result;
    }

    template
    <
        std::size_t Index1,
        std::size_t Index2,
        typename Segment1,
        typename Segment2
    >
    static inline
    bool segment_equals(Segment1 const& a, Segment2 const& b)
    {
        typedef typename geometry::coordinate_type<Segment1>::type coordinate_type;

        coordinate_type const xa = geometry::get<Index1, 0>(a);
        coordinate_type const ya = geometry::get<Index1, 1>(a);

        coordinate_type const xb = geometry::get<Index2, 0>(b);
        coordinate_type const yb = geometry::get<Index2, 1>(b);

        return xa == xb && ya == yb;
    }

    template
    <
        typename Segment1,
        typename Segment2
    >
    static inline
    bool have_common_endpoints(Segment1 const& a, Segment2 const& b,
                               double& fraction_a, double& fraction_b)
    {
        // TODO: we can use the current values of imprecise fractions here
        // to have only one comparison
        boost::array<double, 4> dist;
        dist[0] = comparable_distance_segseg<0, 0>(a, b);
        dist[1] = comparable_distance_segseg<0, 1>(a, b);
        dist[2] = comparable_distance_segseg<1, 0>(a, b);
        dist[3] = comparable_distance_segseg<1, 1>(a, b);

        int index = -1;
        if (dist[0] <= dist[1] && dist[0] <= dist[2] && dist[0] <= dist[3])
        {
            index = 0;
        }
        else if (dist[1] <= dist[0] && dist[1] <= dist[2] && dist[1] <= dist[3])
        {
            index = 1;
        }
        else if (dist[2] <= dist[0] && dist[2] <= dist[1] && dist[2] <= dist[3])
        {
            index = 2;
        }
        else
        {
            index = 3;
        }

        if (index == 0 && segment_equals<0, 0>(a, b))
        {
            fraction_a = 0.0;
            fraction_b = 0.0;
            return true;
        }
        if (index <= 1 && segment_equals<0, 1>(a, b))
        {
            fraction_a = 0.0;
            fraction_b = 1.0;
            return true;
        }
        if (index <= 2 && segment_equals<1, 0>(a, b))
        {
            fraction_a = 1.0;
            fraction_b = 0.0;
            return true;
        }
        if (index <= 3 && segment_equals<1, 1>(a, b))
        {
            fraction_a = 1.0;
            fraction_b = 1.0;
            return true;
        }

        return false;
    }


    template
    <
        int Dimension,
        typename Segment1,
        typename Segment2
    >
    static inline
    bool disjoint_in_dimension(Segment1 const& a, Segment2 const& b)
    {
        typedef typename geometry::coordinate_type<Segment1>::type ct;
        ct a0 = geometry::get<0, Dimension>(a);
        ct a1 = geometry::get<1, Dimension>(a);
        ct b0 = geometry::get<0, Dimension>(b);
        ct b1 = geometry::get<1, Dimension>(b);

        if (a0 > a1) { std::swap(a0, a1); }
        if (b0 > b1) { std::swap(b0, b1); }

        return (a0 < b0 && a1 < b0)
            || (a0 > b1 && a1 > b1);
    }

    template
    <
        typename Segment1,
        typename Segment2,
        typename T
    >
    static inline
    bool disjoint_by_coordinates(Segment1 const& a, Segment2 const& b,
                      arithmetic::general_form<T> const& f)
    {
        return (arithmetic::has_horizontal_component(f) && disjoint_in_dimension<0>(a, b))
            || (arithmetic::has_vertical_component(f)   && disjoint_in_dimension<1>(a, b));
    }

    template
    <
        typename Data,
            typename T
     >
    static inline
    bool disjoint_by_side(Data const& data, T const& threshold)
    {
        for (int i = 0; i <= 2; i += 2)
        {
            if (data[i].value * data[i + 1].value > 0
                    && data[i].abs_value > threshold
                    && data[i + 1].abs_value > threshold)
            {
                return true;
            }
        }
        return false;
    }


    template
    <
        typename Data,
        typename Segment1,
        typename Segment2,
        typename T
    >
    static inline
    void inspect_sides(Data const& data,
                       bool& consider_as_collinear,
                       bool& consider_as_crossing,
                       bool& consider_as_disjoint,
                       Segment1 const& a, Segment2 const& b,
                       arithmetic::general_form<T> const& gf_a,
                       arithmetic::general_form<T> const& gf_b)
    {
        consider_as_collinear = false;
        consider_as_disjoint = false;
        consider_as_disjoint = false;

        // Calculate distance measures (they will all be about the same
        // for 'nearly' collinear lines, but their sign can differ)
        double const dm_a0 = arithmetic::signed_comparable_distance(gf_b, geometry::get<0, 0>(a), geometry::get<0, 1>(a)); // a[0]
        double const dm_a1 = arithmetic::signed_comparable_distance(gf_b, geometry::get<1, 0>(a), geometry::get<1, 1>(a)); // a[1]

        // If product is positive, then A is located at one side of B (and vv)
        // If product is 0, then one endpoint, or both, are truely collinear
        // If product A and B are negative, then A crosses B
        // If one of product A and B is negative,
        //   then (if one of the measures is 0) A touches B (or vv),
        //   or they are just disjoint but they infinite lines are crosing
        double const side_a = dm_a0 * dm_a1;

        double const dm_b0 = arithmetic::signed_comparable_distance(gf_a, geometry::get<0, 0>(b), geometry::get<0, 1>(b)); // b[0]
        double const dm_b1 = arithmetic::signed_comparable_distance(gf_a, geometry::get<1, 0>(b), geometry::get<1, 1>(b)); // b[1]
        double const side_b = dm_b0 * dm_b1;

        // The lines are nearly collinear.
        // They are not (nearly) collinear disjoint
        //   ( in this sense: -----a---->   ------b------> )
        // If sides consequently indicate one side of the other, only then they
        // are considered as disjoint. But if one of them is zero, then we
        // consider them as collinear.

        if (std::fabs(side_a * side_b) > 0.1) // TODO THRESHOLD
        {
            // Not even close to collinear
            // Caused by imprecision of get_intersection_point
            return;
        }
        consider_as_crossing = side_a < 0 && side_b < 0;
        if (consider_as_crossing)
        {
            // OCCURS IN buffer_polygon.cpp, flower35, and more

            // It still can be considered as collinear too
        }
        // TODO: consider_as_touching

        double const dm_abs_a0 = std::fabs(dm_a0);
        double const dm_abs_a1 = std::fabs(dm_a1);
        double const dm_abs_b0 = std::fabs(dm_b0);
        double const dm_abs_b1 = std::fabs(dm_b1);

        double const dm_max_a = std::max(dm_abs_a0, dm_abs_a1);
        double const dm_max_b = std::max(dm_abs_b0, dm_abs_b1);

        double const dm_separation = std::max(dm_max_a, dm_max_b);

        consider_as_disjoint = dm_separation > 1.0e-12; // TODO - use distance here 'magic' distance
        consider_as_collinear = ! consider_as_disjoint;
    }


    template
    <
        typename Segment1,
        typename Segment2,
        typename Point,
        typename T
    >
    static inline
    bool get_fractions_on_segment(Segment1 const& a, Segment2 const& b,
                                  Point const& p,
                                  arithmetic::general_form<T> const& gfa,
                                  arithmetic::general_form<T> const& gfb,
                                  double& fraction_a, double& fraction_b)
    {
        typedef typename geometry::coordinate_type<Point>::type ct;

        bool doubt_a = false;
        bool doubt_b = false;
        bool on_a = on_segment(p, a, gfa, fraction_a, doubt_a);
        bool on_b = on_segment(p, b, gfb, fraction_b, doubt_b);

        if (on_a && on_b)
        {
            return true;
        }

        if (! doubt_a && ! doubt_b)
        {
            return false;
        }

        if (on_a && doubt_b)
        {

            // IP considered on one of the segment, but not on the other, with
            // doubt. Verify with sides.
            int fb = boost::math::iround(fraction_b);
            if (fb == 0)
            {
                double const value = arithmetic::signed_comparable_distance(gfa, geometry::get<0, 0>(b), geometry::get<0, 1>(b)); // b[0] vs a

                if (std::fabs(value) < general_distance_threshold<ct>::get())
                {
                    fraction_b = 0.0;
                    return true;
                }
            }
            if (fb == 1)
            {
                double const value = arithmetic::signed_comparable_distance(gfa, geometry::get<1, 0>(b), geometry::get<1, 1>(b)); // b[1] vs a
                if (std::fabs(value) < general_distance_threshold<ct>::get())
                {
                    fraction_b = 1.0;
                    return true;
                }
            }
        }
        if (on_b && doubt_a)
        {
            // IP considered on one of the segment, but not on the other, with
            // doubt. Verify with sides.
            int fa = boost::math::iround(fraction_a);
            if (fa == 0)
            {
                double const value = arithmetic::signed_comparable_distance(gfb, geometry::get<0, 0>(a), geometry::get<0, 1>(a)); // a[0] vs b
                if (std::fabs(value) < general_distance_threshold<ct>::get())
                {
                    fraction_a = 0.0;
                    return true;
                }
            }
            if (fa == 1)
            {
                double const value = arithmetic::signed_comparable_distance(gfb, geometry::get<1, 0>(a), geometry::get<1, 1>(a)); // a[1] vs b
                if (std::fabs(value) < general_distance_threshold<ct>::get())
                {
                    fraction_a = 1.0;
                    return true;
                }
            }
        }

        // IP not considered as on segment, but there is doubt if the IP is
        // on the end point of at least one of the segments.
        // Verify if the segments themselves have the same point (which is not
        // verified earlier). If so, correct fractions (which were apparently
        // imprecise earlier)

        return have_common_endpoints(a, b, fraction_a, fraction_b);
    }

    // TODO: this should be cleaned up
    template <bool Floating>
    struct side_assorter
    {

        template
        <
            typename Data,
            typename Segment1,
            typename Segment2,
            typename T
         >
        static inline
        void initialize(Data& data, Segment1 const& a, Segment2 const& b,
                   arithmetic::general_form<T> const& gfa,
                   arithmetic::general_form<T> const& gfb,
                   double const& fa,
                   double const& fb)
        {
            // TODO this should be removed
            for (int i = 0; i < 4; i++)
            {
                data[i].side = data[i].on_end
                        ? 0 : data[i].value == 0 ? 99 : data[i].value > 0.0 ? 1 : -1;
            }
            // It can happen that it was NOT on end and still value was absolutely 0,
            // and therefore 'collinear'.
            // In that case try to make it opposite.
            int other_index = 1; // 1, 0, 3, 2
            for (int i = 0; i < 4; i++)
            {
                if (data[i].side == 99 && data[other_index].side != 99)
                {
                    data[i].side = -data[other_index].side;
                }
                if (other_index == 0)
                {
                    other_index = 3;
                }
                else
                {
                    other_index--;
                }
            }
        }

        template <typename Data>
        static inline
        double get_threshold(Data const& data)
        {
            double threshold = 0;
            for (int i = 0; i < 4; i++)
            {
                if (data[i].on_end)
                {
                    if (data[i].abs_value > threshold)
                    {
                        threshold = data[i].abs_value;
                    }
                }
            }
            return threshold;
        }

        template <typename Data>
        static inline
        bool is_consistent(Data const& data)
        {
            // Get threshold (~ epsilon)
            double const threshold = get_threshold(data);

            for (int i = 0; i < 4; i++)
            {
                if (! data[i].on_end && data[i].abs_value <= threshold)
                {
                    // Not consistent: it is not marked as an end-point
                    // (so fraction != 0 && fraction != 1, so the IP
                    // is somewhere in between of the segment
                    std::cout << "Inconsistent threshold " << std::endl;
                    return false;
                }
            }

            for (int i = 0; i < 4; i += 2)
            {
                if (data[i].side != 0 && data[i + 1].side != 0)
                {
                    if (data[i].side != -data[i + 1].side)
                    {
                        // Both points of A on same side of B. Impossible, because
                        // this code is only executed if there is an intersection point
                        return false;
                    }
                }
            }

            if (   (data[0].side == 0 && data[1].side == 0)
                || (data[2].side == 0 && data[3].side == 0))
            {
                // TODO
                return false;
            }

            return true;
        }


        template
        <
            typename Data,
            typename Segment1,
            typename Segment2,
            typename T
         >
        static inline
        bool do_touch(Data& data, Segment1 const& a, Segment2 const& b,
                      arithmetic::general_form<T> const& gfa,
                      arithmetic::general_form<T> const& gfb,
                      double& fa,
                      double& fb)
        {
            // Get threshold (~ epsilon)
            double const threshold = get_threshold(data);

            std::size_t count = 0;
            std::size_t index = 9;
            std::size_t wrong_count = 0;
            std::size_t wrong_index = 9;
            for (int i = 0; i < 4; i++)
            {
                if (data[i].on_end)
                {
                    count++;
                    index = i;
                }
                else if (data[i].abs_value <= threshold)
                {
                    wrong_count++;
                    wrong_index = i;

                }
            }
            if (count != 1 || wrong_count != 1)
            {
                return false;
            }

            double const eps = 0.01; // TODO
            double const eps_f = 0.1; // TODO

            // There are 8 possibilities
            // TODO, clean this up
            if (index == 0 && wrong_index == 2
                    && std::fabs(fb) < eps_f
                    && comparable_distance_segseg<0, 0>(a, b) < eps)
            {
                fb = 0.0;
                data[wrong_index].on_end = true;
                return true;
            }
            if (index == 0 && wrong_index == 3
                    && std::fabs(fb - 1.0) < eps_f
                    && comparable_distance_segseg<0, 1>(a, b) < eps)
            {
                fb = 1.0;
                data[wrong_index].on_end = true;
                return true;
            }
            if (index == 1 && wrong_index == 2
                    && std::fabs(fb) < eps_f
                    && comparable_distance_segseg<1, 0>(a, b) < eps)
            {
                fb = 0.0;
                data[wrong_index].on_end = true;
                return true;
            }
            if (index == 1 && wrong_index == 3
                    && std::fabs(fb - 1.0) < eps_f
                    && comparable_distance_segseg<1, 1>(a, b) < eps)
            {
                fb = 1.0;
                data[wrong_index].on_end = true;
                return true;
            }

            if (index == 2 && wrong_index == 0
                    && std::fabs(fa) < eps_f
                    && comparable_distance_segseg<0, 0>(a, b) < eps)
            {
                fa = 0.0;
                data[wrong_index].on_end = true;
                return true;
            }
            if (index == 2 && wrong_index == 1
                    && std::fabs(fa - 1.0) < eps_f
                    && comparable_distance_segseg<1, 0>(a, b) < eps)
            {
                fa = 1.0;
                data[wrong_index].on_end = true;
                return true;
            }
            if (index == 3 && wrong_index == 0
                    && std::fabs(fa) < eps_f
                    && comparable_distance_segseg<0, 1>(a, b) < eps)
            {
                fa = 0.0;
                data[wrong_index].on_end = true;
                return true;
            }
            if (index == 3 && wrong_index == 1
                    && std::fabs(fa - 1.0) < eps_f
                    && comparable_distance_segseg<1, 1>(a, b) < eps)
            {
                fa = 1.0;
                data[wrong_index].on_end = true;
                return true;
            }

            return false;
        }

        template
        <
            typename Data,
            typename Segment1,
            typename Segment2,
            typename T
         >
        static inline
        bool apply(side_info& sides,
                   double& fa,
                   double& fb,
                   bool& consider_as_collinear,
                   Data& data,
                   Segment1 const& a, Segment2 const& b,
                   arithmetic::general_form<T> const& gf_a,
                   arithmetic::general_form<T> const& gf_b)
        {
            // Sides w.r.t. intersection point IP (*)
            //
            //          b[1]           b[1]~0        -> v[3]
            // a[0]------*-------a[1]  a[0]=left, a[1]=right -> v[0], v[1]
            //           |
            //           |
            //           b[0]          b[0]=right    -> v[2]

            // If fraction of IP on b is one, this means that IP == b[1],
            // this means that b[1] is collinear with a
            // So, in terms used below, side_value(b[1]) ~ 0
            // -> values[3] ~ 0 and on_end[3] == true

            initialize(data, a, b, gf_a, gf_b, fa, fb);

            bool consistent = is_consistent(data);
            if (!consistent)
            {
                if (!consistent
                    && do_touch(data, a, b, gf_a, gf_b, fa, fb))
                {
                    initialize(data, a, b, gf_a, gf_b, fa, fb);
                    consistent = is_consistent(data);
                }

                if (!consistent)
                {
                    bool consider_as_crossing;
                    bool consider_as_disjoint;

                    inspect_sides(data, consider_as_collinear, consider_as_crossing,
                            consider_as_disjoint, a, b, gf_a, gf_b);

                    if (consider_as_collinear)
                    {
                        return false;
                    }

                    // debug only
                    inspect_sides(data, consider_as_collinear, consider_as_crossing,
                            consider_as_disjoint, a, b, gf_a, gf_b);
                    is_consistent(data);
                    have_common_endpoints(a, b, fa, fb);
                    do_touch(data, a, b, gf_a, gf_b, fa, fb);
                    // end debug

                    return false;
                }

            }

            sides.set<0>(data[0].side, data[1].side);
            sides.set<1>(data[2].side, data[3].side);

            return true;
        }
    };


    //--------------------------------------------------------------------------
    template <typename PointType, typename SegmentRatio>
    struct segment_intersection_info
    {
        template <typename Point, typename Segment1, typename Segment2>
        void calculate(Point& point, Segment1 const& , Segment2 const& ) const
        {
            point = m_point;
        }

        PointType m_point;

        // To be removed:
        SegmentRatio robust_ra;
        SegmentRatio robust_rb;
    };

    // to be removed completely (only to create compatable segment_ratio)
    template <typename D, typename W, typename ResultType>
    static inline void cramers_rule(D const& dx_a, D const& dy_a,
        D const& dx_b, D const& dy_b, W const& wx, W const& wy,
        // out:
        ResultType& d, ResultType& da)
    {
        // Cramers rule
        d = geometry::detail::determinant<ResultType>(dx_a, dy_a, dx_b, dy_b);
        da = geometry::detail::determinant<ResultType>(dx_b, dy_b, wx, wy);
        // Ratio is da/d , collinear if d == 0, intersecting if 0 <= r <= 1
        // IntersectionPoint = (x1 + r * dx_a, y1 + r * dy_a)
    }

    //--------------------------------------------------------------------------
    // TODO, this is copied, and this can go
    template
    <
        typename Segment1,
        typename Segment2,
        typename Policy,
        typename RobustPolicy
    >
    static inline typename Policy::return_type
        apply(Segment1 const& a, Segment2 const& b,
              Policy const& policy, RobustPolicy const& robust_policy)
    {
        typedef typename geometry::point_type<Segment1>::type point_type;

        typedef typename geometry::robust_point_type
            <
                point_type, RobustPolicy
            >::type robust_point_type;

        point_type a0, a1, b0, b1;
        robust_point_type a0_rob, a1_rob, b0_rob, b1_rob;

        detail::assign_point_from_index<0>(a, a0);
        detail::assign_point_from_index<1>(a, a1);
        detail::assign_point_from_index<0>(b, b0);
        detail::assign_point_from_index<1>(b, b1);

        geometry::recalculate(a0_rob, a0, robust_policy);
        geometry::recalculate(a1_rob, a1, robust_policy);
        geometry::recalculate(b0_rob, b0, robust_policy);
        geometry::recalculate(b1_rob, b1, robust_policy);

        return apply(a, b, policy, robust_policy, a0_rob, a1_rob, b0_rob, b1_rob);
    }

#ifdef BOOST_GEOMETRY_USE_COMPLEX_SEGMENT_RATIO
    // TODO: if rescaling is gone, this can go too.
    template
    <
        typename Info,
        typename RobustPoint1,
        typename RobustPoint2
    >
    static inline void assign_segment_intersection_info(Info& sinfo,
              RobustPoint1 const& robust_a1, RobustPoint1 const& robust_a2,
              RobustPoint2 const& robust_b1, RobustPoint2 const& robust_b2)
    {
        typedef typename select_most_precise
            <
                typename geometry::coordinate_type<RobustPoint1>::type,
                typename geometry::coordinate_type<RobustPoint2>::type
            >::type robust_coordinate_type;

        robust_coordinate_type robust_da0, robust_da;
        robust_coordinate_type robust_db0, robust_db;

        robust_coordinate_type const robust_dx_a = get<0>(robust_a2) - get<0>(robust_a1);
        robust_coordinate_type const robust_dx_b = get<0>(robust_b2) - get<0>(robust_b1);
        robust_coordinate_type const robust_dy_a = get<1>(robust_a2) - get<1>(robust_a1);
        robust_coordinate_type const robust_dy_b = get<1>(robust_b2) - get<1>(robust_b1);

        cramers_rule(robust_dx_a, robust_dy_a, robust_dx_b, robust_dy_b,
            get<0>(robust_a1) - get<0>(robust_b1),
            get<1>(robust_a1) - get<1>(robust_b1),
            robust_da0, robust_da);

        cramers_rule(robust_dx_b, robust_dy_b, robust_dx_a, robust_dy_a,
            get<0>(robust_b1) - get<0>(robust_a1),
            get<1>(robust_b1) - get<1>(robust_a1),
            robust_db0, robust_db);

        sinfo.robust_ra.assign(robust_da, robust_da0);
        sinfo.robust_rb.assign(robust_db, robust_db0);
    }
#endif

    //--------------------------------------------------------------------------
    template
    <
        typename Policy,
        typename RatioType,
        typename Segment1,
        typename Segment2,
        typename T,
        typename Info,
        typename RobustPoint1,
        typename RobustPoint2
    >
    static inline
    typename Policy::return_type
    handle_common_endpoints(Segment1 const& a, Segment2 const& b,
                            arithmetic::general_form<T> const& gfa,
                            arithmetic::general_form<T> const& gfb,
                            double fa, double fb,
                            Info sinfo,
                            bool a_is_point, bool b_is_point,
                            RobustPoint1 const& robust_a1, RobustPoint1 const& robust_a2,
                            RobustPoint2 const& robust_b1, RobustPoint2 const& robust_b2)
    {
        typedef typename geometry::coordinate_type<Segment1>::type ct;

        // a-b share a common endpoint, fractions indicate which.
        // They intersect, possibly collinearly
        // ------a------*------b------ in picture, fai == 1 and fbi == 0
        const int fai = boost::math::iround(fa);
        const int fbi = boost::math::iround(fb);

        // Verify sides w.r.t. a
        const double value_a = fai == 0
            ? arithmetic::signed_comparable_distance(gfb, geometry::get<1, 0>(a), geometry::get<1, 1>(a))
            : arithmetic::signed_comparable_distance(gfb, geometry::get<0, 0>(a), geometry::get<0, 1>(a));
        const double value_b = fbi == 0
            ? arithmetic::signed_comparable_distance(gfa, geometry::get<1, 0>(b), geometry::get<1, 1>(b))
            : arithmetic::signed_comparable_distance(gfa, geometry::get<0, 0>(b), geometry::get<0, 1>(b));

         if (   std::fabs(value_a) < general_distance_threshold<ct>::get()
             || std::fabs(value_b) < general_distance_threshold<ct>::get())
        {
            // Collinear (if one of them is collinear, both are considered as such)
            return arithmetic::more_horizontal(gfa)
                    ? relate_collinear<0, Policy, RatioType>(a, b, robust_a1, robust_a2, robust_b1, robust_b2, a_is_point, b_is_point)
                    : relate_collinear<1, Policy, RatioType>(a, b, robust_a1, robust_a2, robust_b1, robust_b2, a_is_point, b_is_point);
        }

        // The segments make an angle w.r.t. each other.
        // Set intersection point
        if (fai == 0)
        {
            geometry::set<0>(sinfo.m_point, geometry::get<0, 0>(a));
            geometry::set<1>(sinfo.m_point, geometry::get<0, 1>(a));
        }
        else if (fai == 1)
        {
            geometry::set<0>(sinfo.m_point, geometry::get<1, 0>(a));
            geometry::set<1>(sinfo.m_point, geometry::get<1, 1>(a));
        }

#ifdef BOOST_GEOMETRY_USE_COMPLEX_SEGMENT_RATIO
        assign_segment_intersection_info(sinfo, robust_a1, robust_a2,
                                         robust_b1, robust_b2);
#else
        sinfo.robust_ra.assign(fa);
        sinfo.robust_rb.assign(fb);
#endif

        side_info sides;
        if (fai == 0)
        {
            sides.set<0>(0, value_a > 0 ? 1 : -1);
        }
        else
        {
            sides.set<0>(value_a > 0 ? 1 : -1, 0);
        }
        if (fbi == 0)
        {
            sides.set<1>(0, value_b > 0 ? 1 : -1);
        }
        else
        {
            sides.set<1>(value_b > 0 ? 1 : -1, 0);
        }

        return Policy::segments_crosses(sides, sinfo, a, b);
    }

    //--------------------------------------------------------------------------

    template
    <
        typename Segment1,
        typename Segment2,
        typename Policy,
        typename RobustPolicy,
        typename RobustPoint1,
        typename RobustPoint2
    >
    static inline typename Policy::return_type
        apply(Segment1 const& a, Segment2 const& b,
              Policy const&, RobustPolicy const& /*robust_policy*/,
              RobustPoint1 const& robust_a1, RobustPoint1 const& robust_a2,
              RobustPoint2 const& robust_b1, RobustPoint2 const& robust_b2)
    {
        static int counter = 0;
        counter++;

        BOOST_CONCEPT_ASSERT( (concepts::ConstSegment<Segment1>) );
        BOOST_CONCEPT_ASSERT( (concepts::ConstSegment<Segment2>) );

        typedef typename select_most_precise
            <
                typename geometry::coordinate_type<Segment1>::type,
                typename geometry::coordinate_type<Segment2>::type
            >::type coordinate_type;

        typedef math::detail::equals_factor_policy<coordinate_type> equal_policy_type;

        equal_policy_type equal_policy(geometry::get<0,0>(a),
                       geometry::get<1,1>(a),
                       geometry::get<0,1>(b),
                       geometry::get<1,0>(b));

#if 1
        // todo still using robust points
        using geometry::detail::equals::equals_point_point;
        bool const a_is_point = equals_point_point(robust_a1, robust_a2);
        bool const b_is_point = equals_point_point(robust_b1, robust_b2);

        if (a_is_point && b_is_point)
        {
            return equals_point_point(robust_a1, robust_b2)
                ? Policy::degenerate(a, true)
                : Policy::disjoint()
                ;
        }
#endif

        typedef typename geometry::point_type<Segment1>::type point_type; // TODO: most precise point?
        typedef typename segment_ratio_type
            <
                point_type, RobustPolicy
            >::type ratio_type;

        arithmetic::general_form<coordinate_type> gf_a = arithmetic::construct_line<coordinate_type>(a);
        arithmetic::general_form<coordinate_type> gf_b = arithmetic::construct_line<coordinate_type>(b);

        segment_intersection_info<point_type, ratio_type> sinfo;

        // Verify if lines (not segments) cross each other. If so, IP can be
        // on or of either of the segments.
        // The intersection_doubt flag is set for nearly collinear lines.
        point_type first_result; // Might go later
        bool intersection_doubt = false;
        bool const crossing = get_intersection(first_result, intersection_doubt, gf_a, gf_b);

        bool consider_as_collinear = false;
        boost::array<side_data, 4> data;

#ifdef BOOST_GEOMETRY_GENERIC_INT_SUPPORT_DOUBT1
        if (! crossing && intersection_doubt)
        {
            {
                double fa = -1, fb = -1;
                if (have_common_endpoints(a, b, fa, fb))
                {
                    std::cout << "#";
                    return handle_common_endpoints<Policy, ratio_type>(a, b, gf_a, gf_b, fa, fb, sinfo, a_is_point, b_is_point, robust_a1, robust_a2, robust_b1, robust_b2);
                }
            }

            if (disjoint_by_coordinates(a, b, gf_a))
            {
                return Policy::disjoint();
            }

            initialize_signed_comparable_distances(data, a, b, gf_a, gf_b);

            if (disjoint_by_side(data,
                                 10.0 * general_distance_threshold<coordinate_type>::get()))
            {
                return Policy::disjoint();
            }

            bool consider_as_crossing = false;
            bool consider_as_disjoint = false;

            inspect_sides(data, consider_as_collinear, consider_as_crossing,
                    consider_as_disjoint, a, b, gf_a, gf_b);
            if (consider_as_disjoint)
            {
                return Policy::disjoint();
            }
        }
#endif // BOOST_GEOMETRY_GENERIC_INT_SUPPORT_DOUBT1

        if (crossing)
        {
            sinfo.m_point = first_result;
            double fa = -1, fb = -1;
            bool consistent = true;

            if (get_fractions_on_segment(a, b, first_result, gf_a, gf_b, fa, fb))
            {
                side_info sides;

                initialize_signed_comparable_distances(data, a, b, gf_a, gf_b);
                initialize_on_end(data, fa, fb);

                int fa_i = -1, fb_i = -1;

                if (suspicious(fa_i, fb_i, fa, fb))
                {
                    fix_two_common_endpoints(data, fa, fb, a, b, fa_i, fb_i);
                    fix_one_common_endpoint(data, fa, fb, a, b, gf_a, gf_b, fa_i, fb_i);
                }

#ifdef BOOST_GEOMETRY_USE_COMPLEX_SEGMENT_RATIO
                assign_segment_intersection_info(sinfo, robust_a1, robust_a2,
                                                 robust_b1, robust_b2);
#else
                sinfo.robust_ra.assign(fa);
                sinfo.robust_rb.assign(fb);
#endif

                // Set sides points from A w.r.t. segment B (=q)
                consistent = side_assorter<true>::apply(sides, fa, fb,
                                                        consider_as_collinear,
                                                        data,
                                                        a, b,
                                           gf_a, gf_b);

                if (consistent)
                {
                    return Policy::segments_crosses(sides, sinfo, a, b);
                }
            }

            if (consistent)
            {
                // The two lines (general form) are intersecting,
                // but intersection point is outside segment(s)
                return Policy::disjoint();
            }
        }

        // The two lines (general form) are not intersecting. This means
        // that the two lines are parallel or collinear

        // It can also happen that there was an inconsistency in sides. Also
        // then lines are about collinear, and that is verified below.

        // TODO: only calculate if not yet considered as collinear
        arithmetic::general_form<double> gf_norm_a = arithmetic::normalize_line<double>(gf_a);
        arithmetic::general_form<double> gf_norm_b = arithmetic::normalize_line<double>(gf_b);

        if (consider_as_collinear
            || arithmetic::lines_collinear(gf_norm_a, gf_norm_b, equal_policy))
        {
            // Lines are collinear. There are 0, 1 or 2 intersections
            return arithmetic::more_horizontal(gf_a)
                    ? relate_collinear<0, Policy, ratio_type>(a, b, robust_a1, robust_a2, robust_b1, robust_b2, a_is_point, b_is_point)
                    : relate_collinear<1, Policy, ratio_type>(a, b, robust_a1, robust_a2, robust_b1, robust_b2, a_is_point, b_is_point);
         }

         // Lines are parallel but not collinear, so: disjoint
         return Policy::disjoint();
    }

    //--------------------------------------------------------------------------
    // todo: we can pick this up from intersection or put it more genericly.

    template
    <
        std::size_t Dimension,
        typename Policy,
        typename RatioType,
        typename Segment1,
        typename Segment2,
        typename RobustPoint1,
        typename RobustPoint2
    >
    static inline typename Policy::return_type
        relate_collinear(Segment1 const& a,
                         Segment2 const& b,
                         RobustPoint1 const& robust_a1, RobustPoint1 const& robust_a2,
                         RobustPoint2 const& robust_b1, RobustPoint2 const& robust_b2,
                         bool a_is_point, bool b_is_point)
    {
#ifdef BOOST_GEOMETRY_USE_COMPLEX_SEGMENT_RATIO
        if (a_is_point)
        {
            return relate_one_degenerate<Policy, RatioType>(a,
                get<Dimension>(robust_a1),
                get<Dimension>(robust_b1), get<Dimension>(robust_b2),
                true);
        }
        if (b_is_point)
        {
            return relate_one_degenerate<Policy, RatioType>(b,
                get<Dimension>(robust_b1),
                get<Dimension>(robust_a1), get<Dimension>(robust_a2),
                false);
        }
#else
        if (a_is_point)
        {
            return relate_one_degenerate<Policy, RatioType>(a,
                get<0, Dimension>(a),
                get<0, Dimension>(b), get<1, Dimension>(b),
                true);
        }
        if (b_is_point)
        {
            return relate_one_degenerate<Policy, RatioType>(b,
                get<0, Dimension>(b),
                get<0, Dimension>(a), get<1, Dimension>(a),
                false);
        }
#endif
        return relate_collinear<Policy, RatioType>(a, b,
                                get<Dimension>(robust_a1),
                                get<Dimension>(robust_a2),
                                get<Dimension>(robust_b1),
                                get<Dimension>(robust_b2));
    }

    /// Relate segments known collinear
    template
    <
        typename Policy,
        typename RatioType,
        typename Segment1,
        typename Segment2,
        typename RobustType1,
        typename RobustType2
    >
    static inline typename Policy::return_type
        relate_collinear(Segment1 const& a, Segment2 const& b,
                         RobustType1 oa_1, RobustType1 oa_2,
                         RobustType2 ob_1, RobustType2 ob_2)
    {
        // Calculate the ratios where a starts in b, b starts in a
        //         a1--------->a2         (2..7)
        //                b1----->b2      (5..8)
        // length_a: 7-2=5
        // length_b: 8-5=3
        // b1 is located w.r.t. a at ratio: (5-2)/5=3/5 (on a)
        // b2 is located w.r.t. a at ratio: (8-2)/5=6/5 (right of a)
        // a1 is located w.r.t. b at ratio: (2-5)/3=-3/3 (left of b)
        // a2 is located w.r.t. b at ratio: (7-5)/3=2/3 (on b)
        // A arrives (a2 on b), B departs (b1 on a)

        // If both are reversed:
        //         a2<---------a1         (7..2)
        //                b2<-----b1      (8..5)
        // length_a: 2-7=-5
        // length_b: 5-8=-3
        // b1 is located w.r.t. a at ratio: (8-7)/-5=-1/5 (before a starts)
        // b2 is located w.r.t. a at ratio: (5-7)/-5=2/5 (on a)
        // a1 is located w.r.t. b at ratio: (7-8)/-3=1/3 (on b)
        // a2 is located w.r.t. b at ratio: (2-8)/-3=6/3 (after b ends)

        // If both one is reversed:
        //         a1--------->a2         (2..7)
        //                b2<-----b1      (8..5)
        // length_a: 7-2=+5
        // length_b: 5-8=-3
        // b1 is located w.r.t. a at ratio: (8-2)/5=6/5 (after a ends)
        // b2 is located w.r.t. a at ratio: (5-2)/5=3/5 (on a)
        // a1 is located w.r.t. b at ratio: (2-8)/-3=6/3 (after b ends)
        // a2 is located w.r.t. b at ratio: (7-8)/-3=1/3 (on b)
        RobustType1 const length_a = oa_2 - oa_1; // no abs, see above
        RobustType2 const length_b = ob_2 - ob_1;

#ifdef BOOST_GEOMETRY_USE_COMPLEX_SEGMENT_RATIO
        RatioType ra_from(oa_1 - ob_1, length_b);
        RatioType ra_to(oa_2 - ob_1, length_b);
        RatioType rb_from(ob_1 - oa_1, length_a);
        RatioType rb_to(ob_2 - oa_1, length_a);
#else
        RatioType ra_from(oa_1 - ob_1, length_b);
        RatioType ra_to(oa_2 - ob_1, length_b);
        RatioType rb_from(ob_1 - oa_1, length_a);
        RatioType rb_to(ob_2 - oa_1, length_a);
#endif

        // use absolute measure to detect endpoints intersection
        // NOTE: it'd be possible to calculate bx_wrt_a using ax_wrt_b values
        int const a1_wrt_b = position_value(oa_1, ob_1, ob_2);
        int const a2_wrt_b = position_value(oa_2, ob_1, ob_2);
        int const b1_wrt_a = position_value(ob_1, oa_1, oa_2);
        int const b2_wrt_a = position_value(ob_2, oa_1, oa_2);
        
        // fix the ratios if necessary
        // CONSIDER: fixing ratios also in other cases, if they're inconsistent
        // e.g. if ratio == 1 or 0 (so IP at the endpoint)
        // but position value indicates that the IP is in the middle of the segment
        // because one of the segments is very long
        // In such case the ratios could be moved into the middle direction
        // by some small value (e.g. EPS+1ULP)
        if (a1_wrt_b == 1)
        {
            ra_from.assign(0, 1);
            rb_from.assign(0, 1);
        }
        else if (a1_wrt_b == 3)
        {
            ra_from.assign(1, 1);
            rb_to.assign(0, 1);
        } 

        if (a2_wrt_b == 1)
        {
            ra_to.assign(0, 1);
            rb_from.assign(1, 1);
        }
        else if (a2_wrt_b == 3)
        {
            ra_to.assign(1, 1);
            rb_to.assign(1, 1);
        }

        if ((a1_wrt_b < 1 && a2_wrt_b < 1) || (a1_wrt_b > 3 && a2_wrt_b > 3))
        //if ((ra_from.left() && ra_to.left()) || (ra_from.right() && ra_to.right()))
        {
            return Policy::disjoint();
        }

        bool const opposite = math::sign(length_a) != math::sign(length_b);

        return Policy::segments_collinear(a, b, opposite,
                                          a1_wrt_b, a2_wrt_b, b1_wrt_a, b2_wrt_a,
                                          ra_from, ra_to, rb_from, rb_to);
    }

#ifdef BOOST_GEOMETRY_USE_COMPLEX_SEGMENT_RATIO
    /// Relate segments where one is degenerate
    template
    <
        typename Policy,
        typename RatioType,
        typename DegenerateSegment,
        typename RobustType1,
        typename RobustType2
    >
    static inline typename Policy::return_type
        relate_one_degenerate(DegenerateSegment const& degenerate_segment,
                              RobustType1 d, RobustType2 s1, RobustType2 s2,
                              bool a_degenerate)
    {
        // Calculate the ratios where ds starts in s
        //         a1--------->a2         (2..6)
        //              b1/b2      (4..4)
        // Ratio: (4-2)/(6-2)
        RatioType const ratio(d - s1, s2 - s1);

        if (!ratio.on_segment())
        {
            return Policy::disjoint();
        }

        return Policy::one_degenerate(degenerate_segment, ratio, a_degenerate);
    }
#else
    /// Relate segments where one is degenerate
    template
    <
        typename Policy,
        typename RatioType,
        typename DegenerateSegment,
        typename PointCoordinateType,
        typename SegmentCoordinateType
    >
    static inline typename Policy::return_type
        relate_one_degenerate(DegenerateSegment const& degenerate_segment,
                              PointCoordinateType const& d,
                              SegmentCoordinateType const& s1,
                              SegmentCoordinateType const& s2,
                              bool a_degenerate)
    {
        RatioType const ratio(d - s1 / (s2 - s1));

        if (!ratio.on_segment())
        {
            return Policy::disjoint();
        }

        return Policy::one_degenerate(degenerate_segment, ratio, a_degenerate);
    }
#endif

    template <typename ProjCoord1, typename ProjCoord2>
    static inline int position_value(ProjCoord1 const& ca1,
                                     ProjCoord2 const& cb1,
                                     ProjCoord2 const& cb2)
    {
        // S1x  0   1    2     3   4
        // S2       |---------->
        return math::equals(ca1, cb1) ? 1
             : math::equals(ca1, cb2) ? 3
             : cb1 < cb2 ?
                ( ca1 < cb1 ? 0
                : ca1 > cb2 ? 4
                : 2 )
              : ( ca1 > cb1 ? 0
                : ca1 < cb2 ? 4
                : 2 );
    }
};

#if ! defined(BOOST_GEOMETRY_USE_KRAMER_RULE)
#ifndef DOXYGEN_NO_STRATEGY_SPECIALIZATIONS
namespace services
{

template <typename CalculationType>
struct default_strategy<cartesian_tag, CalculationType>
{
    typedef cartesian_general_segments<CalculationType> type;
};

} // namespace services
#endif // DOXYGEN_NO_STRATEGY_SPECIALIZATIONS
#endif

}} // namespace strategy::intersection

#if ! defined(BOOST_GEOMETRY_USE_KRAMER_RULE)
namespace strategy
{

namespace within { namespace services
{

template <typename Geometry1, typename Geometry2, typename AnyTag1, typename AnyTag2>
struct default_strategy<Geometry1, Geometry2, AnyTag1, AnyTag2, linear_tag, linear_tag, cartesian_tag, cartesian_tag>
{
    typedef strategy::intersection::cartesian_general_segments<> type;
};

template <typename Geometry1, typename Geometry2, typename AnyTag1, typename AnyTag2>
struct default_strategy<Geometry1, Geometry2, AnyTag1, AnyTag2, linear_tag, polygonal_tag, cartesian_tag, cartesian_tag>
{
    typedef strategy::intersection::cartesian_general_segments<> type;
};

template <typename Geometry1, typename Geometry2, typename AnyTag1, typename AnyTag2>
struct default_strategy<Geometry1, Geometry2, AnyTag1, AnyTag2, polygonal_tag, linear_tag, cartesian_tag, cartesian_tag>
{
    typedef strategy::intersection::cartesian_general_segments<> type;
};

template <typename Geometry1, typename Geometry2, typename AnyTag1, typename AnyTag2>
struct default_strategy<Geometry1, Geometry2, AnyTag1, AnyTag2, polygonal_tag, polygonal_tag, cartesian_tag, cartesian_tag>
{
    typedef strategy::intersection::cartesian_general_segments<> type;
};

}} // within::services

namespace covered_by { namespace services
{

template <typename Geometry1, typename Geometry2, typename AnyTag1, typename AnyTag2>
struct default_strategy<Geometry1, Geometry2, AnyTag1, AnyTag2, linear_tag, linear_tag, cartesian_tag, cartesian_tag>
{
    typedef strategy::intersection::cartesian_general_segments<> type;
};

template <typename Geometry1, typename Geometry2, typename AnyTag1, typename AnyTag2>
struct default_strategy<Geometry1, Geometry2, AnyTag1, AnyTag2, linear_tag, polygonal_tag, cartesian_tag, cartesian_tag>
{
    typedef strategy::intersection::cartesian_general_segments<> type;
};

template <typename Geometry1, typename Geometry2, typename AnyTag1, typename AnyTag2>
struct default_strategy<Geometry1, Geometry2, AnyTag1, AnyTag2, polygonal_tag, linear_tag, cartesian_tag, cartesian_tag>
{
    typedef strategy::intersection::cartesian_general_segments<> type;
};

template <typename Geometry1, typename Geometry2, typename AnyTag1, typename AnyTag2>
struct default_strategy<Geometry1, Geometry2, AnyTag1, AnyTag2, polygonal_tag, polygonal_tag, cartesian_tag, cartesian_tag>
{
    typedef strategy::intersection::cartesian_general_segments<> type;
};

}} // within::services

} // strategy
#endif


}} // namespace boost::geometry


#endif // BOOST_GEOMETRY_STRATEGIES_CARTESIAN_GENERAL_INTERSECTION_HPP
