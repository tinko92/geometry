// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2020 Tinko Bartels, Berlin, Germany.

// Contributed and/or modified by Tinko Bartels,
//   as part of Google Summer of Code 2020 program.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_SIDE_FILTERED_HPP
#define BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_SIDE_FILTERED_HPP

#include <boost/geometry/util/select_coordinate_type.hpp>
#include <boost/geometry/core/access.hpp>

#include <boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/stage_a.hpp>
#include <boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/stage_d.hpp>

namespace boost { namespace geometry
{

namespace strategy { namespace side
{

/*!
\brief Check at which side of a segment a point lies:
    left of segment (> 0), right of segment (< 0), on segment (0)
\ingroup strategies
\tparam CalculationType \tparam_calculation
 */
template <typename CalculationType = void>
class side_filtered
{

public:
    typedef cartesian_tag cs_tag;
    template
    <
        typename CoordinateType,
        typename PromotedType,
        typename P1,
        typename P2,
        typename P
    >
    static inline
    PromotedType side_value(P1 const& p1, P2 const& p2, P const& p)
    {
        using boost::geometry::detail::generic_robust_predicates::difference;
        using boost::geometry::detail::generic_robust_predicates::product;
        using boost::geometry::detail::generic_robust_predicates::leaf;
        using expression = 
                difference<
                        product<
                                difference<
                                        leaf<1>,
                                        leaf<5>
                                >,
                                difference<
                                        leaf<4>,
                                        leaf<6>
                                >
                        >,
                        product<
                                difference<
                                        leaf<2>,
                                        leaf<6>
                                >,
                                difference<
                                        leaf<3>,
                                        leaf<5>
                                >
                        >
                >;
        return boost::geometry::detail::generic_robust_predicates::
            approximate_value<expression, PromotedType>(
                get<0>(p1), get<1>(p1),
                get<0>(p2), get<1>(p2),
                get<0>(p), get<1>(p));
    }

    template <typename P1, typename P2, typename P>
    static inline int apply(P1 const& p1, P2 const& p2, P const& p)
    {
        typedef typename coordinate_type<P1>::type coordinate_type1;
        typedef typename coordinate_type<P2>::type coordinate_type2;
        typedef typename coordinate_type<P>::type coordinate_type3;    
            
        typedef typename boost::mpl::if_c 
            <
                boost::is_void<CalculationType>::type::value,
                typename select_most_precise
                    <
                        typename select_most_precise
                            <
                                coordinate_type1, coordinate_type2
                            >::type,
                        coordinate_type3
                    >::type,
                CalculationType
            >::type coordinate_type;
            
        typedef typename select_most_precise
            <
                coordinate_type,
                double
            >::type promoted_type; 
        using boost::geometry::detail::generic_robust_predicates::difference;
        using boost::geometry::detail::generic_robust_predicates::product;
        using boost::geometry::detail::generic_robust_predicates::leaf;
        using expression = 
                difference<
                        product<
                                difference<
                                        leaf<1>,
                                        leaf<5>
                                >,
                                difference<
                                        leaf<4>,
                                        leaf<6>
                                >
                        >,
                        product<
                                difference<
                                        leaf<2>,
                                        leaf<6>
                                >,
                                difference<
                                        leaf<3>,
                                        leaf<5>
                                >
                        >
                >;
        auto stage_a_result = boost::geometry::detail::generic_robust_predicates::
            stage_a<expression, promoted_type>(
                get<0>(p1), get<1>(p1),
                get<0>(p2), get<1>(p2),
                get<0>(p), get<1>(p));
        if(stage_a_result !=
                boost::geometry::detail::generic_robust_predicates::sign_uncertain)
            return stage_a_result;
        else return boost::geometry::detail::generic_robust_predicates::
            stage_d<expression, promoted_type>(
                get<0>(p1), get<1>(p1),
                get<0>(p2), get<1>(p2),
                get<0>(p), get<1>(p));
    }

};

}} // namespace strategy::side

}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_SIDE_FILTERED_HPP
