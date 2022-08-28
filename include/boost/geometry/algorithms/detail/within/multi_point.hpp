// Boost.Geometry

// Copyright (c) 2017-2020, Oracle and/or its affiliates.
// Contributed and/or modified by Adam Wulkiewicz, on behalf of Oracle

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_ALGORITHMS_DETAIL_WITHIN_MULTI_POINT_HPP
#define BOOST_GEOMETRY_ALGORITHMS_DETAIL_WITHIN_MULTI_POINT_HPP


#include <algorithm>
#include <vector>

#include <boost/range/iterator_range_core.hpp>
#include <boost/range/begin.hpp>
#include <boost/range/end.hpp>
#include <boost/range/size.hpp>
#include <boost/range/value_type.hpp>

#include <boost/geometry/algorithms/detail/disjoint/box_box.hpp>
#include <boost/geometry/algorithms/detail/disjoint/point_box.hpp>
#include <boost/geometry/algorithms/detail/expand_by_epsilon.hpp>
#include <boost/geometry/algorithms/detail/within/point_in_geometry.hpp>
#include <boost/geometry/algorithms/envelope.hpp>
#include <boost/geometry/algorithms/detail/partition.hpp>
#include <boost/geometry/core/tag.hpp>
#include <boost/geometry/core/tag_cast.hpp>
#include <boost/geometry/core/tags.hpp>

#include <boost/geometry/geometries/box.hpp>

#include <boost/geometry/index/rtree.hpp>

#include <boost/geometry/policies/compare.hpp>

#include <boost/geometry/strategies/covered_by.hpp>
#include <boost/geometry/strategies/disjoint.hpp>

#include <boost/geometry/util/type_traits.hpp>


namespace boost { namespace geometry {

#ifndef DOXYGEN_NO_DETAIL
namespace detail { namespace within {

struct multi_point_point
{
    template <typename MultiPoint, typename Point, typename Strategy>
    static inline bool apply(MultiPoint const& multi_point,
                             Point const& point,
                             Strategy const& strategy)
    {
        auto const s = strategy.relate(multi_point, point);
        for ( auto const& p : boost::make_iterator_range(multi_point) )
        {
            if (! s.apply(p, point))
            {
                return false;
            }
        }

        // all points of MultiPoint inside Point
        return true;
    }
};

// NOTE: currently the strategy is ignored, math::equals() is used inside geometry::less<>
struct multi_point_multi_point
{
    template <typename MultiPoint1, typename MultiPoint2, typename Strategy>
    static inline bool apply(MultiPoint1 const& multi_point1,
                             MultiPoint2 const& multi_point2,
                             Strategy const& /*strategy*/)
    {
        using point2_type = typename boost::range_value<MultiPoint2>::type;
        using cs_tag = typename Strategy::cs_tag;
        using less_type = geometry::less<void, -1, cs_tag>;

        less_type const less = less_type();

        std::vector<point2_type> points2(boost::begin(multi_point2), boost::end(multi_point2));
        std::sort(points2.begin(), points2.end(), less);

        bool result = false;

        for ( auto const& point1 : boost::make_iterator_range(multi_point1))
        {
            if (! std::binary_search(points2.begin(), points2.end(), point1, less))
            {
                return false;
            }
            else
            {
                result = true;
            }
        }

        return result;
    }
};


// TODO: the complexity could be lesser
//   the second geometry could be "prepared"/sorted
// For Linear geometries partition could be used
// For Areal geometries point_in_geometry() would have to call the winding
//   strategy differently, currently it linearly calls the strategy for each
//   segment. So the segments would have to be sorted in a way consistent with
//   the strategy and then the strategy called only for the segments in range.
template <bool Within>
struct multi_point_single_geometry
{
    template <typename MultiPoint, typename LinearOrAreal, typename Strategy>
    static inline bool apply(MultiPoint const& multi_point,
                             LinearOrAreal const& linear_or_areal,
                             Strategy const& strategy)
    {
        //typedef typename boost::range_value<MultiPoint>::type point1_type;
        using point2_type = typename point_type<LinearOrAreal>::type;
        using box2_type = model::box<point2_type>;

        // Create envelope of geometry
        box2_type box;
        geometry::envelope(linear_or_areal, box, strategy);
        geometry::detail::expand_by_epsilon(box);

        // Test each Point with envelope and then geometry if needed
        // If in the exterior, break
        bool result = false;

        for ( auto const& point : boost::make_iterator_range(multi_point) )
        {
            using point_in_box_type = decltype(strategy.covered_by(point, box));

            int in_val = 0;
            
            // exterior of box and of geometry
            if (! point_in_box_type::apply(point, box)
                || (in_val = point_in_geometry(point, linear_or_areal, strategy)) < 0)
            {
                result = false;
                break;
            }

            // interior : interior/boundary
            if (Within ? in_val > 0 : in_val >= 0)
            {
                result = true;
            }
        }

        return result;
    }
};


// TODO: same here, probably the complexity could be lesser
template <bool Within>
struct multi_point_multi_geometry
{
    template <typename MultiPoint, typename LinearOrAreal, typename Strategy>
    static inline bool apply(MultiPoint const& multi_point,
                             LinearOrAreal const& linear_or_areal,
                             Strategy const& strategy)
    {
        using point2_type = typename point_type<LinearOrAreal>::type;
        using box2_type = model::box<point2_type>;
        constexpr bool is_linear = util::is_linear<LinearOrAreal>::value;

        // TODO: box pairs could be constructed on the fly, inside the rtree

        // Prepare range of envelopes and ids
        std::size_t count2 = boost::size(linear_or_areal);
        using box_pair_type = std::pair<box2_type, std::size_t>;
        using box_pair_vector = std::vector<box_pair_type>;
        box_pair_vector boxes(count2);
        for (std::size_t i = 0 ; i < count2 ; ++i)
        {
            geometry::envelope(linear_or_areal, boxes[i].first, strategy);
            geometry::detail::expand_by_epsilon(boxes[i].first);
            boxes[i].second = i;
        }

        // Create R-tree
        using index_parameters_type = index::parameters<index::rstar<4>, Strategy>;
        index::rtree<box_pair_type, index_parameters_type>
            rtree(boxes.begin(), boxes.end(),
                  index_parameters_type(index::rstar<4>(), strategy));

        // For each point find overlapping envelopes and test corresponding single geometries
        // If a point is in the exterior break
        bool result = false;

        for ( auto const& point : boost::make_iterator_range(multi_point) )
        {
            // TODO: investigate the possibility of using satisfies
            // TODO: investigate the possibility of using iterative queries (optimization below)
            box_pair_vector inters_boxes;
            rtree.query(index::intersects(point), std::back_inserter(inters_boxes));

            bool found_interior = false;
            bool found_boundary = false;
            int boundaries = 0;

            for (auto const& box : inters_boxes)
            {
                int const in_val = point_in_geometry(point,
                    range::at(linear_or_areal, box.second), strategy);

                if (in_val > 0)
                {
                    found_interior = true;
                }
                else if (in_val == 0)
                {
                    ++boundaries;
                }

                // If the result was set previously (interior or
                // interior/boundary found) the only thing that needs to be
                // done for other points is to make sure they're not
                // overlapping the exterior no need to analyse boundaries.
                if (result && in_val >= 0)
                {
                    break;
                }
            }

            if (boundaries > 0)
            {
                if (is_linear && boundaries % 2 == 0)
                {
                    found_interior = true;
                }
                else
                {
                    found_boundary = true;
                }
            }

            // exterior
            if (! found_interior && ! found_boundary)
            {
                result = false;
                break;
            }

            // interior : interior/boundary
            if (Within ? found_interior : (found_interior || found_boundary))
            {
                result = true;
            }
        }

        return result;
    }
};

}} // namespace detail::within
#endif // DOXYGEN_NO_DETAIL

}} // namespace boost::geometry


#endif // BOOST_GEOMETRY_ALGORITHMS_DETAIL_WITHIN_MULTI_POINT_HPP
