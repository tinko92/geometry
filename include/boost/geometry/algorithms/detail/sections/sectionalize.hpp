// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2007-2015 Barend Gehrels, Amsterdam, the Netherlands.
// Copyright (c) 2008-2015 Bruno Lalande, Paris, France.
// Copyright (c) 2009-2015 Mateusz Loskot, London, UK.
// Copyright (c) 2014-2015 Adam Wulkiewicz, Lodz, Poland.

// This file was modified by Oracle on 2013-2021.
// Modifications copyright (c) 2013-2021 Oracle and/or its affiliates.

// Contributed and/or modified by Adam Wulkiewicz, on behalf of Oracle
// Contributed and/or modified by Menelaos Karavelas, on behalf of Oracle

// Parts of Boost.Geometry are redesigned from Geodan's Geographic Library
// (geolib/GGL), copyright (c) 1995-2010 Geodan, Amsterdam, the Netherlands.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_ALGORITHMS_DETAIL_SECTIONS_SECTIONALIZE_HPP
#define BOOST_GEOMETRY_ALGORITHMS_DETAIL_SECTIONS_SECTIONALIZE_HPP

#include <cstddef>
#include <type_traits>
#include <vector>
#include <array>

#include <boost/concept/requires.hpp>
#include <boost/core/ignore_unused.hpp>
#include <boost/range/begin.hpp>
#include <boost/range/end.hpp>
#include <boost/range/size.hpp>
#include <boost/range/value_type.hpp>

#include <boost/geometry/core/config.hpp>

#include <boost/geometry/algorithms/assign.hpp>
#include <boost/geometry/algorithms/envelope.hpp>
#include <boost/geometry/algorithms/detail/interior_iterator.hpp>
#include <boost/geometry/algorithms/detail/recalculate.hpp>
#include <boost/geometry/algorithms/detail/ring_identifier.hpp>
#include <boost/geometry/algorithms/detail/signed_size_type.hpp>

#include <boost/geometry/core/closure.hpp>
#include <boost/geometry/core/exterior_ring.hpp>
#include <boost/geometry/core/point_order.hpp>
#include <boost/geometry/core/static_assert.hpp>
#include <boost/geometry/core/tags.hpp>

#include <boost/geometry/geometries/concepts/check.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/segment.hpp>
#include <boost/geometry/policies/robustness/no_rescale_policy.hpp>
#include <boost/geometry/policies/robustness/robust_point_type.hpp>
#include <boost/geometry/util/sequence.hpp>
#include <boost/geometry/views/detail/closed_clockwise_view.hpp>

namespace boost { namespace geometry
{


/*!
    \brief Structure containing section information
    \details Section information consists of a bounding box, direction
        information (if it is increasing or decreasing, per dimension),
        index information (begin-end, ring, multi) and the number of
        segments in this section

    \tparam Box box-type
    \tparam DimensionCount number of dimensions for this section
    \ingroup sectionalize
 */
template
<
    typename Box,
    std::size_t DimensionCount
>
struct section
{
    using box_type = Box;
    static std::size_t const dimension_count = DimensionCount;

    std::array<int, DimensionCount> directions;
    ring_identifier ring_id;
    Box bounding_box;

    // NOTE: size_type could be passed as template parameter
    // NOTE: these probably also could be of type std::size_t
    signed_size_type begin_index;
    signed_size_type end_index;
    std::size_t count;
    std::size_t range_count;
    bool duplicate;
    signed_size_type non_duplicate_index;

    bool is_non_duplicate_first;
    bool is_non_duplicate_last;

    inline section()
        : begin_index(-1)
        , end_index(-1)
        , count(0)
        , range_count(0)
        , duplicate(false)
        , non_duplicate_index(-1)
        , is_non_duplicate_first(false)
        , is_non_duplicate_last(false)
    {
        assign_inverse(bounding_box);
        directions.fill(0);
    }
};


/*!
    \brief Structure containing a collection of sections
    \note Derived from a vector, proves to be faster than of deque
    \note vector might be templated in the future
    \ingroup sectionalize
 */
template <typename Box, std::size_t DimensionCount>
struct sections : std::vector<section<Box, DimensionCount> >
{
    using box_type = Box;
    static std::size_t const value = DimensionCount;
};


#ifndef DOXYGEN_NO_DETAIL
namespace detail { namespace sectionalize
{

/// @brief Helper class to create sections of a part of a range, on the fly
template<typename DimensionVector>
struct sectionalize_part
{
    static const std::size_t dimension_count
        = util::sequence_size<DimensionVector>::value;

    template
    <
        typename Iterator,
        typename RobustPolicy,
        typename Sections,
        typename Strategy
    >
    static inline void apply(Sections& sections,
                             Iterator begin, Iterator end,
                             RobustPolicy const& robust_policy,
                             Strategy const& strategy,
                             ring_identifier ring_id,
                             std::size_t max_count)
    {
        boost::ignore_unused(robust_policy);

        using section_type = typename boost::range_value<Sections>::type;
        using box_type = typename section_type::box_type;
        using point_type = typename geometry::point_type<box_type>::type;

        BOOST_STATIC_ASSERT
            (
                section_type::dimension_count
                 == util::sequence_size<DimensionVector>::value
            );

        using robust_point_type = typename geometry::robust_point_type
            <
                point_type,
                RobustPolicy
            >::type;

        std::size_t const count = std::distance(begin, end);
        if (count == 0)
        {
            return;
        }

        signed_size_type index = 0;
        signed_size_type ndi = 0; // non duplicate index
        section_type section;

        bool mark_first_non_duplicated = true;
        std::size_t last_non_duplicate_index = sections.size();

        Iterator it = begin;
        robust_point_type previous_robust_point;
        geometry::recalculate(previous_robust_point, *it, robust_policy);

        for (Iterator previous = it++;
             it != end;
             ++previous, ++it, index++)
        {
            robust_point_type current_robust_point;
            geometry::recalculate(current_robust_point, *it, robust_policy);
            model::referring_segment<robust_point_type> robust_segment(
                    previous_robust_point, current_robust_point);

            auto direction_classes = strategy.directions().template apply<DimensionVector>(robust_segment);

            // if "dir" == 0 for all point-dimensions, it is duplicate.
            // Those sections might be omitted, if wished, lateron
            bool duplicate = false;

            if (direction_classes[0] == 0)
            {
                // Recheck because ALL dimensions should be checked,
                // not only first one.
                // (dimension_count might be < dimension<P>::value)
                if (strategy.degenerate_segment().apply(robust_segment))
                {
                    duplicate = true;

                    // Change direction-info to force new section
                    // Note that wo consecutive duplicate segments will generate
                    // only one duplicate-section.
                    // Actual value is not important as long as it is not -1,0,1
                    direction_classes.fill(-99);
                }
            }

            if (section.count > 0
                && ( direction_classes != section.directions
                    || section.count > max_count)
                )
            {
                if (! section.duplicate)
                {
                    last_non_duplicate_index = sections.size();
                }

                sections.push_back(section);
                section = section_type();
            }

            if (section.count == 0)
            {
                section.begin_index = index;
                section.ring_id = ring_id;
                section.duplicate = duplicate;
                section.non_duplicate_index = ndi;
                section.range_count = count;

                if (mark_first_non_duplicated && ! duplicate)
                {
                    section.is_non_duplicate_first = true;
                    mark_first_non_duplicated = false;
                }

                section.directions = direction_classes;

                // In cartesian this is envelope of previous point expanded with current point
                // in non-cartesian this is envelope of a segment

                geometry::model::referring_segment<robust_point_type const> seg(
                                                                            previous_robust_point,
                                                                            current_robust_point);
                geometry::envelope(seg, section.bounding_box, strategy);
            }
            else
            {
                // In cartesian this is expand with current point
                // in non-cartesian this is expand with a segment
                geometry::model::referring_segment<robust_point_type const> seg(
                                                                            previous_robust_point,
                                                                            current_robust_point);
                strategy.expand(section.bounding_box, seg).template apply<true>(section.bounding_box, seg);
            }

            section.end_index = index + 1;
            section.count++;
            if (! duplicate)
            {
                ndi++;
            }
            previous_robust_point = current_robust_point;
        }

        // Add last section if applicable
        if (section.count > 0)
        {
            if (! section.duplicate)
            {
                last_non_duplicate_index = sections.size();
            }

            sections.push_back(section);
        }

        if (last_non_duplicate_index < sections.size()
            && ! sections[last_non_duplicate_index].duplicate)
        {
            sections[last_non_duplicate_index].is_non_duplicate_last = true;
        }
    }
};


template
<
    closure_selector Closure,
    bool Reverse,
    typename DimensionVector
>
struct sectionalize_range
{
    template
    <
        typename Range,
        typename RobustPolicy,
        typename Sections,
        typename Strategy
    >
    static inline void apply(Range const& range,
                             RobustPolicy const& robust_policy,
                             Sections& sections,
                             Strategy const& strategy,
                             ring_identifier ring_id,
                             std::size_t max_count)
    {
        detail::closed_clockwise_view
            <
                Range const,
                Closure,
                Reverse ? counterclockwise : clockwise
            > const view(range);

        std::size_t const n = boost::size(view);
        if (n == 0)
        {
            // Zero points, no section
            return;
        }

        if (n == 1)
        {
            // Line with one point ==> no sections
            return;
        }

        sectionalize_part<DimensionVector>::apply(sections,
            boost::begin(view), boost::end(view),
            robust_policy, strategy,
            ring_id, max_count);
    }
};

template
<
    bool Reverse,
    typename DimensionVector
>
struct sectionalize_polygon
{
    template
    <
        typename Polygon,
        typename RobustPolicy,
        typename Sections,
        typename Strategy
    >
    static inline void apply(Polygon const& poly,
                RobustPolicy const& robust_policy,
                Sections& sections,
                Strategy const& strategy,
                ring_identifier ring_id,
                std::size_t max_count)
    {
        using sectionalizer = sectionalize_range
            <
                closure<Polygon>::value, Reverse, DimensionVector
            >;

        ring_id.ring_index = -1;
        sectionalizer::apply(exterior_ring(poly), robust_policy, sections,
                         strategy, ring_id, max_count);

        ring_id.ring_index++;
        auto const& rings = interior_rings(poly);
        for (auto it = boost::begin(rings); it != boost::end(rings);
             ++it, ++ring_id.ring_index)
        {
            sectionalizer::apply(*it, robust_policy, sections,
                             strategy, ring_id, max_count);
        }
    }
};

template <typename DimensionVector>
struct sectionalize_box
{
    template
    <
        typename Box,
        typename RobustPolicy,
        typename Sections,
        typename Strategy
    >
    static inline void apply(Box const& box,
                RobustPolicy const& robust_policy,
                Sections& sections,
                Strategy const& strategy,
                ring_identifier const& ring_id, std::size_t max_count)
    {
        using point_type = typename point_type<Box>::type;

        assert_dimension<Box, 2>();

        // Add all four sides of the 2D-box as separate section.
        // Easiest is to convert it to a polygon.
        // However, we don't have the polygon type
        // (or polygon would be a helper-type).
        // Therefore we mimic a linestring/std::vector of 5 points

        // TODO: might be replaced by assign_box_corners_oriented
        // or just "convert"
        point_type ll, lr, ul, ur;
        geometry::detail::assign_box_corners(box, ll, lr, ul, ur);

        std::vector<point_type> points;
        points.push_back(ll);
        points.push_back(ul);
        points.push_back(ur);
        points.push_back(lr);
        points.push_back(ll);

        // NOTE: Use cartesian envelope strategy in all coordinate systems
        //       because edges of a box are not geodesic segments
        sectionalize_range
            <
                closed, false, DimensionVector
            >::apply(points, robust_policy, sections,
                     strategy, ring_id, max_count);
    }
};

template <typename DimensionVector, typename Policy>
struct sectionalize_multi
{
    template
    <
        typename MultiGeometry,
        typename RobustPolicy,
        typename Sections,
        typename Strategy
    >
    static inline void apply(MultiGeometry const& multi,
                RobustPolicy const& robust_policy,
                Sections& sections,
                Strategy const& strategy,
                ring_identifier ring_id,
                std::size_t max_count)
    {
        ring_id.multi_index = 0;
        for (auto it = boost::begin(multi);
            it != boost::end(multi);
            ++it, ++ring_id.multi_index)
        {
            Policy::apply(*it, robust_policy, sections,
                          strategy,
                          ring_id, max_count);
        }
    }
};

}} // namespace detail::sectionalize
#endif // DOXYGEN_NO_DETAIL


#ifndef DOXYGEN_NO_DISPATCH
namespace dispatch
{

template
<
    typename Tag,
    typename Geometry,
    bool Reverse,
    typename DimensionVector
>
struct sectionalize
{
    BOOST_GEOMETRY_STATIC_ASSERT_FALSE(
        "Not or not yet implemented for this Geometry type.",
        Tag, Geometry);
};

template
<
    typename Box,
    bool Reverse,
    typename DimensionVector
>
struct sectionalize<box_tag, Box, Reverse, DimensionVector>
    : detail::sectionalize::sectionalize_box<DimensionVector>
{};

template
<
    typename LineString,
    typename DimensionVector
>
struct sectionalize
    <
        linestring_tag,
        LineString,
        false,
        DimensionVector
    >
    : detail::sectionalize::sectionalize_range<closed, false, DimensionVector>
{};

template
<
    typename Ring,
    bool Reverse,
    typename DimensionVector
>
struct sectionalize<ring_tag, Ring, Reverse, DimensionVector>
    : detail::sectionalize::sectionalize_range
        <
            geometry::closure<Ring>::value, Reverse,
            DimensionVector
        >
{};

template
<
    typename Polygon,
    bool Reverse,
    typename DimensionVector
>
struct sectionalize<polygon_tag, Polygon, Reverse, DimensionVector>
    : detail::sectionalize::sectionalize_polygon
        <
            Reverse, DimensionVector
        >
{};

template
<
    typename MultiPolygon,
    bool Reverse,
    typename DimensionVector
>
struct sectionalize<multi_polygon_tag, MultiPolygon, Reverse, DimensionVector>
    : detail::sectionalize::sectionalize_multi
        <
            DimensionVector,
            detail::sectionalize::sectionalize_polygon
                <
                    Reverse,
                    DimensionVector
                >
        >

{};

template
<
    typename MultiLinestring,
    bool Reverse,
    typename DimensionVector
>
struct sectionalize<multi_linestring_tag, MultiLinestring, Reverse, DimensionVector>
    : detail::sectionalize::sectionalize_multi
        <
            DimensionVector,
            detail::sectionalize::sectionalize_range<closed, false, DimensionVector>
        >

{};

} // namespace dispatch
#endif


/*!
    \brief Split a geometry into monotonic sections
    \ingroup sectionalize
    \tparam Geometry type of geometry to check
    \tparam Sections type of sections to create
    \param geometry geometry to create sections from
    \param robust_policy policy to handle robustness issues
    \param sections structure with sections
    \param strategy strategy for envelope calculation
    \param expand_strategy strategy for partitions
    \param source_index index to assign to the ring_identifiers
    \param max_count maximal number of points per section
        (defaults to 10, this seems to give the fastest results)

 */
template
<
    bool Reverse,
    typename DimensionVector,
    typename Geometry,
    typename Sections,
    typename RobustPolicy,
    typename Strategy
>
inline void sectionalize(Geometry const& geometry,
                RobustPolicy const& robust_policy,
                Sections& sections,
                Strategy const& strategy,
                int source_index = 0,
                std::size_t max_count = 10)
{
    concepts::check<Geometry const>();

    using section_type = typename boost::range_value<Sections>::type;

    // Compiletime check for point type of section boxes
    // and point type related to robust policy
    using ctype1 = typename geometry::coordinate_type
    <
        typename section_type::box_type
    >::type;
    using ctype2 = typename geometry::coordinate_type
    <
        typename geometry::robust_point_type
        <
            typename geometry::point_type<Geometry>::type,
            RobustPolicy
        >::type
    >::type;

    BOOST_STATIC_ASSERT((std::is_same<ctype1, ctype2>::value));


    sections.clear();

    ring_identifier ring_id;
    ring_id.source_index = source_index;

    dispatch::sectionalize
        <
            typename tag<Geometry>::type,
            Geometry,
            Reverse,
            DimensionVector
        >::apply(geometry, robust_policy, sections,
                 strategy,
                 ring_id, max_count);

    for (auto& section : sections)
    {
        strategy.postprocess_section_box().apply(section.bounding_box);
    }
}


template
<
    bool Reverse,
    typename DimensionVector,
    typename Geometry,
    typename Sections,
    typename RobustPolicy
>
inline void sectionalize(Geometry const& geometry,
                         RobustPolicy const& robust_policy,
                         Sections& sections,
                         int source_index = 0,
                         std::size_t max_count = 10)
{
    using box_type = typename Sections::box_type;
    using strategy_type = typename strategies::envelope::services::default_strategy
        <
            Geometry,
            box_type
        >::type;

    boost::geometry::sectionalize
        <
            Reverse, DimensionVector
        >(geometry, robust_policy, sections,
          strategy_type(),
          source_index, max_count);
}

}} // namespace boost::geometry


#endif // BOOST_GEOMETRY_ALGORITHMS_DETAIL_SECTIONS_SECTIONALIZE_HPP
