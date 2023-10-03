// Boost.Geometry

// Copyright (c) 2023 Barend Gehrels, Amsterdam, the Netherlands.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_ALGORITHMS_DETAIL_OVERLAY_COLOCATE_CLUSTERS_HPP
#define BOOST_GEOMETRY_ALGORITHMS_DETAIL_OVERLAY_COLOCATE_CLUSTERS_HPP

#include <boost/iterator/iterator_facade.hpp>
#include <boost/iterator/iterator_categories.hpp>

namespace boost { namespace geometry
{

#ifndef DOXYGEN_NO_DETAIL
namespace detail { namespace overlay
{

template <typename IndexSetIt, typename Turns>
struct cluster_points_iterator
    : public boost::iterator_facade
        <
            cluster_points_iterator<IndexSetIt, Turns>,
            typename Turns::value_type::point_type,
            boost::forward_traversal_tag
        >
{
    cluster_points_iterator(IndexSetIt const& it, Turns& turns) : it(it), turns(turns) {}

private:
    inline auto& dereference() const
    {
        return turns[*it].point;
    }

    inline bool equal(cluster_points_iterator<IndexSetIt, Turns> const& other) const
    {
        return it == other.it;
    }

    inline void increment()
    {
        ++it;
    }

    friend class boost::iterator_core_access;
    IndexSetIt it;
    Turns& turns;
};

// Moves intersection points per cluster such that they are identical.
// Because clusters are intersection close together, and
// handled as one location. Then they should also have one location.
// It is necessary to avoid artefacts and invalidities.
template <typename Clusters, typename Turns, typename Strategy>
inline void colocate_clusters(Clusters const& clusters, Turns& turns, Strategy const& strategy)
{
    for (auto const& pair : clusters)
    {
        auto const& turn_indices = pair.second.turn_indices;
        if (turn_indices.size() < 2)
        {
            // Defensive check
            continue;
        }
        cluster_points_iterator<decltype(turn_indices.cbegin()), Turns>
            begin(turn_indices.cbegin(), turns), end(turn_indices.cend(), turns);
        strategy.cluster_colocate(begin, end).apply(begin, end);
    }
}


}} // namespace detail::overlay
#endif //DOXYGEN_NO_DETAIL


}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_ALGORITHMS_DETAIL_OVERLAY_COLOCATE_CLUSTERS_HPP
