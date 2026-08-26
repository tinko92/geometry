// Boost.Geometry

// Licensed under the Boost Software License version 1.0.
// http://www.boost.org/users/license.html

#ifndef BOOST_GEOMETRY_GEOMETRIES_CONCEPTS_DETAIL_MUTABLE_RANGE_HPP
#define BOOST_GEOMETRY_GEOMETRIES_CONCEPTS_DETAIL_MUTABLE_RANGE_HPP

#include <concepts>
#include <cstddef>
#include <utility>

#include <boost/range/begin.hpp>
#include <boost/range/end.hpp>
#include <boost/range/iterator.hpp>
#include <boost/range/traversal.hpp>

#include <boost/geometry/core/mutable_range.hpp>


namespace boost { namespace geometry { namespace concepts { namespace detail
{

template <typename Range>
concept ConstForwardRange =
    requires(Range const& range)
    {
        boost::begin(range);
        boost::end(range);
        typename boost::range_iterator<Range const>::type;
        typename boost::range_traversal<Range const>::type;
    }
    && std::derived_from
        <typename boost::range_traversal<Range const>::type,
         boost::forward_traversal_tag>;

template <typename Range>
concept ConstRandomAccessRange =
    ConstForwardRange<Range>
    && std::derived_from
        <typename boost::range_traversal<Range const>::type,
         boost::random_access_traversal_tag>;

template <typename Range, typename Value>
concept MutableRange =
    ConstRandomAccessRange<Range>
    && requires(Range& range, Value const& value, Value&& rvalue)
    {
        traits::clear<Range>::apply(range);
        traits::resize<Range>::apply(range, std::size_t{});
        traits::push_back<Range>::apply(range, value);
        traits::push_back<Range>::apply(range, std::move(rvalue));
    };

}}}} // namespace boost::geometry::concepts::detail

#endif // BOOST_GEOMETRY_GEOMETRIES_CONCEPTS_DETAIL_MUTABLE_RANGE_HPP
