// Boost.Geometry

// Licensed under the Boost Software License version 1.0.
// http://www.boost.org/users/license.html

#ifndef BOOST_GEOMETRY_GEOMETRIES_CONCEPTS_DETAIL_MUTABLE_RANGE_HPP
#define BOOST_GEOMETRY_GEOMETRIES_CONCEPTS_DETAIL_MUTABLE_RANGE_HPP

#include <concepts>
#include <cstddef>
#include <ranges>
#include <utility>

#include <boost/geometry/core/mutable_range.hpp>


namespace boost { namespace geometry { namespace concepts { namespace detail
{

template <typename Range>
concept ConstForwardRange = std::ranges::forward_range<Range const>;

template <typename Range>
concept ConstRandomAccessRange = std::ranges::random_access_range<Range const>;

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
