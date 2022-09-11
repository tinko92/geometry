// Boost.Geometry (aka GGL, Generic Geometry Library)

// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_STRATEGY_CARTESIAN_DIRECTIONS_HPP
#define BOOST_GEOMETRY_STRATEGY_CARTESIAN_DIRECTIONS_HPP

#include <array>

#include <boost/geometry/core/access.hpp>
#include <boost/geometry/core/coordinate_dimension.hpp>

#include <boost/geometry/util/sequence.hpp>

namespace boost { namespace geometry
{

namespace strategy { namespace directions
{

namespace detail
{
template
<
    typename DimensionVector,
    std::size_t Index,
    std::size_t Count
>
struct get_direction_loop_cart
{
    using dimension = typename util::sequence_element<Index, DimensionVector>::type;

    template <typename Segment, typename Arr>
    static inline void apply(Segment const& seg,
                             Arr& directions)
    {
        auto const& c0 = geometry::get<0, dimension::value>(seg);
        auto const& c1 = geometry::get<1, dimension::value>(seg);

        directions[Index] = c1 > c0 ? 1 : c1 < c0 ? -1 : 0;

        get_direction_loop_cart
        <
            DimensionVector,
            Index + 1,
            Count
        >::apply(seg, directions);
    }
};

template
<
    typename DimensionVector,
    std::size_t Count
>
struct get_direction_loop_cart<DimensionVector, Count, Count>
{
    template <typename Segment, typename Arr>
    static inline void apply(Segment const&, Arr&)
    {}
};

} // detail

struct cartesian
{
    template <typename DimensionVector, typename Segment>
    static auto apply(Segment const& seg)
    {
        constexpr auto dimensions = util::sequence_size<DimensionVector>::value;
        std::array<int, dimensions> directions;
        detail::get_direction_loop_cart<DimensionVector, 0, dimensions>::apply(seg, directions);
        return directions;
    }
};

}} // namespace strategy::directions

}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_STRATEGY_CARTESIAN_DIRECTIONS_HPP
