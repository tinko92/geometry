// Boost.Geometry (aka GGL, Generic Geometry Library)

// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_STRATEGY_SPHERICAL_DIRECTIONS_HPP
#define BOOST_GEOMETRY_STRATEGY_SPHERICAL_DIRECTIONS_HPP

#include <array>

#include <boost/geometry/core/access.hpp>
#include <boost/geometry/core/coordinate_dimension.hpp>
#include <boost/geometry/core/coordinate_type.hpp>
#include <boost/geometry/core/cs.hpp>

#include <boost/geometry/util/normalize_spheroidal_coordinates.hpp>

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
struct get_direction_loop_sph
{
    using dimension = typename util::sequence_element<Index, DimensionVector>::type;

    template <typename Segment>
    static inline void apply(Segment const& seg,
                std::array<int, Count>& directions)
    {
        auto const& c0 = geometry::get<0, dimension::value>(seg);
        auto const& c1 = geometry::get<1, dimension::value>(seg);

        directions[Index] = c1 > c0 ? 1 : c1 < c0 ? -1 : 0;

        get_direction_loop_sph
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
struct get_direction_loop_sph<DimensionVector, 0, Count>
{
    template <typename Segment>
    static inline void apply(Segment const& seg,
                std::array<int, Count>& directions)
    {
        using coordinate_type = typename coordinate_type<Segment>::type;
        using units_t = typename geometry::detail::cs_angular_units<Segment>::type;

        coordinate_type const diff = math::longitude_distance_signed
                                        <
                                            units_t, coordinate_type
                                        >(geometry::get<0, 0>(seg),
                                          geometry::get<1, 0>(seg));

        coordinate_type zero = coordinate_type();
        directions[0] = diff > zero ? 1 : diff < zero ? -1 : 0;

        get_direction_loop_sph
        <
            DimensionVector,
            1,
            Count
        >::apply(seg, directions);
    }
};

template
<
    typename DimensionVector,
    std::size_t Count
>
struct get_direction_loop_sph<DimensionVector, Count, Count>
{
    template <typename Segment>
    static inline void apply(Segment const&, std::array<int, Count>&)
    {}
};

} // detail

struct spherical
{
    template <typename DimensionVector, typename Segment>
    static auto apply(Segment const& seg)
    {
        constexpr auto dimensions = util::sequence_size<DimensionVector>::value;
        std::array<int, dimensions> directions;
        detail::get_direction_loop_sph<DimensionVector, 0, dimensions>::apply(seg, directions);
        return directions;
    }
};

}} // namespace strategy::directions

}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_STRATEGY_SPHERICAL_DIRECTIONS_HPP
