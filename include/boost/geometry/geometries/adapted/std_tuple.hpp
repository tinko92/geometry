// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2026 Tinko Bartels

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_GEOMETRIES_ADAPTED_STD_TUPLE_HPP
#define BOOST_GEOMETRY_GEOMETRIES_ADAPTED_STD_TUPLE_HPP


#include <cstddef>
#include <tuple>
#include <type_traits>

#include <boost/geometry/core/access.hpp>
#include <boost/geometry/core/coordinate_dimension.hpp>
#include <boost/geometry/core/coordinate_type.hpp>
#include <boost/geometry/core/point_type.hpp>
#include <boost/geometry/core/tags.hpp>


namespace boost { namespace geometry
{


#ifndef DOXYGEN_NO_TRAITS_SPECIALIZATIONS
namespace traits
{


template <typename T, typename ...Ts>
struct tag<std::tuple<T, Ts...>>
{
    using type = point_tag;
};


template <typename T, typename ...Ts>
struct coordinate_type<std::tuple<T, Ts...>>
{
    using type = T;
};


template <typename T, typename ...Ts>
struct dimension<std::tuple<T, Ts...>>
    : std::integral_constant<std::size_t, 1 + sizeof...(Ts)>
{};


template <typename T, typename ...Ts, std::size_t Dimension>
struct access<std::tuple<T, Ts...>, Dimension>
{
    static inline T get(std::tuple<T, Ts...> const& point)
    {
        return std::get<Dimension>(point);
    }

    static inline void set(std::tuple<T, Ts...>& point, T const& value)
    {
        std::get<Dimension>(point) = value;
    }
};


} // namespace traits
#endif // DOXYGEN_NO_TRAITS_SPECIALIZATIONS


}} // namespace boost::geometry


#define BOOST_GEOMETRY_REGISTER_STD_TUPLE_CS(CoordinateSystem) \
    namespace boost { namespace geometry { namespace traits { \
    template <typename T, typename ...Ts> \
    struct coordinate_system<std::tuple<T, Ts...>> \
    { \
        using type = CoordinateSystem; \
    }; \
    }}}


#endif // BOOST_GEOMETRY_GEOMETRIES_ADAPTED_STD_TUPLE_HPP
