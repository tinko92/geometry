// Boost.Geometry

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_STRATEGY_PRECEDING_HPP
#define BOOST_GEOMETRY_STRATEGY_PRECEDING_HPP

#include <boost/geometry/core/static_assert.hpp>

namespace boost { namespace geometry
{


namespace strategy { namespace preceding { namespace services
{

/*!
\brief Traits class binding a default preceding strategy to a coordinate system
\ingroup util
\tparam CSTag tag of coordinate system
\tparam CalculationType \tparam_calculation
*/
template <typename CSTag, typename CalculationType = void>
struct default_strategy
{
    BOOST_GEOMETRY_STATIC_ASSERT_FALSE(
        "Not implemented for this type.",
        CSTag);
};

}}} // namespace strategy::preceding::services


}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_STRATEGY_PRECEDING_HPP
