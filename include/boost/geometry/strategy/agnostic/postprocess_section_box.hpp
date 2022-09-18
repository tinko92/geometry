// Boost.Geometry (aka GGL, Generic Geometry Library)

// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_STRATEGY_AGNOSTIC_POSTPROCESS_SECTION_BOX_HPP
#define BOOST_GEOMETRY_STRATEGY_AGNOSTIC_POSTPROCESS_SECTION_BOX_HPP

#include <boost/geometry/algorithms/detail/expand_by_epsilon.hpp>

#include <boost/geometry/core/coordinate_type.hpp>

#include <boost/geometry/util/math.hpp>

namespace boost { namespace geometry
{

namespace strategy { namespace postprocess_section_box
{

template <int N>
struct expand_by_neps
{
    template <typename Box>
    static void apply(Box& b)
    {
        using ct = typename geometry::coordinate_type<Box>::type;
        static ct const eps = math::scaled_epsilon<ct>(N);
        detail::expand_by_epsilon(b, eps);
    }
};

}} // namespace strategy::postprocess_section_box

}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_STRATEGY_AGNOSTIC_POSTPROCESS_SECTION_BOX_HPP
