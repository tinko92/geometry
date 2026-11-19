// Boost.Geometry
// Copyright (c) 2026
// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_ALGORITHMS_DETAIL_BUFFER_ROUND_POLYGON_BACKEND_HPP
#define BOOST_GEOMETRY_ALGORITHMS_DETAIL_BUFFER_ROUND_POLYGON_BACKEND_HPP

#include <cstddef>
#include <stdexcept>

#include <boost/geometry/algorithms/buffer_polygon_certified.hpp>
#include <boost/geometry/algorithms/erode_convex_polygon_certified.hpp>
#include <boost/geometry/algorithms/is_valid.hpp>
#include <boost/geometry/strategies/relate/cartesian.hpp>

namespace boost { namespace geometry { namespace detail { namespace certified_round_buffer
{

namespace common = certified_point_buffer;
using common::exact_number;
using status = certified_buffer_status;

// Validation for the unsupported domain, not a weaker export certificate.
// Preserve stored orientation/closure and use rational predicates; the
// chord backend's disjoint-boundary precondition is stronger than validity.
template <typename Polygon>
status classify_polygon(Polygon const& input)
{
    using point = model::point<exact_number, 2, cs::cartesian>;
    using polygon = model::polygon<point, geometry::point_order<Polygon>::value == clockwise,
                                         closure<Polygon>::value == closed>;
    polygon lifted;
    std::size_t positions = 0;
    auto append = [&](auto const& from, auto& to)
    {
        for (auto const& p : from)
        {
            if (++positions > 4096) { return false; }
            auto const q = common::read_point(p);
            if (!certified_convex_erosion::bounded(q.x) || !certified_convex_erosion::bounded(q.y))
            { return false; }
            to.emplace_back(q.x, q.y);
        }
        return true;
    };
    try
    {
        if (!append(exterior_ring(input), lifted.outer())) { return status::cannot_certify; }
        if (lifted.outer().empty()) { return status::invalid_input; }
        for (auto const& hole : interior_rings(input))
        {
            if (lifted.inners().size() == 4096) { return status::cannot_certify; }
            lifted.inners().emplace_back();
            if (!append(hole, lifted.inners().back())) { return status::cannot_certify; }
        }
        return geometry::is_valid(lifted, strategies::relate::cartesian<exact_number>())
            ? status::cannot_certify : status::invalid_input;
    }
    catch (std::domain_error const&) { return status::invalid_input; }
}

// Integration hook for generic polygon dilation/erosion. Inputs are already
// numeric-parameter checked. Validate the source in the backend. Respect
// both budgets and the actual output model. On success publish only certified output, with
// empty meaning true regularized empty, not a failed export or an omitted
// face. A negative-distance backend must not use chord-band erosion as an
// exact-topology fallback. No change to the public dispatcher is needed to
// extend this hook to bounded native export or irrational side lengths.
template <typename Polygon, typename MultiPolygon>
status polygon_backend(Polygon const& input, MultiPolygon& output,
                       exact_number const& distance, exact_number const& sagitta,
                       exact_number const& /*max_displacement*/)
{
    if (distance > 0 && sagitta > 0)
    {
        // When it succeeds this backend exports the rational chord union
        // losslessly, satisfying any nonnegative native export budget.
        auto const result = geometry::buffer_polygon_certified(input, output, distance, sagitta);
        if (result != status::invalid_input) { return result; }
    }
    return classify_polygon(input);
}

}}}} // namespace boost::geometry::detail::certified_round_buffer

#endif // BOOST_GEOMETRY_ALGORITHMS_DETAIL_BUFFER_ROUND_POLYGON_BACKEND_HPP
