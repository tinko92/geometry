// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.

#ifndef BOOST_GEOMETRY_STRATEGIES_CARTESIAN_EXACT_OUTPUT_POLICY_HPP
#define BOOST_GEOMETRY_STRATEGIES_CARTESIAN_EXACT_OUTPUT_POLICY_HPP

namespace boost { namespace geometry { namespace strategies { namespace relate
{

// Emit the approximation already stored at the output boundary. Predicates
// remain exact, but emission neither refines nor certifies these coordinates.
struct approximate_output {};

// Round floating output to nearest at emission. Root-bound operations use the
// deferred one-ULP layout; legacy point compatibility keeps stored error bounds.
// Integral output follows BG numeric_cast conversion under every policy.
struct nearest_output {};

}}}} // namespace boost::geometry::strategies::relate

#endif
