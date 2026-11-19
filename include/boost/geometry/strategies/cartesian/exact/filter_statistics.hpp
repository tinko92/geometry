// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.
#ifndef BOOST_GEOMETRY_STRATEGIES_CARTESIAN_EXACT_FILTER_STATISTICS_HPP
#define BOOST_GEOMETRY_STRATEGIES_CARTESIAN_EXACT_FILTER_STATISTICS_HPP

#include <cstdint>

namespace boost { namespace geometry { namespace strategies { namespace detail
{ namespace exact_cartesian
{

struct filter_statistics
{
    std::uint64_t native_exact_values = 0;
    std::uint64_t native_rationals = 0;
    std::uint64_t lazy_exact_values = 0;
    std::uint64_t point_exact_values = 0;
    std::uint64_t point_coordinate_values = 0;
    std::uint64_t homogeneous_values = 0;
    std::uint64_t ratio_values = 0;
    std::uint64_t construction_filters = 0;
    std::uint64_t native_point_constructions = 0;
    std::uint64_t small_integer_measure_comparisons = 0;
    std::uint64_t compact_measure_comparisons = 0;
    std::uint64_t full_measure_comparisons = 0;
    std::uint64_t sterbenz_orientations = 0;
    std::uint64_t tagged_orientation_calls = 0;
};

inline filter_statistics& filtering_statistics()
{
    static thread_local filter_statistics value;
    return value;
}

}}}}} // namespace boost::geometry::strategies::detail::exact_cartesian

#endif
