// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.
// Diagnostic only; never use this instrumented build for timings.
#define BG_BENCH_EXACT
#define BOOST_GEOMETRY_EXACT_FILTER_STATISTICS
#define main benchmark_main
#include "exact_cartesian_benchmark.cpp"
#undef main

int main()
{
    namespace ec = bg::strategies::detail::exact_cartesian;
    options opt;
    polygons out;
    lines linear;
    std::cout << "coordinate,workload,operation,calls,native_lifts,lazy_values,point_coordinates,homogeneous_values,ratio_values,construction_filters\n";
    for (char const* shape : {"disjoint", "convex", "star", "holes", "shared"})
    {
        auto const inputs = make_inputs(shape, opt);
        for (char const* name : operation_names)
        {
            auto const operation = operation_code(name);
            ec::filtering_statistics() = {};
            for (auto const& input : inputs) { sink = execute(operation, input, out, linear); }
            auto const& s = ec::filtering_statistics();
            std::cout << coordinate_name << ',' << shape << ',' << name << ',' << inputs.size()
                      << ',' << s.native_rationals << ',' << s.lazy_exact_values
                      << ',' << s.point_coordinate_values << ',' << s.homogeneous_values
                      << ',' << s.ratio_values << ',' << s.construction_filters << '\n';
        }
    }
}
