// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.
// Diagnostic only: requires the isolated profiling header overlay.
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
    std::cout << "coordinate,inputs,workload,operation,calls,views,validations,owned_points,native_sides,side_filter,side_interval,side_expansion,side_fixed,constructed,construction_forward,construction_secondary,coordinate_native,coordinate_filter,coordinate_refine,coordinate_polynomial,cache_entries,cache_expansions,fixed_values,output_certifications,corpus_digest\n";
    for (bool fractional : {false, true})
    {
        if (fractional && std::is_integral<coordinate>::value) continue;
        for (char const* shape : {"disjoint", "convex", "star", "holes", "shared"})
        {
            auto inputs = make_inputs(shape, opt);
            if (fractional) for (auto& g : inputs)
            {
                auto shift = [](point& p)
                {
                    bg::set<0>(p, coordinate(bg::get<0>(p) * 0.125 + 0.0625));
                    bg::set<1>(p, coordinate(bg::get<1>(p) * 0.125 + 0.0625));
                };
                bg::for_each_point(g.a,shift); bg::for_each_point(g.b,shift);
                bg::for_each_point(g.la,shift); bg::for_each_point(g.lb,shift); shift(g.query);
            }
            for (char const* operation : operation_names)
            {
                ec::filtering_statistics() = {};
                std::uint64_t digest = 0;
                for (auto const& g : inputs) digest = mix(digest,execute(operation_code(operation),g,out,linear));
                auto const& s = ec::filtering_statistics();
                std::cout << coordinate_name << ',' << (fractional ? "dyadic" : "integer") << ','
                    << shape << ',' << operation << ',' << inputs.size() << ',' << s.views << ','
                    << s.validations << ',' << s.owned_points << ',' << s.native_sides << ','
                    << s.side_filter << ',' << s.side_interval << ',' << s.side_expansion << ','
                    << s.side_fixed << ',' << s.constructed << ',' << s.construction_forward << ','
                    << s.construction_secondary << ',' << s.coordinate_native << ',' << s.coordinate_filter << ','
                    << s.coordinate_refine << ',' << s.coordinate_polynomial << ',' << s.cache_entries << ','
                    << s.cache_expansions << ',' << s.native_exact_values << ',' << s.construction_filters << ','
                    << digest << '\n';
            }
        }
    }
}
