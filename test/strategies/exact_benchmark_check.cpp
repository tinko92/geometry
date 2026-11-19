// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.
// The registered test defines BG_BENCH_EXACT. The same file may be built with
// native baseline headers for diagnosis; --wkt prints complete output.
#define main benchmark_main
#include "exact_cartesian_benchmark.cpp"
#undef main
#include <boost/multiprecision/cpp_int.hpp>

int main(int argc, char** argv)
{
    bool const dump = argc == 2 && std::string(argv[1]) == "--wkt";
    unsigned failures = 0;
    options opt;
#ifdef BG_BENCH_EXACT
    if (std::is_same<coordinate, double>::value) { projected_self_test(opt); }
#endif
    auto const inputs = make_inputs("star", opt);
    std::cout << std::setprecision((std::numeric_limits<coordinate>::max_digits10));
    for (std::size_t i = 0; i < inputs.size(); ++i)
    {
#ifdef BG_BENCH_EXACT
        namespace ec = bg::strategies::detail::exact_cartesian;
        ec::point_cache_scope<coordinate> operation;
        using exact_polygon = bg::model::polygon<ec::exact_point<coordinate>>;
        bg::model::multi_polygon<exact_polygon> out;
#else
        polygons out;
#endif
        bg::sym_difference(inputs[i].a, inputs[i].b, out, strategy());
        std::string reason;
        using rational = boost::multiprecision::cpp_rational;
        using rational_point = bg::model::d2::point_xy<rational>;
        bg::model::multi_polygon<bg::model::polygon<rational_point>> lifted;
        bg::convert(out, lifted);
        std::string exact_reason;
        bool const rational_valid = bg::is_valid(lifted, exact_reason,
            bg::strategies::relate::cartesian<rational>());
        if (!rational_valid || out.size() != 80) { ++failures; }
#ifdef BG_BENCH_EXACT
        // Approximate export need not be valid. Check the exact result above,
        // and that conversion does not change its selected ring structure.
        polygons exported;
        strategy::append_output(out, exported);
        if (exported.size() != out.size() || bg::num_points(exported) != bg::num_points(out))
            ++failures;
        for (std::size_t j = 0; j < (std::min)(out.size(), exported.size()); ++j)
            if (out[j].inners().size() != exported[j].inners().size()) ++failures;
#endif
        std::cout << i << ' ' << out.size() << ' ' << bg::num_points(out) << ' '
                  << bg::is_valid(out, reason) << ' ' << reason << "; rational: "
                  << rational_valid << ' ' << exact_reason << '\n';
        if (dump) { std::cout << bg::wkt(out) << '\n'; }
    }
    return failures ? 1 : 0;
}
