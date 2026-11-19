// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.

// Reuse the unchanged corpus and call boundary; profile complete call batches.
#define main compact_benchmark_main
#include "simple_cartesian_benchmark.cpp"
#undef main
#include <valgrind/callgrind.h>

int main(int argc, char** argv)
{
    try
    {
        if (argc != 4) throw std::invalid_argument("Expected category, operation, rounds");
        auto const corpus = make_inputs(argv[1]);
        auto const op = operation_code(argv[2]);
        auto const rounds = std::stoul(argv[3]);
        if (rounds == 0) throw std::invalid_argument("Rounds must be positive");
        polygons out;
        lines linear;
        for (auto const& g : corpus) execute(op, g, out, linear);
        execute(op, corpus.front(), out, linear);
        std::uint64_t checksum = 0;
        CALLGRIND_START_INSTRUMENTATION;
        CALLGRIND_ZERO_STATS;
        for (std::size_t i = 0; i < rounds; ++i)
            for (auto const& g : corpus) checksum += execute(op, g, out, linear);
        CALLGRIND_STOP_INSTRUMENTATION;
        std::cout << "category,operation,calls,checksum\n"
                  << argv[1] << ',' << argv[2] << ',' << rounds * corpus.size()
                  << ',' << checksum << '\n';
        return 0;
    }
    catch (std::exception const& e)
    {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
