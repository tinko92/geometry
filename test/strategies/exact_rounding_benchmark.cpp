// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.

// Standalone numerical microbenchmark, not a whole-algorithm benchmark.
// clang++ -std=c++17 -O3 -DNDEBUG -ffp-contract=off -Iinclude -I../..
//   test/strategies/exact_rounding_benchmark.cpp -o /tmp/bg-rounding-bench
// /tmp/bg-rounding-bench --seconds 0.1 --rounds 3 --count 128
// /tmp/bg-rounding-bench --check-only 1 --rounds 1 --count 128
#include "exact_rounding_policy_helpers.hpp"

#include <cassert>
#include <cfenv>
#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>

namespace rp = rounding_policies;
using clock_type = std::chrono::steady_clock;
volatile std::uint64_t benchmark_sink;

struct scenario
{
    char const* name;
    unsigned uses, emit_every;
    rp::predicate kind = rp::predicate::coordinate;
};
std::array<scenario, 16> const scenarios{{
    {"construct_never_emitted", 0, 0}, {"one_use_emitted", 1, 1},
    {"many_uses_emitted", 32, 1}, {"many_uses_never_emitted", 32, 0},
    {"one_use_quarter_emitted", 1, 4},
    {"nni_one_use_emitted", 1, 1, rp::predicate::nni},
    {"nni_many_uses_emitted", 32, 1, rp::predicate::nni},
    {"nni_many_uses_never_emitted", 32, 0, rp::predicate::nni},
    {"nii_one_use_emitted", 1, 1, rp::predicate::nii},
    {"nii_many_uses_emitted", 32, 1, rp::predicate::nii},
    {"nii_many_uses_never_emitted", 32, 0, rp::predicate::nii},
    {"iii_one_use_emitted", 1, 1, rp::predicate::iii},
    {"iii_many_uses_emitted", 32, 1, rp::predicate::iii},
    {"iii_many_uses_never_emitted", 32, 0, rp::predicate::iii},
    {"nni_on_support_emitted", 32, 1, rp::predicate::nni_support},
    {"nii_source_anchor_emitted", 32, 1, rp::predicate::nii_source}
}};

struct options
{
    double seconds = 0.1;
    unsigned rounds = 3;
    std::size_t count = 128;
    std::uint64_t seed = 0x524f554e44494e47;
    std::string shape = "all", scenario_name = "all";
    bool check_only = false;
};

options parse(int argc, char** argv)
{
    options o;
    for (int i = 1; i < argc; ++i)
    {
        std::string const key = argv[i];
        if (i + 1 == argc) throw std::invalid_argument("missing value for " + key);
        std::string const value = argv[++i];
        std::size_t end = 0;
        if (key == "--shape") o.shape = value;
        else if (key == "--scenario") o.scenario_name = value;
        else if (key == "--check-only")
        {
            if (value != "0" && value != "1") throw std::invalid_argument("check-only requires 0 or 1");
            o.check_only = value == "1";
        }
        else if (key == "--seconds")
        {
            o.seconds = std::stod(value, &end);
            if (end != value.size()) throw std::invalid_argument("invalid seconds");
        }
        else if (key == "--rounds" || key == "--count" || key == "--seed")
        {
            auto const n = std::stoull(value, &end, 0);
            if (end != value.size() || value[0] == '-') throw std::invalid_argument("invalid integer option");
            if (key == "--rounds")
            {
                if (n > 1000) throw std::invalid_argument("rounds must be <=1000");
                o.rounds = unsigned(n);
            }
            else if (key == "--count") o.count = n;
            else o.seed = n;
        }
        else throw std::invalid_argument("unknown option " + key);
    }
    if (!(o.seconds > 0 && o.seconds <= 60) || o.count == 0 || o.count > 100000 || o.rounds == 0)
        throw std::invalid_argument("seconds in (0,60], count in [1,100000], rounds in [1,1000] required");
    return o;
}

inline void mix(std::uint64_t& digest, std::uint64_t value)
{ digest ^= value + 0x9e3779b97f4a7c15ULL + (digest << 6) + (digest >> 2); }

template <rp::policy Policy, bool Instrument>
std::uint64_t batch(rp::operation<Policy, Instrument>& op, std::vector<rp::point>& points,
                    std::vector<rp::sample> const& samples, scenario const& s)
{
    op.reset();
    for (std::size_t i = 0; i < samples.size(); ++i) points[i] = op.construct(samples[i].input);
    std::uint64_t digest = 0;
    for (std::size_t i = 0; i < samples.size(); ++i)
    {
        for (unsigned j = 0; j < s.uses; ++j)
        {
            mix(digest, std::uint64_t(rp::evaluate_query(op, points, samples, i, j, s.kind) + 1));
        }
        if (s.emit_every != 0 && i % s.emit_every == 0)
        {
            auto const out = op.emit(points[i]);
            mix(digest, rp::bits(out[0])); mix(digest, rp::bits(out[1]));
        }
        // Consume never-emitted values too. Policies need not have the same
        // internal approximation, so their checksums need not match here.
        else { mix(digest, rp::bits(points[i].approx[0])); mix(digest, rp::bits(points[i].approx[1])); }
    }
    return digest;
}

template <rp::policy Policy>
void measure(options const& o, std::string const& shape, std::vector<rp::sample> const& samples,
             scenario const& s, unsigned round)
{
    // Full rational-oracle checking, allocation and diagnostics are OUTSIDE
    // timing. Instrument=false compiles all counter updates out of timed code.
    rp::verify<Policy>(samples, s.uses, s.emit_every, s.kind);
    rp::operation<Policy, true> diagnostic(samples.size());
    std::vector<rp::point> points(samples.size());
    std::uint64_t digest = batch(diagnostic, points, samples, s);
    benchmark_sink = digest;
    auto const& c = diagnostic.counts;
    std::uint64_t batches = 0;
    double elapsed = 0;
    if (!o.check_only)
    {
        rp::operation<Policy> op(samples.size());
        benchmark_sink = batch(op, points, samples, s); // warmup, including page faults
        auto const start = clock_type::now();
        digest = 0;
        do
        {
            mix(digest, batch(op, points, samples, s));
            ++batches;
            elapsed = std::chrono::duration<double>(clock_type::now() - start).count();
        }
        while (elapsed < o.seconds);
    }
    benchmark_sink = digest;
    double const ns = batches == 0 ? std::numeric_limits<double>::quiet_NaN()
        : elapsed * 1e9 / double(batches * samples.size());
    std::cout << rp::name(Policy) << ',' << shape << ',' << s.name << ',' << round << ','
              << samples.size() << ',' << s.uses << ',' << s.emit_every << ',' << batches << ','
              << elapsed << ',' << ns << ','
              << c.zero_denominators << ',' << c.expansion_entries << ',' << c.fixed_entries << ','
              << c.ulp_certifications << ',' << c.nearest_certifications << ',' << c.boundary_tests << ','
              << c.candidate_steps << ',' << c.estimate_refinements << ',' << c.fixed_rounds << ','
              << c.initial_nearest << ',' << c.emission_nearest << ',' << c.predicate_filters << ','
              << c.predicate_exact << ',' << c.emissions << ',' << digest << ',' << rp::name(s.kind) << ','
              << c.orientation_filters[0] << ',' << c.orientation_exact[0] << ','
              << c.orientation_filters[1] << ',' << c.orientation_exact[1] << ','
              << c.orientation_filters[2] << ',' << c.orientation_exact[2] << ','
              << c.constructions << ',' << (o.check_only ? "check_only" : "timed") << '\n';
}

int main(int argc, char** argv)
{
    try
    {
        assert(std::fegetround() == FE_TONEAREST);
        auto const o = parse(argc, argv);
        std::cerr << "Numerical binary64 policies only; cache capacity reserved before timing.\n"
                  << "Counters are from ONE separate diagnostic batch, not the timed batches.\n"
                  << "ULP predicate filters; existing expansion construction certificates, not whole-overlay timing.\n"
                  << "All constructed points get cache entries, even shared native endpoints; no native reuse fast path.\n"
                  << "General queries and on-support/source-anchor queries are reported separately.\n"
                  << "predicate_filters/exact count coordinate queries; orientation counters are separate.\n"
                  << "point_bytes=" << sizeof(rp::point) << " cache_entry_bytes=" << sizeof(rp::precision_entry)
                  << " seed=" << o.seed << '\n';
        std::cout << std::setprecision(10)
                  << "policy,shape,scenario,round,corpus,uses,emit_every,batches,seconds,ns_per_point,"
                     "zero_denominators,expansion_entries,fixed_entries,ulp_certificates,nearest_certificates,"
                     "boundary_tests,candidate_steps,estimate_refinements,fixed_rounds,initial_nearest,"
                     "emission_nearest,predicate_filters,predicate_exact,emissions,digest,query_family,"
                     "nni_filters,nni_exact,nii_filters,nii_exact,iii_filters,iii_exact,cache_entries,mode\n";
        bool selected_shape = false, selected_scenario = false;
        for (char const* shape : {"ordinary", "near_parallel", "shared_boundary"})
        {
            if (o.shape != "all" && o.shape != shape) continue;
            selected_shape = true;
            auto const samples = rp::corpus(shape, o.count, o.seed);
            for (unsigned round = 0; round < o.rounds; ++round)
            for (auto const& s : scenarios)
            {
                if (o.scenario_name != "all" && o.scenario_name != s.name) continue;
                selected_scenario = true;
                // Rotate order to avoid always favouring one thermal/cache state.
                if ((round & 1) == 0)
                {
                    measure<rp::policy::eager>(o, shape, samples, s, round);
                    measure<rp::policy::deferred>(o, shape, samples, s, round);
                }
                else
                {
                    measure<rp::policy::deferred>(o, shape, samples, s, round);
                    measure<rp::policy::eager>(o, shape, samples, s, round);
                }
            }
        }
        if (!selected_shape || !selected_scenario) throw std::invalid_argument("unknown shape or scenario");
    }
    catch (std::exception const& e)
    {
        std::cerr << e.what() << '\n';
        return EXIT_FAILURE;
    }
}
