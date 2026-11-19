// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.

// Standalone: clang++ -std=c++17 -O3 -DNDEBUG -ffp-contract=off -Iinclude
//   -I../.. test/strategies/exact_rounding_policies.cpp -o /tmp/bg-rounding-test
#include "exact_rounding_policy_helpers.hpp"
#include <boost/core/lightweight_test.hpp>

#include <cfenv>
#include <cstdlib>
#include <iostream>
#include <new>

namespace rp = rounding_policies;
namespace
{
bool count_allocations = false;
std::size_t allocations = 0;
}

void* operator new(std::size_t size)
{
    if (count_allocations) ++allocations;
    if (void* p = std::malloc(size == 0 ? 1 : size)) return p;
    throw std::bad_alloc();
}
void operator delete(void* p) noexcept { std::free(p); }
void* operator new[](std::size_t n) { return ::operator new(n); }
void operator delete[](void* p) noexcept { std::free(p); }
#if defined(__cpp_sized_deallocation)
void operator delete(void* p, std::size_t) noexcept { std::free(p); }
void operator delete[](void* p, std::size_t) noexcept { std::free(p); }
#endif

template <rp::policy Policy>
void check(std::vector<rp::sample> const& samples)
{
    rp::verify<Policy>(samples, 32, 1);
    auto const none = rp::verify<Policy>(samples, 32, 0);
    BOOST_TEST_EQ(none.emissions, 0u);
    BOOST_TEST_EQ(none.emission_nearest, 0u);
    rp::verify<Policy>(samples, 1, 4);
    for (auto kind : {rp::predicate::nni, rp::predicate::nii, rp::predicate::iii,
                      rp::predicate::nni_support, rp::predicate::nii_source})
    {
        rp::verify<Policy>(samples, 3, 1, kind);
        auto const unexported = rp::verify<Policy>(samples, 3, 0, kind);
        BOOST_TEST_EQ(unexported.emission_nearest, 0u);
    }

    rp::operation<Policy, true> op(samples.size());
    std::vector<rp::point> points(samples.size());
    // Only the cache vector grows. Once reserved, construction, fixed/expansion
    // arithmetic, predicates and output must not allocate numerical storage.
    allocations = 0;
    count_allocations = true;
    for (std::size_t i = 0; i < samples.size(); ++i) points[i] = op.construct(samples[i].input);
    for (std::size_t i = 0; i < samples.size(); ++i)
    {
        auto const& s = samples[i];
        auto& p = points[i];
        auto const& q = s.queries[1];
        op.compare(p, q.dimension, q.target);
        for (auto kind : {rp::predicate::nni, rp::predicate::nii, rp::predicate::iii,
                          rp::predicate::nni_support, rp::predicate::nii_source})
            rp::evaluate_query(op, points, samples, i, 0, kind);
        op.emit(p);
    }
    count_allocations = false;
    BOOST_TEST_EQ(allocations, 0u);

    op.reset();
    auto p = op.construct(samples.front().input);
    p.approx = {{17.0, -19.0}}; // Deliberately distinguish output from refinement.
    auto const before = op.counts.boundary_tests + op.counts.fixed_rounds;
    auto const out = op.emit_approximate(p);
    BOOST_TEST_EQ(out[0], 17.0);
    BOOST_TEST_EQ(out[1], -19.0);
    BOOST_TEST_EQ(before, op.counts.boundary_tests + op.counts.fixed_rounds);
}

std::vector<rp::sample> boundaries()
{
    std::vector<rp::sample> samples;
    for (int exponent : {-120, -1, 0, 1, 120})
    for (double sign : {-1.0, 1.0})
    {
        double const scale = std::ldexp(1.0, exponent);
        double const p = std::nextafter(1.0, 0.0), n = std::nextafter(1.0, 2.0);
        for (double lo : {p, 1.0, n})
        for (double height : {0.5, 1.0, 1.5})
        {
            double const hi = std::nextafter(lo, 2.0);
            rp::endpoints v{{sign * lo, 0, sign * hi, 2, 0, height, 3, height}};
            for (auto& x : v) x *= scale;
            // Both segment reversals and source exchange must preserve rounding.
            for (unsigned mask = 0; mask < 8; ++mask)
            {
                auto w = v;
                if (mask & 1) { std::swap(w[0], w[2]); std::swap(w[1], w[3]); }
                if (mask & 2) { std::swap(w[4], w[6]); std::swap(w[5], w[7]); }
                if (mask & 4) for (unsigned i = 0; i < 4; ++i) std::swap(w[i], w[i + 4]);
                samples.push_back(rp::prepare(w));
            }
        }
    }
    samples.push_back(rp::prepare({{0, -1, 0, 1, -1, 0, 1, 0}}));
    // Outside the conservative expansion range, but without FP overflow or
    // underflow. Exercise the format-derived fixed fallback in both directions.
    for (int exponent : {-200, 200})
    {
        double const s = std::ldexp(1.0, exponent);
        samples.push_back(rp::prepare({{0, 0, s, s, 0, s, s, 0}}));
    }
    rp::prepare_orientation_queries(samples);
    return samples;
}

template <rp::policy Policy>
void query_paths()
{
    auto const ordinary = rp::corpus("ordinary", 128, 31);
    for (unsigned k = 0; k < 3; ++k)
    {
        auto const kind = k == 0 ? rp::predicate::nni : k == 1 ? rp::predicate::nii : rp::predicate::iii;
        auto const general = rp::verify<Policy>(ordinary, 3, 0, kind);
        BOOST_TEST(general.orientation_filters[k] > 0);
    }
    auto const on_support = rp::verify<Policy>(ordinary, 3, 0, rp::predicate::nni_support);
    BOOST_TEST_EQ(on_support.orientation_exact[0], 3 * ordinary.size());
    auto const shared = rp::corpus("shared_boundary", 128, 37);
    auto const source_anchor = rp::verify<Policy>(shared, 3, 0, rp::predicate::nii_source);
    BOOST_TEST_EQ(source_anchor.orientation_exact[1], 3 * shared.size());
    // The experiment forces a cache entry even for an intersection equal to
    // a native endpoint. Do not report this as production's native fast path.
    BOOST_TEST_EQ(source_anchor.expansion_entries, shared.size());
    BOOST_TEST(source_anchor.fixed_entries > 0);
}

void relocation()
{
    auto const samples = rp::corpus("ordinary", 64, 7);
    rp::operation<rp::policy::deferred> op(1);
    std::vector<rp::point> points;
    for (auto const& s : samples) points.push_back(op.construct(s.input));
    // Both the point vector and cache have moved, but borrowed input addresses
    // and precision indices remain stable for the lifetime of this operation.
    for (std::size_t i = 0; i < points.size(); ++i)
    {
        auto copy = points[i];
        auto const a = op.emit(copy), b = op.emit(points[i]);
        BOOST_TEST_EQ(rp::bits(a[0]), rp::bits(b[0]));
        BOOST_TEST_EQ(rp::bits(a[1]), rp::bits(b[1]));
        BOOST_TEST_EQ(rp::bits(a[0]), rp::bits(samples[i].rounded[0]));
    }
}

void deferred_is_distinct()
{
    auto const samples = rp::corpus("ordinary", 1024, 0x554c504449464645);
    rp::operation<rp::policy::deferred> op(samples.size());
    std::size_t not_nearest = 0;
    for (auto const& s : samples)
    {
        auto p = op.construct(s.input);
        for (unsigned d = 0; d < 2; ++d)
        {
            if (rp::bits(p.approx[d]) != rp::bits(s.rounded[d])) ++not_nearest;
            BOOST_TEST(rp::one_ulp(p.approx[d], s.exact[d]));
        }
        auto const emitted = op.emit(p);
        BOOST_TEST_EQ(rp::bits(emitted[0]), rp::bits(s.rounded[0]));
        BOOST_TEST_EQ(rp::bits(emitted[1]), rp::bits(s.rounded[1]));
    }
    BOOST_TEST(not_nearest > 0);
    std::cout << not_nearest << " deferred coordinates were initially faithful but not nearest.\n";
}

int main()
{
    try
    {
        BOOST_TEST_EQ(std::fegetround(), FE_TONEAREST);
        std::size_t total = 0;
        for (char const* shape : {"ordinary", "near_parallel", "shared_boundary"})
        {
            auto const samples = rp::corpus(shape, 1024, 0x524f554e44494e47);
            check<rp::policy::eager>(samples);
            check<rp::policy::deferred>(samples);
            if (std::string(shape) == "near_parallel")
            {
                auto const c = rp::verify<rp::policy::deferred>(samples, 1, 0);
                BOOST_TEST(c.zero_denominators > 0);
                BOOST_TEST(c.estimate_refinements > 0);
            }
            total += samples.size();
        }
        auto const edge_cases = boundaries();
        check<rp::policy::eager>(edge_cases);
        check<rp::policy::deferred>(edge_cases);
        auto const counts = rp::verify<rp::policy::deferred>(edge_cases, 32, 1);
        BOOST_TEST(counts.fixed_entries >= 2);
        relocation();
        deferred_is_distinct();
        query_paths<rp::policy::eager>();
        query_paths<rp::policy::deferred>();
        std::cout << "Verified " << total + edge_cases.size()
                  << " intersections under both policies; rational oracle, zero allocations after reserve.\n";
    }
    catch (std::exception const& e)
    {
        count_allocations = false;
        std::cerr << e.what() << '\n';
        BOOST_ERROR("rounding-policy test threw");
    }
    return boost::report_errors();
}
