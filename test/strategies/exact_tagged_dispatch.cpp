// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.

#include <boost/core/lightweight_test.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/strategies/relate/cartesian_exact.hpp>

#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

namespace bg = boost::geometry;
namespace ec = bg::strategies::detail::exact_cartesian;
using point = bg::model::d2::point_xy<double>;
using polygon = bg::model::polygon<point>;
using multipolygon = bg::model::multi_polygon<polygon>;
using line = bg::model::linestring<point>;
using policy = bg::strategies::relate::approximate_output;
using strategy = bg::strategies::relate::cartesian_exact<double, policy>;

enum class operation { intersection, difference, unite, symmetric, relation };
std::array<operation, 5> const operations{{operation::intersection, operation::difference,
    operation::unite, operation::symmetric, operation::relation}};
char const* name(operation op)
{
    return std::array<char const*, 5>{{"intersection", "difference", "union", "symmetric", "relation"}}
        [static_cast<std::size_t>(op)];
}

struct fixture
{
    char const* name;
    polygon a, b;
};
using fixture_set = std::array<fixture, 6>;

polygon rectangle(double x0, double y0, double x1, double y1)
{
    polygon p;
    p.outer() = {{x0, y0}, {x0, y1}, {x1, y1}, {x1, y0}, {x0, y0}};
    return p;
}

fixture cross_strips(char const* label, double origin, double width)
{
    return {label,
        rectangle(origin - 16, origin, origin + 16, origin + width),
        rectangle(origin, origin - 16, origin + width, origin + 16)};
}

fixture near_parallel()
{
    double const d = 0x1p-46;
    auto const strip = [](double left, double right, double width)
    {
        polygon p;
        p.outer() = {{-16, left}, {-16, left + width}, {16, right + width},
                     {16, right}, {-16, left}};
        return p;
    };
    // In coordinates (x, u=y-x-64), the strips are 0<=u<=d and
    // d*(x-4)<=u<=d*(x-2). All four crossings have different x and y,
    // but the turn determinant is tiny compared with its products.
    return {"near-parallel", strip(48, 80, d), strip(48 - 20 * d, 80 + 12 * d, 2 * d)};
}

polygon radial(double x, double y, double offset, bool star)
{
    polygon result;
    for (int i = 0; i != 48; ++i)
    {
        double const angle = offset - double(i) * (2 * std::acos(-1.0) / 48);
        double const radius = star && i % 2 ? 512 : 1024;
        result.outer().emplace_back(x + std::round(radius * std::cos(angle)),
                                    y + std::round(radius * std::sin(angle)));
    }
    result.outer().push_back(result.outer().front());
    return result;
}

fixture_set fixtures()
{
    return {{cross_strips("one-ulp-cross", 0x1p50, 0.25),
        cross_strips("ordinary-cross", 0, 4),
        {"convex", radial(0, 0, 0.137, false), radial(615, 184, 0.213, false)},
        {"star", radial(0, 0, 0.137, true), radial(615, 184, 0.213, true)},
        {"shared-edge", rectangle(0, 0, 8, 8), rectangle(8, 0, 16, 8)},
        near_parallel()}};
}

template <typename Polygon, typename Multi, typename Strategy>
void overlay(operation op, Polygon const& a, Polygon const& b, Multi& out, Strategy const& s)
{
    switch (op)
    {
    case operation::intersection: bg::intersection(a, b, out, s); break;
    case operation::difference: bg::difference(a, b, out, s); break;
    case operation::unite: bg::union_(a, b, out, s); break;
    case operation::symmetric: bg::sym_difference(a, b, out, s); break;
    default: throw std::logic_error("not an overlay operation");
    }
}

BOOST_NOINLINE
std::size_t execute(operation op, fixture const& f, multipolygon& out)
{
    if (op == operation::relation)
    {
        auto const m = bg::relation(f.a, f.b, strategy()).str();
        std::size_t hash = 0;
        for (auto c : m) hash = (hash * 33) ^ std::size_t(c);
        return hash;
    }
    out.clear();
    overlay(op, f.a, f.b, out, strategy());
    return bg::num_points(out) + 33 * out.size();
}

// Keep equality checks on full, round-trip native outputs outside the timer.
std::string snapshot(operation op, fixture const& f)
{
    if (op == operation::relation) return bg::relation(f.a, f.b, strategy()).str();
    multipolygon out;
    overlay(op, f.a, f.b, out, strategy());
    std::ostringstream stream;
    stream << std::setprecision(17) << bg::wkt(out);
    return stream.str();
}

void benchmark(fixture_set const& inputs)
{
    using clock = std::chrono::steady_clock;
    std::cout << "fixture\toperation\tns_per_call\tchecksum\toutput\n";
    for (auto const& f : inputs)
        for (auto op : operations)
        {
            auto const output = snapshot(op, f);
            multipolygon out;
            std::size_t checksum = execute(op, f, out);
            auto const start = clock::now();
            constexpr int repeats = 1000;
            for (int i = 0; i != repeats; ++i) checksum += execute(op, f, out);
            auto const ns = std::chrono::duration<double, std::nano>(clock::now() - start).count();
            std::cout << f.name << '\t' << name(op) << '\t' << ns / repeats << '\t' << checksum
                      << '\t' << output << '\n';
        }
}

#ifdef BOOST_GEOMETRY_EXACT_FILTER_STATISTICS
void traffic(fixture_set const& inputs)
{
    using root = ec::borrowed_root<double, policy, polygon, polygon>;
    std::cout << "borrowed_record_bytes\t" << sizeof(ec::borrowed_record<root>) << '\n';
    std::cout << "fixture\toperation\ttagged_calls\n";
    for (auto const& f : inputs)
        for (auto op : operations)
        {
            ec::filtering_statistics() = {};
            multipolygon out;
            execute(op, f, out);
            std::cout << f.name << '\t' << name(op) << '\t'
                      << ec::filtering_statistics().tagged_orientation_calls << '\n';
        }
}
#endif

#ifndef BG_TAGGED_TIMINGS_ONLY
using rational = ec::exact_number;
using rational_point = bg::model::d2::point_xy<rational>;
using rational_polygon = bg::model::polygon<rational_point>;
using rational_multi = bg::model::multi_polygon<rational_polygon>;

rational_polygon lift(polygon const& p)
{
    rational_polygon r;
    for (auto const& v : p.outer()) r.outer().emplace_back(rational(v.x()), rational(v.y()));
    return r;
}

void oracle(fixture_set const& inputs)
{
    bg::strategies::relate::cartesian<rational> const reference;
    for (auto const& f : inputs)
    {
        rational_polygon const a = lift(f.a), b = lift(f.b);
        BOOST_TEST(bg::is_valid(a, reference));
        BOOST_TEST(bg::is_valid(b, reference));
        BOOST_TEST_EQ(bg::relation(f.a, f.b, strategy()).str(), bg::relation(a, b, reference).str());
        for (auto op : operations)
        {
            if (op == operation::relation) continue;
            rational_multi actual, expected;
            overlay(op, f.a, f.b, actual, strategy());
            overlay(op, a, b, expected, reference);
            BOOST_TEST(bg::is_valid(actual, reference));
            BOOST_TEST_EQ(bg::area(actual, reference), bg::area(expected, reference));
            if (expected.empty()) { BOOST_TEST(actual.empty()); }
            else BOOST_TEST_EQ(bg::relation(actual, expected, reference).str(), "2FFF1FFF2");
            // Independent area/count oracle for the thin cross family.
            if (&f == &inputs.front())
            {
                std::array<rational, 4> const areas{{rational(1, 16), rational(127, 16),
                    rational(255, 16), rational(127, 8)}};
                std::array<std::size_t, 4> const counts{{1, 2, 1, 4}};
                auto const i = static_cast<std::size_t>(op);
                BOOST_TEST_EQ(bg::area(actual, reference), areas[i]);
                BOOST_TEST_EQ(actual.size(), counts[i]);
            }
            if (&f == &inputs.back())
            {
                // Horizontal slices in (x,u) have length 2 for 0<=u<=d.
                // This change of coordinates has determinant 1.
                std::array<int, 4> const factors{{2, 30, 94, 92}};
                rational const area = factors[static_cast<std::size_t>(op)] * rational(0x1p-46);
                BOOST_TEST_EQ(bg::area(actual, reference), area);
            }
        }
    }
}

#ifdef BOOST_GEOMETRY_EXACT_FILTER_STATISTICS
void route_controls()
{
    using root = ec::borrowed_root<double, policy, line, line>;
    using borrowed = ec::borrowed_point<root>;
    using record = ec::borrowed_record<root>;
    struct crossing
    {
        line horizontal, vertical;
        ec::borrowed_stage<root> stage;
        std::array<root::operation_type, 2> operations;
        record turn;
        borrowed p;

        crossing(ec::borrowed_context<root>& context, double x, double y)
            : horizontal{{x - 16, y}, {x + 16, y}}, vertical{{x, y - 16}, {x, y + 16}},
              stage(context, horizontal, vertical), p(turn)
        {
            operations[0].seg_id = bg::segment_identifier(0, -1, -1, 0);
            operations[1].seg_id = bg::segment_identifier(1, -1, -1, 0);
            turn.stage = &stage;
            turn.operations = {{&operations[0], &operations[1]}};
            turn.numeric.init(p.producers());
        }
    };
    ec::borrowed_context<root> context;
    double const base = 0x1p50;
    crossing a(context, base, base), b(context, base + 0.25, base), c(context, base, base + 0.25);
    using cache = ec::point_cache_detail::native_cache<double>;
    auto const ready = [](crossing const& x) { return cache::has_exact(x.turn.numeric.cache_handle()); };
    // Prime only the existing refinement cache, not the fixed values. The
    // additional tagged sign computation must not add or populate entries.
    for (auto const* x : {&a, &b, &c})
    {
        x->p.refine(0);
        x->p.refine(1);
        BOOST_TEST(!ready(*x));
    }
    auto const entries = cache::size();
    ec::filtering_statistics() = {};
    auto& tagged = ec::filtering_statistics().tagged_orientation_calls;
    BOOST_TEST_EQ(ec::root_orientation<root>(a.p, b.p, c.p), 1);
#if BOOST_GEOMETRY_EXACT_TAGGED_ORIENTATION
    BOOST_TEST_EQ(tagged, 1u);
    BOOST_TEST(!ready(a) && !ready(b) && !ready(c));
#else
    BOOST_TEST_EQ(tagged, 0u);
    BOOST_TEST(ready(a) && ready(b) && ready(c));
#endif
    BOOST_TEST_EQ(cache::size(), entries);
    auto const cold_calls = tagged;

    BOOST_TEST_EQ(ec::root_orientation<root>(a.p, a.p, c.p), 0);
    record alias = a.turn;
    BOOST_TEST_EQ(ec::root_orientation<root>(a.p, borrowed(alias), c.p), 0);
    BOOST_TEST_EQ(tagged, cold_calls);
    a.p.homogeneous();
    BOOST_TEST(ready(a));
    BOOST_TEST_EQ(ec::root_orientation<root>(a.p, c.p, b.p), -1);
    BOOST_TEST_EQ(tagged, cold_calls);

    // Cold constructed/constructed/native and query mixes do not enter
    // the new route even if the initial sign is ambiguous.
    crossing d(context, base, base), e(context, base + 0.25, base);
    point const native(base, base + 0.25);
    BOOST_TEST_EQ(ec::root_orientation<root>(d.p, e.p, native), 1);
    auto const middle = ec::query_point<root>::midpoint(d.p, e.p);
    BOOST_TEST_EQ(ec::root_orientation<root>(d.p, middle, e.p), 0);
    BOOST_TEST_EQ(tagged, cold_calls);
}

void public_routes(fixture_set const& inputs)
{
    auto& tagged = ec::filtering_statistics().tagged_orientation_calls;
    multipolygon out;
    tagged = 0;
    execute(operation::intersection, inputs.back(), out);
#if BOOST_GEOMETRY_EXACT_TAGGED_ORIENTATION
    BOOST_TEST(tagged > 0u);
#else
    BOOST_TEST_EQ(tagged, 0u);
#endif
}
#endif
#endif

int main(int argc, char** argv)
{
    auto const inputs = fixtures();
    if (argc == 2 && std::string(argv[1]) == "--benchmark") benchmark(inputs);
    else if (argc == 2 && std::string(argv[1]) == "--check")
    {
        for (auto const& f : inputs)
            for (auto op : operations)
                std::cout << f.name << '\t' << name(op) << '\t' << snapshot(op, f) << '\n';
    }
#ifdef BOOST_GEOMETRY_EXACT_FILTER_STATISTICS
    else if (argc == 2 && std::string(argv[1]) == "--traffic") traffic(inputs);
#endif
#ifndef BG_TAGGED_TIMINGS_ONLY
    else
    {
#ifdef BOOST_GEOMETRY_EXACT_FILTER_STATISTICS
        route_controls();
        public_routes(inputs);
#endif
        oracle(inputs);
    }
#endif
    return boost::report_errors();
}
