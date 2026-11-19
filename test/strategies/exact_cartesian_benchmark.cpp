// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Licensed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt

// Standalone benchmark, deliberately not registered as a unit test.
#include <boost/config.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/strategies/relate/cartesian.hpp>
#if defined(BG_BENCH_EXACT)
#include <boost/geometry/strategies/relate/cartesian_exact.hpp>
#if defined(BG_BENCH_EAGER) || defined(BG_BENCH_DEFERRED)
#include <boost/geometry/strategies/cartesian/exact/construction_rounding.hpp>
#endif
#endif

#if (defined(BG_BENCH_NEAREST) + defined(BG_BENCH_EAGER) + defined(BG_BENCH_DEFERRED)) > 1
#error Select only one exact output policy
#endif
#if !defined(BG_BENCH_EXACT) && (defined(BG_BENCH_NEAREST) || defined(BG_BENCH_EAGER) || defined(BG_BENCH_DEFERRED))
#error Exact output policies require BG_BENCH_EXACT
#endif

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <limits>
#include <new>
#include <random>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>
#if defined(_WIN32)
#include <process.h>
#else
#include <unistd.h>
#endif

#if defined(BG_BENCH_TRACK_ALLOCATIONS)
namespace allocation_counts
{
thread_local bool enabled = false;
thread_local std::uint64_t calls = 0, bytes = 0;
}
void* operator new(std::size_t n)
{
    if (void* p = std::malloc(n == 0 ? 1 : n))
    {
        if (allocation_counts::enabled)
        {
            ++allocation_counts::calls;
            allocation_counts::bytes += n;
        }
        return p;
    }
    throw std::bad_alloc();
}
void* operator new[](std::size_t n) { return ::operator new(n); }
void operator delete(void* p) noexcept { std::free(p); }
void operator delete[](void* p) noexcept { std::free(p); }
void operator delete(void* p, std::size_t) noexcept { std::free(p); }
void operator delete[](void* p, std::size_t) noexcept { std::free(p); }
#endif

namespace bg = boost::geometry;
#if defined(BG_BENCH_INT64)
using coordinate = std::int64_t;
char const* const coordinate_name = "int64";
#else
using coordinate = double;
char const* const coordinate_name = "double";
#endif
using point = bg::model::d2::point_xy<coordinate>;
using polygon = bg::model::polygon<point>;
using polygons = bg::model::multi_polygon<polygon>;
using line = bg::model::linestring<point>;
using lines = bg::model::multi_linestring<line>;
#if defined(BG_BENCH_EXACT)
#if defined(BG_BENCH_EAGER)
using strategy = bg::strategies::relate::cartesian_exact<coordinate,
    bg::strategies::relate::eager_nearest_output>;
#elif defined(BG_BENCH_DEFERRED)
using strategy = bg::strategies::relate::cartesian_exact<coordinate,
    bg::strategies::relate::deferred_nearest_output>;
#elif defined(BG_BENCH_NEAREST)
using strategy = bg::strategies::relate::cartesian_exact<coordinate,
    bg::strategies::relate::nearest_output>;
#else
using strategy = bg::strategies::relate::cartesian_exact<coordinate>;
#endif
#else
using strategy = bg::strategies::relate::cartesian<>;
#endif

struct input
{
    polygon a, b;
    lines la, lb;
    point query;
};

struct options
{
    std::string label = "unspecified", shape = "all", operation = "all";
    double seconds = 0.15;
    std::size_t vertices = 64, corpus = 8, repetitions = 3;
    bool digest_self_test = false, projected_self_test = false;
};

std::array<char const*, 5> const legacy_shapes{{"disjoint", "convex", "star", "holes", "shared"}};
std::array<char const*, 2> const projected_shapes{{"projected", "projected-local"}};
#include "projected_benchmark_fixture.hpp"

point rounded_point(double x, double y)
{
    // Identical exactly representable native inputs for both coordinate builds.
    return point(coordinate(std::llround(x)), coordinate(std::llround(y)));
}

template <typename Ring>
void radial(Ring& ring, std::size_t n, double x, double y, double radius,
            double phase, bool star = false, bool reverse = false)
{
    double const pi = std::acos(-1.0);
    for (std::size_t i = 0; i < n; ++i)
    {
        double const angle = phase + (reverse ? 1 : -1) * 2 * pi * double(i) / double(n);
        double const r = radius * (star && i % 2 ? 0.45 : 1.0);
        ring.push_back(rounded_point(x + r * std::cos(angle), y + r * std::sin(angle)));
    }
    ring.push_back(ring.front());
}

template <typename Ring>
void rectangle(Ring& ring, std::size_t n, double x, double y, double radius)
{
    n = (std::max)(std::size_t(1), n / 4);
    for (std::size_t side = 0; side < 4; ++side)
    {
        for (std::size_t i = 0; i < n; ++i)
        {
            double const t = double(i) / double(n);
            double const px[] = {-1, -1 + 2 * t, 1, 1 - 2 * t};
            double const py[] = {-1 + 2 * t, 1, 1 - 2 * t, -1};
            ring.push_back(rounded_point(x + radius * px[side], y + radius * py[side]));
        }
    }
    ring.push_back(ring.front());
}


input legacy_input(std::string const& shape, options const& opt, std::size_t c,
                   std::mt19937_64& random)
{
    input g;
    double const x = double(std::int64_t(random() % 4000001) - 2000000);
    double const y = double(std::int64_t(random() % 4000001) - 2000000);
    double const r = 1000000 + double(random() % 1000);
    double const dx = shape == "disjoint" ? 3.2 * r : shape == "shared" ? r : 0.61 * r;
    double const dy = shape == "star" ? 0.13 * r : shape == "shared" ? 0 : 0.17 * r;
    if (shape == "shared")
    {
        rectangle(g.a.outer(), opt.vertices, x, y, r);
        rectangle(g.b.outer(), opt.vertices, x + dx, y, r);
    }
    else
    {
        radial(g.a.outer(), opt.vertices, x, y, r, 0.137, shape == "star");
        radial(g.b.outer(), opt.vertices, x + dx, y + dy, r, 0.213, shape == "star");
    }
    if (shape == "holes")
    {
        g.a.inners().resize(1);
        g.b.inners().resize(1);
        radial(g.a.inners().front(), opt.vertices / 2, x, y, 0.3 * r, 0.089, false, true);
        radial(g.b.inners().front(), opt.vertices / 2, x + dx, y + dy, 0.3 * r, 0.271, false, true);
    }
    for (std::size_t i = 0; i < (std::max)(std::size_t(2), opt.vertices / 8); ++i)
    {
        double const offset = r * (double(i) + 0.17) / double((std::max)(std::size_t(2), opt.vertices / 8));
        g.la.push_back(line{rounded_point(x - r, y - r + offset),
                            rounded_point(x + r, y + 0.7 * r - 0.21 * offset)});
        g.lb.push_back(line{rounded_point(x - r + offset, y + r),
                            rounded_point(x + r - 0.15 * offset, y - r)});
    }
    g.query = rounded_point(x + (c % 3 == 0 ? 1.1 : 0.41) * r, y + 0.071 * r);
#if defined(BG_BENCH_FRACTIONAL_INPUTS)
    static_assert(std::is_same<coordinate, double>::value,
                  "The fractional control requires binary64 inputs");
    for_each_input_point(g, [](point& p)
    {
        bg::set<0>(p, bg::get<0>(p) / 8 + 1.0 / 16);
        bg::set<1>(p, bg::get<1>(p) / 16 - 3.0 / 32);
    });
#endif
    return g;
}

std::vector<input> make_inputs(std::string const& shape, options const& opt)
{
    bool const projected = std::find(projected_shapes.begin(), projected_shapes.end(), shape)
                        != projected_shapes.end();
    if (projected && !std::is_same<coordinate, double>::value)
    {
        throw std::invalid_argument("Projected workloads require double coordinates");
    }
    std::mt19937_64 random(0x455841435442454e);
    std::vector<input> result(opt.corpus);
    for (std::size_t c = 0; c < result.size(); ++c)
    {
        auto& g = result[c];
        g = projected ? projected_input(opt.vertices, c, random, shape == projected_shapes.back())
                      : legacy_input(shape, opt, c, random);
        std::string reason;
        if (!bg::is_valid(g.a, reason) || !bg::is_valid(g.b, reason))
        {
            throw std::runtime_error("Invalid benchmark input: " + shape + ": " + reason);
        }
    }
    return result;
}

std::uint64_t mix(std::uint64_t state, std::uint64_t value)
{
    return (state ^ value) * 1099511628211ULL;
}

// Native exported coordinates only, never operation-owned implicit points.
// Quantization at 2^-16 native units suppresses ordinary/nearest ULP noise on
// this corpus. Sorted undirected edges make rings independent of rotation and
// reversal; sorted hole/component hashes preserve their hierarchy. This is a
// spatial regression signal, not an exact geometric-equivalence oracle.
using quantized_point = std::array<std::int64_t, 2>;
quantized_point quantize(point const& p)
{
    quantized_point result;
    for (std::size_t d = 0; d < 2; ++d)
    {
        double const value = (d == 0 ? double(bg::get<0>(p)) : double(bg::get<1>(p))) * 65536;
        double const limit = std::ldexp(1.0, 63);
        if (!std::isfinite(value) || value <= -limit || value >= limit)
            throw std::domain_error("Benchmark output exceeds geometric digest grid");
        result[d] = std::llround(value);
    }
    return result;
}

std::uint64_t point_digest(quantized_point const& point)
{ return mix(mix(0, std::uint64_t(point[0])), std::uint64_t(point[1])); }

std::uint64_t sorted_digest(std::uint64_t tag, std::vector<std::uint64_t> values)
{
    std::sort(values.begin(), values.end());
    auto result = mix(tag, values.size());
    for (auto value : values) result = mix(result, value);
    return result;
}

template <typename Range>
std::uint64_t path_digest(Range const& path, bool closed)
{
    std::vector<quantized_point> points;
    for (auto const& p : path) points.push_back(quantize(p));
    while (closed && points.size() > 1 && points.front() == points.back()) points.pop_back();
    std::vector<std::uint64_t> edges;
    for (std::size_t i = 0; i < points.size(); ++i)
    {
        if (!closed && i + 1 == points.size()) break;
        auto first = points[i], second = points[(i + 1) % points.size()];
        if (second < first) std::swap(first, second);
        edges.push_back(mix(point_digest(first), point_digest(second)));
    }
    auto result = sorted_digest(closed ? 0x52494e47 : 0x4c494e45, std::move(edges));
    if (points.size() == 1) result = mix(result, point_digest(points.front()));
    return result;
}

std::uint64_t geometry_digest(polygons const& output)
{
    std::vector<std::uint64_t> components;
    for (auto const& polygon : output)
    {
        std::vector<std::uint64_t> holes;
        for (auto const& ring : polygon.inners()) holes.push_back(path_digest(ring, true));
        components.push_back(mix(path_digest(polygon.outer(), true), sorted_digest(0x484f4c45, std::move(holes))));
    }
    return sorted_digest(0x504f4c59, std::move(components));
}

std::uint64_t geometry_digest(lines const& output)
{
    std::vector<std::uint64_t> components;
    for (auto const& line : output) components.push_back(path_digest(line, false));
    return sorted_digest(0x4d554c54, std::move(components));
}

void digest_self_test()
{
    polygons first;
    first.resize(1);
    first.front().outer() = {{0, 0}, {0, 8}, {8, 8}, {8, 0}, {0, 0}};
    first.front().inners().resize(2);
    first.front().inners()[0] = {{1, 1}, {2, 1}, {2, 2}, {1, 2}, {1, 1}};
    first.front().inners()[1] = {{4, 4}, {5, 4}, {5, 5}, {4, 5}, {4, 4}};
    auto same = first;
    auto& outer = same.front().outer();
    outer.pop_back();
    std::rotate(outer.begin(), outer.begin() + 1, outer.end());
    outer.push_back(outer.front());
    std::reverse(outer.begin(), outer.end());
    std::swap(same.front().inners()[0], same.front().inners()[1]);
    if (geometry_digest(first) != geometry_digest(same)) throw std::logic_error("Polygon digest is order dependent");
    bg::set<0>(same.front().outer()[1], coordinate(7));
    if (geometry_digest(first) == geometry_digest(same)) throw std::logic_error("Polygon digest missed a coordinate change");
    lines line_output{{{0, 0}, {2, 3}, {4, 5}}};
    auto reverse = line_output;
    std::reverse(reverse.front().begin(), reverse.front().end());
    if (geometry_digest(line_output) != geometry_digest(reverse)) throw std::logic_error("Line digest is order dependent");
    bg::set<1>(reverse.front()[1], coordinate(4));
    if (geometry_digest(line_output) == geometry_digest(reverse)) throw std::logic_error("Line digest missed a coordinate change");
    std::cerr << "Geometric digest self-test passed.\n";
}

volatile std::uint64_t sink = 0;

enum class operation_id
{
    union_, intersection, difference, sym_difference,
    relate_aa, relate_la, relate_ll, covered_by, clip
};

std::array<char const*, 9> const operation_names{{"union", "intersection", "difference",
    "sym_difference", "relate_aa", "relate_la", "relate_ll", "covered_by", "clip"}};

operation_id operation_code(std::string const& name)
{
    auto const found = std::find(operation_names.begin(), operation_names.end(), name);
    if (found == operation_names.end()) { throw std::invalid_argument("Unknown operation: " + name); }
    return static_cast<operation_id>(found - operation_names.begin());
}

// A call boundary and observable result prevent loop-invariant hoisting.
BOOST_NOINLINE
std::uint64_t execute(operation_id operation, input const& g, polygons& out, lines& linear)
{
    strategy const s;
    std::uint64_t value = 0;
    if (operation == operation_id::relate_aa || operation == operation_id::relate_la
        || operation == operation_id::relate_ll)
    {
        auto const matrix = operation == operation_id::relate_aa ? bg::relation(g.a, g.b, s)
            : operation == operation_id::relate_la ? bg::relation(g.la, g.a, s) : bg::relation(g.la, g.lb, s);
        for (char c : matrix.str()) { value = mix(value, std::uint64_t(c)); }
    }
    else if (operation == operation_id::covered_by)
    {
        value = bg::covered_by(g.query, g.a, s);
    }
    else if (operation == operation_id::clip)
    {
        linear.clear();
        bg::intersection(g.la, g.a, linear, s);
        value = mix(linear.size(), bg::num_points(linear));
    }
    else
    {
        out.clear();
        if (operation == operation_id::union_) { bg::union_(g.a, g.b, out, s); }
        else if (operation == operation_id::intersection) { bg::intersection(g.a, g.b, out, s); }
        else if (operation == operation_id::difference) { bg::difference(g.a, g.b, out, s); }
        else { bg::sym_difference(g.a, g.b, out, s); }
        value = mix(out.size(), bg::num_points(out));
        for (auto const& p : out) { value = mix(value, p.inners().size()); }
    }
    sink = value;
    return value;
}

void projected_self_test(options const& opt)
{
    auto const projected = make_inputs(projected_shapes.front(), opt);
    auto local = make_inputs(projected_shapes.back(), opt);
    auto const rebase = [](point& p)
    {
        bg::set<0>(p, bg::get<0>(p) - easting);
        bg::set<1>(p, bg::get<1>(p) - northing);
    };
    for (std::size_t i = 0; i < local.size(); ++i)
    {
        bool west = false, east = false, south = false, north = false;
        for_each_input_point(local[i], [&](point const& p)
        {
            west |= bg::get<0>(p) < 0;
            east |= bg::get<0>(p) > 0;
            south |= bg::get<1>(p) < 0;
            north |= bg::get<1>(p) > 0;
        });
        if (!(west && east && south && north))
            throw std::logic_error("Local fixture must straddle zero in both axes");

        for (auto const& op : operation_names)
        {
            auto const id = operation_code(op);
            polygons shifted_out, local_out;
            lines shifted_lines, local_lines;
            auto const shifted_result = execute(id, projected[i], shifted_out, shifted_lines);
            auto const local_result = execute(id, local[i], local_out, local_lines);
            bool match = shifted_result == local_result;
            if (id <= operation_id::sym_difference)
            {
                bg::for_each_point(shifted_out, rebase);
                match = match && !local_out.empty()
                    && geometry_digest(shifted_out) == geometry_digest(local_out);
            }
            else if (id == operation_id::clip)
            {
                bg::for_each_point(shifted_lines, rebase);
                match = match && !local_lines.empty()
                    && geometry_digest(shifted_lines) == geometry_digest(local_lines);
            }
            else if (id == operation_id::covered_by)
            {
                match = match && local_result == (i % 3 != 0);
            }
            if (!match)
                throw std::logic_error("Projected fixture mismatch: " + std::string(op)
                                       + " case " + std::to_string(i));
        }
    }
    std::cerr << "Projected fixture self-test passed.\n";
}

void measure(std::string const& shape, std::string const& operation, options const& opt)
{
    auto const corpus = make_inputs(shape, opt);
    auto const id = operation_code(operation);
    polygons out;
    lines linear;
    // Only native output containers survive public calls. Each exact operation
    // owns and destroys its own Root records/caches before returning here.
    std::uint64_t digest = 0, spatial_digest = 0;
    for (auto const& g : corpus)
    {
        auto const observed = execute(id, g, out, linear);
        digest = mix(digest, observed);
        auto const spatial = id == operation_id::clip ? geometry_digest(linear)
            : id <= operation_id::sym_difference ? geometry_digest(out) : observed;
        spatial_digest = mix(spatial_digest, spatial);
    }
    for (std::size_t repetition = 0; repetition < opt.repetitions; ++repetition)
    {
        using clock = std::chrono::steady_clock;
        auto const warm_start = clock::now();
        execute(id, corpus.front(), out, linear);
        double const warm_seconds = std::chrono::duration<double>(clock::now() - warm_start).count();
        auto const estimate = std::size_t((std::max)(1.0, (std::min)(64.0, 0.003 / warm_seconds)));
        // Every timed batch visits the complete corpus, including slow cases.
        auto const batch = ((estimate + corpus.size() - 1) / corpus.size()) * corpus.size();
        std::uint64_t calls = 0, checksum = 0;
#if defined(BG_BENCH_TRACK_ALLOCATIONS)
        allocation_counts::calls = allocation_counts::bytes = 0;
        allocation_counts::enabled = true;
#endif
        auto const cpu_start = std::clock();
        auto const start = clock::now();
        auto end = start;
        do
        {
            for (std::size_t i = 0; i < batch; ++i)
            {
                checksum += execute(id, corpus[calls % corpus.size()], out, linear);
                ++calls;
            }
            end = clock::now();
        } while (std::chrono::duration<double>(end - start).count() < opt.seconds);
        auto const cpu_end = std::clock();
#if defined(BG_BENCH_TRACK_ALLOCATIONS)
        allocation_counts::enabled = false;
        auto const allocations = allocation_counts::calls, bytes = allocation_counts::bytes;
#else
        std::uint64_t const allocations = 0, bytes = 0;
#endif
        double const seconds = std::chrono::duration<double>(end - start).count();
        std::cout << opt.label << ',' << coordinate_name << ',' << shape << ',' << opt.vertices
                  << ',' << operation << ',' << repetition << ',' << calls << ',' << seconds
                  << ',' << double(cpu_end - cpu_start) / CLOCKS_PER_SEC
                  << ',' << seconds * 1e9 / double(calls) << ',' << allocations << ',' << bytes
                  << ',' << checksum << ',' << digest << ',' << spatial_digest << '\n' << std::flush;
    }
}

int main(int argc, char** argv)
{
    try
    {
        options opt;
        for (int i = 1; i < argc; i += 2)
        {
            if (i + 1 == argc) { throw std::invalid_argument("Missing option value"); }
            std::string const name(argv[i]), value(argv[i + 1]);
            if (name == "--label") { opt.label = value; }
            else if (name == "--workload") { opt.shape = value; }
            else if (name == "--operation") { opt.operation = value; }
            else if (name == "--seconds") { opt.seconds = std::stod(value); }
            else if (name == "--vertices") { opt.vertices = std::stoul(value); }
            else if (name == "--corpus") { opt.corpus = std::stoul(value); }
            else if (name == "--repeat") { opt.repetitions = std::stoul(value); }
            else if (name == "--digest-self-test" && (value == "0" || value == "1")) { opt.digest_self_test = value == "1"; }
            else if (name == "--projected-self-test" && (value == "0" || value == "1")) { opt.projected_self_test = value == "1"; }
            else { throw std::invalid_argument("Unknown option: " + name); }
        }
        if (!std::isfinite(opt.seconds) || opt.seconds <= 0 || opt.vertices < 8
            || opt.corpus == 0 || opt.repetitions == 0)
        {
            throw std::invalid_argument("Invalid benchmark bounds");
        }
        if (opt.digest_self_test || opt.projected_self_test)
        {
            if (opt.digest_self_test) { digest_self_test(); }
            if (opt.projected_self_test) { projected_self_test(opt); }
            return 0;
        }
        bool const known_shape = std::find(legacy_shapes.begin(), legacy_shapes.end(), opt.shape) != legacy_shapes.end()
                              || std::find(projected_shapes.begin(), projected_shapes.end(), opt.shape) != projected_shapes.end();
        if (!known_shape && opt.shape != "all" && opt.shape != "projected-pair")
            throw std::invalid_argument("Unknown workload: " + opt.shape);
        if (opt.operation != "all") { (void) operation_code(opt.operation); }
#if defined(_WIN32)
        std::cerr << "pid=" << _getpid() << '\n';
#else
        std::cerr << "pid=" << getpid() << '\n';
#endif
        std::cout << std::setprecision(12)
            << "label,coordinate,workload,vertices,operation,repeat,calls,seconds,cpu_seconds,ns_per_call,allocations,allocated_bytes,checksum,corpus_digest,geometry_digest\n";
        auto const run = [&](std::string const& shape)
        {
            for (auto const& op : operation_names)
            {
                if (opt.operation == "all" || opt.operation == op) { measure(shape, op, opt); }
            }
        };
        if (opt.shape == "all") { for (auto shape : legacy_shapes) { run(shape); } }
        else if (opt.shape == "projected-pair") { for (auto shape : projected_shapes) { run(shape); } }
        else { run(opt.shape); }
        return 0;
    }
    catch (std::exception const& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
