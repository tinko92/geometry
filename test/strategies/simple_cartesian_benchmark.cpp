// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Licensed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt

// Standalone benchmark, deliberately not registered as a unit test.
#include <boost/config.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/strategies/relate/cartesian.hpp>
#if defined(BG_BENCH_SIMPLE)
#include <boost/geometry/strategies/relate/simple_cartesian.hpp>
#endif
#if defined(BG_BENCH_SIMPLE) && defined(BG_BENCH_EXACT)
#error Select one benchmark strategy
#endif
#if defined(BG_BENCH_EAGER) && defined(BG_BENCH_DEFERRED)
#error Select one exact output policy
#endif
#if (defined(BG_BENCH_EAGER) || defined(BG_BENCH_DEFERRED)) && !defined(BG_BENCH_EXACT)
#error Exact output policies require BG_BENCH_EXACT
#endif
#if defined(BG_BENCH_EXACT)
#include <boost/geometry/strategies/relate/cartesian_exact.hpp>
#if defined(BG_BENCH_EAGER) || defined(BG_BENCH_DEFERRED)
#include <boost/geometry/strategies/cartesian/exact/construction_rounding.hpp>
#endif
#endif

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#if defined(_WIN32)
#include <process.h>
#else
#include <unistd.h>
#endif

namespace bg = boost::geometry;
using coordinate = double;
char const* const coordinate_name = "double";
using point = bg::model::d2::point_xy<coordinate>;
using polygon = bg::model::polygon<point>;
using polygons = bg::model::multi_polygon<polygon>;
using line = bg::model::linestring<point>;
using lines = bg::model::multi_linestring<line>;
#if defined(BG_BENCH_SIMPLE)
using strategy = bg::strategies::relate::simple_cartesian;
#elif defined(BG_BENCH_EXACT)
#if defined(BG_BENCH_ULPS)
using strategy = bg::strategies::relate::cartesian_exact<coordinate,
    bg::strategies::relate::bounded_ulp_output<BG_BENCH_ULPS>>;
#elif defined(BG_BENCH_EAGER)
using strategy = bg::strategies::relate::cartesian_exact<coordinate,
    bg::strategies::relate::eager_nearest_output>;
#elif defined(BG_BENCH_DEFERRED)
using strategy = bg::strategies::relate::cartesian_exact<coordinate,
    bg::strategies::relate::deferred_nearest_output>;
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
    double seconds = 0.1;
    std::size_t repetitions = 3;
    bool digest_self_test = false, check_only = false;
};

#include "projected_benchmark_fixture.hpp"

struct wkt_input
{
    char const* category;
    char const* a;
    char const* b;
    char const* la;
    char const* lb;
    char const* query;
};

// Three valid native inputs per category. Degeneracy is between geometries,
// not invalid rings. Decimal coordinates are literal thousandths throughout.
std::array<wkt_input, 9> const input_wkts{{
    {"degenerate",
     "POLYGON((0 0,0 6,4 6,4 4,6 4,6 0,0 0))",
     "POLYGON((6 0,6 4,8 4,8 6,12 6,12 0,6 0))",
     "MULTILINESTRING((-1 2,6 2,13 2),(6 -1,6 4,8 4))",
     "MULTILINESTRING((6 -1,6 2,6 7),(-1 4,6 4,13 4))", "POINT(6 2)"},
    {"degenerate",
     "POLYGON((0 0,0 4,2 6,6 6,8 4,8 0,0 0))",
     "POLYGON((6 6,8 10,12 10,14 8,14 6,10 6,6 6))",
     "MULTILINESTRING((-1 6,6 6,15 6),(2 -1,2 6,8 10))",
     "MULTILINESTRING((6 0,6 6,6 12),(0 0,6 6,14 8))", "POINT(6 6)"},
    {"degenerate",
     "POLYGON((0 0,0 6,0 10,6 10,10 10,10 4,10 0,4 0,0 0),(3 3,7 3,7 7,3 7,3 3))",
     "POLYGON((3 3,3 5,4 7,6 7,7 5,7 3,3 3))",
     "MULTILINESTRING((-1 5,5 5,11 5),(3 1,3 7,8 7))",
     "MULTILINESTRING((1 3,7 3,9 3),(5 1,5 5,5 9))", "POINT(4 3)"},
    {"general",
     "POLYGON((0 0,-1 4,2 8,7 7,9 3,6 -1,0 0))",
     "POLYGON((3 -2,1 3,4 9,10 8,12 4,9 -2,3 -2))",
     "MULTILINESTRING((-2 2,4 4,13 3),(0 -3,4 4,8 10))",
     "MULTILINESTRING((-1 7,5 2,13 6),(2 -2,6 5,11 9))", "POINT(3 3)"},
    {"general",
     "POLYGON((0 0,0 10,3 12,9 12,12 9,12 3,9 0,0 0),(3 3,7 3,9 5,7 8,4 8,3 5,3 3))",
     "POLYGON((5 -2,5 4,2 4,2 10,10 10,10 6,14 6,14 -2,5 -2))",
     "MULTILINESTRING((-1 5,6 5,15 7),(2 -1,6 6,10 13))",
     "MULTILINESTRING((-2 9,7 4,15 3),(1 1,8 9,13 11))", "POINT(5 5)"},
    {"general",
     "POLYGON((-12 -4,-14 0,-12 5,-7 5,-5 1,-7 -4,-12 -4))",
     "POLYGON((4 -3,2 1,5 6,10 5,12 0,9 -4,4 -3))",
     "MULTILINESTRING((-16 1,-8 1,14 2),(-10 -6,-9 2,-6 7))",
     "MULTILINESTRING((-15 -2,-6 4,13 0),(-11 7,-8 0,-3 -5))", "POINT(0 0)"},
    {"decimal",
     "POLYGON((0.125 0.375,-1.205 4.125,1.875 8.333,7.125 7.875,9.005 3.125,6.625 -1.375,0.125 0.375))",
     "POLYGON((3.015 -2.125,1.375 3.005,4.625 9.125,10.375 8.005,12.125 4.625,9.875 -2.125,3.015 -2.125))",
     "MULTILINESTRING((-2.005 2.125,4.375 4.005,13.125 3.625),(0.125 -3.005,4.005 4.375,8.625 10.005))",
     "MULTILINESTRING((-1.125 7.005,5.375 2.625,13.005 6.125),(2.625 -2.005,6.375 5.125,11.005 9.625))", "POINT(3.125 3.375)"},
    {"decimal",
     "POLYGON((0.111 0.222,0.111 10.222,3.111 12.222,9.111 12.222,12.111 9.222,12.111 3.222,9.111 0.222,0.111 0.222),(3.111 3.222,7.111 3.222,9.111 5.222,7.111 8.222,4.111 8.222,3.111 5.222,3.111 3.222))",
     "POLYGON((5.333 -1.555,5.333 4.445,2.333 4.445,2.333 10.445,10.333 10.445,10.333 6.445,14.333 6.445,14.333 -1.555,5.333 -1.555))",
     "MULTILINESTRING((-1.125 5.375,6.625 5.005,15.125 7.375),(2.005 -1.125,6.375 6.625,10.005 13.125))",
     "MULTILINESTRING((-2.625 9.125,7.005 4.375,15.625 3.125),(1.005 1.375,8.125 9.625,13.375 11.005))", "POINT(5.111 5.222)"},
    {"decimal",
     "POLYGON((0.111 0.222,0.111 6.222,4.111 6.222,4.111 4.222,6.111 4.222,6.111 0.222,0.111 0.222))",
     "POLYGON((6.111 0.222,6.111 4.222,8.111 4.222,8.111 6.222,12.111 6.222,12.111 0.222,6.111 0.222))",
     "MULTILINESTRING((-1.111 2.222,6.111 2.222,13.111 2.222),(6.111 -1.222,6.111 4.222,8.111 4.222))",
     "MULTILINESTRING((6.111 -1.222,6.111 2.222,6.111 7.222),(-1.111 4.222,6.111 4.222,13.111 4.222))", "POINT(6.111 2.222)"}
}};

std::vector<input> make_inputs(std::string const& category)
{
    std::vector<input> result;
#if defined(BG_BENCH_EXACT) && defined(BG_BENCH_ULPS)
    if (category == "denominator")
    {
        // Two clockwise parallelograms with primary edges crossing at t=1/2.
        // Their exact determinant is -1 while the naive head is -1 or zero.
        auto strip = [](point const& p, point const& q)
        {
            polygon a;
            a.outer() = {p, {p.x(), p.y() + 4}, {q.x(), q.y() + 4}, q, p};
            return a;
        };
        for (int exponent : {26, 30, 36})
        {
            double const n = std::ldexp(1.0, exponent);
            point const a(n, n), b(2*n+1, 2*n), c(n+0.5, n+0.5), d(2*n+0.5, 2*n-0.5);
            input g;
            g.a = strip(a, b); g.b = strip(c, d);
            g.la.push_back({a, b}); g.lb.push_back({c, d});
            g.query = {n+2, n+3};
            result.push_back(std::move(g));
        }
        return result;
    }
#endif
    if (category == "projected" || category == "projected-local")
    {
        std::mt19937_64 random(0x455841435442454e);
        for (std::size_t i = 0; i < 3; ++i)
        {
            auto g = projected_input(64, i, random, category == "projected-local");
            if (!bg::is_valid(g.a) || !bg::is_valid(g.b))
                throw std::runtime_error("Invalid projected benchmark polygon");
            result.push_back(std::move(g));
        }
        return result;
    }
    for (auto const& text : input_wkts)
    {
        if (category != text.category) continue;
        input g;
        bg::read_wkt(text.a, g.a); bg::read_wkt(text.b, g.b);
        bg::read_wkt(text.la, g.la); bg::read_wkt(text.lb, g.lb);
        bg::read_wkt(text.query, g.query);
        bg::correct(g.a); bg::correct(g.b);
        std::string reason;
        if (!bg::is_valid(g.a, reason) || !bg::is_valid(g.b, reason))
            throw std::runtime_error("Invalid WKT corpus: " + category + ": " + reason);
        result.push_back(std::move(g));
    }
    if (result.size() != 3) throw std::invalid_argument("Unknown category: " + category);
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

void measure(std::string const& shape, std::string const& operation, options const& opt,
             std::vector<input> const& corpus)
{
    auto const id = operation_code(operation);
    polygons out;
    lines linear;
    // Only native output containers survive public calls. Each exact operation
    // owns and destroys its own Root records/caches before returning here.
    std::uint64_t digest = 0, spatial_digest = 0, components = 0, holes = 0, points = 0;
    for (auto const& g : corpus)
    {
        auto const observed = execute(id, g, out, linear);
        digest = mix(digest, observed);
        auto const spatial = id == operation_id::clip ? geometry_digest(linear)
            : id <= operation_id::sym_difference ? geometry_digest(out) : observed;
        spatial_digest = mix(spatial_digest, spatial);
        if (id == operation_id::clip) { components += linear.size(); points += bg::num_points(linear); }
        else if (id <= operation_id::sym_difference)
        {
            components += out.size(); points += bg::num_points(out);
            for (auto const& p : out) holes += p.inners().size();
        }
    }
    if (opt.check_only)
    {
        std::cout << opt.label << ',' << coordinate_name << ',' << shape
                  << ',' << operation << ",0,0,0,0,0,0," << digest << ',' << spatial_digest
                  << ',' << corpus.size() << ',' << components << ',' << holes << ',' << points << '\n';
        return;
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
        double const seconds = std::chrono::duration<double>(end - start).count();
        std::cout << opt.label << ',' << coordinate_name << ',' << shape
                  << ',' << operation << ',' << repetition << ',' << calls << ',' << seconds
                  << ',' << double(cpu_end - cpu_start) / CLOCKS_PER_SEC
                  << ',' << seconds * 1e9 / double(calls)
                  << ',' << checksum << ',' << digest << ',' << spatial_digest
                  << ',' << corpus.size() << ',' << components << ',' << holes << ',' << points << '\n' << std::flush;
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
            else if (name == "--repeat") { opt.repetitions = std::stoul(value); }
            else if (name == "--digest-self-test" && (value == "0" || value == "1")) { opt.digest_self_test = value == "1"; }
            else if (name == "--check" && (value == "0" || value == "1")) { opt.check_only = value == "1"; }
            else { throw std::invalid_argument("Unknown option: " + name); }
        }
        if (!std::isfinite(opt.seconds) || opt.seconds <= 0 || opt.repetitions == 0)
        {
            throw std::invalid_argument("Invalid benchmark bounds");
        }
        if (opt.digest_self_test) { digest_self_test(); return 0; }
#if defined(_WIN32)
        std::cerr << "pid=" << _getpid() << '\n';
#else
        std::cerr << "pid=" << getpid() << '\n';
#endif
        std::cout << std::setprecision(12)
            << "label,coordinate,workload,operation,repeat,calls,seconds,cpu_seconds,ns_per_call,checksum,corpus_digest,geometry_digest,corpus_size,output_components,output_holes,output_points\n";
        std::array<char const*, 5> const shapes{{"degenerate", "general", "decimal", "projected", "projected-local"}};
#if defined(BG_BENCH_EXACT) && defined(BG_BENCH_ULPS)
        if (opt.shape == "denominator")
        {
            if (opt.operation != "all") (void)operation_code(opt.operation);
            auto const corpus = make_inputs(opt.shape);
            for (auto const& op : operation_names)
                if (opt.operation == "all" || opt.operation == op) measure(opt.shape, op, opt, corpus);
            return 0;
        }
#endif
        if (opt.shape != "all" && std::find(shapes.begin(), shapes.end(), opt.shape) == shapes.end())
        {
            throw std::invalid_argument("Unknown workload: " + opt.shape);
        }
        if (opt.operation != "all") { (void) operation_code(opt.operation); }
        for (auto const& shape : shapes)
        {
            if (opt.shape != "all" && opt.shape != shape) { continue; }
            auto const corpus = make_inputs(shape); // Parse once, outside every operation timer.
            for (auto const& op : operation_names)
            {
                if (opt.operation == "all" || opt.operation == op) { measure(shape, op, opt, corpus); }
            }
        }
        return 0;
    }
    catch (std::exception const& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}

