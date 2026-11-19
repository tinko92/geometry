// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.

#ifndef BOOST_GEOMETRY_EXACT_FILTER_STATISTICS
#define BOOST_GEOMETRY_EXACT_FILTER_STATISTICS
#endif

#include <boost/core/lightweight_test.hpp>
#include <boost/geometry/algorithms/intersection.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/strategies/relate/cartesian_exact.hpp>
#include <boost/multiprecision/cpp_int.hpp>

#include <array>
#include <cfenv>
#include <cmath>
#include <iterator>
#include <limits>
#include <random>
#include <vector>

namespace bg = boost::geometry;
namespace ec = bg::strategies::detail::exact_cartesian;
namespace sr = bg::strategies::relate;
namespace st = ec::sterbenz;

using point = bg::model::d2::point_xy<double>;
using line = bg::model::linestring<point>;
using segment = bg::model::segment<point>;

int oracle(point const& a, point const& b, point const& p)
{
    using rational = boost::multiprecision::cpp_rational;
    rational const bx = rational(b.x()) - rational(a.x()), by = rational(b.y()) - rational(a.y());
    rational const px = rational(p.x()) - rational(a.x()), py = rational(p.y()) - rational(a.y());
    rational const left = bx * py, right = by * px;
    return (left > right) - (left < right);
}

void certification()
{
    line const positive{{1, 2}, {2, 4}};
    line const negative{{-1, -2}, {-2, -4}};
    line const mixed_axes{{1, -2}, {2, -4}};
    BOOST_TEST(st::certify_sources(positive, positive));
    BOOST_TEST(st::certify_sources(negative, negative));
    BOOST_TEST(st::certify_sources(mixed_axes, mixed_axes));
    BOOST_TEST(!st::certify_sources(positive, negative));
    BOOST_TEST(!st::certify_sources(positive, line{{1, 2}, {std::nextafter(2., 3.), 4}}));
    BOOST_TEST(!st::certify_sources(positive, line{{1, -0.0}, {2, 4}}));
    BOOST_TEST(!st::certify_sources(positive, line{{1, std::numeric_limits<double>::infinity()}}));
    BOOST_TEST(!st::certify_sources(positive, line{{std::numeric_limits<double>::quiet_NaN(), 4}}));
    BOOST_TEST(!st::certify_sources(line{}, line{}));
    for (double x : {0x1p151, 0x1p-151})
        BOOST_TEST(!st::certify_sources(line{{x, 2}, {x, 4}}, line{{x, 2}, {x, 4}}));

    using polygon = bg::model::polygon<point>;
    polygon with_hole;
    with_hole.outer() = {{1, 2}, {1, 4}, {2, 4}, {2, 2}, {1, 2}};
    with_hole.inners().resize(1);
    with_hole.inners()[0] = {{1.25, 2.5}, {1.75, 2.5}, {1.5, 3.5}, {1.25, 2.5}};
    BOOST_TEST(st::certify_sources(with_hole, positive));
    with_hole.inners()[0][1].x(std::nextafter(2., 3.));
    BOOST_TEST(!st::certify_sources(with_hole, positive));
    static_assert(!st::source<bg::model::linestring<bg::model::d2::point_xy<float>>>, "Double only");
}

void exact_signs()
{
    // The two rounded products tie at +/-0.25. Exact residuals separate
    // them by 2^-104. Reversing the endpoint order must reverse the sign.
    double const e = 0x1p-52;
    for (double y : {1., 2.})
    {
        double const sign = y == 1 ? 1 : -1;
        point const a(1, y), b(1.5, y + sign * (0.5 - e)), p(1.5 + e, y + sign * 0.5);
        BOOST_TEST_EQ(ec::arithmetic_detail::rounded_mul(b.x() - a.x(), p.y() - a.y()),
                      ec::arithmetic_detail::rounded_mul(b.y() - a.y(), p.x() - a.x()));
        BOOST_TEST_EQ(st::orientation(a, b, p), oracle(a, b, p));
        BOOST_TEST_EQ(st::orientation(a, b, p), int(sign));
        BOOST_TEST_EQ(st::orientation(a, p, b), -int(sign));
    }
    BOOST_TEST_EQ(st::orientation(point(1, 1), point(2, 2), point(1.5, 1.5)), 0);
    BOOST_TEST_EQ(st::orientation(point(1, 1), point(1, 2), point(1, 1.5)), 0);

    std::mt19937_64 random(0x5354455242454e5a);
    std::array<int, 4> const exponents{{-150, -1, 26, 149}};
    for (std::size_t i = 0; i < 1024; ++i)
    {
        auto coordinate = [&](std::size_t index)
        { return std::ldexp(1 + double(random() >> 12) * 0x1p-52, exponents[index]); };
        auto make = [&]
        { return point(i % 2 ? -coordinate(i % 4) : coordinate(i % 4), coordinate((i / 4) % 4)); };
        point const a = make(), b = make(), p = make();
        BOOST_TEST(st::certify_sources(line{a, b}, line{p}));
        int const expected = oracle(a, b, p);
        BOOST_TEST_EQ(st::orientation(a, b, p), expected);
        BOOST_TEST_EQ(ec::native_orientation<double>(a, b, p), expected);
    }
}

void public_intersection()
{
    auto& stats = ec::filtering_statistics();
    sr::cartesian_exact<double> const strategy;
    line const a{{1, 1}, {2, 2}}, b{{1, 2}, {2, 1}};
    std::vector<point> output;
    auto const before = stats.sterbenz_orientations;
    bg::intersection(a, b, output, strategy);
    BOOST_TEST(stats.sterbenz_orientations > before);
    BOOST_TEST_EQ(output.size(), 1u);
    if (!output.empty())
    {
        BOOST_TEST_EQ(output.front().x(), 1.5);
        BOOST_TEST_EQ(output.front().y(), 1.5);
    }

    line const c{{-1, -1}, {1, 1}}, d{{-1, 1}, {1, -1}};
    auto const after = stats.sterbenz_orientations;
    output.clear();
    bg::intersection(c, d, output, strategy);
    BOOST_TEST_EQ(stats.sterbenz_orientations, after);
    BOOST_TEST_EQ(output.size(), 1u);
    if (!output.empty())
    {
        BOOST_TEST_EQ(output.front().x(), 0.);
        BOOST_TEST_EQ(output.front().y(), 0.);
    }
}

void operand_provenance()
{
    auto& stats = ec::filtering_statistics();
    using root = ec::borrowed_root<double, sr::approximate_output, segment, segment>;
    using borrowed = ec::borrowed_point<root>;
    point const base(1, 1), right(1.25, 1);
    segment const a(base, point(2, std::nextafter(1., 2.))), b(right, point(1.25, 2));
    ec::borrowed_context<root> context(a, b);
    BOOST_TEST(ec::borrowed_stage<root>(context, b, a).native_differences_exact);

    std::vector<borrowed> output;
    auto emit = [](borrowed const& p) { return p; };
    ec::borrowed_intersect_segments(context, a, b, std::back_inserter(output), emit);
    BOOST_TEST_EQ(output.size(), 1u);
    if (output.empty()) return;
    BOOST_TEST(!output.front().is_native());
    BOOST_TEST_EQ(output.front().approximate(1), 1.);
    // The exact intersection sits above y=1; its binary64 approximation
    // lies ON y=1. Even an enabled native certificate must not be used.
    auto const before = stats.sterbenz_orientations;
    BOOST_TEST_EQ(ec::root_orientation<root>(base, right, output.front(), true), 1);
    BOOST_TEST_EQ(stats.sterbenz_orientations, before);

    segment const c(point(-1, -1), point(1, 1)), d(point(-1, 1), point(1, -1));
    BOOST_TEST(!ec::borrowed_stage<root>(context, c, d).native_differences_exact);
    output.clear();
    ec::borrowed_intersect_segments(context, c, d, std::back_inserter(output), emit);
    BOOST_TEST_EQ(output.size(), 1u);
    BOOST_TEST_EQ(stats.sterbenz_orientations, before);
}

int main()
{
    BOOST_TEST_EQ(std::fegetround(), FE_TONEAREST);
    certification();
    exact_signs();
    public_intersection();
    operand_provenance();
    return boost::report_errors();
}
