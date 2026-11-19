// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Licensed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt

#include <boost/core/lightweight_test.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/strategies/relate/cartesian.hpp>
#include <boost/geometry/strategies/relate/cartesian_exact.hpp>

#include "exact_additional_issue_cases.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace bg = boost::geometry;
namespace cases = exact_additional_issue_cases;
namespace ec = bg::strategies::detail::exact_cartesian;
using rational = ec::exact_number;
using rational_point = bg::model::d2::point_xy<rational>;
using rational_line = bg::model::linestring<rational_point>;
using rational_lines = bg::model::multi_linestring<rational_line>;
using rational_polygon = bg::model::polygon<rational_point>;
using rational_polygons = bg::model::multi_polygon<rational_polygon>;
using rational_points = bg::model::multi_point<rational_point>;

cases::typed_pair const& fixture(char const* name)
{
    auto const& items = cases::pairs();
    auto const found = std::find_if(items.begin(), items.end(), [&](cases::typed_pair const& item)
    {
        return std::string(item.data.name) == name;
    });
    if (found == items.end())
    {
        throw std::logic_error(std::string("Missing fixture: ") + name);
    }
    return *found;
}

bool same_point(rational_point const& a, rational_point const& b)
{
    return bg::get<0>(a) == bg::get<0>(b) && bg::get<1>(a) == bg::get<1>(b);
}

template <typename Point>
rational_point lift(Point const& p)
{
    return rational_point(ec::exact_value(bg::get<0>(p)), ec::exact_value(bg::get<1>(p)));
}

rational cross(rational_point const& a, rational_point const& b)
{
    return bg::get<0>(a) * bg::get<1>(b) - bg::get<1>(a) * bg::get<0>(b);
}

rational_point subtract(rational_point const& a, rational_point const& b)
{
    return rational_point(bg::get<0>(a) - bg::get<0>(b), bg::get<1>(a) - bg::get<1>(b));
}

bool on_segment(rational_point const& p, rational_point const& a, rational_point const& b)
{
    if (cross(subtract(p, a), subtract(b, a)) != 0)
    {
        return false;
    }
    return bg::get<0>(p) >= (std::min)(bg::get<0>(a), bg::get<0>(b))
        && bg::get<0>(p) <= (std::max)(bg::get<0>(a), bg::get<0>(b))
        && bg::get<1>(p) >= (std::min)(bg::get<1>(a), bg::get<1>(b))
        && bg::get<1>(p) <= (std::max)(bg::get<1>(a), bg::get<1>(b));
}

void append_unique(std::vector<rational_point>& points, rational_point const& point)
{
    if (std::none_of(points.begin(), points.end(), [&](rational_point const& p)
        { return same_point(p, point); }))
    {
        points.push_back(point);
    }
}

// Independent rational 2x2 solve, including collinear endpoints. No BG
// predicates/constructions are used to produce this segment-pair oracle.
void segment_oracle(rational_point const& a, rational_point const& b,
                    rational_point const& c, rational_point const& d,
                    std::vector<rational_point>& result)
{
    auto const r = subtract(b, a);
    auto const s = subtract(d, c);
    rational const denominator = cross(r, s);
    if (denominator == 0)
    {
        for (auto const& p : {a, b, c, d})
        {
            if (on_segment(p, a, b) && on_segment(p, c, d))
            {
                append_unique(result, p);
            }
        }
        return;
    }
    auto const delta = subtract(c, a);
    rational const t = cross(delta, s) / denominator;
    rational const u = cross(delta, r) / denominator;
    if (t >= 0 && t <= 1 && u >= 0 && u <= 1)
    {
        append_unique(result, rational_point(bg::get<0>(a) + t * bg::get<0>(r),
                                             bg::get<1>(a) + t * bg::get<1>(r)));
    }
}

double binary64(std::uint64_t bits)
{
    static_assert(sizeof(double) == sizeof(bits) && std::numeric_limits<double>::is_iec559
                  && std::numeric_limits<double>::digits == 53,
                  "Issue 528 requires IEEE-754 binary64");
    double result;
    std::memcpy(&result, &bits, sizeof(result));
    BOOST_TEST(std::isfinite(result));
    return result;
}

void issue_528()
{
    using point = bg::model::d2::point_xy<double>;
    bg::model::polygon<point> polygon;
    bg::model::linestring<point> line;
    for (auto const& p : cases::issue_528_polygon())
    {
        polygon.outer().emplace_back(binary64(p.x), binary64(p.y));
    }
    for (auto const& p : cases::issue_528_line())
    {
        line.emplace_back(binary64(p.x), binary64(p.y));
    }
    std::vector<rational_point> expected;
    for (std::size_t i = 1; i < polygon.outer().size(); ++i)
    {
        segment_oracle(lift(polygon.outer()[i - 1]), lift(polygon.outer()[i]),
                       lift(line[0]), lift(line[1]), expected);
    }
    bg::strategies::relate::cartesian_exact<double> strategy;
    bg::correct(polygon, strategy);
    rational_points actual;
    bg::intersection(polygon, line, actual, strategy);
    BOOST_TEST_EQ(actual.size(), expected.size());
    for (auto const& p : actual)
    {
        BOOST_TEST(on_segment(p, lift(line[0]), lift(line[1])));
        BOOST_TEST(std::any_of(expected.begin(), expected.end(), [&](rational_point const& q)
            { return same_point(p, q); }));
    }
    for (auto const& p : expected)
    {
        BOOST_TEST(std::any_of(actual.begin(), actual.end(), [&](rational_point const& q)
            { return same_point(p, q); }));
    }
    std::cout << "528 rational boundary oracle points=" << expected.size() << '\n';
}

void issue_440()
{
    using point = bg::model::d2::point_xy<std::int32_t>;
    using segment = bg::model::segment<point>;
    std::array<point, 2> const a = {{point(4710043, -2665034), point(4685409, -2666385)}};
    std::array<point, 2> const b = {{point(4700000, -2500000), point(4700000, -3500000)}};
    rational const y = rational(-2665034) - rational(1351) * 10043 / 24634;
    bg::strategies::relate::cartesian_exact<std::int32_t> strategy;
    for (unsigned reverse_a = 0; reverse_a < 2; ++reverse_a)
    {
        for (unsigned reverse_b = 0; reverse_b < 2; ++reverse_b)
        {
            segment const first(a[reverse_a], a[1 - reverse_a]);
            segment const second(b[reverse_b], b[1 - reverse_b]);
            for (bool swap : {false, true})
            {
                rational_points result;
                bg::intersection(swap ? second : first, swap ? first : second, result, strategy);
                BOOST_TEST_EQ(result.size(), 1u);
                if (!result.empty())
                {
                    BOOST_TEST(bg::get<0>(result.front()) == 4700000);
                    BOOST_TEST(bg::get<1>(result.front()) == y);
                }
            }
        }
    }
}

void issue_831()
{
    using point = bg::model::d2::point_xy<float>;
    bg::model::linestring<point> line;
    bg::model::polygon<point> polygon;
    auto const& item = fixture("issue_831_crosses").data;
    bg::read_wkt(item.first, line);
    bg::read_wkt(item.second, polygon);
    bg::strategies::relate::cartesian_exact<float> strategy;
    bg::correct(polygon, strategy);
    auto const matrix = bg::relation(line, polygon, strategy).str();
    BOOST_TEST(bg::crosses(line, polygon, strategy));
    bg::reverse(line);
    BOOST_TEST(bg::crosses(line, polygon, strategy));
    BOOST_TEST_EQ(bg::relation(line, polygon, strategy).str(), matrix);

    rational_polygon exact_polygon;
    rational_line exact_line;
    for (auto const& p : polygon.outer()) exact_polygon.outer().push_back(lift(p));
    for (auto const& p : line) exact_line.push_back(lift(p));
    bg::strategies::relate::cartesian<rational> validation;
    BOOST_TEST_EQ(bg::relation(exact_line, exact_polygon, validation).str(), matrix);
    std::cout << "831 exact matrix=" << matrix << '\n';
}

void issue_656()
{
    using point = bg::model::d2::point_xy<double>;
    bg::model::multi_linestring<bg::model::linestring<point>> lines;
    bg::model::multi_polygon<bg::model::polygon<point>> polygons;
    auto const& item = fixture("issue_656_linear_intersection");
    bg::read_wkt(item.data.first, lines);
    bg::read_wkt(item.data.second, polygons);
    bg::strategies::relate::cartesian_exact<double> strategy;
    bg::correct(polygons, strategy);
    rational_lines expected;
    bg::read_wkt(item.expected_set, expected);
    bg::strategies::relate::cartesian<rational> validation;
    for (bool reverse : {false, true})
    {
        if (reverse) bg::reverse(lines);
        rational_lines result;
        bg::intersection(lines, polygons, result, strategy);
        BOOST_TEST_EQ(result.size(), 1u);
        BOOST_TEST(bg::equals(result, expected, validation));
        BOOST_TEST(bg::is_valid(result, validation));
    }
}

void issue_1232()
{
    using point = bg::model::d2::point_xy<double>;
    bg::model::linestring<point> line;
    bg::model::multi_linestring<bg::model::linestring<point>> other;
    auto const& item = fixture("issue_1232_line_difference").data;
    bg::read_wkt(item.first, line);
    bg::read_wkt(item.second, other);
    bg::strategies::relate::cartesian_exact<double> strategy;
    bg::strategies::relate::cartesian<rational> validation;
    rational_lines forward, reverse;
    bg::difference(line, other, forward, strategy);
    bg::reverse(line);
    bg::difference(line, other, reverse, strategy);
    BOOST_TEST(bg::covered_by(rational_point(7, 4), forward, validation));
    BOOST_TEST(bg::covered_by(rational_point(5, -8), forward, validation));
    BOOST_TEST(bg::equals(forward, reverse, validation));
    BOOST_TEST(bg::is_valid(forward, validation));
}

void issue_1235()
{
    using point = bg::model::d2::point_xy<double>;
    bg::model::polygon<point> polygon;
    bg::model::multi_polygon<bg::model::polygon<point>> polygons;
    auto const& item = fixture("issue_1235_tuple").data;
    bg::read_wkt(item.first, polygon);
    bg::read_wkt(item.second, polygons);
    bg::strategies::relate::cartesian_exact<double> strategy;
    bg::correct(polygon, strategy);
    bg::correct(polygons, strategy);
    std::tuple<rational_points, rational_lines, rational_polygons> result;
    bg::intersection(polygon, polygons, result, strategy);
    BOOST_TEST(std::get<0>(result).empty());
    BOOST_TEST(std::get<1>(result).empty());
    BOOST_TEST_EQ(std::get<2>(result).size(), 2u);
    bg::strategies::relate::cartesian<rational> validation;
    BOOST_TEST(bg::is_valid(std::get<2>(result), validation));
    rational_polygons areal;
    bg::intersection(polygon, polygons, areal, strategy);
    BOOST_TEST(bg::relate(areal, std::get<2>(result), bg::de9im::mask("T*F**FFF*"), validation));
}

void issue_858()
{
    using point = bg::model::d2::point_xy<std::int32_t>;
    bg::strategies::relate::cartesian_exact<std::int32_t> strategy;
    for (int y = 0; y < 8; ++y)
    {
        bg::model::linestring<point> const line{{-1, y}, {17, y}};
        rational_point const shared(rational(9) - rational(2) * y / 9, rational(y));
        for (auto const* wkt : {cases::data::issue_858_0, cases::data::issue_858_1})
        {
            bg::model::polygon<point> polygon;
            bg::read_wkt(wkt, polygon);
            bg::correct(polygon, strategy);
            rational_lines result;
            bg::intersection(line, polygon, result, strategy);
            BOOST_TEST_EQ(result.size(), 1u);
            if (!result.empty() && !result.front().empty())
            {
                BOOST_TEST(same_point(result.front().front(), shared)
                           || same_point(result.front().back(), shared));
            }
        }
    }
}

template <typename Coordinate>
void collapsed_clip(bg::model::linestring<bg::model::d2::point_xy<Coordinate>> const& input,
                    bg::model::box<bg::model::d2::point_xy<Coordinate>> const& box,
                    rational_point const& first, rational_point const& second)
{
    using point = bg::model::d2::point_xy<Coordinate>;
    bg::strategies::relate::cartesian_exact<Coordinate> strategy;
    rational_lines exact;
    bg::intersection(input, box, exact, strategy);
    BOOST_TEST_EQ(exact.size(), 1u);
    if (!exact.empty())
    {
        BOOST_TEST_EQ(exact.front().size(), 2u);
        if (exact.front().size() == 2)
        {
            BOOST_TEST(same_point(exact.front()[0], first));
            BOOST_TEST(same_point(exact.front()[1], second));
            BOOST_TEST(!same_point(first, second));
        }
    }
    bg::model::multi_linestring<bg::model::linestring<point>> rounded;
    bg::intersection(input, box, rounded, strategy);
    BOOST_TEST_EQ(rounded.size(), 1u);
    if (!rounded.empty())
    {
        BOOST_TEST_EQ(rounded.front().size(), 2u);
        if (rounded.front().size() == 2)
        {
            BOOST_TEST(bg::get<0>(rounded.front()[0]) == bg::get<0>(rounded.front()[1]));
            BOOST_TEST(bg::get<1>(rounded.front()[0]) == bg::get<1>(rounded.front()[1]));
        }
    }
}

void issue_1448()
{
    using point = bg::model::d2::point_xy<std::int32_t>;
    using line = bg::model::linestring<point>;
    bg::model::box<point> const box(point(0, 0), point(4096, 4096));
    bg::strategies::relate::cartesian_exact<std::int32_t> strategy;
    line input{{4096, 4000}, {4100, 4010}};
    for (bool reverse : {false, true})
    {
        if (reverse) bg::reverse(input);
        rational_lines exact;
        bg::intersection(input, box, exact, strategy);
        bg::model::multi_linestring<line> rounded;
        bg::intersection(input, box, rounded, strategy);
#ifdef BOOST_GEOMETRY_TEST_FAILURES
        BOOST_TEST(exact.empty());
        BOOST_TEST(rounded.empty());
#endif
        std::cout << "1448 dimensional advisory exact_lines=" << exact.size()
                  << " exact_points=" << bg::num_points(exact)
                  << "; legacy singleton contact retained, not rounding\n";
    }

    // Cardinality is decided before native rounding. These two exact endpoints
    // round to one integer point, but the emitted line must retain both slots.
    collapsed_clip<std::int32_t>(line{{-1, 2}, {2, -5}},
        bg::model::box<point>(point(-10, -10), point(0, 0)),
        rational_point(-rational(1) / 7, 0), rational_point(0, -rational(1) / 3));

    using double_point = bg::model::d2::point_xy<double>;
    double const base = 9007199254740992.0;
    collapsed_clip<double>(bg::model::linestring<double_point>{
        {base - 2, base + 6}, {base + 4, base - 13}},
        bg::model::box<double_point>({base - 20, base - 20}, {base, base}),
        rational_point(rational(9007199254740992LL) - rational(2) / 19, rational(9007199254740992LL)),
        rational_point(rational(9007199254740992LL), rational(9007199254740992LL) - rational(1) / 3));
}

void literal_orientation(unsigned issue)
{
    using point = bg::model::d2::point_xy<double>;
    for (auto const& item : cases::foundation_orientations())
    {
        if (item.issue != issue) continue;
        auto const& v = item.coordinates;
        point const a(v[0], v[1]), b(v[2], v[3]), c(v[4], v[5]);
        rational const determinant = cross(subtract(lift(b), lift(a)), subtract(lift(c), lift(a)));
        int const expected = determinant > 0 ? 1 : determinant < 0 ? -1 : 0;
        BOOST_TEST_EQ(expected, item.expected_sign);
        auto const x = ec::make_point(a), y = ec::make_point(b), z = ec::make_point(c);
        BOOST_TEST_EQ(ec::orientation(x, y, z), expected);
        BOOST_TEST_EQ(ec::orientation(x, z, y), -expected);
        BOOST_TEST_EQ(ec::orientation(y, z, x), expected);
        BOOST_TEST_EQ(ec::orientation(z, x, y), expected);
        std::cout << item.name << " exact_sign=" << expected
                  << " reported_old_sign=" << item.reported_old_sign << '\n';
    }
}

void issue_702() { literal_orientation(702); }
void issue_918() { literal_orientation(918); }

template <typename Coordinate>
void epsilon_ordering()
{
    ec::point_cache_scope<Coordinate> operation;
    Coordinate const epsilon = std::numeric_limits<Coordinate>::epsilon();
    std::array<Coordinate, 3> const native{{1, Coordinate(1 - epsilon), Coordinate(1 - epsilon / 2)}};
    std::array<rational, 3> const values{{ec::exact_value(native[0]),
        ec::exact_value(native[1]), ec::exact_value(native[2])}};
    bg::strategies::relate::cartesian_exact<Coordinate> strategy;
    auto const topology = strategy.ring_topology();
    using rational_point = bg::model::d2::point_xy<rational>;
    // Width ordering is coordinate ordering, not an area-magnitude policy.
    auto const less = [&](rational const& a, rational const& b)
    {
        return ec::static_compare_coordinate<0, Coordinate>(rational_point(a, 0), rational_point(b, 0)) < 0;
    };
    BOOST_TEST(less(values[1], values[2]));
    BOOST_TEST(less(values[2], values[0]));
    rational const width = values[2] - values[0];
    bg::model::ring<rational_point> const reversed_thin{{0, 0}, {0, 1}, {width, 1}, {width, 0}, {0, 0}};
    auto const analysis = topology.analyze(reversed_thin);
    BOOST_TEST(analysis.valid);
    BOOST_TEST_EQ(analysis.direction, -1);
    for (auto const& a : values)
    {
        BOOST_TEST(!less(a, a));
        for (auto const& b : values)
        {
            BOOST_TEST(less(a, b) == (a < b));
            BOOST_TEST(!(less(a, b) && less(b, a)));
            for (auto const& c : values)
                BOOST_TEST(!(less(a, b) && less(b, c)) || less(a, c));
        }
    }
    std::array<unsigned, 3> order{{0, 1, 2}};
    do
    {
        std::array<rational, 3> const permutation{{values[order[0]], values[order[1]], values[order[2]]}};
        auto const extremes = std::minmax_element(permutation.begin(), permutation.end(), less);
        BOOST_TEST(*extremes.first == values[1]);
        BOOST_TEST(*extremes.second == values[0]);
    }
    while (std::next_permutation(order.begin(), order.end()));
}

void issue_1148()
{
    // Literal issue examples around 1. Legacy math::smaller is not changed or
    // expected to become an exact strict-order comparator by this opt-in path.
    epsilon_ordering<float>();
    epsilon_ordering<double>();
}

void issue_1253()
{
    using point = bg::model::d2::point_xy<double>;
    bg::model::linestring<point> line;
    bg::model::polygon<point> polygon;
    bg::read_wkt(cases::data::issue_1253_0, line);
    bg::read_wkt(cases::data::issue_1253_1, polygon);
    std::vector<rational_point> expected;
    for (std::size_t i = 1; i < polygon.outer().size(); ++i)
        segment_oracle(lift(polygon.outer()[i - 1]), lift(polygon.outer()[i]),
            lift(line[0]), lift(line[1]), expected);
    BOOST_TEST(expected.empty());
    // Both segment endpoints lie outside the triangle's x extent. With no
    // boundary crossing, a wholly-interior segment is impossible here.
    for (auto const& p : polygon.outer())
    {
        BOOST_TEST(bg::get<0>(line.front()) < bg::get<0>(p));
        BOOST_TEST(bg::get<0>(line.back()) > bg::get<0>(p));
    }
    bg::strategies::relate::cartesian_exact<double> strategy;
    bg::correct(polygon, strategy);
    rational_lines result;
    bg::intersection(line, polygon, result, strategy);
    std::size_t const count = expected.size() == 2 ? 1 : 0;
    BOOST_TEST_EQ(result.size(), count);
    if (expected.size() == 2 && !result.empty())
    {
        BOOST_TEST_EQ(result.front().size(), 2u);
        if (result.front().size() == 2)
            for (auto const& p : result.front())
                BOOST_TEST(same_point(p, expected[0]) || same_point(p, expected[1]));
        BOOST_TEST(!same_point(result.front().front(), result.front().back()));
    }
    bg::model::multi_linestring<bg::model::linestring<point>> rounded;
    bg::intersection(line, polygon, rounded, strategy);
    BOOST_TEST_EQ(rounded.size(), count);
    if (!rounded.empty()) BOOST_TEST_EQ(rounded.front().size(), 2u);
    std::cout << "1253 rational_boundary_points=" << expected.size()
              << " exact_lines=" << result.size() << '\n';
}

rational top_chain_area(bg::model::polygon<bg::model::d2::point_xy<double>> const& polygon)
{
    rational const xmin = ec::exact_value(-5.8), xmax = ec::exact_value(-5.5);
    rational const bottom = ec::exact_value(bg::get<1>(polygon.outer().front()));
    rational area = 0;
    for (std::size_t i = 1; i < polygon.outer().size(); ++i)
    {
        auto const a = lift(polygon.outer()[i - 1]), b = lift(polygon.outer()[i]);
        if (bg::get<0>(a) >= bg::get<0>(b)) continue;
        rational const lo = (std::max)(xmin, bg::get<0>(a));
        rational const hi = (std::min)(xmax, bg::get<0>(b));
        if (lo >= hi) continue;
        rational const slope = (bg::get<1>(b) - bg::get<1>(a)) / (bg::get<0>(b) - bg::get<0>(a));
        rational const ylo = bg::get<1>(a) + (lo - bg::get<0>(a)) * slope;
        rational const yhi = bg::get<1>(a) + (hi - bg::get<0>(a)) * slope;
        area += (hi - lo) * ((ylo + yhi) / 2 - bottom);
    }
    return area;
}

void issue_1044()
{
    using point = bg::model::d2::point_xy<double>;
    using polygon = bg::model::polygon<point>;
    bg::strategies::relate::cartesian_exact<double> strategy;
    bg::strategies::relate::cartesian<rational> validation;
    for (auto const* second : {cases::data::issue_1044_1, cases::data::issue_1044_control_1})
    {
        for (bool shifted : {false, true})
        {
            polygon a, b;
            bg::read_wkt(cases::data::issue_1044_0, a);
            bg::read_wkt(second, b);
            if (shifted)
                for (auto& p : b.outer()) bg::set<1>(p, bg::get<1>(p) + 1e-13);
            // The downloaded arrays have an increasing-x upper chain. Integrate
            // its clipped trapezoids before orientation normalization; this
            // oracle uses no BG overlay and retains the tiny top-edge notch.
            rational const expected_area = top_chain_area(b);
            rational const bottom = ec::exact_value(bg::get<1>(b.outer().front()));
            double const top = shifted ? 5.3 + 1e-13 : 5.3;
            bg::correct(a, strategy);
            bg::correct(b, strategy);
            for (bool swap : {false, true})
            {
                rational_polygons result;
                bg::intersection(swap ? b : a, swap ? a : b, result, strategy);
                BOOST_TEST_EQ(result.size(), 1u);
                BOOST_TEST(bg::is_valid(result, validation));
                BOOST_TEST(bg::area(result, validation) == expected_area);
                bg::model::box<rational_point> bounds;
                bg::envelope(result, bounds);
                BOOST_TEST((bg::get<bg::min_corner, 0>(bounds) == ec::exact_value(-5.8)));
                BOOST_TEST((bg::get<bg::max_corner, 0>(bounds) == ec::exact_value(-5.5)));
                BOOST_TEST((bg::get<bg::min_corner, 1>(bounds) == bottom));
                BOOST_TEST((bg::get<bg::max_corner, 1>(bounds) == ec::exact_value(top)));
            }
        }
    }
}

template <typename Multi>
void rotate_starts(Multi& geometry)
{
    auto rotate = [](auto& ring)
    {
        if (ring.size() <= 3) return;
        ring.pop_back();
        std::rotate(ring.begin(), ring.begin() + 1, ring.end());
        ring.push_back(ring.front());
    };
    for (auto& polygon : geometry)
    {
        rotate(polygon.outer());
        for (auto& ring : polygon.inners()) rotate(ring);
    }
}

void issue_669()
{
    ec::point_cache_scope operation;
    using point = bg::model::d2::point_xy<double>;
    using polygon = bg::model::polygon<point, false>;
    using multi = bg::model::multi_polygon<polygon>;
    using working = bg::model::multi_polygon<bg::model::polygon<ec::exact_point<double>, false>>;
    multi input, subtractor;
    bg::read_wkt(cases::data::issue_669_0, input);
    bg::read_wkt(cases::data::issue_669_1, subtractor);
    bg::strategies::relate::cartesian_exact<double> strategy;
    bg::strategies::relate::cartesian<rational> validation;
    bg::correct(input, strategy);
    bg::correct(subtractor, strategy);
    auto native_valid = [&](multi const& source)
    {
        ec::point_cache_scope operation;
        working exact;
        for (auto const& polygon : source)
        {
            exact.emplace_back();
            for (auto const& p : polygon.outer()) exact.back().outer().push_back(ec::make_point(p));
            exact.back().inners().resize(polygon.inners().size());
            for (std::size_t i = 0; i < polygon.inners().size(); ++i)
                for (auto const& p : polygon.inners()[i])
                    exact.back().inners()[i].push_back(ec::make_point(p));
        }
        return bg::is_valid(exact, validation);
    };
    BOOST_TEST(native_valid(input));
    BOOST_TEST(native_valid(subtractor));
    auto check_difference = [&](auto source, char const* variant)
    {
        working previous;
        for (bool rotate : {false, true})
        {
            if (rotate) rotate_starts(source);
            working result, overlap;
            bg::difference(source, subtractor, result, strategy);
            bg::intersection(source, subtractor, overlap, strategy);
            BOOST_TEST(bg::is_valid(result, validation));
            BOOST_TEST(bg::area(result, validation) + bg::area(overlap, validation) == bg::area(source, validation));
            BOOST_TEST(bg::covered_by(result, source, strategy));
            BOOST_TEST(bg::relate(result, subtractor, bg::de9im::mask("F********"), strategy));
            if (rotate)
                BOOST_TEST(bg::relate(previous, result, bg::de9im::mask("T*F**FFF*"), validation));
            else previous = result;
            std::cout << "669 " << variant << " rotate=" << rotate
                      << " components=" << result.size() << " holes=" << bg::num_interior_rings(result)
                      << " overlap_area=" << static_cast<double>(bg::area(overlap, validation)) << '\n';
        }
    };
    check_difference(input, "full");
    bg::model::box<point> box;
    bg::envelope(subtractor, box);
    bg::set<bg::min_corner, 0>(box, bg::get<bg::min_corner, 0>(box) - 1);
    bg::set<bg::min_corner, 1>(box, bg::get<bg::min_corner, 1>(box) - 1);
    bg::set<bg::max_corner, 0>(box, bg::get<bg::max_corner, 0>(box) + 1);
    bg::set<bg::max_corner, 1>(box, bg::get<bg::max_corner, 1>(box) + 1);
    working cropped;
    bg::intersection(box, input, cropped, strategy);
    BOOST_TEST(bg::is_valid(cropped, validation));
    check_difference(cropped, "exact_crop");
    multi rounded;
    bg::intersection(box, input, rounded, strategy);
    // Native-feedback validity is its own precondition, not a rational-output
    // oracle. Do not call predicates on an invalid rounded geometry.
    if (native_valid(rounded)) check_difference(rounded, "rounded_exact_crop");
    else std::cout << "669 rounded_exact_crop invalid; native-feedback branch not a valid exact oracle\n";
    multi legacy_crop;
    bg::intersection(box, input, legacy_crop);
    if (native_valid(legacy_crop)) check_difference(legacy_crop, "source_default_crop");
    else std::cout << "669 source_default_crop invalid; original feedback branch not replayed\n";
}

int main(int argc, char** argv)
{
    // An optional issue number permits a focused replay while other cases are
    // being integrated. This unit intentionally has no stock-baseline claim.
    struct test_case { unsigned issue; void (*run)(); };
    test_case const tests[] = {{528, issue_528}, {440, issue_440}, {831, issue_831},
        {656, issue_656}, {1232, issue_1232}, {1235, issue_1235}, {858, issue_858},
        {1448, issue_1448}, {702, issue_702}, {918, issue_918}, {1148, issue_1148},
        {1253, issue_1253}, {1044, issue_1044}, {669, issue_669}};
    bool matched = argc == 1;
    for (auto const& test : tests)
    {
        if (argc == 1 || std::string(argv[1]) == std::to_string(test.issue))
        {
            matched = true;
            std::cout << "Replaying issue " << test.issue << std::endl;
            test.run();
        }
    }
    BOOST_TEST(matched);
    return boost::report_errors();
}
