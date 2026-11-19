// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.

#include <boost/core/lightweight_test.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/index/equal_to.hpp>

#include <cmath>
#include <deque>
#include <limits>
#include <sstream>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace bg = boost::geometry;
template <typename C> using point = bg::model::d2::point_xy<C>;
using polygon = bg::model::polygon<point<double>>;
using polygons = bg::model::multi_polygon<polygon>;
using variant = boost::variant<polygon>;
using static_disjoint = bg::de9im::static_mask<'F', 'F', '*', 'F', 'F', '*', '*', '*', '*'>;
using static_equals = bg::de9im::static_mask<'T', '*', 'F', '*', '*', 'F', 'F', 'F', '*'>;
using static_edge_touch = bg::de9im::static_mask<'F', '*', '*', '*', '1'>;
auto const legacy = bg::default_strategy{};

struct observed_legacy : bg::strategies::relate::cartesian<>
{
    using base = bg::strategies::relate::cartesian<>;

    explicit observed_legacy(unsigned& calls) : calls(calls) {}

    template <typename A, typename B, typename Apply>
    decltype(auto) with_sources(A const&, B const&, Apply&& apply) const
    {
        ++calls;
        return std::forward<Apply>(apply)(static_cast<base const&>(*this));
    }

    unsigned& calls;
};

#if defined(BOOST_GEOMETRY_TEST_EXPECT_EXACT_DEFAULT_PILOT) \
    && !defined(BOOST_GEOMETRY_DETAIL_CARTESIAN_EXACT_DEFAULT_PILOT)
#error The strict C++17 test must activate the pilot.
#endif
#if defined(BOOST_GEOMETRY_TEST_EXPECT_LEGACY_DEFAULT_PILOT) \
    && defined(BOOST_GEOMETRY_STRATEGIES_RELATE_CARTESIAN_EXACT_HPP)
#error The legacy build must not see the exact strategy header.
#endif

#ifdef BOOST_GEOMETRY_DETAIL_CARTESIAN_EXACT_DEFAULT_PILOT
namespace pilot = bg::detail::cartesian_default_pilot;
using approximate = bg::strategies::relate::cartesian_exact<double, bg::strategies::relate::approximate_output>;
static_assert(std::is_same<pilot::strategy, approximate>::value, "Keep approximate output");
static_assert(pilot::admits_pair<polygon, polygon>, "Only this concrete input model is enabled");
static_assert(pilot::admits_output<polygon, polygon, polygons>, "Use the canonical multi output");
static_assert(!pilot::admits_pair<variant, polygon>, "Reject before visiting dynamic types");
static_assert(!pilot::admits_pair<polygons, polygon>, "No multipart inputs");
static_assert(!pilot::admits_pair<bg::model::polygon<point<float>>, polygon>, "No mixed precision");
static_assert(!pilot::admits_pair<bg::model::polygon<point<double>, false>, polygon>, "No other winding");
static_assert(!pilot::admits_pair<bg::model::polygon<point<double>, true, false>, polygon>, "No open rings");
static_assert(!pilot::admits_pair<bg::model::polygon<point<double>, true, true, std::deque>, polygon>,
              "No other storage models");
static_assert(!pilot::admits_output<polygon, polygon, std::vector<polygon>>, "No other destinations");
static_assert(!pilot::admits_output<polygon, polygon, polygons const>, "No const destinations");
static_assert(!pilot::admits_output<polygon, polygon, std::tuple<polygons>>, "No tuple destinations");
static_assert(std::is_same<bg::strategies::relate::services::default_strategy<polygon, polygon>::type,
                           bg::strategies::relate::cartesian<>>::value, "Do not replace the umbrella");
auto const reference = approximate{};
#else
auto const reference = legacy;
#endif

template <typename Geometry = polygon>
Geometry read(char const* wkt)
{
    Geometry result;
    bg::read_wkt(wkt, result);
    return result;
}

template <typename Geometry>
std::string signature(Geometry const& g)
{
    std::ostringstream out;
    out << std::hexfloat << bg::wkt(g);
    return out.str();
}

void test_relation()
{
    auto const a = read("POLYGON((0 0,0 1,1 1,1 0,0 0))");
    auto b = read("POLYGON((1 0,1 1,2 1,2 0,1 0))");
    double const n = std::nextafter(1.0, 2.0);
    for (auto& p : b.outer())
        if (bg::get<0>(p) == 1)
            bg::set<0>(p, n);
    bg::de9im::mask const disjoint("FF*FF****");

    BOOST_TEST_EQ(bg::relation(a, b, legacy).str(), "FF2F11212");
    BOOST_TEST(!bg::relate(a, b, disjoint, legacy));
    BOOST_TEST_EQ(bg::relation(a, b).str(), bg::relation(a, b, reference).str());
    BOOST_TEST_EQ(bg::relate(a, b, disjoint), bg::relate(a, b, disjoint, reference));
#ifdef BOOST_GEOMETRY_DETAIL_CARTESIAN_EXACT_DEFAULT_PILOT
    BOOST_TEST_EQ(bg::relation(a, b).str(), "FF2FF1212");
    BOOST_TEST(pilot::admits_values(a, b));
    BOOST_TEST(bg::relate(a, b, disjoint));
#else
    BOOST_TEST_EQ(bg::relation(a, b).str(), "FF2F11212");
    BOOST_TEST(!bg::relate(a, b, disjoint));
#endif
    // These public and internal users of the shared umbrella are not opted in.
    BOOST_TEST(bg::intersects(a, b));
    BOOST_TEST(bg::touches(a, b));
    BOOST_TEST_EQ(bg::relate(a, b, static_disjoint{}), bg::relate(a, b, disjoint));
    BOOST_TEST(!bg::relate(a, b, static_disjoint{}, legacy));
    BOOST_TEST(!bg::relate(a, b, static_disjoint{} || static_equals{}));
    BOOST_TEST(!bg::relate(a, b, disjoint || bg::de9im::mask("T*F**FFF*")));
    BOOST_TEST_EQ(bg::relation(variant(a), b).str(), bg::relation(a, b, legacy).str());
    BOOST_TEST(!bg::relate(a, variant(b), disjoint));
    BOOST_TEST(!bg::relate(a, variant(b), static_disjoint{}));
    BOOST_TEST_EQ(bg::relation(variant(b), variant(a)).str(), bg::relation(b, a, legacy).str());

    unsigned calls = 0;
    observed_legacy const custom(calls);
    BOOST_TEST_EQ(bg::relation(a, b, custom).str(), bg::relation(a, b, legacy).str());
    BOOST_TEST(!bg::relate(a, b, disjoint, custom));
    BOOST_TEST(!bg::relate(a, b, static_disjoint{}, custom));
    BOOST_TEST(bg::touches(a, b, custom));
    BOOST_TEST(bg::intersects(a, b, custom));
    BOOST_TEST_EQ(calls, 5u);

    for (auto& p : b.outer())
        if (bg::get<1>(p) == 0)
            bg::set<1>(p, std::ldexp(1.0, -65));
    BOOST_TEST_EQ(bg::relation(a, b).str(), bg::relation(a, b, legacy).str());
    BOOST_TEST_EQ(bg::relate(a, b, static_disjoint{}), bg::relate(a, b, static_disjoint{}, legacy));
}

template <char... Cells>
void check_mask(polygon const& a, polygon const& b, std::string const& matrix,
                bg::de9im::static_mask<Cells...> mask)
{
    char const pattern[] = {Cells..., '\0'};
    bool expected = true;
    for (std::size_t i = 0; i < matrix.size(); ++i)
        expected = expected && (pattern[i] == '*' || pattern[i] == matrix[i]
                                  || (pattern[i] == 'T' && matrix[i] != 'F'));
    BOOST_TEST_EQ(bg::relate(a, b, mask), expected);
    BOOST_TEST_EQ(bg::relate(a, b, bg::de9im::mask(pattern)), expected);
}

void test_static_masks()
{
    auto const a = read("POLYGON((0 0,0 1,1 1,1 0,0 0))");
    auto const edge = read("POLYGON((1 0,1 1,2 1,2 0,1 0))");
    auto const corner = read("POLYGON((1 1,1 2,2 2,2 1,1 1))");
    auto const overlap = read("POLYGON((0.5 0.5,0.5 2,2 2,2 0.5,0.5 0.5))");
    for (auto const* b : {&a, &edge, &corner, &overlap})
    {
#ifdef BOOST_GEOMETRY_DETAIL_CARTESIAN_EXACT_DEFAULT_PILOT
        BOOST_TEST(pilot::admits_values(a, *b));
#endif
        auto const matrix = bg::relation(a, *b).str();
        check_mask(a, *b, matrix, static_equals{});
        check_mask(a, *b, matrix, static_disjoint{});
        check_mask(a, *b, matrix, static_edge_touch{});
        check_mask(a, *b, matrix, bg::de9im::static_mask<'F', '*', '*', '*', '0'>{});
        check_mask(a, *b, matrix, bg::de9im::static_mask<'2', '1', '2', '1', '0', '1', '2', '1', '2'>{});
    }
}

void test_named_equality_boundary()
{
    auto const a = read("POLYGON((0 0,0 1,1 1,1 0,0 0))");
    auto b = a;
    for (auto& p : b.outer())
        if (bg::get<0>(p) == 1)
            bg::set<0>(p, std::nextafter(1.0, 2.0));
    BOOST_TEST(bg::is_valid(a));
    BOOST_TEST(bg::is_valid(b));
#ifdef BOOST_GEOMETRY_DETAIL_CARTESIAN_EXACT_DEFAULT_PILOT
    BOOST_TEST(pilot::admits_values(a, b));
    // The explicit exact umbrella still uses the optimized equals algorithm.
    BOOST_TEST(bg::equals(a, b, approximate{}));
    BOOST_TEST(!bg::relate(a, b, static_equals{}, approximate{}));
#endif
    BOOST_TEST_EQ(bg::relation(a, b).str(), "2FF11F212");
    BOOST_TEST(!bg::relate(a, b, static_equals{}));
    BOOST_TEST(bg::equals(a, b));
    BOOST_TEST(bg::equals(a, b, legacy));

    // Even an explicitly legacy index comparator calls no-strategy equals
    // for polygon payloads. Changing that overload would leak into it.
    using value = std::pair<point<double>, polygon>;
    bg::index::equal_to<value> const equal;
    BOOST_TEST(equal(value{{0, 0}, a}, value{{0, 0}, b}, bg::strategies::relate::cartesian<>{}));
}

auto const intersect = [](auto const& a, auto const& b, auto& out, auto const&... strategy)
{ return bg::intersection(a, b, out, strategy...); };
auto const unite = [](auto const& a, auto const& b, auto& out, auto const&... strategy)
{ bg::union_(a, b, out, strategy...); };
auto const subtract = [](auto const& a, auto const& b, auto& out, auto const&... strategy)
{ bg::difference(a, b, out, strategy...); };
auto const xor_ = [](auto const& a, auto const& b, auto& out, auto const&... strategy)
{ bg::sym_difference(a, b, out, strategy...); };

template <typename Call>
void compare_overlay(Call const& operation, polygon const& a, polygon const& b, polygon const& seed)
{
    polygons actual{seed}, expected{seed};
    operation(a, b, actual);
    operation(a, b, expected, reference);
    BOOST_TEST(actual.size() > 1u);
    BOOST_TEST_EQ(signature(actual), signature(expected));
    BOOST_TEST_EQ(signature(actual.front()), signature(seed));
}

#ifdef BOOST_GEOMETRY_DETAIL_CARTESIAN_EXACT_DEFAULT_PILOT
template <typename Call>
void compare_alias(Call const& operation, polygon const& a, polygon const& b)
{
    polygons actual{a, b};
    while (actual.size() < actual.capacity())
        actual.push_back(a);
    polygons expected = actual;
    auto const before = actual.size();
    operation(a, b, expected, reference);
    operation(actual[0], actual[1], actual);
    BOOST_TEST(actual.size() > before);
    BOOST_TEST_EQ(signature(actual), signature(expected));
}
#endif

void test_overlay()
{
    auto const a = read("POLYGON((0 0,0 4,4 4,4 0,0 0))");
    auto const b = read("POLYGON((1 -1,1 5,2 5,2 -1,1 -1))");
    auto const seed = read("POLYGON((10 10,10 11,11 11,11 10,10 10))");
    compare_overlay(intersect, a, b, seed);
    compare_overlay(unite, b, a, seed);
    compare_overlay(subtract, a, b, seed);
    compare_overlay(xor_, b, a, seed);
#ifdef BOOST_GEOMETRY_DETAIL_CARTESIAN_EXACT_DEFAULT_PILOT
    compare_alias(intersect, a, b);
    compare_alias(unite, a, seed); // Multiple results initially borrow both inputs.
    compare_alias(subtract, a, b);
    compare_alias(xor_, a, b);
#endif

    // Different output container and dynamic input both bypass selection.
    std::vector<polygon> vector_out, vector_expected;
    BOOST_TEST_EQ(intersect(a, b, vector_out), intersect(a, b, vector_expected, legacy));
    BOOST_TEST_EQ(vector_out.size(), 1u);
    BOOST_TEST_EQ(vector_expected.size(), 1u);
    if (vector_out.size() == 1 && vector_expected.size() == 1)
        BOOST_TEST_EQ(signature(vector_out.front()), signature(vector_expected.front()));
    polygons dynamic, expected;
    BOOST_TEST_EQ(intersect(variant(a), b, dynamic), intersect(a, b, expected, legacy));
    BOOST_TEST_EQ(signature(dynamic), signature(expected));

    using integer_polygon = bg::model::polygon<point<int>>;
    auto const i = read<integer_polygon>("POLYGON((0 0,0 4,4 4,4 0,0 0))");
    auto const j = read<integer_polygon>("POLYGON((1 -1,1 5,2 5,2 -1,1 -1))");
    bg::model::multi_polygon<integer_polygon> integer, integer_expected;
    BOOST_TEST_EQ(intersect(i, j, integer), intersect(i, j, integer_expected, legacy));
    BOOST_TEST_EQ(signature(integer), signature(integer_expected));
}

#ifdef BOOST_GEOMETRY_DETAIL_CARTESIAN_EXACT_DEFAULT_PILOT
void test_emission()
{
    auto const a = read("POLYGON((43 7,-24 81,40 75,43 7))");
    auto const b = read("POLYGON((-78 22,63 95,30 30,-78 22))");
    polygons actual, approx, nearest, ordinary, dynamic;
    BOOST_TEST(bg::intersection(a, b, actual));
    BOOST_TEST(bg::intersection(a, b, approx, reference));
    using nearest_strategy = bg::strategies::relate::cartesian_exact<double, bg::strategies::relate::nearest_output>;
    BOOST_TEST(bg::intersection(a, b, nearest, nearest_strategy{}));
    BOOST_TEST(bg::intersection(a, b, ordinary, legacy));
    BOOST_TEST(bg::intersection(variant(a), b, dynamic));
    BOOST_TEST_EQ(signature(actual), signature(approx));
    BOOST_TEST(signature(actual) != signature(nearest));
    BOOST_TEST(signature(actual) != signature(ordinary));
    BOOST_TEST_EQ(signature(dynamic), signature(ordinary));

    // #663: a nonempty exact intersection can collapse on native emission.
    auto ca = read("POLYGON((140.10975222826562 135.99989776999999,140.10975222826562 126.97244753050222,"
                   "140.16255693232961 135.99989776999999,140.10975222826562 135.99989776999999))");
    auto cb = read("POLYGON((151.99995129999999 123.99989777000000,139.99994018999999 126.99990054000000,"
                   "129.99995129999999 124.99989777000000,151.99995129999999 123.99989777000000))");
    bg::correct(ca);
    bg::correct(cb);
    polygons collapsed;
    BOOST_TEST_EQ(bg::relation(ca, cb).str(), "212101212");
    BOOST_TEST(bg::intersection(ca, cb, collapsed));
    BOOST_TEST_EQ(collapsed.size(), 1u);
    BOOST_TEST_EQ(bg::area(collapsed), 0.0);
    BOOST_TEST(!bg::is_valid(collapsed));
}

void test_domain()
{
    BOOST_TEST(pilot::admits_coordinate(-0.0));
    BOOST_TEST(pilot::admits_coordinate(-0x1p-64));
    BOOST_TEST(pilot::admits_coordinate(0x1p16));
    BOOST_TEST(pilot::admits_coordinate(std::nextafter(0x1p16, 0.0)));
    BOOST_TEST(!pilot::admits_coordinate(0x1p-65));
    BOOST_TEST(!pilot::admits_coordinate(std::nextafter(0x1p-64, 1.0)));
    BOOST_TEST(!pilot::admits_coordinate(std::nextafter(0x1p16, 0x1p17)));
    BOOST_TEST(!pilot::admits_coordinate(std::numeric_limits<double>::denorm_min()));
    BOOST_TEST(!pilot::admits_coordinate(std::numeric_limits<double>::infinity()));
    BOOST_TEST(!pilot::admits_coordinate(std::numeric_limits<double>::quiet_NaN()));

    auto const a = read("POLYGON((0 0,0 1,1 1,1 0,0 0))");
    auto b = read("POLYGON((1 0,1 1,2 1,2 0,1 0))");
    for (auto& p : b.outer())
        if (bg::get<0>(p) == 1)
            bg::set<0>(p, std::nextafter(1.0, 2.0));
    BOOST_TEST(pilot::admits_values(a, b));

    // Scan holes as well as exteriors. Nothing is snapped onto the grid.
    auto const off_grid = read("POLYGON((0 0,0 4,4 4,4 0,0 0),"
                               "(1 1,2 0.00000001,2 2,1 2,1 1))");
    BOOST_TEST(!pilot::admits_values(a, off_grid));
}
#endif

int main()
{
    test_relation();
    test_static_masks();
    test_named_equality_boundary();
    test_overlay();
#ifdef BOOST_GEOMETRY_DETAIL_CARTESIAN_EXACT_DEFAULT_PILOT
    test_emission();
    test_domain();
#endif
    return boost::report_errors();
}
