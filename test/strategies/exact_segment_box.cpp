// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.

#define BOOST_GEOMETRY_EXACT_FILTER_STATISTICS
#include <boost/core/lightweight_test.hpp>
#include <boost/geometry/algorithms/disjoint.hpp>
#include <boost/geometry/geometries/pointing_segment.hpp>
#include <boost/geometry/strategies/relate/cartesian_exact.hpp>
#include <boost/multiprecision/cpp_int.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <memory>
#include <new>
#include <random>
#include <type_traits>

namespace bg = boost::geometry;
namespace ec = bg::strategies::detail::exact_cartesian;
using rational = boost::multiprecision::cpp_rational;
template <typename C> using point = bg::model::d2::point_xy<C>;
template <typename C> using box = bg::model::box<point<C>>;

namespace
{
bool counting = false;
std::size_t allocations = 0;
}
void* operator new(std::size_t size)
{
    if (counting) { ++allocations; }
    if (void* result = std::malloc(size ? size : 1)) { return result; }
    throw std::bad_alloc();
}
void operator delete(void* value) noexcept { std::free(value); }
void* operator new[](std::size_t size) { return ::operator new(size); }
void operator delete[](void* value) noexcept { ::operator delete(value); }
#if defined(__cpp_sized_deallocation)
void operator delete(void* value, std::size_t) noexcept { std::free(value); }
void operator delete[](void* value, std::size_t) noexcept { std::free(value); }
#endif

template <typename C, typename Tag>
struct indexed_geometry
{
    std::array<C, 4> const& values;
    unsigned& reads;
    indexed_geometry(std::array<C, 4> const& v, unsigned& r) : values(v), reads(r) {}
    indexed_geometry(indexed_geometry const&) = delete;
    indexed_geometry& operator=(indexed_geometry const&) = delete;
};
namespace boost { namespace geometry { namespace traits
{
template <typename C, typename Tag> struct tag<indexed_geometry<C, Tag>> { using type = Tag; };
template <typename C, typename Tag> struct point_type<indexed_geometry<C, Tag>> { using type = ::point<C>; };
template <typename C, typename Tag, std::size_t I, std::size_t D>
struct indexed_access<indexed_geometry<C, Tag>, I, D>
{
    static C get(indexed_geometry<C, Tag> const& value)
    { ++value.reads; return value.values[2 * I + D]; }
};
}}}

// Independent exact slab clipping. Unlike the tested separating-axis predicate,
// this intersects the permissible parameter intervals on the segment.
template <typename Segment, typename Box>
bool oracle_disjoint(Segment const& segment, Box const& bounds)
{
    std::array<rational, 2> const a{{rational(bg::get<0, 0>(segment)), rational(bg::get<0, 1>(segment))}};
    std::array<rational, 2> const b{{rational(bg::get<1, 0>(segment)), rational(bg::get<1, 1>(segment))}};
    std::array<rational, 2> const low{{rational(bg::get<0, 0>(bounds)), rational(bg::get<0, 1>(bounds))}};
    std::array<rational, 2> const high{{rational(bg::get<1, 0>(bounds)), rational(bg::get<1, 1>(bounds))}};
    rational first = 0, last = 1;
    for (std::size_t d = 0; d != 2; ++d)
    {
        if (low[d] > high[d]) { return true; }
        rational const delta = b[d] - a[d];
        if (delta == 0)
        {
            if (a[d] < low[d] || a[d] > high[d]) { return true; }
            continue;
        }
        rational enter = (low[d] - a[d]) / delta, leave = (high[d] - a[d]) / delta;
        if (enter > leave) { std::swap(enter, leave); }
        if (enter > first) { first = enter; }
        if (leave < last) { last = leave; }
        if (first > last) { return true; }
    }
    return false;
}

template <typename C, typename Segment, typename Box>
void check(Segment const& segment, Box const& bounds)
{
    bool const expected = oracle_disjoint(segment, bounds);
    auto const direct = [&] { return ec::segment_box<C>::apply(segment, bounds); };
    BOOST_TEST_EQ(direct(), expected); // Prime any permitted indexed cache growth.
    ec::filtering_statistics() = {};
    allocations = 0;
    counting = true;
    bool const actual = direct();
    counting = false;
    BOOST_TEST_EQ(actual, expected);
    BOOST_TEST_EQ(allocations, 0u);
    BOOST_TEST_EQ(ec::filtering_statistics().construction_filters, 0u);
    BOOST_TEST_EQ(ec::filtering_statistics().native_rationals, 0u);
    BOOST_TEST_EQ(ec::filtering_statistics().lazy_exact_values, 0u);
    BOOST_TEST_EQ(ec::filtering_statistics().point_exact_values, 0u);
    BOOST_TEST_EQ(ec::filtering_statistics().point_coordinate_values, 0u);
    BOOST_TEST_EQ(ec::filtering_statistics().homogeneous_values, 0u);
    // Mixed internal bounds are not a request to extend the public input domain.
    if constexpr (!ec::static_box_detail::is_search_box<Box>::value
                  && ec::fixed_native_supported<bg::coordinate_type_t<Segment>>::value
                  && ec::fixed_native_supported<bg::coordinate_type_t<Box>>::value
                  && ec::lossless_coordinate_conversion<bg::coordinate_type_t<Segment>, C>::value
                  && ec::lossless_coordinate_conversion<bg::coordinate_type_t<Box>, C>::value)
    {
        bg::strategies::relate::cartesian_exact<C> const strategy;
        BOOST_TEST_EQ(bg::disjoint(segment, bounds, strategy), expected);
        BOOST_TEST_EQ(bg::disjoint(bounds, segment, strategy), expected);
    }
}

template <typename C>
void test_native()
{
    ec::point_cache_scope<C> operation;
    using segment = bg::model::segment<point<C>>;
    ec::search_box const unit({0, 0}, {1, 1});
    box<C> const native_unit({0, 0}, {1, 1});
    std::array<std::array<int, 4>, 12> const cases{{
        {{-1, 0, 2, 0}}, {{0, 0, 1, 1}}, {{0, 3, 3, 0}}, {{1, 2, 2, 1}},
        {{-1, 1, 1, -1}}, {{0, 2, 2, 0}}, {{0, 1, 1, 0}}, {{2, 0, 2, 1}},
        {{0, 0, 0, 0}}, {{1, 1, 1, 1}}, {{2, 0, 2, 0}}, {{-1, -1, 2, 2}}
    }};
    std::array<ec::search_box, 5> const degenerate{{
        {{0.5, 0}, {0.5, 1}}, {{0, 0.5}, {1, 0.5}}, {{0.5, 0.5}, {0.5, 0.5}},
        {{0.5, 0.75}, {0.5, 0.75}}, {{2, 2}, {1, 1}}
    }};
    for (auto const& v : cases)
    {
        segment const s{{C(v[0]), C(v[1])}, {C(v[2]), C(v[3])}};
        check<C>(s, unit);
        check<C>(s, native_unit);
        check<C>(bg::model::referring_segment<point<C> const>(s.second, s.first), unit);
        for (auto const& bounds : degenerate) { check<C>(s, bounds); }
    }
    std::mt19937_64 random(20261115);
    for (int i = 0; i != 1024; ++i)
    {
        auto const coordinate = [&] { return C(int(random() % 33) - 16); };
        point<C> const a(coordinate(), coordinate()), b(coordinate(), coordinate());
        double const x0 = (int(random() % 65) - 32) / 4.0, x1 = (int(random() % 65) - 32) / 4.0;
        double const y0 = (int(random() % 65) - 32) / 4.0, y1 = (int(random() % 65) - 32) / 4.0;
        check<C>(segment(a, b), ec::search_box({(std::min)(x0, x1), (std::min)(y0, y1)},
                                             {(std::max)(x0, x1), (std::max)(y0, y1)}));
    }
}

void test_mixed_extremes()
{
    using C = std::int64_t;
    ec::point_cache_scope<C> operation;
    ec::point_cache_scope<double> double_operation;
    using segment = bg::model::segment<point<C>>;
    C const n = 9007199254740992LL;
    C const low = (std::numeric_limits<C>::min)(), high = (std::numeric_limits<C>::max)();
    segment const above(point<C>(n + 1, n + 1), point<C>(n + 1, n + 1));
    ec::search_box const rounded({double(n), double(n)}, {double(n), double(n)});
    BOOST_TEST(ec::segment_box<C>::apply(above, rounded));
    check<C>(above, rounded);
    check<C>(above, ec::search_box({double(n), double(n)}, {double(n + 2), double(n + 2)}));
    segment const endpoint(point<C>(high, high), point<C>(high, high));
    double const outside = std::ldexp(1.0, 63), below = std::nextafter(outside, 0.0);
    check<C>(endpoint, ec::search_box({outside, outside}, {outside, outside}));
    check<C>(endpoint, ec::search_box({below, below}, {outside, outside}));
    ec::search_box const origin({0, 0}, {0, 0});
    check<C>(segment(point<C>(low, low), point<C>(high, high - 1)), origin);
    check<C>(segment(point<C>(low, low), point<C>(high, high)), origin);
    box<C> const integer_box(point<C>(high - 4, high - 4), point<C>(high, high));
    check<double>(bg::model::segment<point<double>>({outside, outside}, {outside, outside}), integer_box);

    ec::point_cache_scope<float> float_operation;
    double const next = 1.0 + std::ldexp(1.0, -40);
    check<float>(bg::model::segment<point<float>>({1, 0}, {1, 1}), ec::search_box({next, 0}, {2, 1}));
    std::mt19937_64 random(515);
    for (int i = 0; i != 256; ++i)
    {
        auto const draw = [&] { return static_cast<C>(random()); };
        segment const s(point<C>(draw(), draw()), point<C>(draw(), draw()));
        double x0 = double(draw()), x1 = double(draw()), y0 = double(draw()), y1 = double(draw());
        check<C>(s, ec::search_box({(std::min)(x0, x1), (std::min)(y0, y1)},
                                  {(std::max)(x0, x1), (std::max)(y0, y1)}));
    }
}

template <typename C>
struct fractional_point
{
    using implicit = ec::exact_point<C>;
    std::array<implicit, 4> support;
    implicit value;
    fractional_point(C x, C y)
        : support{{implicit(x, y), implicit(x + 1, y + 1), implicit(x, y + 1), implicit(x + 1, y - 1)}}
        , value(ec::construct_intersection(support[0], support[1], support[2], support[3]))
    {}
};

template <typename C>
void test_implicit()
{
    ec::point_cache_scope<C> operation;
    using implicit = ec::exact_point<C>;
    fractional_point<C> const a(0, 0), b(1, 1), c(0, 1);
    BOOST_TEST(bg::get<0>(a.value) == rational(1) / 3);
    BOOST_TEST(bg::get<1>(a.value) == rational(1) / 3);
    bg::model::referring_segment<implicit const> const diagonal(a.value, b.value), single(a.value, a.value);
    ec::search_box const crossing({0.5, 0.5}, {1, 1}), miss({0.5, 0.75}, {0.5, 0.75});
    check<C>(diagonal, crossing);
    check<C>(diagonal, miss);
    check<C>(single, crossing);
    check<C>(single, ec::search_box({1.0 / 3, 1.0 / 3}, {1.0 / 3, 1.0 / 3}));
    check<C>(bg::model::referring_segment<implicit const>(a.value, c.value), crossing);
    // Tuple/GC cleanup traverses boundary_view with segment_iterator, whose
    // value type is pointing_segment rather than referring_segment.
    bg::model::pointing_segment<implicit const> const iterated(a.value, b.value);
    check<C>(iterated, crossing);
    check<C>(iterated, miss);
    check<C>(bg::model::pointing_segment<implicit const>(a.value, a.value), crossing);
    bg::model::box<implicit> const exact(a.value, b.value);
    check<C>(diagonal, exact);
    check<C>(iterated, exact);
    check<C>(single, exact);
    check<C>(bg::model::segment<point<C>>({0, 2}, {2, 0}), exact);
    check<C>(bg::model::segment<point<C>>({0, 3}, {3, 0}), exact);
    for (int x = -3; x <= 3; ++x)
    {
        fractional_point<C> const p(C(x), 0), q(C(x + 1), 1);
        check<C>(bg::model::referring_segment<implicit const>(p.value, q.value), crossing);
        check<C>(bg::model::referring_segment<implicit const>(p.value, q.value), exact);
    }
}

void test_bounded_literal_operands()
{
    ec::point_cache_scope<double> operation;
    rational const third = rational(1) / 3;
    bg::model::segment<point<rational>> const extracted({third, third}, {1 + third, 1 + third});
    std::array<ec::search_box, 3> const bounds{{
        {{0.5, 0.5}, {1, 1}}, {{0.5, 0.75}, {0.5, 0.75}}, {{-1, -1}, {0, 0}}
    }};
    // Internal indexed adapters may expose already exact rational coordinates.
    // Materialization uses the known kernel capacity, never exact_point<rational>.
    for (auto const& b : bounds)
        BOOST_TEST_EQ(ec::segment_box<double>::apply(extracted, b), oracle_disjoint(extracted, b));
    box<rational> const rational_bounds({third, third}, {1 + third, 1 + third});
    BOOST_TEST_EQ(ec::segment_box<double>::apply(extracted, rational_bounds),
                  oracle_disjoint(extracted, rational_bounds));
    bg::model::segment<point<double>> const native({0, 3}, {3, 0});
    BOOST_TEST_EQ(ec::segment_box<double>::apply(native, rational_bounds),
                  oracle_disjoint(native, rational_bounds));

    ec::point_cache_scope<std::int64_t> integer_operation;
    using implicit = ec::exact_point<std::int64_t>;
    std::int64_t const high = (std::numeric_limits<std::int64_t>::max)();
    fractional_point<std::int64_t> const near_high(high - 4, high - 4);
    bg::model::pointing_segment<implicit const> const single(near_high.value, near_high.value);
    double const outside = std::ldexp(1.0, 63), below = std::nextafter(outside, 0.0);
    check<std::int64_t>(single, ec::search_box({outside, outside}, {outside, outside}));
    check<std::int64_t>(single, ec::search_box({below, below}, {outside, outside}));
}

void test_borrowed_adaptors()
{
    std::array<int, 4> const segment_values{{0, 3, 3, 0}};
    std::array<double, 4> const bounds_values{{0, 0, 1, 1}};
    unsigned segment_reads = 0, box_reads = 0;
    indexed_geometry<int, bg::segment_tag> const segment(segment_values, segment_reads);
    indexed_geometry<double, bg::box_tag> const bounds(bounds_values, box_reads);
    auto const a = ec::static_box_detail::borrowed_corner<0>(segment);
    ec::static_box_detail::indexed_point<decltype(bounds), 1, 0> const corner{bounds};
    static_assert(sizeof(corner) == sizeof(void const*), "Corners only borrow their indexed geometry");
    BOOST_TEST(std::addressof(a.source) == std::addressof(segment));
    BOOST_TEST(std::addressof(corner.source) == std::addressof(bounds));
    BOOST_TEST_EQ(segment_reads + box_reads, 0u);
    BOOST_TEST_EQ(bg::get<0>(corner), 1.0);
    BOOST_TEST_EQ(box_reads, 1u);
    check<int>(segment, bounds);
}

void test_unbounded_enclosures()
{
    double const infinity = (std::numeric_limits<double>::infinity)();
    std::array<point<int>, 9> const points{{{-3, -3}, {-3, 0}, {-3, 3},
        {0, -3}, {0, 0}, {0, 3}, {3, -3}, {3, 0}, {3, 3}}};
    for (double xlow : {-infinity, 0.0})
    for (double xhigh : {0.0, infinity})
    for (double ylow : {-infinity, 0.0})
    for (double yhigh : {0.0, infinity})
    {
        ec::search_box const bounds({xlow, ylow}, {xhigh, yhigh});
        // Every segment lies inside [-4,4]^2, so clipping an unbounded box
        // there gives an independent finite rational reference for the query.
        box<int> const finite({std::isinf(xlow) ? -4 : 0, std::isinf(ylow) ? -4 : 0},
                              {std::isinf(xhigh) ? 4 : 0, std::isinf(yhigh) ? 4 : 0});
        for (auto const& a : points)
        for (auto const& b : points)
        {
            bg::model::referring_segment<point<int> const> const segment(a, b);
            BOOST_TEST_EQ(ec::segment_box<int>::apply(segment, bounds), oracle_disjoint(segment, finite));
        }
    }
}

int main()
{
    test_native<float>();
    test_native<double>();
    test_native<int>();
    test_native<std::int64_t>();
    test_mixed_extremes();
    test_implicit<float>();
    test_implicit<double>();
    test_implicit<std::int64_t>();
    test_borrowed_adaptors();
    test_bounded_literal_operands();
    test_unbounded_enclosures();
    return boost::report_errors();
}
