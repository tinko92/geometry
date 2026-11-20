// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.

#define BOOST_GEOMETRY_EXACT_FILTER_STATISTICS
#include <boost/core/lightweight_test.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/strategies/detail/operation_binding.hpp>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <new>
#include <type_traits>
#include <utility>

#if __cplusplus >= 201703L
#include <boost/geometry/strategies/relate/cartesian_exact.hpp>
#endif

namespace bg = boost::geometry;
namespace sd = bg::strategies::detail;
namespace
{
bool count_allocations = false;
std::size_t allocations = 0;
}
void* operator new(std::size_t size)
{
    if (count_allocations) ++allocations;
    if (void* p = std::malloc(size ? size : 1)) return p;
    throw std::bad_alloc();
}
void operator delete(void* p) noexcept { std::free(p); }
void* operator new[](std::size_t size) { return ::operator new(size); }
void operator delete[](void* p) noexcept { std::free(p); }

#ifdef BOOST_GEOMETRY_BORROWED_TEST_INT64
using coordinate = std::int64_t;
#else
using coordinate = double;
#endif
using point = bg::model::d2::point_xy<coordinate>;
using line = bg::model::linestring<point>;
using polygon = bg::model::polygon<point>;
using multi_point = bg::model::multi_point<point>;
using segment = bg::model::segment<point>;

struct counters { unsigned entries = 0; };
struct bound_probe : bg::strategies::relate::cartesian<>
{
    template <typename A, typename B, typename Apply>
    decltype(auto) with_sources(A const&, B const&, Apply&& apply) const
    { return std::forward<Apply>(apply)(*this); }
};
struct umbrella_probe : bg::strategies::relate::cartesian<>
{
    counters& stats;
    explicit umbrella_probe(counters& stats) : stats(stats) {}
    template <typename A, typename B, typename Apply>
    decltype(auto) with_sources(A const&, B const&, Apply&& apply) const
    {
        ++stats.entries;
        return std::forward<Apply>(apply)(bound_probe{});
    }
};

template <typename Apply>
void check_call(counters& stats, bool expected, Apply apply)
{
    stats.entries = 0;
    BOOST_TEST_EQ(apply(), expected);
    BOOST_TEST_EQ(stats.entries, 1u);
}

template <typename A, typename B, typename Strategy>
void check_all(A const& a, B const& b, Strategy const& strategy, counters& stats,
               bool disjoint, bool touches, bool within, bool covered, bool equals, bool overlaps, bool crosses)
{
    check_call(stats, disjoint, [&]{ return bg::disjoint(a, b, strategy); });
    check_call(stats, !disjoint, [&]{ return bg::intersects(a, b, strategy); });
    check_call(stats, touches, [&]{ return bg::touches(a, b, strategy); });
    check_call(stats, within, [&]{ return bg::within(a, b, strategy); });
    check_call(stats, covered, [&]{ return bg::covered_by(a, b, strategy); });
    check_call(stats, equals, [&]{ return bg::equals(a, b, strategy); });
    check_call(stats, overlaps, [&]{ return bg::overlaps(a, b, strategy); });
    check_call(stats, crosses, [&]{ return bg::crosses(a, b, strategy); });
}

template <typename Strategy>
void matrix(Strategy const& strategy, counters& stats)
{
    line const a{{0, 0}, {3, 3}}, b{{0, 3}, {3, 0}};
    polygon p, q;
    p.outer() = {{0, 0}, {0, 4}, {4, 4}, {4, 0}, {0, 0}};
    q.outer() = {{2, -1}, {2, 3}, {6, 3}, {6, -1}, {2, -1}};
    check_all(a, b, strategy, stats, false, false, false, false, false, false, true);
    check_all(p, q, strategy, stats, false, false, false, false, false, true, false);
    line const crossing{{-1, 2}, {5, 2}};
    check_all(crossing, p, strategy, stats, false, false, false, false, false, false, true);
    multi_point const mp1{{0, 0}, {1, 1}}, mp2{{1, 1}, {2, 2}};
    check_all(mp1, mp2, strategy, stats, false, false, false, false, false, true, false);
    multi_point const mixed{{1, 1}, {10, 10}};
    check_all(mixed, p, strategy, stats, false, false, false, false, false, false, true);
    check_all(p.outer(), q.outer(), strategy, stats, false, false, false, false, false, true, false);
    segment const first(a[0], a[1]), second(b[0], b[1]);
    check_call(stats, false, [&]{ return bg::disjoint(first, second, strategy); });
    check_call(stats, true, [&]{ return bg::intersects(first, second, strategy); });
    check_call(stats, false, [&]{ return bg::equals(first, second, strategy); });
    segment const cut(crossing[0], crossing[1]);
    check_call(stats, false, [&]{ return bg::disjoint(cut, p.outer(), strategy); });
    check_call(stats, true, [&]{ return bg::intersects(cut, p.outer(), strategy); });
    point const inside(1, 1);
    check_call(stats, true, [&]{ return bg::within(inside, p, strategy); });
    check_call(stats, true, [&]{ return bg::covered_by(inside, p, strategy); });
    bg::model::box<point> const box(point(0, 0), point(4, 4));
    check_call(stats, true, [&]{ return bg::covered_by(a, box, strategy); });
    check_call(stats, false, [&]{ return bg::disjoint(first, box, strategy); });
}

#if __cplusplus >= 201703L
namespace ec = sd::exact_cartesian;
#ifndef BOOST_GEOMETRY_BORROWED_TEST_POLICY
#define BOOST_GEOMETRY_BORROWED_TEST_POLICY 0
#endif
#if BOOST_GEOMETRY_BORROWED_TEST_POLICY == 1
using policy = bg::strategies::relate::nearest_output;
#else
using policy = bg::strategies::relate::approximate_output;
#endif
struct root_probe : bg::strategies::relate::cartesian_exact<coordinate, policy>
{
    using base = bg::strategies::relate::cartesian_exact<coordinate, policy>;
    using base::relate;
    counters& stats;
    explicit root_probe(counters& stats) : stats(stats) {}
    void relate() const = delete;
    template <typename A, typename B, typename Apply>
    decltype(auto) with_sources(A const& a, B const& b, Apply&& apply) const
    {
        ++stats.entries;
        return base::with_sources(a, b, [&](auto const& bound)
        {
            using bound_type = std::decay_t<decltype(bound)>;
            static_assert(std::is_same<typename bound_type::segment_position_type,
                ec::borrowed_segment_position_tag>::value, "Native predicates must bind the borrowed Root");
            return std::forward<Apply>(apply)(bound);
        });
    }
};

template <typename Apply>
std::size_t count_new(Apply apply)
{
    allocations = 0;
    count_allocations = true;
    try { apply(); }
    catch (...) { count_allocations = false; throw; }
    count_allocations = false;
    return allocations;
}

void test_root()
{
    counters stats;
    root_probe strategy(stats);
    ec::filtering_statistics() = {};
    matrix(strategy, stats);
    BOOST_TEST_EQ(ec::filtering_statistics().native_point_constructions, 0u);

    line const a{{0, 0}, {3, 3}}, b{{0, 3}, {3, 0}};
    stats.entries = 0;
    strategy.with_sources(a, b, [&](auto const& bound)
    {
        auto* const cache = ec::operation_detail::active<coordinate>();
        sd::with_operation_sources(bound, a, b, [&](auto const& reused)
        {
            BOOST_TEST(&reused == &bound);
            BOOST_TEST(ec::operation_detail::active<coordinate>() == cache);
            BOOST_TEST(bg::intersects(a, b, reused));
            BOOST_TEST(ec::operation_detail::active<coordinate>() == cache);
        });
    });
    BOOST_TEST_EQ(stats.entries, 1u);

    point const p(1, 1), q(2, 2);
    bg::model::box<point> const box(point(0, 0), point(4, 4));
    segment const first(point(0, 0), point(3, 3)), second(point(0, 3), point(3, 0));
    std::size_t const pure = count_new([&]
    {
        BOOST_TEST(!bg::equals(p, q, strategy));
        BOOST_TEST(bg::within(p, box, strategy));
        BOOST_TEST(bg::covered_by(p, first, strategy));
        BOOST_TEST(!bg::disjoint(first, second, strategy));
    });
    BOOST_TEST_EQ(pure, 0u);
    line const left{{0, 0}, {1, 1}}, right{{10, 10}, {11, 11}};
    std::size_t const ordinary_work = count_new([&]
    { BOOST_TEST(bg::disjoint(left, right, bg::strategies::relate::cartesian<>())); });
    std::size_t const borrowed_work = count_new([&]
    { BOOST_TEST(bg::disjoint(left, right, strategy)); });
    std::printf("pure_native_allocations=%zu ordinary_disjoint_work=%zu borrowed_disjoint_work=%zu\n",
                pure, ordinary_work, borrowed_work);
    BOOST_TEST(borrowed_work <= ordinary_work);
    BOOST_TEST_EQ(ec::filtering_statistics().native_point_constructions, 0u);
}
#endif

int main()
{
    counters stats;
    matrix(umbrella_probe(stats), stats);
#if __cplusplus >= 201703L
    test_root();
#endif
    return boost::report_errors();
}
