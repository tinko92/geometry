// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.

// clang++ -std=c++17 -O1 -DNDEBUG -ffp-contract=off -Iinclude -I../..
//   test/strategies/exact_borrowed_predicates.cpp -o /tmp/bg-borrowed-predicates
#define BOOST_GEOMETRY_EXACT_FILTER_STATISTICS
#include <boost/geometry/strategies/cartesian/exact/borrowed_predicates.hpp>
#include <boost/geometry/strategies/cartesian/exact/borrowed_query.hpp>
#include <boost/geometry/strategies/cartesian/exact/static_predicates.hpp>
#include <boost/geometry/strategies/relate/cartesian_exact.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/views/detail/boundary_view.hpp>
#include <boost/core/lightweight_test.hpp>

#include <algorithm>
#include <cfenv>
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <new>
#include <random>
#include <vector>

namespace bg = boost::geometry;
namespace ec = bg::strategies::detail::exact_cartesian;
namespace bp = ec::borrowed_predicates_detail;
namespace policy = bg::strategies::relate;
using q = boost::multiprecision::cpp_rational;
using qp = std::array<q, 2>;

namespace
{
bool counting = false;
std::size_t allocations = 0, triples = 0, comparisons = 0, allocation_checks = 0, total_allocations = 0;
}
void* operator new(std::size_t n)
{
    if (counting) ++allocations;
    if (void* p = std::malloc(n == 0 ? 1 : n)) return p;
    throw std::bad_alloc();
}
void operator delete(void* p) noexcept { std::free(p); }
void* operator new[](std::size_t n) { return ::operator new(n); }
void operator delete[](void* p) noexcept { std::free(p); }
#if defined(__cpp_sized_deallocation)
void operator delete(void* p, std::size_t) noexcept { std::free(p); }
void operator delete[](void* p, std::size_t) noexcept { std::free(p); }
#endif

template <typename F> void no_allocations(F function)
{
    allocations = 0; counting = true;
    try { function(); } catch (...) { counting = false; throw; }
    counting = false;
    ++allocation_checks; total_allocations += allocations;
    BOOST_TEST_EQ(allocations, 0u);
}
int sign(q const& value) { return int(value > 0) - int(value < 0); }
q absolute(q const& value) { return value < 0 ? q(-value) : value; }
q det(qp const& a, qp const& b, qp const& p)
{ return (b[0] - a[0]) * (p[1] - a[1]) - (b[1] - a[1]) * (p[0] - a[0]); }
q dot(qp const& a, qp const& b, qp const& p)
{ return (b[0] - a[0]) * (p[0] - b[0]) + (b[1] - a[1]) * (p[1] - b[1]); }
q distance(qp const& a, qp const& b)
{ q const x = b[0] - a[0], y = b[1] - a[1]; return x * x + y * y; }

template <typename C> qp intersection(std::array<C, 8> const& v)
{
    std::array<q, 8> p;
    for (std::size_t i = 0; i < 8; ++i) p[i] = q(v[i]);
    q const a = p[1] - p[3], b = p[2] - p[0], c = p[0] * p[3] - p[1] * p[2];
    q const d = p[5] - p[7], e = p[6] - p[4], f = p[4] * p[7] - p[5] * p[6];
    q const w = a * e - b * d;
    return {{(b * f - c * e) / w, (c * d - a * f) / w}};
}

template <typename C, typename Policy>
struct fixture
{
    using native = bg::model::d2::point_xy<C>;
    using line = bg::model::linestring<native>;
    using root = ec::borrowed_root<C, Policy, line, line>;
    using point = ec::borrowed_point<root>;
    line a, b;
    ec::borrowed_context<root> context;
    ec::borrowed_stage<root> stage;
    std::array<typename root::operation_type, 2> operations;
    ec::borrowed_record<root> record;
    point p;
    qp exact;

    explicit fixture(std::array<C, 8> const& v)
        : a{{v[0], v[1]}, {v[2], v[3]}}, b{{v[4], v[5]}, {v[6], v[7]}},
          stage(context, a, b), p(record), exact(intersection(v))
    {
        operations[0].seg_id = bg::segment_identifier(0, -1, -1, 0);
        operations[1].seg_id = bg::segment_identifier(1, -1, -1, 0);
        record.stage = &stage;
        record.operations = {{&operations[0], &operations[1]}};
        record.numeric.init(p.producers());
    }
};

template <typename Root, typename A, typename B, typename P>
void check(A const& a, B const& b, P const& p, qp const& x, qp const& y, qp const& z)
{
    auto const side = policy::cartesian_exact<typename Root::coordinate_type, typename Root::output_policy>{}.side();
    BOOST_TEST_EQ(side.apply(a, b, p), sign(det(x, y, z)));
    BOOST_TEST_EQ(side.direction_code(a, b, p), sign(dot(x, y, z)));
    BOOST_TEST_EQ((ec::root_compare_coordinate<0, Root>(a, b)), sign(x[0] - y[0]));
    BOOST_TEST_EQ((ec::root_compare_coordinate<1, Root>(a, b)), sign(x[1] - y[1]));
    BOOST_TEST_EQ(ec::root_compare_coordinate<Root>(b, a, 0), sign(y[0] - x[0]));
    BOOST_TEST_EQ(ec::root_equal<Root>(a, b), x == y);
    BOOST_TEST(ec::root_equal<Root>(p, p));
    ++triples;
}

struct counted_native_iterator
{
    using iterator_category = std::forward_iterator_tag;
    using value_type = bg::model::d2::point_xy<double>;
    using difference_type = std::ptrdiff_t;
    using pointer = value_type const*;
    using reference = value_type const&;

    pointer point;
    unsigned* reads;

    reference operator*() const { ++*reads; return *point; }
};

template <typename Policy>
void stored_filter_reads_once()
{
    using F = fixture<double, Policy>;
    using P = typename F::native;
    auto const side = policy::cartesian_exact<double, Policy>{}.side();
    // An exact x after refinement, but y remains uncertain (7/3). The filter
    // must use its stored error rather than treating the record as native.
    F f({{0, 0, 0, 5, -2, 1, 1, 3}});
    f.p.refine(0);
    BOOST_TEST_EQ(f.p.approximation(0).error, 0.0);
    BOOST_TEST(f.p.approximation(1).error > 0.0);
    unsigned reads_a = 0, reads_b = 0;
    P const pa(-2, 0), pb(-1, 0);
    ec::iterator_point_view<counted_native_iterator> const
        a(counted_native_iterator{&pa, &reads_a}), b(counted_native_iterator{&pb, &reads_b});
    no_allocations([&]
    {
        BOOST_TEST_EQ(side.apply(a, b, f.p), 1);
        BOOST_TEST_EQ(reads_a, 1u);
        BOOST_TEST_EQ(reads_b, 1u);
        reads_a = reads_b = 0;
        BOOST_TEST_EQ(side.direction_code(a, b, f.p), 1);
        BOOST_TEST_EQ(reads_a, 1u);
        BOOST_TEST_EQ(reads_b, 1u);
    });

    double const rounded_y = 7.0 / 3.0;
    P const c(-2, rounded_y), d(0, rounded_y), e(0, rounded_y - 1);
    BOOST_TEST_EQ(side.apply(c, d, f.p), sign(det({{-2, q(rounded_y)}}, {{0, q(rounded_y)}}, f.exact)));
    BOOST_TEST_EQ(side.direction_code(e, d, f.p), sign(dot({{0, q(rounded_y - 1)}}, {{0, q(rounded_y)}}, f.exact)));

    // Exactly (2, 0). Full zero/cancellation cases must still agree with the
    // underlying topology even if they enter a different numerical bound.
    F exact({{2, -4, 2, 4, -4, 0, 4, 0}});
    exact.p.refine(0); exact.p.refine(1);
    BOOST_TEST_EQ(side.apply(pa, pb, exact.p), 0);
    BOOST_TEST_EQ(side.direction_code(pa, pb, exact.p), 1);
    BOOST_TEST_EQ(side.direction_code(P(2, 1), P(2, 0), exact.p), 0);

    // Translate a crossing to the origin to exercise two certified zero axes.
    F origin({{0, -4, 0, 4, -4, 0, 4, 0}});
    origin.p.refine(0); origin.p.refine(1);
    BOOST_TEST_EQ(origin.p.approximation(0).error, 0.0);
    BOOST_TEST_EQ(origin.p.approximation(1).error, 0.0);
    BOOST_TEST_EQ(side.apply(P(1, 1), P(2, 1), origin.p), -1);
    BOOST_TEST_EQ(side.direction_code(P(1, 1), P(2, 1), origin.p), -1);
}

template <typename Root>
void permutations(std::array<ec::borrowed_point<Root>, 3> const& points, std::array<qp, 3> const& exact)
{
    std::array<unsigned, 3> order{{0, 1, 2}};
    do { check<Root>(points[order[0]], points[order[1]], points[order[2]],
                    exact[order[0]], exact[order[1]], exact[order[2]]); }
    while (std::next_permutation(order.begin(), order.end()));
}

template <typename Value> void bound(Value const& value, q const& exact)
{
    auto const a = value.bounds();
    if (std::isfinite(a.error)) BOOST_TEST(absolute(q(a.approximate) - exact) <= q(a.error));
}

template <typename A, typename B> void compare_measure(A const& a, B const& b, q const& x, q const& y)
{
    int const expected = sign(x - y);
    BOOST_TEST_EQ(ec::root_compare_distance(a, b), expected);
    BOOST_TEST_EQ(ec::compare_distance(a, b), expected);
    BOOST_TEST_EQ(a < b, expected < 0); BOOST_TEST_EQ(a <= b, expected <= 0);
    BOOST_TEST_EQ(a > b, expected > 0); BOOST_TEST_EQ(a >= b, expected >= 0);
    BOOST_TEST_EQ(a == b, expected == 0); BOOST_TEST_EQ(a != b, expected != 0);
    ++comparisons;
}

template <typename C, typename Policy> void cases(C base)
{
    using F = fixture<C, Policy>;
    using Root = typename F::root;
    using Native = typename F::native;
    using Point = typename F::point;
    F first({{base, base, C(base + 3), C(base + 3), base, C(base + 2), C(base + 3), base}});
    F second({{base, base, C(base + 3), C(base + 3), base, C(base + 4), C(base + 3), C(base + 1)}});
    F third({{base, C(base + 3), C(base + 4), C(base + 3), C(base + 1), base, C(base + 1), C(base + 5)}});
    Native const a(base, base), b(C(base + 3), C(base + 3)), c(base, C(base + 4));
    Native const right(C(base + 4), C(base + 2));
    qp const qa{{q(base), q(base)}}, qb{{q(base) + 3, q(base) + 3}}, qc{{q(base), q(base) + 4}};
    auto const pa = Point::borrow(a), pb = Point::borrow(b), pc = Point::borrow(c);
    permutations<Root>({{pa, pb, pc}}, {{qa, qb, qc}});
    permutations<Root>({{pa, pb, first.p}}, {{qa, qb, first.exact}}); // Original support.
    permutations<Root>({{pa, pb, third.p}}, {{qa, qb, third.exact}}); // Also when wide N casts coincide.
    permutations<Root>({{pa, first.p, second.p}}, {{qa, first.exact, second.exact}});
    permutations<Root>({{first.p, second.p, third.p}}, {{first.exact, second.exact, third.exact}});
    check<Root>(a, b, first.p, qa, qb, first.exact); // Raw natives too.
    check<Root>(a, b, right, qa, qb, {{q(base) + 4, q(base) + 2}}); // Noncollinear, dot == 0.
    auto const iterator = Point::from_iterator(first.a.cbegin());
    check<Root>(iterator, b, first.p, qa, qb, first.exact);
    bg::model::d2::point_xy<std::int32_t> const small(0, 1);
    check<Root>(small, first.p, second.p, {{0, 1}}, first.exact, second.exact);

    auto const d1 = ec::root_squared_distance<Root>(first.p, a);
    auto const d2 = ec::root_squared_distance<Root>(first.p, second.p);
    auto const s1 = ec::root_signed_determinant<Root>(a, b, third.p);
    auto const s2 = ec::root_distance_measure<Root>(first.p, second.p, third.p);
    static_assert(std::is_same<decltype(d1.a), Point>::value, "Descriptors store borrowed references");
    static_assert(std::is_same<decltype(s1.p), Point>::value, "No owning exact_point descriptor");
    q const ed1 = distance(first.exact, qa), ed2 = distance(first.exact, second.exact);
    q const es1 = det(qa, qb, third.exact), es2 = det(first.exact, second.exact, third.exact);
    bound(d1, ed1); bound(d2, ed2); bound(s1, es1); bound(s2, es2);
    compare_measure(d1, d2, ed1, ed2); compare_measure(d1, d1, ed1, ed1);
    compare_measure(s1, s2, es1, es2); compare_measure(s1, -s1, es1, -es1);
    compare_measure(d1, s2, ed1, es2); compare_measure(s1, d2, es1, ed2);
    auto negative = -s1;
    bound(negative, -es1); negative.bounds(true); bound(negative, -es1);
    BOOST_TEST_EQ(s1.is_negative(), es1 < 0); BOOST_TEST_EQ(s1.is_positive(), es1 > 0);
    BOOST_TEST_EQ(s2.is_zero(), es2 == 0);
    BOOST_TEST_EQ(d1 > 0, ed1 > 0); BOOST_TEST_EQ(0 <= d2, q(0) <= ed2);
    ec::RootSquaredDistanceValue<Root> zero, unset(-1);
    ec::RootSignedMeasure<Root> signed_zero;
    compare_measure(zero, signed_zero, 0, 0); compare_measure(unset, d1, -1, ed1);
    BOOST_TEST(zero == 0); BOOST_TEST(signed_zero == 0); BOOST_TEST(unset < 0);

    // Query expressions have no provenance/record and use their actual larger
    // homogeneous format. Also cover the legacy auxiliary midpoint fallback.
    auto const midpoint = ec::query_point<Root>::midpoint(first.p, second.p);
    qp const qm{{(first.exact[0] + second.exact[0]) / 2, (first.exact[1] + second.exact[1]) / 2}};
    check<Root>(first.p, midpoint, second.p, first.exact, qm, second.exact);
    check<Root>(a, midpoint, third.p, qa, qm, third.exact);
    auto const legacy_a = ec::make_point<C>(a), legacy_b = ec::make_point<C>(b);
    auto const legacy = ec::exact_point<C>::from_homogeneous(ec::fixed_midpoint(legacy_a.fixed_homogeneous(), legacy_b.fixed_homogeneous()));
    qp const qlegacy{{(qa[0] + qb[0]) / 2, (qa[1] + qb[1]) / 2}};
    check<Root>(legacy, first.p, midpoint, qlegacy, first.exact, qm);
    // Existing strategy bridge should dispatch before accessing diagnostic get.
    BOOST_TEST_EQ(ec::static_orientation(a, midpoint, third.p), sign(det(qa, qm, third.exact)));
    BOOST_TEST_EQ(ec::static_direction(first.p, midpoint, second.p), sign(dot(first.exact, qm, second.exact)));
    BOOST_TEST_EQ(ec::static_equal_points(first.p, first.p), true);

    // All operation caches and source registrations have been warmed. Predicates
    // and descriptors must allocate neither owning points nor numerical storage.
    for (std::size_t d = 0; d < 2; ++d) { first.p.refine(d); second.p.refine(d); third.p.refine(d); }
    first.p.homogeneous(); second.p.homogeneous(); third.p.homogeneous();
    d1.fixed_pair(); d2.fixed_pair(); s1.fixed_pair(); s2.fixed_pair();
    int const expected_side = sign(es2), expected_measure = sign(es1 - es2);
    no_allocations([&]
    {
        for (unsigned i = 0; i < 8; ++i)
        {
            BOOST_TEST_EQ(ec::root_orientation<Root>(pa, pb, first.p), 0);
            BOOST_TEST_EQ(ec::root_orientation<Root>(first.p, second.p, third.p), expected_side);
            BOOST_TEST_EQ((ec::root_compare_coordinate<0, Root>(first.p, first.p)), 0);
            BOOST_TEST(ec::root_equal<Root>(first.p, first.p));
            BOOST_TEST_EQ(ec::root_direction<Root>(pa, first.p, second.p), 1);
            BOOST_TEST_EQ(ec::root_compare_distance(d1, d1), 0);
            BOOST_TEST_EQ(ec::root_compare_distance(s1, s2), expected_measure);
            auto const borrowed = ec::root_squared_distance<Root>(a, first.p);
            BOOST_TEST_EQ(ec::root_compare_distance(borrowed, d1), 0);
        }
    });
}

template <typename Policy> void random_cases()
{
    using F = fixture<double, Policy>;
    using Root = typename F::root;
    std::mt19937_64 random(0x424f52524f574544);
    for (unsigned i = 0; i < 32; ++i)
    {
        double const x = double(random() % 10000) / 64, y = double(random() % 10000) / 64;
        double const r = 1 + double(random() % 64) / 8;
        F a({{x - r, y - r, x + r, y + r + 0.25, x - r, y + r, x + r, y - r - 0.5}});
        F b({{x - r, y + 1, x + r, y + r, x, y - r, x + 1, y + r}});
        F c({{x, y - r, x, y + r, x - r, y, x + r, y}});
        permutations<Root>({{a.p, b.p, c.p}}, {{a.exact, b.exact, c.exact}});
        auto const d1 = ec::root_squared_distance<Root>(a.p, b.p), d2 = ec::root_squared_distance<Root>(a.p, c.p);
        auto const s = ec::root_signed_determinant<Root>(a.p, b.p, c.p);
        q const es = det(a.exact, b.exact, c.exact);
        bound(d1, distance(a.exact, b.exact)); bound(s, es);
        compare_measure(d1, d2, distance(a.exact, b.exact), distance(a.exact, c.exact));
        BOOST_TEST_EQ(ec::root_compare_distance(s, 0), sign(es));
    }
}

void indexed_native()
{
    using P = bg::model::d2::point_xy<double>;
    using Box = bg::model::box<P>;
    using Segment = bg::model::segment<P>;
    using Root = ec::borrowed_root<double, policy::deferred_nearest_output, Box, Segment>;
    Box const box(P(0, 0), P(3, 5));
    Segment const segment(P(0, 5), P(3, 0));
    auto const a = ec::borrowed_point<Root>::borrow(ec::borrowed_detail::indexed_vertex<Box>{&box, 0});
    auto const b = ec::borrowed_point<Root>::borrow(ec::borrowed_detail::indexed_vertex<Box>{&box, 2});
    auto const p = ec::borrowed_point<Root>::borrow(ec::borrowed_detail::indexed_vertex<Segment>{&segment, 0});
    no_allocations([&] { BOOST_TEST_EQ(ec::root_orientation<Root>(a, b, p), 1); });
    BOOST_TEST(ec::operation_detail::current<double>() == nullptr);
}

// A query's approximate() can be used without touching its interval accessor.
// Deliberately no Point traits, record(), producers(), or native-coordinate API.
template <typename Root> struct probe_query
{
    using root_type = Root;
    using C = typename Root::coordinate_type;
    decltype(ec::fixed_input_homogeneous(C(0), C(0))) value;
    mutable unsigned interval_reads = 0, refinement_reads = 0;
    double approximate(std::size_t d) const { return d == 0 ? 2 : 3; }
    ec::interval approximation(std::size_t d) const { ++interval_reads; return {approximate(d), 0}; }
    ec::interval refine(std::size_t d) const { ++refinement_reads; return {approximate(d), 0}; }
    auto homogeneous() const { return value; }
};

void lean_does_not_reconstruct_errors()
{
    using F = fixture<double, policy::deferred_nearest_output>;
    using Root = F::root;
    F::native const a(0, 0), b(1, 0);
    F::native const aligned(2, 0);
    probe_query<Root> query{ec::fixed_input_homogeneous(2.0, 3.0)};
    no_allocations([&]
    {
        BOOST_TEST_EQ(ec::root_orientation<Root>(a, b, query), 1);
        BOOST_TEST_EQ((ec::root_compare_coordinate<0, Root>(a, query)), -1);
        BOOST_TEST_EQ(ec::root_direction<Root>(a, b, query), 1);
        BOOST_TEST(!ec::root_equal<Root>(aligned, query));
    });
    BOOST_TEST_EQ(query.interval_reads, 0u);
    BOOST_TEST_EQ(query.refinement_reads, 0u);
}

void coordinate_filter_reads_one_axis()
{
    using F = fixture<double, policy::approximate_output>;
    using Root = F::root;
    F::native const a(0, 0);
    probe_query<Root> query{ec::fixed_input_homogeneous(2.0, 3.0)};
    no_allocations([&]
    {
        BOOST_TEST_EQ((ec::root_compare_coordinate<0, Root>(a, query)), -1);
        BOOST_TEST_EQ(query.interval_reads, 1u);
        BOOST_TEST_EQ((ec::root_compare_coordinate<1, Root>(a, query)), -1);
        BOOST_TEST_EQ(query.interval_reads, 2u);
    });
    BOOST_TEST_EQ(query.refinement_reads, 0u);
}

template <typename Policy> void polygon_line_audit(char const* label)
{
    using Native = bg::model::d2::point_xy<double>;
    using Polygon = bg::model::polygon<Native, true, false>;
    using Line = bg::model::linestring<Native>;
    using Root = ec::borrowed_root<double, Policy, Polygon, Line>;
    using Point = ec::borrowed_point<Root>;
    Polygon polygon;
    polygon.outer() = {{0, 0}, {4, 0}, {4, 4}, {0, 4}};
    Line const line{{1, 1}, {2, 3}};
    auto const native = Point::from_iterator(line.cbegin());
    auto const descriptor = [&]
    {
        using Ring = std::decay_t<decltype(polygon.outer())>;
        bg::detail::boundary_view<Ring const> boundary(polygon.outer());
        ec::iterator_point_view<decltype(boundary.begin())> view(boundary.begin());
        static_assert(!boost::mp11::mp_contains<typename Root::native_access, decltype(view)>::value,
                      "Exercise an unregistered legacy boundary accessor");
        return ec::RootSquaredDistanceValue<Root>(view, native);
    }(); // The temporary boundary and iterator accessor have expired.
    BOOST_TEST(descriptor.fixed_pair().N * 1 == 2 * descriptor.fixed_pair().D);
    auto const reverse = ec::iterator_point_view<decltype(polygon.outer().crbegin())>(polygon.outer().crbegin());
    ec::RootSignedMeasure<Root> const signed_value(reverse, polygon.outer().front(), native);
    BOOST_TEST_EQ(signed_value.sign(), 1);
    std::vector<Point> const working{native, Point::borrow(polygon.outer().front())};
    auto const working_view = ec::iterator_point_view<decltype(working.cbegin())>(working.cbegin());
    ec::RootSquaredDistanceValue<Root> const working_distance(working_view, working.back());
    no_allocations([&]
    {
        BOOST_TEST_EQ(ec::root_compare_distance(descriptor, working_distance), 0);
        auto const retained = ec::root_squared_distance<Root>(reverse, working_view);
        BOOST_TEST_EQ(retained.sign(), 1);
        BOOST_TEST_EQ(ec::root_orientation<Root>(working_view, polygon.outer().front(), polygon.outer()[1]), 1);
    });
    std::cout << label << " polygon+line: numeric=" << sizeof(typename Root::numeric_state)
              << " record=" << sizeof(ec::borrowed_record<Root>) << " point=" << sizeof(Point)
              << " distance=" << sizeof(ec::RootSquaredDistanceValue<Root>)
              << " signed=" << sizeof(ec::RootSignedMeasure<Root>) << '\n';
}

template <typename Policy> void cache_free_support_zero()
{
    using F = fixture<double, Policy>;
    using Root = typename F::root;
    F sample({{1048576, 1048576, 1048592, 1048608,
               1048576, 1048600, 1048592, 1048584}});
    BOOST_TEST_EQ(ec::point_cache_detail::native_cache<double>::size(), 0u);
    auto const a = F::point::from_iterator(sample.a.cbegin());
    auto const b = F::point::from_iterator(sample.a.cbegin() + 1);
    typename F::native const extended_a(1048560, 1048544), extended_b(1048608, 1048640);
    no_allocations([&]
    {
        for (unsigned i = 0; i < 32; ++i)
        {
            BOOST_TEST_EQ(ec::root_orientation<Root>(a, b, sample.p), 0);
            BOOST_TEST_EQ(ec::root_orientation<Root>(b, sample.p, a), 0);
            BOOST_TEST_EQ(ec::root_orientation<Root>(sample.p, a, b), 0);
            BOOST_TEST_EQ(ec::root_orientation<Root>(extended_a, extended_b, sample.p), 0);
        }
    });
    BOOST_TEST_EQ(ec::point_cache_detail::native_cache<double>::size(), 0u);
    for (auto& operation : sample.operations)
    {
        auto const edge = sample.stage.support(operation.seg_id);
        operation.seg_id = bg::segment_identifier(edge.range, -1, -1, edge.segment);
    }
    sample.stage.finalized = true;
    sample.stage.sources = {};
    no_allocations([&] { BOOST_TEST_EQ(ec::root_orientation<Root>(extended_a, extended_b, sample.p), 0); });
    BOOST_TEST_EQ(ec::point_cache_detail::native_cache<double>::size(), 0u);
}

template <typename C, typename Policy>
void native_measure_fallbacks()
{
    using F = fixture<C, Policy>;
    using Root = typename F::root;
    using P = typename F::native;
    C const limit = C(1073741824);
    P const low(-limit, -limit), high(limit, limit);
    auto const distance = ec::root_squared_distance<Root>(low, high);
    auto const distance_copy = distance;
    ec::filtering_statistics() = {};
    no_allocations([&] { BOOST_TEST_EQ(ec::root_compare_distance(distance, distance_copy), 0); });
    BOOST_TEST_EQ(ec::filtering_statistics().small_integer_measure_comparisons, 1u);
    BOOST_TEST_EQ(ec::filtering_statistics().full_measure_comparisons, 0u);
    std::uint64_t exact;
    BOOST_TEST(bp::small_integer_distance(distance, exact));
    BOOST_TEST_EQ(exact, std::uint64_t(1) << 63);
    P const outside(C(limit + 1), C(limit + 1));
    BOOST_TEST(!bp::small_integer_distance(ec::root_squared_distance<Root>(low, outside), exact));

    // A fractional or large integral native is still handled with exact compact
    // homogeneous arithmetic, never by rounding it into the integer shortcut.
    C const x = std::is_integral<C>::value ? C(std::int64_t(1) << 54) : C(0.25);
    P const a(x, x), b(C(x + 2), C(x + 4));
    auto const compact = ec::root_squared_distance<Root>(a, b);
    auto const compact_copy = compact;
    ec::filtering_statistics() = {};
    no_allocations([&] { BOOST_TEST_EQ(ec::root_compare_distance(compact, compact_copy), 0); });
    BOOST_TEST_EQ(ec::filtering_statistics().compact_measure_comparisons, 1u);
    BOOST_TEST_EQ(ec::filtering_statistics().full_measure_comparisons, 0u);
    typename decltype(compact)::compact_pair_type pair;
    BOOST_TEST(compact.compact_pair(pair));
    BOOST_TEST_EQ(ec::fixed_pair_compare(pair, compact.fixed_pair()), 0);
    P const c(C(x + 1), C(x - 3));
    auto const signed_value = ec::root_signed_determinant<Root>(a, b, c);
    typename decltype(signed_value)::compact_pair_type signed_pair;
    BOOST_TEST(signed_value.compact_pair(signed_pair));
    BOOST_TEST_EQ(ec::fixed_pair_compare(signed_pair, signed_value.fixed_pair()), 0);
    if constexpr (std::is_floating_point<C>::value)
    {
        P const tiny(C(std::ldexp(1.0, -200)), C(std::ldexp(1.0, 200)));
        auto const full = ec::root_squared_distance<Root>(a, tiny);
        auto const full_copy = full;
        BOOST_TEST(!full.compact_pair(pair));
        ec::filtering_statistics() = {};
        no_allocations([&] { BOOST_TEST_EQ(ec::root_compare_distance(full, full_copy), 0); });
        BOOST_TEST_EQ(ec::filtering_statistics().full_measure_comparisons, 1u);
    }
}

void unsigned_measure_limits()
{
    using F = fixture<std::uint64_t, policy::approximate_output>;
    using Root = F::root;
    F::native const origin(0, 0), small(std::uint64_t(1) << 30, std::uint64_t(1) << 30);
    F::native const maximum((std::numeric_limits<std::uint64_t>::max)(), 0);
    auto const distance = ec::root_squared_distance<Root>(origin, small);
    std::uint64_t exact;
    BOOST_TEST(bp::small_integer_distance(distance, exact));
    BOOST_TEST_EQ(exact, std::uint64_t(1) << 61);
    auto const large = ec::root_squared_distance<Root>(origin, maximum);
    BOOST_TEST(!bp::small_integer_distance(large, exact));
    typename decltype(large)::compact_pair_type compact;
    BOOST_TEST(large.compact_pair(compact));
    BOOST_TEST_EQ(ec::fixed_pair_compare(compact, large.fixed_pair()), 0);
    no_allocations([&] { BOOST_TEST_EQ(ec::root_compare_distance(large, large), 0); });
}

void shared_intersection_provenance()
{
    using F = fixture<double, policy::approximate_output>;
    using Root = F::root;
    F first({{0, 0, 3, 3, 0, 2, 3, 0}});
    ec::borrowed_stage<Root> stage(first.context, first.a, first.b);
    auto operations = first.operations;
    std::swap(operations[0], operations[1]);
    ec::borrowed_record<Root> record;
    record.stage = &stage;
    record.operations = {{&operations[0], &operations[1]}};
    F::point const second(record);
    record.numeric.init(second.producers());
    // Equal local source/segment IDs in another stage need not denote the
    // same original supports, even within one Root context.
    F::line const other_a{{10, 10}, {13, 13}}, other_b{{10, 12}, {13, 10}};
    ec::borrowed_stage<Root> other_stage(first.context, other_a, other_b);
    auto other_operations = first.operations;
    ec::borrowed_record<Root> other_record;
    other_record.stage = &other_stage;
    other_record.operations = {{&other_operations[0], &other_operations[1]}};
    F::point const other(other_record);
    other_record.numeric.init(other.producers());
    auto const check = [&]
    {
        no_allocations([&]
        {
            BOOST_TEST(bp::same_intersection(first.p, second));
            BOOST_TEST(!bp::same_intersection(first.p, other));
            BOOST_TEST(ec::root_equal<Root>(first.p, second));
            BOOST_TEST_EQ((ec::root_compare_coordinate<0, Root>(first.p, second)), 0);
            BOOST_TEST_EQ(ec::root_orientation<Root>(first.p, second, first.a.front()), 0);
            BOOST_TEST_EQ(ec::root_direction<Root>(first.a.front(), first.p, second), 0);
        });
        BOOST_TEST_EQ(ec::point_cache_detail::native_cache<double>::size(), 0u);
    };
    check();
    auto finalize = [&](auto& source, auto& ops)
    {
        for (auto& operation : ops)
        {
            auto const edge = source.support(operation.seg_id);
            operation.seg_id = bg::segment_identifier(edge.range, -1, -1, edge.segment);
        }
        source.finalized = true;
        source.sources = {};
    };
    finalize(first.stage, first.operations);
    check();
    finalize(stage, operations);
    check();
}

int main()
{
    try
    {
        BOOST_TEST_EQ(std::fegetround(), FE_TONEAREST);
        cases<double, policy::approximate_output>(0);
        cases<double, policy::eager_nearest_output>(0);
        cases<double, policy::deferred_nearest_output>(0);
        cases<double, policy::nearest_output>(1048576);
        cases<std::int64_t, policy::approximate_output>((std::int64_t(1) << 54) + 1);
        cases<std::int64_t, policy::eager_nearest_output>((std::int64_t(1) << 54) + 1);
        cases<std::int64_t, policy::deferred_nearest_output>((std::numeric_limits<std::int64_t>::max)() - 64);
        cases<std::int64_t, policy::deferred_nearest_output>((std::numeric_limits<std::int64_t>::min)() + 64);
        random_cases<policy::approximate_output>();
        random_cases<policy::eager_nearest_output>();
        random_cases<policy::deferred_nearest_output>();
        stored_filter_reads_once<policy::approximate_output>();
        stored_filter_reads_once<policy::bounded_ulp_output<4>>();
        indexed_native();
        lean_does_not_reconstruct_errors();
        coordinate_filter_reads_one_axis();
        cache_free_support_zero<policy::approximate_output>();
        cache_free_support_zero<policy::eager_nearest_output>();
        cache_free_support_zero<policy::deferred_nearest_output>();
        native_measure_fallbacks<double, policy::approximate_output>();
        native_measure_fallbacks<double, policy::eager_nearest_output>();
        native_measure_fallbacks<double, policy::deferred_nearest_output>();
        native_measure_fallbacks<std::int64_t, policy::approximate_output>();
        unsigned_measure_limits();
        shared_intersection_provenance();
        polygon_line_audit<policy::approximate_output>("ordinary");
        polygon_line_audit<policy::eager_nearest_output>("eager");
        polygon_line_audit<policy::deferred_nearest_output>("deferred");
        std::cout << triples << " rational predicate triples; " << comparisons << " measure comparisons; "
                  << allocation_checks << " allocation checks, " << total_allocations << " allocations.\n";
    }
    catch (std::exception const& error)
    {
        counting = false;
        std::cerr << error.what() << '\n';
        BOOST_ERROR("Borrowed predicate test threw");
    }
    return boost::report_errors();
}
