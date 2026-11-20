// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.

#include <boost/core/bit.hpp>
#include <boost/core/lightweight_test.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/strategies/relate/cartesian_exact.hpp>
#include <boost/multiprecision/cpp_int.hpp>

#include <algorithm>
#include <array>
#include <cfenv>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <new>
#include <random>
#include <type_traits>
#include <utility>
#include <vector>

namespace bg = boost::geometry;
namespace ec = bg::strategies::detail::exact_cartesian;
namespace sr = bg::strategies::relate;

namespace
{
struct counts { std::size_t calls = 0, bytes = 0; };
counts* counter = nullptr;
}

void* operator new(std::size_t n)
{
    if (counter) { ++counter->calls; counter->bytes += n; }
    if (void* p = std::malloc(n ? n : 1)) return p;
    throw std::bad_alloc();
}
void* operator new[](std::size_t n) { return ::operator new(n); }
void operator delete(void* p) noexcept { std::free(p); }
void operator delete[](void* p) noexcept { ::operator delete(p); }
#if defined(__cpp_sized_deallocation)
void operator delete(void* p, std::size_t) noexcept { std::free(p); }
void operator delete[](void* p, std::size_t) noexcept { std::free(p); }
#endif

namespace
{
using q = boost::multiprecision::cpp_rational;
using xy = std::array<q, 2>;
using endpoints = std::array<double, 8>;
using point = bg::model::d2::point_xy<double>;
using line = bg::model::linestring<point>;
using polygon = bg::model::polygon<point>;
template <std::size_t N> using lean = sr::lean_bounded_ulp_output<N>;

template <typename F> counts count_new(F f)
{
    counts result;
    counter = &result;
    try { f(); } catch (...) { counter = nullptr; throw; }
    counter = nullptr;
    return result;
}

int signum(q const& r) { return (r > 0) - (r < 0); }

xy crossing(endpoints const& v)
{
    q const ax(v[0]), ay(v[1]), bx(v[2]), by(v[3]);
    q const cx(v[4]), cy(v[5]), dx(v[6]), dy(v[7]);
    q const w = (bx-ax)*(dy-cy) - (by-ay)*(dx-cx);
    q const t = ((cx-ax)*(dy-cy) - (cy-ay)*(dx-cx)) / w;
    return {{ax+t*(bx-ax), ay+t*(by-ay)}};
}

// Find the two adjacent doubles by exact comparisons. No production rounding
// helper or rational-to-floating conversion, even for subnormal halfway ties.
double nearest(q const& r)
{
    if (r < 0) return -nearest(-r);
    using word = std::uint64_t;
    word lo = 0, hi = boost::core::bit_cast<word>((std::numeric_limits<double>::max)());
    while (hi - lo > 1)
    {
        word const m = lo + (hi - lo) / 2;
        if (q(boost::core::bit_cast<double>(m)) <= r) lo = m;
        else hi = m;
    }
    double const a = boost::core::bit_cast<double>(lo), b = boost::core::bit_cast<double>(hi);
    int const side = signum(2*r - q(a) - q(b));
    return side < 0 || (side == 0 && !(lo & 1)) ? a : b;
}

std::uint64_t rank64(double x)
{
    constexpr auto sign = std::uint64_t(1) << 63;
    auto const bits = boost::core::bit_cast<std::uint64_t>(x);
    return bits & sign ? sign - (bits & ~sign) : sign + bits;
}

constexpr std::size_t budget_cap = 1u << 24;

// This is a real-number bound, NOT a lattice radius around RN(X). Form
// even the subnormal floor product in Q, independently of point_state::bounds.
template <std::size_t N> void producer_bound(double z, q const& exact)
{
    BOOST_TEST(std::isfinite(z));
    if (!std::isfinite(z)) return;
    constexpr std::size_t k = N < budget_cap ? N : budget_cap;
    q const height((std::max)(std::abs(z), (std::numeric_limits<double>::min)()));
    q const limit = height * (k + 2) / (std::uint64_t(1) << 53);
    if (abs(q(z) - exact) > limit)
        std::fprintf(stderr, "producer bound: N=%zu K=%zu z=%a\n", N, k, z);
    BOOST_TEST(abs(q(z) - exact) <= limit);
}

template <std::size_t N> void check_coordinate(double x, q const& exact)
{
    auto const a = rank64(x), b = rank64(nearest(exact));
    auto const steps = a < b ? b-a : a-b;
    if (!std::isfinite(x) || steps > N)
        std::fprintf(stderr, "N=%zu emitted=%a nearest=%a\n", N, x, nearest(exact));
    BOOST_TEST(std::isfinite(x) && steps <= N);
}

void encloses(ec::interval const& i, q const& exact)
{
    BOOST_TEST(i.error >= 0);
    if (std::isfinite(i.error))
    {
        BOOST_TEST(std::isfinite(i.approximate));
        if (std::isfinite(i.approximate))
            BOOST_TEST(abs(exact - q(i.approximate)) <= q(i.error));
    }
}

template <std::size_t N> void construction(endpoints const& v)
{
    using state = ec::construction_rounding::point_state<double, lean<N>>;
    ec::point_cache_scope<double> scope;
    state p;
    p.init(v);
    xy const exact = crossing(v);
    for (std::size_t d = 0; d < 2; ++d)
    {
        double out;
        auto const alloc = count_new([&]{ out = p.template emit<double>(d, v); });
        BOOST_TEST_EQ(alloc.calls, 0u);
        producer_bound<N>(p.approximate(d), exact[d]);
        check_coordinate<N>(out, exact[d]);
        encloses(p.bounds(d), exact[d]);
        auto const refined = p.refine(v, d);
        BOOST_TEST_EQ(refined.approximate, nearest(exact[d]));
        producer_bound<N>(p.approximate(d), exact[d]);
        encloses(refined, exact[d]);
    }
}

endpoints const witness{{-32.487, 5.108, 26.755, 42.581, -24.619, 48.479, 11.438, -6.773}};

// Proper crossing at t=2/7 on the first segment. The ordinary determinant
// bounds reject x even at N=100; tightening certifies its original FMA even
// at N=8. That FMA is not RN(3/7), so exact fallback cannot mask a wiring bug.
endpoints const det_witness{{-1, 1000, 4, 1004, -1, 999, 1, 1002}};

template <std::size_t N, bool Store, bool DetRetry>
ec::forward_error::point_approximation det_case(endpoints const& v, xy const& exact,
                                                std::size_t missed, bool retry_passes)
{
    using policy = std::conditional_t<Store, sr::bounded_ulp_output<N, false, DetRetry>,
                                      sr::lean_bounded_ulp_output<N, false, DetRetry>>;
    using state = ec::construction_rounding::point_state<double, policy>;
    static_assert(state::stores_errors == Store, "Exercise the requested layout");
    ec::forward_error::point_approximation first;
    first.error.fill(-1.);
    auto const pass = ec::construction_rounding::detail::bounded_intersection
        <N, false, true, true, Store, DetRetry>(v, first);
    ec::point_cache_scope<double> scope;
    state p;
    p.init(v);
    for (std::size_t d = 0; d < 2; ++d)
    {
        BOOST_TEST_EQ(pass[d], d != missed || (DetRetry && retry_passes));
        double const rn = nearest(exact[d]);
        double const emitted = p.template emit<double>(d, v);
        if constexpr (!Store) producer_bound<N>(emitted, exact[d]);
        if (pass[d])
        {
            check_coordinate<N>(first.value[d], exact[d]);
            if (d == missed)
            {
                double const candidate = std::fma(nearest(q(2, 7)), 5., v[d]);
                BOOST_TEST_NE(candidate, rn);
                BOOST_TEST_EQ(first.value[d], candidate);
            }
            BOOST_TEST_EQ(emitted, first.value[d]);
            if constexpr (Store)
            {
                BOOST_TEST(std::isfinite(first.error[d]));
                encloses({first.value[d], first.error[d]}, exact[d]);
            }
        }
        else BOOST_TEST_EQ(emitted, rn);
        if constexpr (!Store) BOOST_TEST_EQ(first.error[d], -1.);
        BOOST_TEST(std::isfinite(p.bounds(d).error));
        encloses(p.bounds(d), exact[d]);
    }
    return first;
}

template <std::size_t N, bool Store> void det_cases()
{
    auto check = [](endpoints const& v, std::size_t missed, bool retry_passes)
    {
        auto const exact = crossing(v);
        auto const off = det_case<N, Store, false>(v, exact, missed, retry_passes);
        auto const on = det_case<N, Store, true>(v, exact, missed, retry_passes);
        std::size_t const kept = 1 - missed;
        BOOST_TEST_EQ(off.value[kept], on.value[kept]);
        if constexpr (Store) BOOST_TEST_EQ(off.error[kept], on.error[kept]);
    };
    check(det_witness, 0, true);
    endpoints transposed = det_witness, shifted = det_witness;
    for (std::size_t i = 0; i < 8; i += 2)
    {
        std::swap(transposed[i], transposed[i+1]);
        shifted[i] -= 7./16;
    }
    check(transposed, 1, true);
    // Now x=-1/112 and |t*ux|/|x|=160: even the optimistic retry fails.
    check(shifted, 0, false);
}

// Limit public algorithm instantiations to one budget and layout.
template <bool DetRetry> void public_det_case()
{
    sr::cartesian_exact<double, sr::lean_bounded_ulp_output<8, false, DetRetry>> strategy;
    auto const& v = det_witness;
    line const a{{v[0], v[1]}, {v[2], v[3]}}, b{{v[4], v[5]}, {v[6], v[7]}};
    auto const exact = crossing(v);
    std::vector<point> points;
    bg::intersection(a, b, points, strategy);
    BOOST_TEST_EQ(points.size(), 1u);
    if (points.size() == 1)
    {
        check_coordinate<8>(bg::get<0>(points[0]), exact[0]);
        check_coordinate<8>(bg::get<1>(points[0]), exact[1]);
    }
    BOOST_TEST_EQ(bg::relation(a, b, strategy).str(), "0F1FF0102");
    if constexpr (DetRetry)
    {
        std::vector<bg::model::d2::point_xy<q>> rational;
        bg::intersection(a, b, rational, strategy);
        BOOST_TEST_EQ(rational.size(), 1u);
        if (rational.size() == 1)
        {
            BOOST_TEST(bg::get<0>(rational[0]) == exact[0]);
            BOOST_TEST(bg::get<1>(rational[0]) == exact[1]);
        }
    }
}

template <std::size_t N> void construction_cases()
{
    construction<N>(witness);
    construction<N>({{-1, -1, 1, 1, -1, 1, 1, -1}});
    double const tiny = (std::numeric_limits<double>::denorm_min)();
    for (double s : {-1., 1.})
    {
        construction<N>({{0, 0, s*tiny, 1, -1, 0.5, 1, 0.5}});
        for (int e : {-1022, -700, -64, 0, 64, 700, 1022})
        {
            double const x = s*std::ldexp(1., e);
            double const a = std::nextafter(x, 0.), b = std::nextafter(x, 2*x);
            construction<N>({{a, -2, b, 1, a, 0, b, 0}});
        }
    }
    // Large nearly parallel products, exact denominator -1.
    double const k = 0x1p27;
    construction<N>({{0, 0, k+1, k, 0, 1, k, k}});
    std::mt19937 rng(0x20261120);
    for (unsigned i = 0; i < 48; ++i)
    {
        // Opposite diagonals across the same strip: always a proper crossing.
        auto offset = [&] { return double(1 + rng()%1000)/17; };
        double const base = i%2 ? 0x1p20 : -0x1p-12;
        construction<N>({{base, base-offset(), base+2, base+offset(),
                          base, base+offset(), base+2, base-offset()}});
    }
}

// Pin the path by its acceptance mask; a silently substituted nearest
// fallback could otherwise pass every producer-bound assertion.
template <std::size_t N, bool Retry, bool DetRetry, bool Rescue = true, bool Cross = true>
std::array<double, 2> lean_path(endpoints const& v, std::array<bool, 2> expected)
{
    auto const exact = crossing(v);
    ec::forward_error::point_approximation candidate;
    auto const accepted = ec::construction_rounding::detail::bounded_intersection
        <N, Retry, Rescue, Cross, false, DetRetry>(v, candidate);
    BOOST_TEST(accepted == expected);
    for (std::size_t d = 0; d != 2; ++d)
        if (accepted[d]) producer_bound<N>(candidate.value[d], exact[d]);
    if constexpr (Rescue && Cross)
    {
        using policy = sr::lean_bounded_ulp_output<N, Retry, DetRetry>;
        ec::point_cache_scope<double> scope;
        ec::construction_rounding::point_state<double, policy> p;
        p.init(v);
        for (std::size_t d = 0; d != 2; ++d)
        {
            BOOST_TEST_EQ(p.approximate(d), accepted[d] ? candidate.value[d] : nearest(exact[d]));
            producer_bound<N>(p.approximate(d), exact[d]);
            p.refine(v, d);
            BOOST_TEST_EQ(p.approximate(d), nearest(exact[d]));
            producer_bound<N>(p.approximate(d), exact[d]);
        }
    }
    return candidate.value;
}

void retry_paths()
{
    // Nonzero residual of the rebuilt 1/3 quotient; the ordinary rounded
    // denominator is zero. Disabling rescue must reject both coordinates.
    double const m = 0x1p30;
    endpoints const rescue{{m, m, 2*m+1, 2*m, m+1, m+1, 2*m-1, 2*m-2}};
    lean_path<1, false, false, false>(rescue, {{false, false}});
    lean_path<1, false, false>(rescue, {{true, true}});

    // Translate the same D=-3, T=-1 geometry so the rebuilt FMA is
    // millions of steps from RN(x). The y coordinate still passes at N=1.
    double const n = 0x1p29, a = 32 - (n+1)/3;
    endpoints const translated{{a, n, a+n+1, 2*n, a+1, n+1, a+n-1, 2*n-2}};
    lean_path<budget_cap, false, false, false>(translated, {{false, false}});
    lean_path<1, false, false>(translated, {{false, true}});
    auto const rebuilt = lean_path<budget_cap, false, false>(translated, {{true, true}});
    BOOST_TEST(crossing(translated)[0] == 32);
    BOOST_TEST(rank64(32.) - rank64(rebuilt[0]) > 1000000);

    endpoints const division{{1.6875, 9, 2.6875, 10, 1.6875, 10, 3.6875, 8}};
    lean_path<4, false, false>(division, {{false, true}});
    lean_path<4, true, false>(division, {{true, true}});

    // Nonzero subtraction tail on the shared vertical direction.
    endpoints const tails{{4, 8, 5, 8, 4.5, -0x1p-10, 4.5, 0x1p54}};
    lean_path<1, false, true, true, false>(tails, {{false, true}});
    lean_path<1, false, true>(tails, {{true, true}});
}

void wide_integer_fallback()
{
    using integer = std::int64_t;
    integer const a = (integer(1) << 58) + 1;
    // The double endpoints would collapse both x directions. For the
    // original integers the proper crossing has x=a+2/3 and y=2/3.
    std::array<integer, 8> const v{{a, 0, a+2, 2, a, 1, a+2, 0}};
    xy const exact{{q(a) + q(2, 3), q(2, 3)}};
    ec::point_cache_scope<integer> scope;
    ec::construction_rounding::point_state<integer, lean<1>> p;
    p.init(v);
    for (std::size_t d = 0; d != 2; ++d)
    {
        BOOST_TEST_EQ(p.approximate(d), nearest(exact[d]));
        producer_bound<1>(p.approximate(d), exact[d]);
    }
}

template <std::size_t N> void wide_budget(std::size_t scale)
{
    // Exact t=1/3, x=delta. The original FMA misses by exactly 2^-54:
    // at scale=2^24, delta=2/scale the accepted head misses by 2^22
    // relative roundoff units. A nearest-only corpus cannot exercise this.
    for (double sign : {-1., 1.})
    {
        auto input = [&](double delta) -> endpoints
        { return {{sign*(-1+delta), 2, sign*(2+delta), 5, -sign, 3, 2*sign, 3}}; };
        auto v = input(32./scale);
        auto z = lean_path<N, false, false>(v, {{true, true}});
        BOOST_TEST(abs(q(z[0]) - crossing(v)[0]) == q(0x1p-54));
        v = input(2./scale);
        lean_path<N, false, false>(v, {{false, true}});
        z = lean_path<N, false, true>(v, {{true, true}});
        BOOST_TEST(abs(q(z[0]) - crossing(v)[0]) == q(0x1p-54));
        // The direction contribution alone prevents passing at the cap,
        // even when the nominal policy would have ample room.
        lean_path<N, false, true>(input(1./scale), {{false, true}});
    }
}

template <std::size_t N> struct fixture
{
    using root = ec::borrowed_root<double, lean<N>, line, line>;
    line a, b;
    ec::borrowed_context<root> context;
    ec::borrowed_stage<root> stage;
    std::array<typename root::operation_type, 2> ops;
    ec::borrowed_record<root> record;
    ec::borrowed_point<root> p;

    explicit fixture(endpoints const& v)
        : a{{v[0], v[1]}, {v[2], v[3]}}, b{{v[4], v[5]}, {v[6], v[7]}},
          stage(context, a, b), p(record)
    {
        ops[0].seg_id = bg::segment_identifier(0, -1, -1, 0);
        ops[1].seg_id = bg::segment_identifier(1, -1, -1, 0);
        record.stage = &stage;
        record.operations = {{&ops[0], &ops[1]}};
        reset();
    }
    void reset() { record.numeric.init(p.producers()); }
};

void operand_heads()
{
    namespace bp = ec::borrowed_predicates_detail;
    using root = fixture<100>::root;
    using query = ec::query_point<root>;
    using wide = bg::model::d2::point_xy<std::int64_t>;
    using foreign = fixture<1>::root;
    static_assert(!bp::step_operand<root, wide> && !bp::step_operand<root, ec::exact_point<double>>
               && !bp::step_operand<root, ec::borrowed_point<foreign>>
               && !bp::step_operand<root, ec::query_point<foreign>>, "No borrowed certificate");
    fixture<100> f(witness);
    auto check = [](auto const& operand, xy const& exact)
    {
        bp::visit<root>(operand, [&](auto const& selected)
        {
            using P = std::decay_t<decltype(selected)>;
            static_assert(bp::step_filter<root, P>, "Exercise admitted operands");
            auto const head = bp::bounded_coordinates(selected);
            for (std::size_t d = 0; d != 2; ++d)
            {
                producer_bound<100>(head[d], exact[d]);
                if constexpr (bp::exact_double<P>::value) BOOST_TEST(q(head[d]) == exact[d]);
            }
        });
    };
    xy const exact = crossing(witness), a{{q(witness[0]), q(witness[1])}};
    ec::iterator_point_view<decltype(&f.p)> const record_view(&f.p);
    check(record_view, exact);
    auto native = ec::borrowed_context<root>::native_vertex(f.a, 0);
    check(ec::iterator_point_view<decltype(&native)>(&native), a);
    check(ec::iterator_point_view<line::const_iterator>(f.a.cbegin()), a);

    auto const midpoint = query::midpoint(f.p, native);
    auto const splice = query::splice_coordinates(f.p, 1, native, 0);
    check(ec::iterator_point_view<decltype(&midpoint)>(&midpoint),
          {{(exact[0]+a[0])/2, (exact[1]+a[1])/2}});
    check(splice, {{exact[1], a[0]}});
    check(query::midpoint(native, native), a);
    check(query::splice_coordinates(native, 1, native, 0), {{a[1], a[0]}});
    double const tiny = (std::numeric_limits<double>::denorm_min)();
    check(query::midpoint(point(0, -tiny), point(tiny, -2*tiny)), {{q(tiny)/2, -3*q(tiny)/2}});
    f.p.refine(0);
    check(record_view, exact);
}

template <std::size_t N> void ambiguous(endpoints const& v)
{
    using root = typename fixture<N>::root;
    fixture<N> f(v);
    xy const z = crossing(v);
    double const x = nearest(z[0]), y = nearest(z[1]);
    // Expand from the current head (not the oracle). Neither endpoint of the
    // new box may prune this crossing, even before its first predicate.
    ec::search_point const head(f.p.approximate(0), f.p.approximate(1));
    ec::search_box bounds(head, head);
    ec::root_expand<root>::apply(bounds, f.p);
    BOOST_TEST((q(bg::get<bg::min_corner, 0>(bounds)) <= z[0]));
    BOOST_TEST((q(bg::get<bg::min_corner, 1>(bounds)) <= z[1]));
    BOOST_TEST((q(bg::get<bg::max_corner, 0>(bounds)) >= z[0]));
    BOOST_TEST((q(bg::get<bg::max_corner, 1>(bounds)) >= z[1]));
    // Probe both sides of the true crossing and the rounding cell containing
    // it. Reset between predicates so a prior exact fallback cannot mask an
    // unsound initial filter. The wide-budget witness's FMA lies 5 steps away.
    for (double t : {std::nextafter(x, x-1), x, std::nextafter(x, x+1)})
    {
        point const center(t, y), above(t, y+1), left(t-1, y), right(t+1, y);
        int const cmp = signum(z[0] - q(t));
        f.reset();
        BOOST_TEST_EQ((ec::root_compare_coordinate<0, root>(f.p, center)), cmp);
        f.reset();
        BOOST_TEST_EQ(ec::root_equal<root>(f.p, center), z[0] == q(t) && z[1] == q(y));
        f.reset();
        BOOST_TEST_EQ(ec::root_orientation<root>(center, above, f.p), -cmp);
        f.reset();
        BOOST_TEST_EQ(ec::root_direction<root>(left, center, f.p), cmp);
        f.reset();
        auto const dl = ec::root_squared_distance<root>(f.p, left);
        auto const dr = ec::root_squared_distance<root>(f.p, right);
        auto square = [](q a) -> q { return a*a; };
        q const dy2 = square(z[1]-q(y));
        q const l = square(z[0]-q(bg::get<0>(left)))+dy2;
        q const r = square(z[0]-q(bg::get<0>(right)))+dy2;
        encloses(dl.bounds(), l);
        encloses(dr.bounds(), r);
        BOOST_TEST_EQ(ec::root_compare_distance(dl, dr), signum(l-r));
        f.reset();
        auto const det = ec::root_signed_determinant<root>(center, above, f.p);
        q const d = (q(y+1)-q(y))*(q(t)-z[0]);
        encloses(det.bounds(), d);
        BOOST_TEST_EQ(ec::root_compare_distance(det, -det), signum(d));
    }
    f.reset();
    auto const m = ec::query_point<root>::midpoint(f.p, f.a.front());
    auto const corner = ec::query_point<root>::splice_coordinates(f.p, 0, f.b.front(), 1);
    auto const mh = m.homogeneous(), ch = corner.homogeneous();
    BOOST_TEST(q(mh.X)/q(mh.W) == (z[0]+q(v[0]))/2);
    BOOST_TEST(q(mh.Y)/q(mh.W) == (z[1]+q(v[1]))/2);
    BOOST_TEST(q(ch.X)/q(ch.W) == z[0]);
    BOOST_TEST(q(ch.Y)/q(ch.W) == q(v[5]));
    BOOST_TEST_EQ(ec::root_orientation<root>(f.a.front(), m, f.p), 0);
    BOOST_TEST_EQ((ec::root_compare_coordinate<0, root>(f.p, corner)), 0);
}

template <std::size_t N> void owning_point_admission()
{
    using root = typename fixture<N>::root;
    using owning = ec::exact_point<double>;
    point const native(0.5, 1.5);
    double const h = 0x1p160;
    // The second support fixes x=1. Cancellation on the first support can
    // produce x=0. An owning point does not inherit the root's N-step proof.
    for (bool reverse : {false, true})
    {
        ec::point_cache_scope<double> scope;
        auto const legacy = ec::construct_intersection(owning(h, 0), owning(-h, 3),
                                                       owning(1, 1), owning(1, 2));
        int const cmp = reverse ? ec::root_compare_coordinate<0, root>(native, legacy)
                                : ec::root_compare_coordinate<0, root>(legacy, native);
        BOOST_TEST_EQ(cmp, reverse ? -1 : 1);
    }
}

template <typename Policy> counts public_allocations()
{
    sr::cartesian_exact<double, Policy> strategy;
    point const p(11, 11), q0(12, 12);
    bg::model::box<point> const box(point(10, 10), point(20, 20));
    auto const native = count_new([&] {
        BOOST_TEST(!bg::equals(p, q0, strategy));
        BOOST_TEST(bg::within(p, box, strategy));
    });
    BOOST_TEST_EQ(native.calls, 0u);
    line const a{{100, 100}, {103, 103}}, b{{100, 102}, {103, 100}};
    std::vector<point> points;
    points.reserve(1);
    auto const alloc = count_new([&]{ bg::intersection(a, b, points, strategy); });
    BOOST_TEST_EQ(points.size(), 1u);
    if (points.size() == 1)
    {
        check_coordinate<Policy::steps>(bg::get<0>(points[0]), q(506, 5));
        check_coordinate<Policy::steps>(bg::get<1>(points[0]), q(506, 5));
    }
    return alloc;
}

template <std::size_t N> void topology()
{
    sr::cartesian_exact<double, lean<N>> strategy;
    double const b = std::nextafter(1., 2.);
    polygon triangle, clip;
    triangle.outer() = {{1,-1}, {1,1}, {b,-1}, {1,-1}};
    clip.outer() = {{0,0}, {0,2}, {2,2}, {2,0}, {0,0}};
    strategy.with_sources(triangle, clip, [](auto const& bound)
    {
        using selected = typename std::decay_t<decltype(bound)>::template intersection_point_type<point>;
        using root = ec::borrowed_root<double, lean<N>, polygon, polygon>;
        static_assert(std::is_same<selected, ec::borrowed_point<root>>::value,
                      "Native public polygons must bind the tested lean record type");
    });
    // x=b is attained only at the triangle's lower right vertex. Moving the
    // rectangle's left edge one step to the right opens a strictly positive gap.
    polygon right;
    for (double edge : {b, std::nextafter(b, 2.)})
    {
        right.outer() = {{edge,-2}, {edge,2}, {2,2}, {2,-2}, {edge,-2}};
        BOOST_TEST_EQ(bg::relation(triangle, right, strategy).str(),
                      edge == b ? "FF2F01212" : "FF2FF1212");
    }
    // Cutting at y=0 leaves the top half of the triangle. Its two base
    // vertices coincide under nearest-double emission but are exactly distinct.
    q const h = (1+q(b))/2;
    using qp = bg::model::d2::point_xy<q>;
    std::vector<bg::model::polygon<qp>> result;
    BOOST_TEST_EQ(bg::relation(triangle, clip, strategy).str(), "212101212");
    bg::intersection(triangle, clip, result, strategy);
    BOOST_TEST_EQ(result.size(), 1u);
    if (result.size() != 1) return;
    BOOST_TEST(result[0].inners().empty());
    auto const& ring = result[0].outer();
    BOOST_TEST_EQ(ring.size(), 4u);
    if (ring.size() != 4) return;
    std::array<xy, 3> const clockwise{{{{1,0}}, {{1,1}}, {{h,0}}}};
    for (unsigned start = 0; start < 3; ++start)
    {
        bool matches = true;
        for (unsigned i = 0; i < 4; ++i)
        {
            auto const& e = clockwise[(start+i)%3];
            matches &= bg::get<0>(ring[i]) == e[0] && bg::get<1>(ring[i]) == e[1];
        }
        if (matches) return;
    }
    BOOST_ERROR("Exact closed boundary differs from the analytically clipped triangle");
}

template <std::size_t N> void run()
{
    using stored = sr::bounded_ulp_output<N>;
    using root = ec::borrowed_root<double, lean<N>, polygon, polygon>;
    using other = ec::borrowed_root<double, stored, polygon, polygon>;
    static_assert(sizeof(typename root::numeric_state) + 2*sizeof(double)
                  == sizeof(typename other::numeric_state), "No per-axis error fields");
    static_assert(sizeof(ec::borrowed_record<root>) < sizeof(ec::borrowed_record<other>),
                  "The smaller state must reach the turn record");
    construction_cases<N>();
    owning_point_admission<N>();
    ambiguous<N>(witness);
    endpoints transposed = witness;
    for (std::size_t i = 0; i < 8; i += 2) std::swap(transposed[i], transposed[i+1]);
    ambiguous<N>(transposed);
    topology<N>();
    auto const a = public_allocations<lean<N>>(), b = public_allocations<stored>();
    std::printf("N=%zu public crossing: lean %zu calls/%zu bytes, stored %zu calls/%zu bytes\n",
                N, a.calls, a.bytes, b.calls, b.bytes);
    // Total allocator bytes include the algorithm's other work buffers.
    // The layout assertions above apply specifically to constructed records.
    BOOST_TEST(a.calls <= b.calls);
}
} // namespace

int main()
{
    BOOST_TEST_EQ(std::fegetround(), FE_TONEAREST);
    run<1>(); run<2>(); run<4>(); run<8>(); run<100>();
    retry_paths();
    wide_integer_fallback();
    wide_budget<1024>(1024);
    wide_budget<budget_cap-1>(budget_cap);
    wide_budget<budget_cap>(budget_cap);
    wide_budget<4*budget_cap>(budget_cap);
    operand_heads();
    det_cases<8, true>(); det_cases<8, false>();
    det_cases<100, true>(); det_cases<100, false>();
    public_det_case<false>(); public_det_case<true>();
    return boost::report_errors();
}
