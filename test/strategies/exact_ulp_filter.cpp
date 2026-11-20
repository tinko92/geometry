// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.

// From this worktree:
// clang++ -std=c++17 -O3 -DNDEBUG -ffp-contract=off -Iinclude -I../..
//   test/strategies/exact_ulp_filter.cpp -o /tmp/bg-ulp-filter
// /tmp/bg-ulp-filter
#include <boost/geometry/strategies/cartesian/exact/ulp_filter.hpp>
#include <boost/geometry/strategies/cartesian/exact/borrowed_predicates.hpp>
#include <boost/geometry/strategies/cartesian/exact/borrowed_query.hpp>
#include <boost/geometry/geometries/linestring.hpp>
#include <boost/geometry/strategy/cartesian/side_robust.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/core/lightweight_test.hpp>

#include <array>
#include <cfenv>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <cstdlib>
#include <new>
#include <random>
#include <utility>
#include <type_traits>

namespace ec = boost::geometry::strategies::detail::exact_cartesian;
namespace uf = ec::ulp_filter;
using rational = ec::exact_number;
using coordinates = uf::coordinates;
using exact_coordinates = std::array<rational, 2>;

namespace
{
bool allocation_guard = false;
std::size_t allocations = 0;
}

void* operator new(std::size_t size)
{
    if (allocation_guard) ++allocations;
    if (void* p = std::malloc(size == 0 ? 1 : size)) return p;
    throw std::bad_alloc();
}
void operator delete(void* p) noexcept { std::free(p); }
void* operator new[](std::size_t n) { return ::operator new(n); }
void operator delete[](void* p) noexcept { std::free(p); }
#if defined(__cpp_sized_deallocation)
void operator delete(void* p, std::size_t) noexcept { std::free(p); }
void operator delete[](void* p, std::size_t) noexcept { std::free(p); }
#endif

double previous(double a) { return std::nextafter(a, -std::numeric_limits<double>::infinity()); }
double next(double a) { return std::nextafter(a, std::numeric_limits<double>::infinity()); }
int sign(rational const& q) { return int(q > 0) - int(q < 0); }

int orientation_reference(exact_coordinates const& a, exact_coordinates const& b, exact_coordinates const& c)
{
    // Independent arbitrary-precision rational arithmetic, never the filter's
    // floating magnitude expression or its error coefficient.
    rational const determinant = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
    return sign(determinant);
}

exact_coordinates native(coordinates const& p)
{ return {{ec::exact_value(p[0]), ec::exact_value(p[1])}}; }

void proof_constants()
{
    static_assert(std::is_same<uf::step_certificate<0>, uf::nearest_certificate>::value);
    static_assert(std::is_same<uf::producer_certificate<0>, uf::nearest_certificate>::value);
    rational const u = ec::exact_value(uf::detail::unit_roundoff);
    rational const one_minus_u = 1 - u;
    rational const square = one_minus_u * one_minus_u;
    rational const eta = ec::exact_value((std::numeric_limits<double>::denorm_min)());
    BOOST_TEST_EQ(2 * u * ec::exact_value(uf::detail::minimum_normal), eta);
    rational const k = 1u << 24, c = 2*k + 3;
    BOOST_TEST((2*k+2)*2*k*u <= 1);
    BOOST_TEST(2*(1+u)*c + c*c*u + (2+u)*(1+u)*(1+u) < 2*c + rational(5)/2);
    BOOST_TEST((1+8*u)*square*square >= 1+2*u);
    BOOST_TEST(rational(13)/2*(1+8*u)+4 < 12*(1-u)-1);
    BOOST_TEST((4*k+9)*(1+8*u)+4 < (4*k+16)*(1-u)-1);
    BOOST_TEST((1-18*u)*(1+u)*(1+u) < square*square);
    BOOST_TEST((2*k+9)*(1+8*u)+4 < (2*k+16)*(1-u)-1);
}

// Independently generate the extreme points allowed by the public certificate:
// move N binary64 steps, then almost half the next rounding cell. In particular
// these errors need not satisfy the tighter bound from our FMA construction.
template <unsigned N>
rational step_extreme(double head, bool up)
{
    constexpr unsigned k = N < (1u << 24) ? N : (1u << 24);
    if constexpr (k <= 100)
        for (unsigned i = 0; i < k; ++i) head = up ? next(head) : previous(head);
    else
    {
        // Jump in magnitude space, reflecting at zero (test inputs stay finite).
        double magnitude = std::abs(head);
        std::uint64_t bits;
        std::memcpy(&bits, &magnitude, sizeof(bits));
        bool const away = up == !std::signbit(head);
        bool const negative = away || bits >= k ? std::signbit(head) : !std::signbit(head);
        bits = away ? bits + k : bits >= k ? bits - k : k - bits;
        std::memcpy(&magnitude, &bits, sizeof(bits));
        head = negative ? -magnitude : magnitude;
    }
    rational const r = ec::exact_value(head);
    return r + (ec::exact_value(up ? next(head) : previous(head)) - r) * 511 / 1024;
}

// Monotonicity orders even adjacent nearest heads. Equal heads can hide
// either exact order, also across signed zero, and must not certify equality.
template <bool AI, bool BI> void nearest_order_case(double a, double b)
{
    for (bool up : {false, true})
    {
        rational const x = AI ? step_extreme<0>(a, up) : rational(a);
        rational const y = BI ? step_extreme<0>(b, !up) : rational(b);
        int result = 17;
        bool const pass = uf::try_compare_n<0, AI, BI>(a, b, result);
        BOOST_TEST_EQ(pass, a != b || !(AI || BI));
        BOOST_TEST_EQ(result, pass ? sign(x-y) : 17);
    }
}

void nearest_order()
{
    double const tiny = (std::numeric_limits<double>::denorm_min)();
    for (double a : {-1., -tiny, -0., 0., tiny, uf::detail::minimum_normal, previous(1.), 1.})
    for (double b : {previous(a), a, next(a), -a})
    {
        nearest_order_case<false,false>(a,b);
        nearest_order_case<false,true>(a,b);
        nearest_order_case<true,false>(a,b);
        nearest_order_case<true,true>(a,b);
    }
}

template <unsigned N, bool Producer>
rational ordinate_extreme(double head, bool up)
{
    rational const loose = step_extreme<N>(head, up);
    if constexpr (!Producer) return loose;
    // Intersect the step interval with the producer's absolute enclosure so
    // the same exact coordinates also exercise the scalar-order check.
    constexpr unsigned k = N < (1u << 24) ? N : (1u << 24);
    rational const radius = (k+2) * ec::exact_value(uf::detail::unit_roundoff)
        * ec::exact_value((std::max)(std::abs(head), uf::detail::minimum_normal));
    rational const edge = ec::exact_value(head) + (up ? radius : -radius);
    return up ? (std::min)(loose, edge) : (std::max)(loose, edge);
}

template <unsigned N, unsigned Mask, bool Producer = false>
void bounded_triple(coordinates a, coordinates b, coordinates p)
{
    constexpr bool ai = (Mask & 1) != 0, bi = (Mask & 2) != 0, pi = (Mask & 4) != 0;
    auto make = [](coordinates z, bool uncertain, bool up) -> exact_coordinates {
        if (!uncertain) return native(z);
        return {{ordinate_extreme<N,Producer>(z[0], up), ordinate_extreme<N,Producer>(z[1], !up)}};
    };
    using cert = std::conditional_t<Producer, uf::producer_certificate<N>, uf::step_certificate<N>>;
    constexpr unsigned k = N < (1u << 24) ? N : 1u << 24;
    static_assert(cert::rho == (N == 0 ? 1 : (Producer ? k : 2*k) + 2));
    static_assert(cert::bilinear_units == (N == 0 ? 12 : (Producer ? 2*k : 4*k) + 16));
    auto const x = make(a, ai, true), y = make(b, bi, false), z = make(p, pi, true);
    rational const dx = y[0]-x[0], dy = y[1]-x[1];
    rational const det = dx*(z[1]-x[1]) - dy*(z[0]-x[0]);
    rational const dot = dx*(z[0]-y[0]) + dy*(z[1]-y[1]);
    rational const dist = dx*dx + dy*dy;
    ec::forward_error::scalar_approximation d{}, s{};
    int side = 2, direction = 2, cmp = 2;
    allocation_guard = true;
    bool const bounded_d = uf::try_signed_determinant<cert,ai,bi,pi>(a,b,p,d);
    bool const bounded_s = uf::try_squared_distance<cert,ai,bi>(a,b,s);
    bool const cs = uf::try_orientation<cert,ai,bi,pi>(a,b,p,side);
    bool const cd = uf::try_dot<cert,ai,bi,pi>(a,b,p,direction);
    auto const dot_bound = uf::detail::n_bilinear<std::conditional_t<ai || bi || pi, cert, uf::nearest_certificate>,true>(
        uf::detail::n_difference<ai,bi>(a,b), uf::detail::n_difference<bi,pi>(b,p));
    bool const cc = uf::try_compare_n<N,ai,bi>(a[0],b[0],cmp);
    allocation_guard = false;
    BOOST_TEST(bounded_d && bounded_s);
    BOOST_TEST(abs(det - ec::exact_value(d.value)) <= ec::exact_value(d.error));
    BOOST_TEST(abs(dist - ec::exact_value(s.value)) <= ec::exact_value(s.error));
    BOOST_TEST(abs(dot - ec::exact_value(dot_bound.value)) <= ec::exact_value(dot_bound.error));
    BOOST_TEST_EQ(side, cs ? sign(det) : 2);
    BOOST_TEST_EQ(direction, cd ? sign(dot) : 2);
    BOOST_TEST_EQ(cmp, cc ? sign(x[0] - y[0]) : 2);
}

void randomized()
{
    std::mt19937_64 random(0x554c5046494c5445);
    auto value = [&]()
    {
        double const significand = 1 + std::ldexp(double(random() >> 11), -53);
        return std::ldexp(random() & 1 ? significand : -significand, int(random() % 241) - 120);
    };
    for (unsigned i = 0; i < 12000; ++i)
    {
        coordinates a{{value(), value()}}, b{{value(), value()}}, c{{value(), value()}};
        if (i % 7 == 0) { a[0] = 0; b[1] = 0; }
        if (i % 11 == 0)
        {
            // Nearly collinear, translated triangles whose tiny differences
            // must not hide the uncertainty of their large coordinates.
            a = {{value(), value()}};
            b = {{next(a[0]), next(a[1])}};
            c = {{next(b[0]), next(b[1])}};
        }
        bounded_triple<0,4>(a,b,c);
        bounded_triple<0,6>(a,b,c);
        bounded_triple<0,7>(a,b,c);
        bounded_triple<0,7>(c,a,b);

        // NNN must continue to use BG's existing adaptive Shewchuk predicate.
        using point = boost::geometry::model::d2::point_xy<double>;
        using robust = boost::geometry::strategy::side::side_robust<double,
            boost::geometry::strategy::side::fp_equals_policy, 3>;
        int const actual = robust::apply(point(a[0], a[1]), point(b[0], b[1]), point(c[0], c[1]));
        BOOST_TEST_EQ(actual, orientation_reference(native(a), native(b), native(c)));
    }
}

template <unsigned N, unsigned... M>
void bounded_filters(std::integer_sequence<unsigned,M...>)
{
    double const tiny = (std::numeric_limits<double>::denorm_min)();
    for (double t : {0., -0., tiny, -tiny, 0x1p-1022, -0x1p-1022, 0x1p-501,
                     .5, -1., previous(2.), next(-2.), 0x1p499})
    {
        coordinates const a{{t,-t}}, b{{next(t),previous(-t)}};
        (bounded_triple<N,M>(a,b,{{previous(t),next(-t)}}), ...);
        (bounded_triple<N,M>(a,{{1.,t}},{{-t,.25}}), ...);
        if constexpr (N == 1 || N == 100 || N == (1u << 24))
            (bounded_triple<N,M,true>(a,{{-t,t}},{{next(t),t}}), ...);
    }
    // Worst admitted opposite-sign sums, inexact tiny products and a
    // cancelling dot whose two large terms have opposite signs.
    (bounded_triple<N,M>({{0x1p500,-0x1p500}},{{-0x1p500,0x1p500}},{{0,0}}), ...);
    (bounded_triple<N,M>({{0,0}},{{0x1p-537,0x1p-538}},{{0x1p-538,0x1p-537}}), ...);
    (bounded_triple<N,M>({{-1,-1}},{{0,0}},{{1,next(-1)}}), ...);
    if constexpr (N <= 100)
    {
        for (double a : {-1., -tiny, -0., 0., tiny, previous(1.), 1.})
        {
            double b = a;
            for (unsigned k = 0; k <= 2*N+1; ++k, b = next(b))
            {
                int r = 2;
                bool const ok = uf::try_compare_n<N,true,true>(a,b,r);
                BOOST_TEST_EQ(ok, k > 2*N);
                BOOST_TEST_EQ(r, ok ? -1 : 2);
            }
        }
    }
}

void producer_witnesses()
{
    coordinates const a{{0,0}}, b{{1,1}}, p{{1,1+0x1p-44}};
    int r = 2;
    BOOST_TEST((!uf::try_orientation<uf::step_certificate<100>,false,false,true>(a,b,p,r)));
    BOOST_TEST_EQ(r, 2);
    BOOST_TEST((uf::try_orientation<uf::producer_certificate<100>,false,false,true>(a,b,p,r)));
    BOOST_TEST_EQ(r, 1); // exact det = 512u

    // Same heads, two possible contracts. The producer call is valid for
    // exact (1,0); do NOT apply it to the displaced generic coordinate.
    coordinates const unit{{1,0}};
    ec::forward_error::scalar_approximation raw{}, tight{};
    BOOST_TEST((uf::try_squared_distance<uf::step_certificate<100>,false,true>(a,unit,raw)));
    BOOST_TEST((uf::try_squared_distance<uf::producer_certificate<100>,false,true>(a,unit,tight)));
    rational const u = ec::exact_value(uf::detail::unit_roundoff);
    rational const loose = 1+200*u, err = loose*loose-1;
    BOOST_TEST_EQ(raw.value, 1.);
    BOOST_TEST_EQ(tight.value, 1.);
    BOOST_TEST_EQ(ec::exact_value(raw.error), 416*u);
    BOOST_TEST_EQ(ec::exact_value(tight.error), 216*u);
    BOOST_TEST(err <= ec::exact_value(raw.error));
    BOOST_TEST(err > ec::exact_value(tight.error));
}

void step_edges()
{
    constexpr unsigned cap = 1u << 24;
    double const limit = 1.0 + std::ldexp(double(2 * cap), -52);
    int r = 2;
    BOOST_TEST((!uf::try_compare_n<cap+1,true,true>(1, limit, r)));
    BOOST_TEST_EQ(r, 2);
    BOOST_TEST((uf::try_compare_n<cap+1,true,true>(1, next(limit), r)));
    BOOST_TEST_EQ(r, -1);
    for (double bad : {next(0x1p500), -next(0x1p500), std::numeric_limits<double>::infinity(),
                       std::numeric_limits<double>::quiet_NaN()})
    {
        coordinates const a{{0,0}}, b{{1,1}}, p{{0,bad}};
        ec::forward_error::scalar_approximation s{17,19};
        r = 2;
        BOOST_TEST((!uf::try_signed_determinant<uf::nearest_certificate,false,false,true>(a,b,p,s)));
        BOOST_TEST((!uf::try_squared_distance<uf::step_certificate<cap>,true,true>(p,b,s)));
        BOOST_TEST((!uf::try_dot<uf::nearest_certificate,false,false,true>(a,b,p,r)));
        BOOST_TEST_EQ(s.value, 17);
        BOOST_TEST_EQ(s.error, 19);
        BOOST_TEST_EQ(r, 2);
        if (!std::isfinite(bad)) BOOST_TEST((!uf::try_compare_n<0,true,true>(bad,0,r)));
    }
    coordinates const a{{1,0}}, b{{next(1),0}};
    auto mixed = uf::detail::n_difference<false,true>(a,b);
    BOOST_TEST_EQ(mixed.scale[0], 2.0);
    BOOST_TEST_EQ(mixed.scale[1], 2 * uf::detail::minimum_normal);
    BOOST_TEST_EQ((uf::detail::n_difference<false,false>(a,b).scale[0]), next(1)-1);
}

void cold_nearest()
{
    namespace bg = boost::geometry;
    namespace bp = ec::borrowed_predicates_detail;
    namespace policy = bg::strategies::relate;
    using point = bg::model::d2::point_xy<double>;
    using line = bg::model::linestring<point>;
    using root = ec::borrowed_root<double, policy::nearest_output, line, line>;
    using stored = ec::borrowed_root<double, policy::bounded_ulp_output<1>, line, line>;
    static_assert(bp::step_filter<root, point, ec::borrowed_point<root>>);
    static_assert(!bp::step_filter<stored, ec::borrowed_point<stored>>);
    static_assert(!bp::step_filter<root, ec::borrowed_point<stored>>);
    static_assert(!bp::step_filter<root, ec::exact_point<double>>);
    static_assert(!bp::step_filter<root, bg::model::d2::point_xy<std::int64_t>>);

    ec::point_cache_scope<double> scope;
    line const first{{1048576,1048576}, {1048592,1048608}};
    line const second{{1048576,1048600}, {1048592,1048584}};
    ec::borrowed_context<root> context;
    ec::borrowed_stage<root> stage(context, first, second);
    std::array<root::operation_type,2> operations;
    operations[0].seg_id = bg::segment_identifier(0,-1,-1,0);
    operations[1].seg_id = bg::segment_identifier(1,-1,-1,0);
    ec::borrowed_record<root> record;
    record.stage = &stage;
    record.operations = {{&operations[0], &operations[1]}};
    ec::borrowed_point<root> const p(record);
    record.numeric.init(p.producers());
    point const adjacent(next(1048584),1048592), a(0,0), b(1,0);
    BOOST_TEST_EQ(ec::point_cache_detail::native_cache<double>::size(), 0u);
    allocation_guard = true;
    BOOST_TEST_EQ((ec::root_compare_coordinate<0,root>(p, adjacent)), -1);
    BOOST_TEST(!ec::root_equal<root>(adjacent,p));
    BOOST_TEST_EQ(ec::root_orientation<root>(a,b,p), 1);
    BOOST_TEST_EQ(ec::root_direction<root>(a,b,p), 1);
    auto const det = ec::root_signed_determinant<root>(a,b,p).bounds();
    auto const dist = ec::root_squared_distance<root>(a,p).bounds();
    allocation_guard = false;
    BOOST_TEST_EQ(det.approximate, 1048592.);
    BOOST_TEST_EQ(dist.approximate, 1048584.*1048584. + 1048592.*1048592.);
    BOOST_TEST(std::isfinite(det.error) && std::isfinite(dist.error));
    BOOST_TEST_EQ(ec::point_cache_detail::native_cache<double>::size(), 0u);

    // Both midpoints are exact rationals; ties round to the even heads.
    // The left query and the native value 1 have equal heads but differ.
    point const even(1,0), odd(next(1),0), even2(next(next(1)),0);
    using query = ec::query_point<root>;
    auto const left = query::midpoint(even,odd), right = query::midpoint(odd,even2);
    BOOST_TEST_EQ(left.approximate(0), 1);
    BOOST_TEST_EQ(right.approximate(0), bg::get<0>(even2));
    int result = 2;
    BOOST_TEST((!bp::coordinate_filter<0,root>(left,even,result)));
    BOOST_TEST_EQ(result, 2);
    BOOST_TEST((bp::coordinate_filter<0,root>(left,odd,result)));
    BOOST_TEST_EQ(result, -1);
    BOOST_TEST((bp::coordinate_filter<0,root>(right,odd,result)));
    BOOST_TEST_EQ(result, 1);
    // Adjacent uncertain heads prove strict order. A tied pair needs the
    // exact ordinates: the midpoints on either side of 1 both round to 1.
    auto const below = query::midpoint(point(previous(1),0), even);
    auto const adjacent_query = query::midpoint(odd, odd);
    BOOST_TEST_EQ(below.approximate(0), left.approximate(0));
    BOOST_TEST((!bp::coordinate_filter<0,root>(below,left,result)));
    BOOST_TEST((bp::coordinate_filter<0,root>(left,adjacent_query,result)));
    BOOST_TEST_EQ(result, -1);
    auto ordered = [](auto const& x, auto const& y, rational const& delta)
    {
        BOOST_TEST_EQ((ec::root_compare_coordinate<0,root>(x,y)), sign(delta));
        BOOST_TEST_EQ((ec::root_compare_coordinate<0,root>(y,x)), -sign(delta));
        BOOST_TEST_EQ(ec::root_equal<root>(x,y), delta == 0);
    };
    rational const half = (rational(1)+rational(next(1)))/2;
    ordered(below, left, (rational(previous(1))+1)/2 - half);
    ordered(left, adjacent_query, half-rational(next(1)));
    ordered(left, bg::model::d2::point_xy<float>(1,0), half-1);
    ordered(left, bg::model::d2::point_xy<int>(1,0), half-1);

    // At 2^53 both h and the integer h+1 convert to the query's nearest
    // double head. The wide native must retain its original integer value.
    std::int64_t const h = std::int64_t(1) << 53;
    auto const wide_mid = query::midpoint(point(double(h),0), point(double(h+2),0));
    BOOST_TEST_EQ(wide_mid.approximate(0), double(h));
    for (std::int64_t i : {h, h+1, h+2})
        ordered(wide_mid, bg::model::d2::point_xy<std::int64_t>(i,0), rational(h)+1-i);
}

int main()
{
    BOOST_TEST_EQ(std::fegetround(), FE_TONEAREST);
    static_assert(sizeof(coordinates) == 2 * sizeof(double), "no stored errors or hidden point state");
    proof_constants();
    nearest_order();
    randomized();
    using kinds = std::make_integer_sequence<unsigned,8>;
    producer_witnesses();
    cold_nearest();
    step_edges();
    bounded_filters<0>(kinds{});
    bounded_filters<1>(kinds{}); bounded_filters<2>(kinds{});
    bounded_filters<4>(kinds{}); bounded_filters<8>(kinds{}); bounded_filters<100>(kinds{});
    bounded_filters<(1u<<24)>(kinds{});
    bounded_filters<(1u<<24)+1>(kinds{});
    BOOST_TEST_EQ(allocations, 0u);
    return boost::report_errors();
}
