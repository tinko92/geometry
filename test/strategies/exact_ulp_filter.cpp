// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.

// From this worktree:
// clang++ -std=c++17 -O3 -DNDEBUG -ffp-contract=off -Iinclude -I../..
//   test/strategies/exact_ulp_filter.cpp -o /tmp/bg-ulp-filter
// /tmp/bg-ulp-filter
#include <boost/geometry/strategies/cartesian/exact/ulp_filter.hpp>
#include <boost/geometry/strategy/cartesian/side_robust.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/core/lightweight_test.hpp>

#include <array>
#include <cfenv>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <new>
#include <random>
#include <utility>
#include <vector>

namespace ec = boost::geometry::strategies::detail::exact_cartesian;
namespace uf = ec::ulp_filter;
using rational = ec::exact_number;
using coordinates = uf::coordinates;
using exact_coordinates = std::array<rational, 2>;

namespace
{
bool allocation_guard = false;
std::size_t allocations = 0;
std::array<std::size_t, 3> certified{}, uncertain{}, zeros{};
std::size_t comparisons = 0, coordinate_certificates = 0;
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

// Oracle data has independent exact errors on both sides, INCLUDING a=0.
// It is not restricted to nearest-rounded values or symmetric binade spacing.
rational implicit_value(double a, int numerator, unsigned denominator = 4)
{
    rational const q = ec::exact_value(a);
    rational const step = numerator < 0 ? q - ec::exact_value(previous(a))
                                        : ec::exact_value(next(a)) - q;
    return q + step * numerator / denominator;
}

void coordinate_check(double a, double b, rational const& exact_a, rational const& exact_b, bool first_native)
{
    int result = 17;
    bool equal = true;
    allocation_guard = true;
    bool const certain = first_native ? uf::try_compare_ni(a, b, result) : uf::try_compare_ii(a, b, result);
    bool const equality_certain = first_native ? uf::try_equal_ni(a, b, equal) : uf::try_equal_ii(a, b, equal);
    allocation_guard = false;
    int const expected = sign(exact_a - exact_b);
    ++comparisons;
    if (certain)
    {
        ++coordinate_certificates;
        BOOST_TEST_EQ(result, expected);
        BOOST_TEST(result != 0);
    }
    else BOOST_TEST_EQ(result, 17);
    BOOST_TEST_EQ(certain, equality_certain);
    if (equality_certain) { BOOST_TEST(!equal); BOOST_TEST(expected != 0); }
    else BOOST_TEST(equal);
}

int orientation_reference(exact_coordinates const& a, exact_coordinates const& b, exact_coordinates const& c)
{
    // Independent arbitrary-precision rational arithmetic, never the filter's
    // floating magnitude expression or its error coefficient.
    rational const determinant = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
    return sign(determinant);
}

template <unsigned Kind>
void orientation_check(coordinates const& a, coordinates const& b, coordinates const& c,
                       exact_coordinates const& ea, exact_coordinates const& eb, exact_coordinates const& ec_)
{
    int result = 17;
    allocation_guard = true;
    bool certain;
    if constexpr (Kind == 0) certain = uf::try_orientation_nni(a, b, c, result);
    else if constexpr (Kind == 1) certain = uf::try_orientation_nii(a, b, c, result);
    else certain = uf::try_orientation_iii(a, b, c, result);
    allocation_guard = false;
    if (certain)
    {
        ++certified[Kind];
        BOOST_TEST_EQ(result, orientation_reference(ea, eb, ec_));
        if (result == 0) ++zeros[Kind];
    }
    else { ++uncertain[Kind]; BOOST_TEST_EQ(result, 17); }
}

exact_coordinates native(coordinates const& p)
{ return {{ec::exact_value(p[0]), ec::exact_value(p[1])}}; }

exact_coordinates implicit(coordinates const& p, int x, int y)
{ return {{implicit_value(p[0], x), implicit_value(p[1], y)}}; }

void scalars()
{
    std::vector<double> values{0, -0.0, 1, -1,
        (std::numeric_limits<double>::denorm_min)(), -(std::numeric_limits<double>::denorm_min)(),
        (std::numeric_limits<double>::min)(), -(std::numeric_limits<double>::min)()};
    for (int e : {-120, -1, 0, 1, 120})
    {
        double const x = std::ldexp(1.0, e);
        for (double a : {previous(x), x, next(x)}) { values.push_back(a); values.push_back(-a); }
    }
    for (double a : values)
    {
        for (int offset : {-8, -1, 0, 1, 8})
        {
            double b = a;
            for (int i = 0; i < std::abs(offset); ++i) b = offset < 0 ? previous(b) : next(b);
            int comparison = 17;
            bool equal = false;
            BOOST_TEST(uf::try_compare_nn(a, b, comparison));
            BOOST_TEST_EQ(comparison, sign(ec::exact_value(a) - ec::exact_value(b)));
            BOOST_TEST(uf::try_equal_nn(a, b, equal));
            BOOST_TEST_EQ(equal, a == b);
            for (int eb : {-4, -1, 0, 1, 4})
            {
                coordinate_check(a, b, ec::exact_value(a), implicit_value(b, eb), true);
                for (int ea : {-4, -1, 0, 1, 4})
                    coordinate_check(a, b, implicit_value(a, ea), implicit_value(b, eb), false);
            }
        }
    }
    // Overlapping uncertainty must never become exact equality, notably at zero.
    int result = 17;
    BOOST_TEST(!uf::try_compare_ni(0, 0, result));
    BOOST_TEST_EQ(result, 17);
    BOOST_TEST(!uf::try_compare_ii(0, 0, result));
    BOOST_TEST_EQ(result, 17);
    BOOST_TEST(!uf::try_compare_ii(1, 1, result));
    double const maximum = (std::numeric_limits<double>::max)();
    BOOST_TEST(!uf::try_compare_ni(0, maximum, result));
    BOOST_TEST_EQ(result, 17);
    BOOST_TEST(!uf::try_compare_ii(maximum, previous(maximum), result));
}

void adversarial()
{
    for (int exponent : {-120, -1, 0, 1, 120})
    {
        double const s = std::ldexp(1.0, exponent);
        std::array<std::array<coordinates, 3>, 6> const triples{{
            {{{{0, 0}}, {{s, 0}}, {{0, s}}}},
            {{{{0, -s}}, {{0, s}}, {{0, 0}}}},
            {{{{-s, -s}}, {{s, s}}, {{0, 0}}}},
            {{{{s, s}}, {{next(s), next(s)}}, {{next(next(s)), next(next(s))}}}},
            {{{{s, s}}, {{s, s}}, {{previous(s), next(s)}}}},
            {{{{0, 0}}, {{0, 0}}, {{0, 0}}}}
        }};
        for (auto const& t : triples)
        for (int ex : {-4, 0, 4}) for (int ey : {-4, 0, 4})
        {
            auto const na = native(t[0]), nb = native(t[1]);
            auto const ia = implicit(t[0], ex, ey), ib = implicit(t[1], ey, -ex);
            auto const ic = implicit(t[2], -ex, -ey);
            orientation_check<0>(t[0], t[1], t[2], na, nb, ic);
            orientation_check<0>(t[1], t[0], t[2], nb, na, ic);
            orientation_check<1>(t[0], t[1], t[2], na, ib, ic);
            orientation_check<1>(t[0], t[2], t[1], na, ic, ib);
            orientation_check<2>(t[0], t[1], t[2], ia, ib, ic);
            orientation_check<2>(t[1], t[2], t[0], ib, ic, ia);
            orientation_check<2>(t[2], t[0], t[1], ic, ia, ib);
            orientation_check<2>(t[1], t[0], t[2], ib, ia, ic);
        }
    }
    // Equal approximate zero is not exact zero: the oracle realizes both signs.
    coordinates const a{{0, -1}}, b{{0, 1}}, c{{0, 0}};
    for (int ex : {-4, 4})
    {
        auto const exact_c = implicit(c, ex, 0);
        BOOST_TEST(orientation_reference(native(a), native(b), exact_c) != 0);
        int result = 17;
        BOOST_TEST(!uf::try_orientation_nni(a, b, c, result));
        BOOST_TEST_EQ(result, 17);
    }

    // An implicit y=0 can contribute a NORMAL determinant through a huge
    // native x. Treating its uncertainty as eps*abs(0)=0 gets the sign wrong.
    coordinates const large{{std::ldexp(1.0, 1000), 1}}, origin{{0, 0}};
    coordinates const small{{std::ldexp(1.0, -75), 0}};
    auto const eo = native(origin), el = native(large);
    for (int error : {-4, 4})
    {
        auto const es = implicit(small, 0, error);
        BOOST_TEST_EQ(orientation_reference(eo, el, es), error < 0 ? -1 : 1);
        int result = 17;
        BOOST_TEST(!uf::try_orientation_nni(origin, large, small, result));
        BOOST_TEST_EQ(result, 17);
        orientation_check<1>(origin, large, small, eo, implicit(large, 0, 0), es);
        orientation_check<2>(origin, large, small, implicit(origin, 0, 0), implicit(large, 0, 0), es);
    }

    // A subnormal magnitude bound is rounded, not silently flushed to zero.
    // The value product here is exact; no inexact underflow is required.
    coordinates const unit{{1, 0}};
    coordinates const tiny{{0, 8 * (std::numeric_limits<double>::denorm_min)()}};
    orientation_check<0>(origin, unit, tiny, eo, native(unit), implicit(tiny, 4, -4));
}

void proof_constants()
{
    rational const u = ec::exact_value(uf::detail::unit_roundoff);
    rational const one_minus_u = 1 - u;
    rational const square = one_minus_u * one_minus_u;
    BOOST_TEST(3 * square > 2);
    BOOST_TEST(7 * square * square * square > (6 + 5 * u) * (1 + 2 * u));
    rational const eta = ec::exact_value((std::numeric_limits<double>::denorm_min)());
    BOOST_TEST_EQ(2 * u * ec::exact_value(uf::detail::minimum_normal), eta);
}

void randomized()
{
    std::mt19937_64 random(0x554c5046494c5445);
    auto value = [&]()
    {
        double const significand = 1 + std::ldexp(double(random() >> 11), -53);
        return std::ldexp(random() & 1 ? significand : -significand, int(random() % 241) - 120);
    };
    auto error = [&]() { return int(random() % 9) - 4; };
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
        auto const na = native(a), nb = native(b);
        auto const ia = implicit(a, error(), error());
        auto const ib = implicit(b, error(), error());
        auto const ic = implicit(c, error(), error());
        orientation_check<0>(a, b, c, na, nb, ic);
        orientation_check<1>(a, b, c, na, ib, ic);
        orientation_check<2>(a, b, c, ia, ib, ic);
        orientation_check<2>(c, a, b, ic, ia, ib);
        orientation_check<2>(a, c, b, ia, ic, ib);
        coordinate_check(a[0], c[0], na[0], ic[0], true);
        coordinate_check(a[1], c[1], ia[1], ic[1], false);

        // NNN must continue to use BG's existing adaptive Shewchuk predicate.
        using point = boost::geometry::model::d2::point_xy<double>;
        using robust = boost::geometry::strategy::side::side_robust<double,
            boost::geometry::strategy::side::fp_equals_policy, 3>;
        int const actual = robust::apply(point(a[0], a[1]), point(b[0], b[1]), point(c[0], c[1]));
        BOOST_TEST_EQ(actual, orientation_reference(na, nb, native(c)));
    }
}

// Independently generate the extreme points allowed by the public certificate:
// move N binary64 steps, then almost half the next rounding cell. In particular
// these errors need not satisfy the tighter bound from our FMA construction.
template <unsigned N>
rational step_extreme(double head, bool up)
{
    for (unsigned i = 0; i < N; ++i) head = up ? next(head) : previous(head);
    rational const r = ec::exact_value(head);
    return r + (ec::exact_value(up ? next(head) : previous(head)) - r) * 511 / 1024;
}

template <unsigned N, unsigned Mask>
void bounded_triple(coordinates a, coordinates b, coordinates p)
{
    constexpr bool ai = (Mask & 1) != 0, bi = (Mask & 2) != 0, pi = (Mask & 4) != 0;
    auto make = [](coordinates z, bool uncertain, bool up) -> exact_coordinates {
        if (!uncertain) return native(z);
        return {{step_extreme<N>(z[0], up), step_extreme<N>(z[1], !up)}};
    };
    auto const x = make(a, ai, true), y = make(b, bi, false), z = make(p, pi, true);
    rational const dx = y[0]-x[0], dy = y[1]-x[1];
    rational const det = dx*(z[1]-x[1]) - dy*(z[0]-x[0]);
    rational const dot = dx*(z[0]-y[0]) + dy*(z[1]-y[1]);
    rational const dist = dx*dx + dy*dy;
    ec::forward_error::scalar_approximation d{}, s{};
    int side = 2, direction = 2, cmp = 2;
    allocation_guard = true;
    bool const bounded_d = uf::try_signed_determinant_n<N,ai,bi,pi>(a,b,p,d);
    bool const bounded_s = uf::try_squared_distance_n<N,ai,bi>(a,b,s);
    bool const cs = uf::try_orientation_n<N,ai,bi,pi>(a,b,p,side);
    bool const cd = uf::try_dot_n<N,ai,bi,pi>(a,b,p,direction);
    bool const cc = uf::try_compare_n<N,ai,bi>(a[0],b[0],cmp);
    allocation_guard = false;
    BOOST_TEST(bounded_d && bounded_s);
    BOOST_TEST(abs(det - ec::exact_value(d.value)) <= ec::exact_value(d.error));
    BOOST_TEST(abs(dist - ec::exact_value(s.value)) <= ec::exact_value(s.error));
    BOOST_TEST_EQ(side, cs ? sign(det) : 2);
    BOOST_TEST_EQ(direction, cd ? sign(dot) : 2);
    BOOST_TEST_EQ(cmp, cc ? sign(x[0] - y[0]) : 2);
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
    }
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

int main()
{
    BOOST_TEST_EQ(std::fegetround(), FE_TONEAREST);
    static_assert(sizeof(coordinates) == 2 * sizeof(double), "no stored errors or hidden point state");
    proof_constants();
    scalars();
    adversarial();
    randomized();
    using kinds = std::make_integer_sequence<unsigned,8>;
    bounded_filters<1>(kinds{}); bounded_filters<2>(kinds{});
    bounded_filters<4>(kinds{}); bounded_filters<8>(kinds{}); bounded_filters<100>(kinds{});
    BOOST_TEST_EQ(allocations, 0u);
    BOOST_TEST(coordinate_certificates > 0);
    BOOST_TEST(zeros[0] > 0);
    BOOST_TEST_EQ(zeros[1], 0u);
    BOOST_TEST_EQ(zeros[2], 0u);
    std::cout << comparisons << " coordinate cases; " << coordinate_certificates << " certified.\n";
    for (unsigned i = 0; i < 3; ++i)
    {
        BOOST_TEST(certified[i] > 0);
        BOOST_TEST(uncertain[i] > 0);
        std::cout << (i == 0 ? "NNI" : i == 1 ? "NII" : "III") << ": "
                  << certified[i] << " certified, " << uncertain[i] << " uncertain, "
                  << zeros[i] << " structural zeros.\n";
    }
    return boost::report_errors();
}
