// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.

// clang++ -std=c++17 -O3 -DNDEBUG -ffp-contract=off -Iinclude -I../..
//   test/strategies/exact_construction_rounding.cpp -o /tmp/bg-construction-rounding
#include <boost/geometry/strategies/cartesian/exact/construction_rounding.hpp>
#include <boost/geometry/strategies/cartesian/exact/ulp_filter.hpp>
#include <boost/geometry/strategy/cartesian/side_robust.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/core/lightweight_test.hpp>

#include <cfenv>
#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <new>
#include <random>
#include <string>
#include <utility>
#include <vector>

namespace bg = boost::geometry;
namespace ec = bg::strategies::detail::exact_cartesian;
namespace cr = ec::construction_rounding;
namespace uf = ec::ulp_filter;
namespace policy = bg::strategies::relate;
using rational = boost::multiprecision::cpp_rational;
using exact_point = std::array<rational, 2>;
using producers = std::array<double, 8>;
using ordinary = cr::point_state<double, policy::approximate_output>;
using nearest_state = cr::point_state<double, policy::nearest_output>;
static_assert(std::is_same<policy::eager_nearest_output, policy::nearest_output>::value,
              "The compatibility spelling is the canonical policy");

template <typename T, typename = void> struct has_error : std::false_type {};
template <typename T> struct has_error<T, std::void_t<decltype(std::declval<T>().error)>> : std::true_type {};
static_assert(has_error<ordinary::storage_type>::value, "Ordinary bounds must be retained");
static_assert(!has_error<nearest_state::storage_type>::value, "Nearest must not store errors");
constexpr std::size_t lean_payload = 2 * sizeof(double) + sizeof(ec::point_cache_detail::native_cache_handle);
constexpr auto state_alignment = alignof(nearest_state);
static_assert(sizeof(nearest_state) == (lean_payload + state_alignment - 1) / state_alignment * state_alignment,
              "Lean state must contain only two approximations and a plain cache index");
static_assert(sizeof(ordinary) == sizeof(nearest_state) + 2 * sizeof(double), "Only ordinary stores errors");
static_assert(std::is_trivially_copyable<nearest_state>::value && std::is_trivially_copyable<ordinary>::value,
              "No producer ownership or reference counting in numerical state");

namespace
{
bool counting = false;
std::size_t allocations = 0, checked = 0;
std::size_t forward_certificates = 0;
std::size_t small_integer_certificates = 0;
std::size_t ordinary_certificates = 0;
}

void* operator new(std::size_t size)
{
    if (counting) ++allocations;
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

template <typename F> void without_allocations(F function)
{
    allocations = 0;
    counting = true;
    try { function(); }
    catch (...) { counting = false; throw; }
    counting = false;
    BOOST_TEST_EQ(allocations, 0u);
}

template <typename T> using word = typename std::conditional<sizeof(T) == 4, std::uint32_t, std::uint64_t>::type;
template <typename T> word<T> bits(T value)
{
    word<T> result;
    std::memcpy(&result, &value, sizeof(value));
    return result;
}
template <typename T> T previous(T value) { return std::nextafter(value, -(std::numeric_limits<T>::infinity)()); }
template <typename T> T next(T value) { return std::nextafter(value, (std::numeric_limits<T>::infinity)()); }
rational absolute(rational const& value) { return value < 0 ? rational(-value) : value; }

// The oracle crosses rational line equations, independently of the production
// segment-parameter forward bound, floating expansions and fixed integer code.
template <typename C> bool oracle(std::array<C, 8> const& v, exact_point& output)
{
    std::array<rational, 8> p;
    for (std::size_t i = 0; i < 8; ++i) p[i] = rational(v[i]);
    rational const a = p[1] - p[3], b = p[2] - p[0], c = p[0] * p[3] - p[1] * p[2];
    rational const d = p[5] - p[7], e = p[6] - p[4], f = p[4] * p[7] - p[5] * p[6];
    rational const w = a * e - b * d;
    if (w == 0) return false;
    output = {{(b * f - c * e) / w, (c * d - a * f) / w}};
    return true;
}

// convert_to only supplies a candidate. Rational midpoint comparisons, not
// that converter's rounding, determine the reference answer and tie ownership.
template <typename T> T nearest(rational const& exact)
{
    if (exact == 0) return T(0);
    T candidate = exact.convert_to<T>();
    for (unsigned i = 0; i < 16; ++i)
    {
        rational const lower = (rational(previous(candidate)) + rational(candidate)) / 2;
        rational const upper = (rational(candidate) + rational(next(candidate))) / 2;
        bool const even = (bits(candidate) & 1) == 0;
        if (exact < lower || (exact == lower && !even)) candidate = previous(candidate);
        else if (exact > upper || (exact == upper && !even)) candidate = next(candidate);
        else return candidate;
    }
    BOOST_ERROR("Rational rounding oracle did not converge");
    return candidate;
}

template <typename C, typename Policy>
void check_policy(std::array<C, 8> const& v, exact_point const& expected)
{
    using state = cr::point_state<C, Policy>;
    using traits = cr::policy_traits<Policy>;
    ec::point_cache_scope<C> scope;
    state point;
    ec::forward_error::point_approximation native;
    bool native_ready = false;
    if constexpr (traits::store_errors)
        native_ready = ec::forward_error::try_intersection(v, native);
    if (native_ready)
    {
        without_allocations([&] { point.init(v); });
        BOOST_TEST(ec::operation_detail::current<C>() == nullptr);
        ++ordinary_certificates;
        for (std::size_t d = 0; d < 2; ++d)
        {
            double const value = native.value[d] == 0 && native.error[d] == 0 ? 0 : native.value[d];
            BOOST_TEST_EQ(bits(point.approximate(d)), bits(value));
            BOOST_TEST_EQ(bits(point.bounds(d).error), bits(native.error[d]));
        }
        ec::forward_error::point_approximation tight;
        if (cr::detail::forward_bound(v, tight))
            for (std::size_t d = 0; d < 2; ++d)
                BOOST_TEST_EQ(bits(point.approximate(d)), bits(tight.value[d]));
    }
    else point.init(v);
    for (std::size_t d = 0; d < 2; ++d)
    {
        double const first = point.approximation(d).approximate;
        if constexpr (traits::store_errors)
        {
            auto const bound = point.bounds(d);
            BOOST_TEST(absolute(rational(first) - expected[d]) <= rational(bound.error));
            BOOST_TEST_EQ(bits(point.template emit<double>(d, v)), bits(first));
        }
        else
        {
            BOOST_TEST_EQ(bits(first), bits(nearest<double>(expected[d])));
        }
        // Integral exports must not refine, allocate a precision entry, or
        // apply the exact rational's nearest-integer rule.
        if (std::abs(first) < 1e9)
        {
            auto const before = ec::point_cache_detail::native_cache<C>::size();
            BOOST_TEST_EQ(point.template emit<int>(d, v), bg::util::numeric_cast<int>(first));
            BOOST_TEST_EQ(bits(point.approximation(d).approximate), bits(first));
            BOOST_TEST_EQ(ec::point_cache_detail::native_cache<C>::size(), before);
        }
        state const& lazy = point;
        double const refined = lazy.refine(v, d).approximate;
        BOOST_TEST_EQ(bits(refined), bits(nearest<double>(expected[d])));
        BOOST_TEST_EQ(bits(point.approximation(d).approximate), bits(refined));
        BOOST_TEST_EQ(bits(point.template emit<double>(d, v)), bits(refined));
        auto const bound = point.bounds(d);
        BOOST_TEST(absolute(rational(bound.approximate) - expected[d]) <= rational(bound.error));
        if constexpr (traits::nearest_output)
        {
            BOOST_TEST_EQ(bits(lazy.template emit<double>(d, v)), bits(nearest<double>(expected[d])));
            if (expected[d] == 0 || (absolute(expected[d]) >= rational(2 * (std::numeric_limits<float>::min)())
                && absolute(expected[d]) < rational((std::numeric_limits<float>::max)()) / 2))
                BOOST_TEST_EQ(bits(lazy.template emit<float>(d, v)), bits(nearest<float>(expected[d])));
        }
    }
    ++checked;
}

template <typename C> void check(std::array<C, 8> const& v)
{
    exact_point expected;
    if (!oracle(v, expected)) return;
    ec::forward_error::point_approximation bound;
    if (cr::detail::forward_bound(v, bound))
    {
        for (std::size_t d = 0; d < 2; ++d)
        {
            BOOST_TEST(absolute(rational(bound.value[d]) - expected[d]) <= rational(bound.error[d]));
            if (cr::try_certify_nearest({bound.value[d], bound.error[d]}))
            {
                BOOST_TEST_EQ(bits(bound.value[d]), bits(nearest<double>(expected[d])));
                ++forward_certificates;
            }
        }
    }
    check_policy<C, policy::approximate_output>(v, expected);
    check_policy<C, policy::nearest_output>(v, expected);
}

void certificates()
{
    BOOST_TEST(cr::try_certify_nearest({0, 0}));
    BOOST_TEST(!cr::try_certify_nearest({0, (std::numeric_limits<double>::min)()}));
    BOOST_TEST(!cr::try_certify_nearest({1, -1}));
    BOOST_TEST(!cr::try_certify_nearest({(std::numeric_limits<double>::infinity)(), 0}));
    BOOST_TEST(!cr::try_certify_nearest({1, (std::numeric_limits<double>::quiet_NaN)()}));
    for (double value : {-2.0, -next(1.0), -1.0, previous(1.0), 1.0, next(1.0), 2.0})
    {
        double const step = (std::min)(value - previous(value), next(value) - value);
        BOOST_TEST_EQ(cr::try_certify_nearest({value, step / 2}), (bits(value) & 1) == 0);
        BOOST_TEST(!cr::try_certify_nearest({value, next(step / 2)}));
    }
    // The narrower side of the binade must control the nearest certificate.
    BOOST_TEST(!cr::try_certify_nearest({1, (next(1.0) - 1) / 2}));
}

void boundaries()
{
    for (int exponent : {-120, -1, 0, 1, 120})
    for (double sign : {-1.0, 1.0})
    for (double lo : {previous(1.0), 1.0, next(1.0)})
    for (double height : {0.5, 1.0, 1.5})
    {
        producers v{{sign * lo, 0, sign * next(lo), 2, 0, height, 3, height}};
        for (double& x : v) x = std::ldexp(x, exponent);
        for (unsigned mask = 0; mask < 8; ++mask)
        {
            auto w = v;
            if (mask & 1) { std::swap(w[0], w[2]); std::swap(w[1], w[3]); }
            if (mask & 2) { std::swap(w[4], w[6]); std::swap(w[5], w[7]); }
            if (mask & 4) for (std::size_t i = 0; i < 4; ++i) std::swap(w[i], w[i + 4]);
            check(w);
        }
    }
    check(producers{{0, -1, 0, 1, -1, 0, 1, 0}});
    check(producers{{-0.0, 0, -0.0, 1, -1, -1, 1, -1}});
    check(producers{{-2, -1, 2, 1, -1, 1, 1, -1}});
    for (int exponent : {-200, 200})
    {
        double const s = std::ldexp(1.0, exponent);
        check(producers{{0, 0, s, s, 0, s, s, 0}});
    }
    double const n = std::ldexp(1.0, 27);
    check(producers{{0, 0, n + 1, n, 0, 1, n, n}}); // Rounded D=0, exact D=-1.
    {
        // T=1/6 rounds down. Its product with a three-ULP direction puts
        // the FMA just below a midpoint whose correct ties-even owner is UP.
        // Check the constructed head before any emission or predicate can refine it.
        double const odd = next(1.0);
        producers const v{{odd, 0, next(next(next(odd))), 3, 0, 0.5, 2, 0.5}};
        nearest_state point;
        point.init(v);
        BOOST_TEST_EQ(bits(std::fma(1. / 6, v[2] - odd, odd)), bits(odd));
        BOOST_TEST_EQ(bits(point.approximate(0)), bits(next(odd)));
        BOOST_TEST_EQ(bits(point.emit<double>(0, v)), bits(next(odd)));
        check(v);
    }
    for (double sign : {-1.0, 1.0})
    for (double direction : {-1.0, 1.0})
    {
        // The exact x is just to one side of a float midpoint, but rounds to
        // that midpoint in binary64. A second rounding would choose the wrong
        // float in two of these four cases.
        double const midpoint = sign * (1 + std::ldexp(1.0, -24));
        double const adjacent = direction < 0 ? previous(midpoint) : next(midpoint);
        producers const v{{midpoint, 0, adjacent, 1, 0, 0.25, 2, 0.25}};
        exact_point expected;
        BOOST_TEST(oracle(v, expected));
        check(v);
        nearest_state point;
        point.init(v);
        BOOST_TEST_EQ(bits(point.emit<float>(0, v)), bits(nearest<float>(expected[0])));
        BOOST_TEST_EQ(bits(point.approximation(0).approximate), bits(nearest<double>(expected[0])));
        if ((sign > 0) == (direction > 0))
            BOOST_TEST(bits(static_cast<float>(nearest<double>(expected[0]))) != bits(nearest<float>(expected[0])));
    }
}

void random_cases()
{
    std::mt19937_64 random(0x434f4e5354525543);
    auto fraction = [&] { return std::ldexp(double(random() >> 11), -53); };
    for (unsigned i = 0; i < 512; ++i)
    {
        double const x = double(std::int64_t(random() % 2048) - 1024) + fraction();
        double const y = double(std::int64_t(random() % 2048) - 1024) + fraction();
        double const radius = 1 + double(random() % 64) + fraction();
        producers v{{x - radius, y - radius, x + radius, y + radius + 0.25,
                     x - radius, y + radius, x + radius, y - radius - 0.5}};
        if (i & 1) for (double& coordinate : v) coordinate += (i & 2 ? 4096 : 1048576);
        check(v);
        if (i < 64)
        {
            double const e = std::ldexp(1.0, -48 - int(random() % 4));
            check(producers{{0, 0, 1, 1, 0, e, 1, 1 - e}});
        }
    }
    BOOST_TEST(forward_certificates > 100);
}

std::uint64_t ulp_distance(double a, double b)
{
    auto rank = [](double value)
    {
        constexpr std::uint64_t sign = std::uint64_t(1) << 63;
        auto const word = bits(value), magnitude = word & ~sign;
        return word & sign ? sign - magnitude : sign + magnitude;
    };
    auto const x = rank(a), y = rank(b);
    return x < y ? y - x : x - y;
}

struct bounded_counts
{
    std::size_t accepted = 0, rescued = 0, different = 0, beyond_one = 0;
};

std::vector<producers> bounded_corpus()
{
    std::mt19937_64 random(0x424f554e444544);
    std::vector<producers> inputs(512);
    for (std::size_t i = 0; i < inputs.size(); ++i)
    {
        for (double& x : inputs[i])
            x = double(std::int64_t(random() % 100001) - 50000) / 1000;
        if (i % 3 == 0) for (double& x : inputs[i]) x += 0x1p20;
    }
    return inputs;
}

template <std::size_t N, bool Retry, bool DeterminantRetry = policy::bounded_ulp_output<N, Retry>::determinant_retry>
void check_bounded(producers const& input, exact_point const& exact,
                   std::array<double, 2> const& expected, bounded_counts& counts)
{
    using state = cr::point_state<double, policy::bounded_ulp_output<N, Retry, DeterminantRetry>>;
    using cache = ec::point_cache_detail::native_cache<double>;
    static_assert(sizeof(state) == sizeof(ordinary), "Reuse the stored-error layout");
    static_assert(std::is_trivially_copyable<state>::value, "No producer ownership");
    ec::forward_error::point_approximation first, tighter;
    auto const plain = cr::detail::bounded_intersection<N, false, true, true, true, DeterminantRetry>(input, first);
    auto const accepted = cr::detail::bounded_intersection<N, Retry, true, true, true, DeterminantRetry>(input, tighter);
    ec::point_cache_scope<double> scope;
    state point;
    auto const init = [&] { point.init(input); };
    if (accepted[0] && accepted[1])
    {
        without_allocations(init);
        BOOST_TEST(ec::operation_detail::current<double>() == nullptr);
    }
    else
    {
        init();
        BOOST_TEST_EQ(cache::size(), 1u);
    }
    auto const entries = cache::size();
    for (std::size_t d = 0; d < 2; ++d)
    {
        auto const bound = point.bounds(d);
        auto const emitted = point.template emit<double>(d, input);
        BOOST_TEST_EQ(bits(emitted), bits(bound.approximate));
        auto const steps = ulp_distance(emitted, expected[d]);
        BOOST_TEST(steps <= N);
        BOOST_TEST(absolute(rational(emitted) - exact[d]) <= rational(bound.error));
        if (plain[d]) BOOST_TEST(accepted[d] && bits(first.value[d]) == bits(tighter.value[d]));
        if (accepted[d])
        {
            BOOST_TEST_EQ(bits(emitted), bits(tighter.value[d]));
            ++counts.accepted;
            counts.rescued += !plain[d];
            counts.different += steps != 0;
            counts.beyond_one += steps > 1;
        }
        else BOOST_TEST_EQ(bits(emitted), bits(expected[d]));
    }
    BOOST_TEST_EQ(cache::size(), entries);
    for (std::size_t d = 0; d < 2; ++d)
    {
        auto const refined = point.refine(input, d);
        BOOST_TEST_EQ(bits(refined.approximate), bits(expected[d]));
        BOOST_TEST(absolute(rational(refined.approximate) - exact[d]) <= rational(refined.error));
    }
}

void bounded_cases()
{
    std::array<bounded_counts, 8> counts{};
    std::size_t cases = 0;
    auto check_input = [&](producers const& input)
    {
        exact_point exact;
        if (!oracle(input, exact)) return;
        std::array<double, 2> const expected{{nearest<double>(exact[0]), nearest<double>(exact[1])}};
        check_bounded<1, false>(input, exact, expected, counts[0]);
        check_bounded<1, true>(input, exact, expected, counts[1]);
        check_bounded<4, false>(input, exact, expected, counts[2]);
        check_bounded<4, true>(input, exact, expected, counts[3]);
        check_bounded<100, false>(input, exact, expected, counts[4]);
        check_bounded<100, true>(input, exact, expected, counts[5]);
        check_bounded<100, false, false>(input, exact, expected, counts[6]);
        check_bounded<100, true, false>(input, exact, expected, counts[7]);
        ++cases;
    };
    // Exact t=1/2. The default selective retry now certifies x as well as y.
    producers const rescue{{1.6875, 9, 2.6875, 10, 1.6875, 10, 3.6875, 8}};
    ec::forward_error::point_approximation first;
    auto const plain = cr::detail::bounded_intersection<4, false>(rescue, first);
    auto const retry = cr::detail::bounded_intersection<4, true>(rescue, first);
    BOOST_TEST(plain[0] && plain[1] && retry == plain);
    check_input(rescue);

    for (double sign : {-1.0, 1.0})
    for (double a : {previous(1.0), 1.0, next(1.0)})
    {
        // Place x at a rounding boundary on both sides of a binade.
        check_input({{sign * a, 0, sign * next(next(next(a))), 3, 0, 0.5, 2, 0.5}});
    }
    check_input({{0, -1, 0, 1, -1, 0, 1, 0}}); // Zero candidates.
    double const n = 0x1p27;
    check_input({{0, 0, n + 1, n, 0, 1, n, n}}); // Rounded denominator zero.
    producers const near_parallel{{0, 0, n + 1, n, 0, 1, n - 2, n - 2}};
    auto const heads = ec::forward_error::detail::intersection_heads(near_parallel, false);
    BOOST_TEST(heads.denominator != 0);
    auto const rejected = cr::detail::bounded_intersection<100, true, false>(near_parallel, first);
    BOOST_TEST(!rejected[0] && !rejected[1]);
    check_input(near_parallel);
    for (int exponent : {-700, -64, 60, 700})
    {
        auto input = rescue;
        for (double& x : input) x = std::ldexp(x, exponent);
        check_input(input);
    }
    auto const tiny = (std::numeric_limits<double>::denorm_min)();
    check_input({{0, 0, tiny, 0, 0, -tiny, tiny, 3 * tiny}});

    for (auto const& input : bounded_corpus()) check_input(input);
    BOOST_TEST(counts[0].accepted > 100 && counts[0].different > 0);
    BOOST_TEST_EQ(counts[3].accepted, counts[2].accepted);
    BOOST_TEST(counts[4].accepted > counts[0].accepted && counts[4].beyond_one > 0);
    BOOST_TEST(counts[7].rescued > 0); // Residual-only retry remains independently selectable.
    for (std::size_t i = 0; i < counts.size(); ++i)
        std::cout << "bounded " << (i < 2 ? 1 : i < 4 ? 4 : 100)
                  << (i == 6 ? " no retry: " : i == 7 ? " residual only: "
                      : i % 2 ? " determinant (residual flag set): " : " determinant: ")
                  << counts[i].accepted << '/' << 2 * cases << " axes accepted, "
                  << counts[i].rescued << " rescued, " << counts[i].different
                  << " non-nearest, " << counts[i].beyond_one << " beyond one step\n";
}

void bounded_other_formats()
{
    using bounded = policy::bounded_ulp_output<100, true>;
    // A double midpoint cannot be used as a proxy for exact float rounding.
    double const midpoint = 1 + 0x1p-24;
    producers const input{{midpoint, 0, next(midpoint), 1, 0, 0.25, 2, 0.25}};
    cr::point_state<double, bounded> a;
    a.init(input);
    BOOST_TEST_EQ(a.emit<double>(0, input), midpoint);
    BOOST_TEST_EQ(a.emit<float>(0, input), next(1.0f));
    BOOST_TEST_EQ(a.emit<int>(0, input), 1);

    std::array<float, 8> const f{{0.25f, -1.5f, 1.75f, 2.25f, -2.5f, 0.75f, 3.25f, -1.25f}};
    std::int64_t const base = (std::int64_t(1) << 54) + 1;
    std::array<std::int64_t, 8> const j{{base, base, base + 17, base + 31,
                                        base, base + 25, base + 17, base + 2}};
    auto check_input = [](auto const& values)
    {
        using C = typename std::decay_t<decltype(values)>::value_type;
        exact_point exact;
        BOOST_TEST(oracle(values, exact));
        cr::point_state<C, bounded> p;
        p.init(values);
        for (std::size_t d = 0; d < 2; ++d)
        {
            BOOST_TEST_EQ(p.template emit<double>(d, values), nearest<double>(exact[d]));
            BOOST_TEST_EQ(p.template emit<float>(d, values), nearest<float>(exact[d]));
        }
    };
    check_input(f);
    check_input(j);
}

void determinant_budgets()
{
    namespace ar = ec::arithmetic_detail;
    double const u = ec::forward_error::detail::u;
    auto check = [&](double ax, double ay, double bx, double by)
    {
        double const k = ar::rounded_sub(ar::rounded_mul(ax, by), ar::rounded_mul(ay, bx));
        rational const exact = rational(ax) * rational(by) - rational(ay) * rational(bx);
        double const budget = cr::detail::determinant_budget(ax, ay, bx, by, k);
        BOOST_TEST(absolute(exact - rational(k)) <= rational(u) * rational(budget));
        auto const rebuilt = cr::detail::rebuild_determinant(ax, ay, bx, by, k);
        BOOST_TEST(absolute(exact - rational(rebuilt.value)) <= rational(rebuilt.error));
        return budget;
    };
    double const n = 0x1p27;
    // Equal rounded products can conceal a nonzero determinant.
    BOOST_TEST(check(n + 1, n, n, n - 1) > 0);
    BOOST_TEST_EQ(check(1.25, 2, 4, 5.5), 0.);
    std::mt19937_64 random(0x4445544255444745);
    auto value = [&]
    {
        double const x = 1 + double(random() >> 12) * 0x1p-52;
        int const exponent = int(random() % 128) - 64;
        return std::ldexp((random() & 1) ? -x : x, exponent);
    };
    for (unsigned i = 0; i < 512; ++i) check(value(), value(), value(), value());
}

void selective_one_ulp_retry()
{
    namespace ar = ec::arithmetic_detail;
    namespace fe = ec::forward_error;
    constexpr double u = fe::detail::u;
    auto check = [](producers const& v, std::array<bool, 2> expected,
                    bool numerator_tightens, bool denominator_tightens)
    {
        auto const heads = fe::detail::intersection_heads(v, false);
        double const et = 8 * heads.mt, ed = 8 * heads.md;
        double const lower = (1 - 2 * u) * (std::abs(heads.denominator) - u * ed);
        BOOST_TEST(lower > 0);
        double const t = ar::rounded_div(heads.numerator, heads.denominator), at = std::abs(t);
        double const threshold = (1 - 18 * u) * lower;
        std::array<double, 2> const z{{std::fma(t, heads.ux, v[0]), std::fma(t, heads.uy, v[1])}};
        double const h = std::fma(at, lower, std::fma(8, std::fma(at, heads.md, heads.mt),
                                                       std::abs(heads.numerator)));
        std::array<bool, 2> const initial{{cr::detail::bounded_axis_passes(heads.ux, z[0], h, threshold),
                                          cr::detail::bounded_axis_passes(heads.uy, z[1], h, threshold)}};
        BOOST_TEST(!initial[0]);
        // In these fixtures every tightened determinant is exact. Checking
        // the returned height distinguishes doing one, both or no retries.
        double const residual = cr::detail::scaled_division_residual(t, heads.denominator, heads.numerator);
        double const want = std::fma(at, lower, ar::rounded_add(residual,
            std::fma(at, denominator_tightens ? 0. : ed, numerator_tightens ? 0. : et)));
        BOOST_TEST_EQ(cr::detail::retry_determinants<false>(v, heads, t, lower, threshold, z, initial), want);

        fe::point_approximation result{{{1, 2}}, {{3, 4}}};
        auto const passed = cr::detail::bounded_intersection<1, false, true, false>(v, result);
        BOOST_TEST(passed == expected);
        exact_point exact;
        BOOST_TEST(oracle(v, exact));
        for (std::size_t d = 0; d < 2; ++d)
        {
            if (!passed[d])
            {
                BOOST_TEST_EQ(result.value[d], double(d + 1));
                BOOST_TEST_EQ(result.error[d], double(d + 3));
                continue;
            }
            BOOST_TEST_EQ(bits(result.value[d]), bits(z[d])); // No head/candidate replacement.
            BOOST_TEST(ulp_distance(result.value[d], nearest<double>(exact[d])) <= 1);
            BOOST_TEST(absolute(rational(result.value[d]) - exact[d]) <= rational(result.error[d]));
        }
    };
    // Both budgets block x; y already passes.
    check({{1.6875, 9, 2.6875, 10, 1.6875, 10, 3.6875, 8}}, {{true, true}}, true, true);
    // Only the dominant budget needs tightening. Retain the other one.
    check({{5.5, 9, 6.5, 10, 5.5, 10, 7.5, 8}}, {{true, true}}, true, false);
    check({{5, 10, 6, 11, 5, 10.25, 6, 10.75}}, {{true, true}}, false, true);
    // Using the exact division residual is sufficient; no determinant retry.
    check({{8.25, 9, 9.25, 10, 8.25, 10, 10.25, 8}}, {{true, true}}, false, false);
    // Direction error alone blocks x, or the remaining division residual does.
    check({{-0.25, 9, 0.75, 10, -0.25, 10, 1.75, 8}}, {{false, true}}, false, false);
    check({{0.25, 2, 3.25, 2, 1.25, 1, 1.25, 3}}, {{false, true}}, false, false);
    // The shared endpoint subtraction is not exact; keep both old budgets.
    check({{4, 8, 5, 8, 4.5, -0x1p-10, 4.5, 0x1p54}}, {{false, true}}, false, false);
    // No certified positive lower bound: keep the original exact fallback.
    double const n = 0x1p27;
    fe::point_approximation result;
    auto const rejected = cr::detail::bounded_intersection<1, false, false>
        ({{0, 0, n + 1, n, 0, 1, n - 2, n - 2}}, result);
    BOOST_TEST(!rejected[0] && !rejected[1]);
}

void denominator_rescue()
{
    namespace fe = ec::forward_error;
    auto check = [](auto steps, producers const& p)
    {
        constexpr auto n = decltype(steps)::value;
        fe::point_approximation original, rescued, composed;
        // Cross terms cannot replace denominator rescue, or change its
        // rebuilt candidate and bound when both optimizations are enabled.
        auto const before = cr::detail::bounded_intersection<n, false, false, true>(p, original);
        BOOST_TEST(!before[0] && !before[1]);
        auto const after = cr::detail::bounded_intersection<n, false, true, false>(p, rescued);
        BOOST_TEST(after[0] && after[1]);
        auto const both = cr::detail::bounded_intersection<n, false, true, true>(p, composed);
        BOOST_TEST(both == after);
        exact_point reference;
        BOOST_TEST(oracle(p, reference));
        for (std::size_t d = 0; d < 2; ++d)
        {
            if (after[d] && both[d])
            {
                BOOST_TEST_EQ(bits(composed.value[d]), bits(rescued.value[d]));
                BOOST_TEST_EQ(bits(composed.error[d]), bits(rescued.error[d]));
            }
            BOOST_TEST(ulp_distance(rescued.value[d], nearest<double>(reference[d])) <= n);
            BOOST_TEST(absolute(rational(rescued.value[d]) - reference[d]) <= rational(rescued.error[d]));
        }
    };
    for (int exponent : {26, 30, 49})
    {
        double const n = std::ldexp(1.0, exponent);
        producers const p{{n, n, 2*n+1, 2*n, n+0.5, n+0.5, 2*n+0.5, 2*n-0.5}};
        auto const h = fe::detail::intersection_heads(p, false);
        BOOST_TEST_EQ(h.denominator == 0, exponent != 26);
        BOOST_TEST(std::abs(h.denominator) <= 8 * fe::detail::u * h.md);
        check(std::integral_constant<std::size_t, 1>{}, p);
        check(std::integral_constant<std::size_t, 4>{}, p);
        check(std::integral_constant<std::size_t, 100>{}, p);
        auto reversed = p;
        for (std::size_t d = 0; d < 2; ++d) std::swap(reversed[d], reversed[2+d]);
        check(std::integral_constant<std::size_t, 1>{}, reversed);
    }
    {
        // Exact D=-3 and T=-1: rebuild the quotient and its nonzero residual.
        double const n = 0x1p30;
        check(std::integral_constant<std::size_t, 1>{},
              {{n, n, 2*n+1, 2*n, n+1, n+1, 2*n-1, 2*n-2}});
    }
    for (producers const p : {producers{{1, 1, 3, 3, 1, 2, 3, 4}}, // exact parallel
                              producers{{1, 1, 0x1p55, 0x1p55, 2, 2, 0x1p55, 0x1p55}}})
    {
        fe::point_approximation sentinel{{{5, 6}}, {{1, 2}}};
        auto const accepted = cr::detail::bounded_intersection<100, true>(p, sentinel);
        BOOST_TEST(!accepted[0] && !accepted[1]);
        BOOST_TEST(sentinel.value[0] == 5 && sentinel.value[1] == 6);
    }
}

struct cross_term_counts
{
    std::size_t before = 0, after = 0, cache_before = 0, cache_after = 0;
};

template <std::size_t N>
std::array<bool, 2> check_cross_term_retry(producers const& input, exact_point const& exact,
                                           std::array<double, 2> const& nearest_coords, cross_term_counts& counts)
{
    using ec::forward_error::point_approximation;
    point_approximation plain, corrected{{{1, 2}}, {{3, 4}}};
    auto const before = cr::detail::bounded_intersection<N, false, true, false>(input, plain);
    auto const after = cr::detail::bounded_intersection<N, false, true, true>(input, corrected);
    counts.cache_before += !(before[0] && before[1]);
    counts.cache_after += !(after[0] && after[1]);
    for (std::size_t i = 0; i < 2; ++i)
    {
        counts.before += before[i];
        counts.after += after[i];
        if (before[i])
        {
            BOOST_TEST(after[i]);
            BOOST_TEST_EQ(bits(corrected.value[i]), bits(plain.value[i]));
            BOOST_TEST_EQ(bits(corrected.error[i]), bits(plain.error[i]));
        }
        if (after[i])
        {
            BOOST_TEST(ulp_distance(corrected.value[i], nearest_coords[i]) <= N);
            BOOST_TEST(absolute(rational(corrected.value[i]) - exact[i]) <= rational(corrected.error[i]));
        }
        else
        {
            BOOST_TEST_EQ(corrected.value[i], double(i + 1));
            BOOST_TEST_EQ(corrected.error[i], double(i + 3));
        }
    }
    return {{after[0] && !before[0], after[1] && !before[1]}};
}

void nonzero_cross_terms()
{
    namespace ar = ec::arithmetic_detail;
    namespace fe = ec::forward_error;
    using bg::detail::precise_math::two_diff_tail;

    auto check_budget = [](producers const& p, unsigned i)
    {
        double const x = ar::rounded_sub(p[i], p[0]), y = ar::rounded_sub(p[i + 1], p[1]);
        double const vx = ar::rounded_sub(p[6], p[4]), vy = ar::rounded_sub(p[7], p[5]);
        double const xt = two_diff_tail(p[i], p[0], x), yt = two_diff_tail(p[i + 1], p[1], y);
        double const vxt = two_diff_tail(p[6], p[4], vx), vyt = two_diff_tail(p[7], p[5], vy);
        auto const h = fe::detail::intersection_heads(p, false);
        double const k = i == 2 ? h.denominator : h.numerator, m = i == 2 ? h.md : h.mt;
        double const budget = cr::detail::cross_term_budget(x, y, vx, vy, k, xt, yt, vxt, vyt, m);
        rational const exact = (rational(p[i]) - rational(p[0])) * (rational(p[7]) - rational(p[5]))
                             - (rational(p[i + 1]) - rational(p[1])) * (rational(p[6]) - rational(p[4]));
        BOOST_TEST(absolute(exact - rational(k)) <= rational(fe::detail::u) * rational(budget));
        return budget;
    };
    // Both rounded products are 1. The linear terms cancel; the residual
    // -e^2 must be covered even though the head and signed correction are 0.
    double const e = 0x1p-54;
    producers const quadratic{{-e, 0, 1, 1, 0, e, 1, 1}};
    BOOST_TEST(check_budget(quadratic, 2) > 0);
    fe::point_approximation output;
    auto const rejected = cr::detail::bounded_intersection<1, false, true, true>(quadratic, output);
    BOOST_TEST(!rejected[0] && !rejected[1]); // No positive L; rescue rejects nonzero tails.

    std::array<cross_term_counts, 4> counts{};
    std::size_t checked = 0;
    auto check = [&](producers const& p)
    {
        check_budget(p, 2);
        check_budget(p, 4);
        exact_point exact;
        if (!oracle(p, exact)) return std::array<bool, 2>{{false, false}};
        std::array<double, 2> const nearest_coords{{nearest<double>(exact[0]), nearest<double>(exact[1])}};
        ++checked;
        auto const rescue = check_cross_term_retry<1>(p, exact, nearest_coords, counts[0]);
        check_cross_term_retry<2>(p, exact, nearest_coords, counts[1]);
        check_cross_term_retry<3>(p, exact, nearest_coords, counts[2]);
        check_cross_term_retry<4>(p, exact, nearest_coords, counts[3]);
        return rescue;
    };
    // Proper crossing. The shared vertical difference rounds away 2^-10;
    // exactly x was refused by the zero-tail retry.
    auto const rescue = check({{4, 8, 5, 8, 4.5, -0x1p-10, 4.5, 0x1p54}});
    BOOST_TEST(rescue[0] && !rescue[1]);
    for (auto const& p : bounded_corpus()) check(p);
    std::mt19937_64 random(0x54574f44494646);
    for (std::size_t i = 0; i < 2048; ++i)
    {
        producers p;
        for (double& x : p)
        {
            double const m = 1 + double(random() >> 12) * 0x1p-52;
            bool const negative = (random() & 1) != 0;
            int const exponent = int(random() % 128) - 64;
            x = std::ldexp(negative ? -m : m, exponent);
        }
        check(p);
    }
    for (std::size_t i = 0; i < counts.size(); ++i)
    {
        auto const& count = counts[i];
        std::cout << "TwoDiff N=" << i + 1 << ": " << count.before << " -> " << count.after
                  << '/' << 2 * checked << " axes; " << count.cache_before << " -> " << count.cache_after
                  << " constructions with a fallback\n";
    }
}

void integer_inputs()
{
    for (std::int64_t base : {std::int64_t(0), (std::int64_t(1) << 54) + 1,
         (std::numeric_limits<std::int64_t>::max)() - 128,
         (std::numeric_limits<std::int64_t>::min)() + 64})
    {
        check(std::array<std::int64_t, 8>{{base, base, base + 17, base + 31,
              base, base + 25, base + 17, base + 2}});
    }
    std::uint64_t const base = (std::numeric_limits<std::uint64_t>::max)() - 128;
    check(std::array<std::uint64_t, 8>{{base, base, base + 17, base + 31,
          base, base + 25, base + 17, base + 2}});
    check(std::array<float, 8>{{0.25f, -1.5f, 1.75f, 2.25f, -2.5f, 0.75f, 3.25f, -1.25f}});
}

void ordinary_first_stage()
{
    // Cancellation, binade neighbours, and near-parallel supports exercise the
    // inexpensive native bound independently of the lean EFT certificate.
    std::array<producers, 4> const inputs{{
        {{-1, -1, 1, 1, -1, 1, 1, -1}},
        {{previous(1.0), -1, next(1.0), 1, -2, 0, 2, 0}},
        {{-3.25, 0.5, 7.75, -1.5, -2.5, -4.25, 4.5, 6.75}},
        {{0, 0, 1, 1, 0, std::ldexp(1.0, -45), 1, 1 - std::ldexp(1.0, -45)}}
    }};
    for (auto input : inputs)
    for (int exponent : {-200, -53, -1, 0, 1, 53, 200})
    for (bool translated : {false, true})
    {
        auto v = input;
        for (double& coordinate : v)
            coordinate = std::ldexp(coordinate + (translated ? 1024 : 0), exponent);
        for (unsigned mask = 0; mask < 8; ++mask)
        {
            auto w = v;
            if (mask & 1) { std::swap(w[0], w[2]); std::swap(w[1], w[3]); }
            if (mask & 2) { std::swap(w[4], w[6]); std::swap(w[5], w[7]); }
            if (mask & 4) for (std::size_t i = 0; i < 4; ++i) std::swap(w[i], w[i + 4]);
            check(w);
        }
    }
    BOOST_TEST(ordinary_certificates > 1000);
}

template <typename C, typename = void> struct integer_fast_path : std::false_type {};
template <typename C> struct integer_fast_path<C, std::void_t<decltype(cr::detail::small_integer(
    std::declval<std::array<C, 8> const&>(), std::declval<std::array<double, 2>&>()))>> : std::true_type {};
static_assert(!integer_fast_path<float>::value && !integer_fast_path<double>::value,
              "No floating-value integer detection");

template <typename C, typename Policy>
void small_integer_state(std::array<C, 8> const& input, exact_point const& exact)
{
    ec::point_cache_scope<C> scope;
    cr::point_state<C, Policy> value;
    BOOST_TEST(ec::operation_detail::current<C>() == nullptr);
    without_allocations([&] { value.init(input); });
    BOOST_TEST(ec::operation_detail::current<C>() == nullptr);
    for (std::size_t d = 0; d < 2; ++d)
        BOOST_TEST_EQ(bits(value.approximate(d)), bits(nearest<double>(exact[d])));
    // A nearest-double approximation is not sufficient at float midpoints.
    // First float export must still use the target-format fixed rounding cache.
    for (std::size_t d = 0; d < 2; ++d)
        BOOST_TEST_EQ(bits(value.template emit<float>(d, input)), bits(nearest<float>(exact[d])));
    auto const* owner = ec::operation_detail::current<C>();
    BOOST_TEST(owner != nullptr);
    if (owner)
    {
        BOOST_TEST_EQ(owner->storage().native.size(), 1u);
        BOOST_TEST(owner->storage().native.front().exact_ready);
    }
}

template <typename C>
void small_integer_case(std::array<C, 8> const& input, bool in_range = true)
{
    std::array<double, 2> fast{{17, 19}};
    exact_point exact;
    bool const unique = oracle(input, exact);
    bool const accepted = cr::detail::small_integer(input, fast);
    BOOST_TEST_EQ(accepted, in_range && unique);
    if (accepted)
    {
        ++small_integer_certificates;
        auto const fixed = ec::direct_native_homogeneous(input);
        for (std::size_t d = 0; d < 2; ++d)
        {
            double expected;
            BOOST_TEST(ec::try_fixed_round(fixed, d, expected));
            BOOST_TEST_EQ(bits(fast[d]), bits(expected));
            BOOST_TEST_EQ(bits(fast[d]), bits(nearest<double>(exact[d])));
        }
        small_integer_state<C, policy::nearest_output>(input, exact);
        // Ordinary output must retain the original naive forward approximation,
        // even when the small-integral quotient could instead be rounded now.
        ec::forward_error::point_approximation initial;
        if (ec::forward_error::try_intersection(input, initial))
        {
            cr::point_state<C, policy::approximate_output> ordinary;
            ordinary.init(input);
            for (std::size_t d = 0; d < 2; ++d)
            {
                double const value = initial.value[d] == 0 && initial.error[d] == 0 ? 0 : initial.value[d];
                BOOST_TEST_EQ(bits(ordinary.approximate(d)), bits(value));
                BOOST_TEST_EQ(bits(ordinary.bounds(d).error), bits(initial.error[d]));
            }
        }
    }
    else { BOOST_TEST_EQ(fast[0], 17); BOOST_TEST_EQ(fast[1], 19); }
    check(input);
}

template <typename C>
void small_integer_boundaries(std::vector<C> const& limits)
{
    for (C boundary : limits)
    for (std::size_t i = 0; i < 8; ++i)
    {
        std::array<C, 8> input{{0, 0, 3, 3, 0, 3, 3, 0}};
        input[i] = boundary;
        bool const in_range = rational(boundary) > -65536 && rational(boundary) < 65536;
        small_integer_case(input, in_range);
    }
    std::mt19937_64 random(0x4231364341504143);
    for (unsigned i = 0; i < 512; ++i)
    {
        std::array<C, 8> input;
        for (C& value : input)
        {
            if constexpr (!std::is_signed<C>::value) value = C(random() % 65536);
            else if constexpr (std::numeric_limits<C>::digits <= 15) value = C(std::int64_t(random() % 65536) - 32768);
            else value = C(std::int64_t(random() % 131071) - 65535);
        }
        small_integer_case(input);
    }
}

void small_integers()
{
    small_integer_boundaries<std::int16_t>({-32768, -32767, -1, 0, 1, 32766, 32767});
    small_integer_boundaries<std::uint16_t>({0, 1, 32767, 32768, 65534, 65535});
    small_integer_boundaries<std::int64_t>({(std::numeric_limits<std::int64_t>::min)(),
        -65537, -65536, -65535, -1, 0, 1, 65535, 65536, 65537, (std::numeric_limits<std::int64_t>::max)()});
    small_integer_boundaries<std::uint64_t>({0, 1, 65535, 65536, 65537, (std::numeric_limits<std::uint64_t>::max)()});
    // Near-parallel supports with exact |W|=1 and large degree-three numerators.
    small_integer_case<std::int16_t>({{-32768, -32768, -16385, -16386, -32768, 0, -16386, 16381}});
    small_integer_case<std::uint16_t>({{0, 0, 32767, 32766, 0, 32768, 32766, 65533}});
    small_integer_case<std::int64_t>({{-65535, -65535, 0, -1, -65534, 1, 0, 65534}});
    // Exact coordinates just ABOVE a float midpoint, yet nearest-double lands
    // ON the midpoint. A cast of the small-integer approximation is incorrect.
    std::array<std::uint16_t, 8> const tie{{1, 0, 2, 32768, 32705, 1, 65473, 2}};
    exact_point exact;
    BOOST_TEST(oracle(tie, exact));
    BOOST_TEST(bits(static_cast<float>(nearest<double>(exact[0]))) != bits(nearest<float>(exact[0])));
    small_integer_case(tie);
    small_integer_case<std::int16_t>({{-32767, -1, -32766, 32767, -32767, 31, 1, 32}});
    small_integer_case<std::int16_t>({{-1, -1, 1, 1, -1, 1, 1, -1}}); // Canonical zero, either W sign.
    small_integer_case<std::int64_t>({{0, 0, 1, 1, 0, 1, 1, 2}}); // Parallel: no false certificate.
    BOOST_TEST(small_integer_certificates > 2000);
}

void cache_lifetime_and_paths()
{
    producers const simple{{1048576, 1048576, 1048592, 1048608,
                           1048576, 1048600, 1048592, 1048584}};
    {
        ec::point_cache_scope<double> scope;
        BOOST_TEST(ec::operation_detail::current<double>() == nullptr);
        ordinary a; nearest_state b;
        without_allocations([&] { a.init(simple); b.init(simple); });
        BOOST_TEST(ec::operation_detail::current<double>() == nullptr);
        BOOST_TEST_EQ(ec::point_cache_detail::native_cache<double>::size(), 0u);
        without_allocations([&] { BOOST_TEST_EQ(b.emit<double>(0, simple), 1048584); });
        BOOST_TEST(ec::operation_detail::current<double>() == nullptr);
    }
    {
        ec::point_cache_scope<double> scope;
        double const n = std::ldexp(1.0, 27);
        producers const v{{0, 0, n + 1, n, 0, 1, n, n}};
        nearest_state point;
        point.init(v);
        auto& cache = ec::operation_detail::require<double>().storage().native;
        BOOST_TEST_EQ(cache.size(), 1u);
        BOOST_TEST(cache.front().expanded[0].attempted);
        BOOST_TEST(!cache.front().exact_ready);
        auto const first = point.emit<double>(0, v);
        auto translated = [&](std::size_t i)
        {
            auto input = v;
            for (std::size_t j = 0; j < 8; j += 2) input[j] += 8 * double(i + 1);
            return input;
        };
        std::vector<nearest_state> relocated;
        for (unsigned i = 0; i < 40; ++i)
        { nearest_state other; other.init(translated(i)); relocated.push_back(other); }
        auto copy = point;
        BOOST_TEST_EQ(copy.refine(v, 0).approximate, first);
        // Emitting double now reads only the head. Refine distinct crossings
        // to exercise the relocated cache indices, not just the copied heads.
        for (std::size_t i = 0; i < relocated.size(); ++i)
        {
            auto const input = translated(i);
            exact_point exact;
            BOOST_TEST(oracle(input, exact));
            BOOST_TEST_EQ(relocated[i].refine(input, 0).approximate, nearest<double>(exact[0]));
        }
        point.init(simple); // Reinitialization must discard the old precision index.
        BOOST_TEST_EQ(point.emit<double>(0, simple), 1048584);
    }
    {
        ec::point_cache_scope<double> scope;
        double const s = std::ldexp(1.0, 200);
        producers const v{{0, 0, s, s, 0, s, s, 0}};
        nearest_state point;
        point.init(v);
        auto const& cache = ec::operation_detail::require<double>().storage().native;
        BOOST_TEST_EQ(cache.size(), 1u);
        BOOST_TEST(cache.front().exact_ready);
        BOOST_TEST_EQ(point.emit<double>(0, v), s / 2);
    }
    // There is no retained producer address, including across a stack lifetime.
    nearest_state point = [&] { auto local = simple; nearest_state p; p.init(local); return p; }();
    auto relocated = simple;
    without_allocations([&] { BOOST_TEST_EQ(point.emit<double>(0, relocated), 1048584); });
}

void warm_filters()
{
    ec::point_cache_scope<double> scope;
    ec::operation_detail::require<double>().storage().native.reserve(4);
    producers const v{{10, 10, 13, 17, 9, 15, 17, 11}};
    nearest_state a, b, c;
    auto v2 = v, v3 = v;
    for (double& x : v2) x += 16;
    for (double& x : v3) x -= 8;
    a.init(v); b.init(v2); c.init(v3);
    for (std::size_t d = 0; d < 2; ++d)
    {
        a.emit<float>(d, v); b.emit<float>(d, v2); c.emit<float>(d, v3);
    }
    using side = bg::strategy::side::side_robust<double, bg::strategy::side::fp_equals_policy, 3>;
    using native_point = bg::model::d2::point_xy<double>;
    native_point const p(0, 0), q(1, 0), r(0, 1);
    without_allocations([&]
    {
        for (unsigned i = 0; i < 256; ++i)
        {
            uf::coordinates const x{{a.approximation(0).approximate, a.approximation(1).approximate}};
            uf::coordinates const y{{b.approximation(0).approximate, b.approximation(1).approximate}};
            uf::coordinates const z{{c.approximation(0).approximate, c.approximation(1).approximate}};
            int result;
            uf::try_compare_n<0, false, true>(0, x[0], result);
            uf::try_compare_n<0, true, true>(x[0], y[0], result);
            uf::try_orientation<uf::nearest_certificate, false, false, true>({{0, 0}}, {{1, 0}}, x, result);
            uf::try_orientation<uf::nearest_certificate, false, true, true>({{0, 0}}, x, y, result);
            uf::try_orientation<uf::nearest_certificate, true, true, true>(x, y, z, result);
            BOOST_TEST_EQ(side::apply(p, q, r), 1);
            for (std::size_t d = 0; d < 2; ++d)
            { a.refine(v, d); a.emit<double>(d, v); a.emit<float>(d, v); a.emit<int>(d, v); }
        }
    });
}

void invalid_inputs()
{
    nearest_state point;
    bool rejected = false;
    try { point.init({{0, 0, 1, 1, 0, 1, 1, 2}}); }
    catch (std::domain_error const&) { rejected = true; }
    BOOST_TEST(rejected);
    rejected = false;
    try { point.approximation(2); }
    catch (std::out_of_range const&) { rejected = true; }
    BOOST_TEST(rejected);
}

template <bool Emit, typename C, typename Policy>
BOOST_NOINLINE std::uint64_t benchmark_sweep(std::vector<std::array<C, 8>> const& inputs,
    std::vector<cr::point_state<C, Policy>>& points, ec::operation_detail::storage_for<C>& storage)
{
    storage.native.clear();
    std::uint64_t checksum = 0;
    for (std::size_t i = 0; i < inputs.size(); ++i)
    {
        points[i].init(inputs[i]);
        if constexpr (Emit)
        {
            checksum += bits(points[i].template emit<double>(0, inputs[i]));
            checksum ^= bits(points[i].template emit<double>(1, inputs[i]));
        }
        else
        {
            checksum += bits(points[i].approximate(0));
            checksum ^= bits(points[i].approximate(1));
        }
    }
    return checksum;
}

volatile std::uint64_t benchmark_sink = 0;

template <typename C, typename Policy, bool Emit = false>
void benchmark_corpus(std::vector<std::array<C, 8>> const& inputs,
                      char const* type, char const* policy_name, char const* corpus)
{
    auto const count = inputs.size();
    ec::point_cache_scope<C> scope;
    auto& storage = ec::operation_detail::require<C>().storage();
    storage.native.reserve(count);
    std::vector<cr::point_state<C, Policy>> points(count);
    auto run = [&](std::size_t sweeps)
    {
        auto const begin = std::chrono::steady_clock::now();
        for (std::size_t i = 0; i < sweeps; ++i)
            benchmark_sink ^= benchmark_sweep<Emit>(inputs, points, storage);
        return std::chrono::duration<double>(std::chrono::steady_clock::now() - begin).count();
    };
    double const pilot = run(64);
    auto const sweeps = std::size_t((std::max)(8.0, (std::min)(262144.0, 64 * 0.03 / pilot)));
    std::array<double, 5> times;
    for (double& t : times) t = run(sweeps) * 1e9 / double(sweeps * count);
    std::sort(times.begin(), times.end());
    std::cout << type << ',' << policy_name << ',' << corpus << ','
              << std::setprecision(10) << times[2] << ',' << times.front() << ',' << times.back()
              << ',' << benchmark_sweep<Emit>(inputs, points, storage) << ',' << storage.native.size() << '\n';
}

template <typename C, typename Policy>
void benchmark_case(char const* type, char const* policy_name, bool outside = false)
{
    std::vector<std::array<C, 8>> inputs;
    std::mt19937_64 random(0x534d414c4c494e54);
    for (std::size_t i = 0; i < 128; ++i)
    {
        std::int64_t const radius = 1 + random() % 128;
        auto position = [&]
        {
            if (outside) return std::int64_t(1048576 + random() % 2001);
            if constexpr (!std::is_signed<C>::value) return std::int64_t(256 + random() % 64000);
            else return std::int64_t(random() % 60001) - 30000;
        };
        auto const x = position(), y = position();
        inputs.push_back({{C(x - radius), C(y - radius), C(x + radius), C(y + radius + 1),
            C(x - radius), C(y + radius), C(x + radius), C(y - radius - 2)}});
    }
    benchmark_corpus<C, Policy>(inputs, type, policy_name, outside ? "outside" : "small");
}

void benchmark_bounded()
{
    auto const inputs = bounded_corpus();
    auto run = [&](auto tag, char const* name)
    { benchmark_corpus<double, decltype(tag), true>(inputs, "double", name, "decimal"); };
    std::cout << "coordinate,policy,corpus,median_ns,min_ns,max_ns,checksum,cache_entries\n";
    run(policy::approximate_output{}, "ordinary");
    run(policy::nearest_output{}, "nearest");
    run(policy::bounded_ulp_output<1>{}, "bounded1");
    run(policy::bounded_ulp_output<1, true>{}, "bounded1-retry");
    run(policy::bounded_ulp_output<4>{}, "bounded4");
    run(policy::bounded_ulp_output<4, true>{}, "bounded4-retry");
    run(policy::bounded_ulp_output<100>{}, "bounded100");
    run(policy::bounded_ulp_output<100, true>{}, "bounded100-retry");
}

void benchmark()
{
    std::cout << "coordinate,policy,corpus,median_ns,min_ns,max_ns,checksum,cache_entries\n";
    benchmark_case<std::int16_t, policy::nearest_output>("int16", "nearest");
    benchmark_case<std::uint16_t, policy::nearest_output>("uint16", "nearest");
    benchmark_case<std::int64_t, policy::nearest_output>("int64", "nearest");
    benchmark_case<std::int64_t, policy::nearest_output>("int64", "nearest", true);
    benchmark_case<std::int64_t, policy::approximate_output>("int64", "ordinary");
}

int main(int argc, char** argv)
{
    try
    {
        BOOST_TEST_EQ(std::fegetround(), FE_TONEAREST);
        if (argc == 2 && std::string(argv[1]) == "--benchmark") { benchmark(); return boost::report_errors(); }
        if (argc == 2 && std::string(argv[1]) == "--benchmark-bounded")
        { benchmark_bounded(); return boost::report_errors(); }
        certificates();
        boundaries();
        random_cases();
        bounded_cases();
        bounded_other_formats();
        determinant_budgets();
        selective_one_ulp_retry();
        denominator_rescue();
        nonzero_cross_terms();
        integer_inputs();
        small_integers();
        ordinary_first_stage();
        cache_lifetime_and_paths();
        warm_filters();
        invalid_inputs();
        std::cout << checked << " policy constructions; " << forward_certificates
                  << " forward nearest certificates.\n";
        std::cout << small_integer_certificates << " exact small-integer certificates.\n";
        std::cout << ordinary_certificates << " allocation-free native ordinary certificates.\n";
    }
    catch (std::exception const& e)
    {
        counting = false;
        std::cerr << e.what() << '\n';
        BOOST_ERROR("Construction-rounding test threw");
    }
    return boost::report_errors();
}
