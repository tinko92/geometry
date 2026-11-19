// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.

#ifndef BOOST_GEOMETRY_TEST_EXACT_ROUNDING_POLICY_HELPERS_HPP
#define BOOST_GEOMETRY_TEST_EXACT_ROUNDING_POLICY_HELPERS_HPP

// Standalone numerical experiment, not a production point representation.
// Contract: binary64, round-to-nearest, no reassociation/contraction, no FP
// overflow/underflow. Compile with -ffp-contract=off and without fast-math.
#include <boost/geometry/strategies/cartesian/exact/construction_filter.hpp>
#include <boost/geometry/strategies/cartesian/exact/fixed_rounding.hpp>
#include <boost/geometry/strategies/cartesian/exact/ulp_filter.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <optional>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace rounding_policies
{
namespace ec = boost::geometry::strategies::detail::exact_cartesian;
namespace cd = ec::construction_detail;
namespace uf = ec::ulp_filter;
using endpoints = std::array<double, 8>;
using coordinates = uf::coordinates;
using rational = ec::exact_number;
using fixed_point = ec::native_fixed_traits<double>::intersection_type;

enum class policy { eager, deferred };
enum class predicate { coordinate, nni, nii, iii, nni_support, nii_source };

inline char const* name(policy p)
{ return p == policy::eager ? "eager_nearest" : "deferred_nearest"; }

inline char const* name(predicate p)
{
    switch (p)
    {
    case predicate::coordinate: return "coordinate";
    case predicate::nni: return "NNI";
    case predicate::nii: return "NII";
    case predicate::iii: return "III";
    case predicate::nni_support: return "NNI_on_support";
    case predicate::nii_source: return "NII_source_anchor";
    }
    throw std::invalid_argument("unknown predicate");
}

inline std::uint64_t bits(double value)
{
    std::uint64_t result;
    std::memcpy(&result, &value, sizeof(result));
    return result;
}

inline bool even(double value) { return (bits(value) & 1) == 0; }
inline double previous(double x)
{ return std::nextafter(x, -std::numeric_limits<double>::infinity()); }
inline double next(double x)
{ return std::nextafter(x, std::numeric_limits<double>::infinity()); }

struct statistics
{
    std::uint64_t constructions = 0, zero_denominators = 0;
    std::uint64_t expansion_entries = 0, fixed_entries = 0;
    std::uint64_t ulp_certifications = 0, nearest_certifications = 0;
    std::uint64_t boundary_tests = 0, candidate_steps = 0;
    std::uint64_t estimate_refinements = 0, fixed_rounds = 0;
    std::uint64_t initial_nearest = 0, emission_nearest = 0;
    std::uint64_t predicate_filters = 0, predicate_exact = 0, emissions = 0;
    std::array<std::uint64_t, 3> orientation_filters{}, orientation_exact{};
};

struct point
{
    std::array<double, 2> approx;
    endpoints const* source;
    std::size_t cache_index;
    unsigned nearest_mask = 0;
};

struct precision_entry
{
    precision_entry() {}
    cd::expansion<32> x, y;
    cd::expansion<16> w;
    std::optional<fixed_point> fixed;
    bool expanded = false;
};

// Bounds used below, including zero/binade boundaries:
//
//   prev(a) <= N/W <= next(a)                           (faithful / one ULP)
//   (prev(a)+a)/2 <= N/W <= (a+next(a))/2               (nearest)
//
// All comparisons use exact expansion signs of 2N-(b+c)W, multiplied
// by sign(W). The nearest cell includes a midpoint iff a's significand is
// even. Thus no assumed error estimate or rounded midpoint enters either
// certificate. The first bound implies |a-N/W| <= the adjacent spacing in
// the direction of N/W, also at a binade boundary. At zero we require N=0
// instead of forming subnormal neighbours (underflow is outside the contract).
// No error field is stored. Failed bounded expansion certification goes to
// checked, coordinate-format-derived fixed arithmetic, never a guessed sign.
//
// The two policies deliberately share this conservative certifier. This
// measures rounding scheduling with existing arithmetic, not the performance
// of future specialized forward-error construction filters.
template <policy Policy, bool Instrument = false>
class operation
{
    std::vector<precision_entry> m_cache;

    void count(std::uint64_t statistics::* member)
    {
        if constexpr (Instrument) { ++(counts.*member); }
    }

    static bool supported(endpoints const& v)
    {
        // Existing expansion capacity/range proof; no release-mode fenv check.
        for (double x : v)
            if (!std::isfinite(x) || !cd::supported_value(x, 150)) return false;
        return true;
    }

    precision_entry& precision(point const& p) { return m_cache[p.cache_index]; }

    fixed_point const& fixed(point const& p)
    {
        auto& e = precision(p);
        if (!e.fixed)
        {
            e.fixed.emplace(ec::direct_native_homogeneous(*p.source));
            if (e.fixed->W == 0) throw std::domain_error("parallel supporting lines");
            count(&statistics::fixed_entries);
        }
        return *e.fixed;
    }

    void expand(point const& p)
    {
        auto& e = precision(p);
        auto const& v = *p.source;
        if (!supported(v)) return;
        // Same degree-three X/Y and degree-two W formulas as
        // homogeneous_coordinate(), sharing the line equations for x and y.
        auto const ax = cd::difference(v[2], v[0]);
        auto const ay = cd::difference(v[1], v[3]);
        auto const bx = cd::difference(v[6], v[4]);
        auto const by = cd::difference(v[5], v[7]);
        cd::expansion<4> a, b;
        e.expanded = a.product(v[0], v[3]) && a.product(v[1], v[2], true)
            && b.product(v[4], v[7]) && b.product(v[5], v[6], true)
            && cd::add_product(e.w, ay, bx) && cd::add_product(e.w, ax, by, true)
            && cd::add_product(e.x, ax, b) && cd::add_product(e.x, a, bx, true)
            && cd::add_product(e.y, a, by) && cd::add_product(e.y, ay, b, true);
        if (e.expanded)
        {
            e.x.compress(); e.y.compress(); e.w.compress();
            if (e.w.sign() == 0) throw std::domain_error("parallel supporting lines");
            count(&statistics::expansion_entries);
        }
    }

    bool boundary(point const& p, unsigned d, double a, double b, int& sign)
    {
        count(&statistics::boundary_tests);
        auto const& e = precision(p);
        return cd::compare_boundary(d == 0 ? e.x : e.y, e.w, a, b, sign);
    }

    bool certify(point const& p, unsigned d, bool nearest, double& candidate)
    {
        auto const& e = precision(p);
        if (!e.expanded) return false;
        auto const& n = d == 0 ? e.x : e.y;
        if (n.sign() == 0)
        {
            candidate = 0;
            count(nearest ? &statistics::nearest_certifications : &statistics::ulp_certifications);
            return true;
        }
        if (!cd::rounding_cell<double>::allowed(candidate)) return false;
        for (unsigned correction = 0; correction < 4; ++correction)
        {
            int lo, hi;
            bool const ok = nearest
                ? boundary(p, d, previous(candidate), candidate, lo)
                    && boundary(p, d, candidate, next(candidate), hi)
                : boundary(p, d, previous(candidate), previous(candidate), lo)
                    && boundary(p, d, next(candidate), next(candidate), hi);
            if (!ok) return false;
            bool const include = !nearest || even(candidate);
            bool const above = lo > 0 || (lo == 0 && include);
            bool const below = hi < 0 || (hi == 0 && include);
            if (above && below)
            {
                count(nearest ? &statistics::nearest_certifications : &statistics::ulp_certifications);
                return true;
            }
            candidate = above ? next(candidate) : previous(candidate);
            count(&statistics::candidate_steps);
        }
        return false;
    }

    void refine(point& p, unsigned d, bool nearest)
    {
        double candidate = p.approx[d];
        if (!certify(p, d, nearest, candidate))
        {
            auto const& e = precision(p);
            bool ready = false;
            if (e.expanded)
            {
                count(&statistics::estimate_refinements);
                candidate = (d == 0 ? e.x : e.y).estimate() / e.w.estimate();
                ready = certify(p, d, nearest, candidate);
            }
            if (!ready)
            {
                if (!ec::try_fixed_round(fixed(p), d, candidate) || !std::isfinite(candidate))
                    throw std::domain_error("unsupported intersection output");
                count(&statistics::fixed_rounds);
                // Fixed fallback already returns nearest, also satisfying one ULP.
                p.nearest_mask |= 1u << d;
            }
        }
        p.approx[d] = candidate;
        if (nearest) p.nearest_mask |= 1u << d;
    }

public:
    statistics counts;

    explicit operation(std::size_t capacity) { m_cache.reserve(capacity); }
    void reset() { m_cache.clear(); }

    point construct(endpoints const& v)
    {
        count(&statistics::constructions);
        point p{{{0, 0}}, &v, m_cache.size()};
        m_cache.emplace_back();
        // Ordinary segment-parameter construction, not precomputed rational
        // coordinates. A zero rounded denominator is refined now, never left
        // as an unusable approximate point until the first predicate.
        double const ux = v[2] - v[0], uy = v[3] - v[1];
        double const vx = v[6] - v[4], vy = v[7] - v[5];
        double const wx = v[4] - v[0], wy = v[5] - v[1];
        double const denominator = ux * vy - uy * vx;
        if (denominator == 0)
        {
            count(&statistics::zero_denominators);
            p.approx.fill(std::numeric_limits<double>::quiet_NaN());
        }
        else
        {
            double const t = (wx * vy - wy * vx) / denominator;
            p.approx = {{std::fma(t, ux, v[0]), std::fma(t, uy, v[1])}};
        }
        expand(p);
        for (unsigned d = 0; d < 2; ++d)
        {
            if constexpr (Policy == policy::eager)
            {
                count(&statistics::initial_nearest);
                refine(p, d, true);
            }
            else { refine(p, d, false); }
        }
        return p;
    }

    int compare(point const& p, unsigned d, double target)
    {
        int sign;
        // Canonical NI compares native-implicit; this method returns the
        // opposite order. The exact expansion/fixed fallback is unchanged.
        if (uf::try_compare_ni(target, p.approx[d], sign))
        {
            count(&statistics::predicate_filters);
            return -sign;
        }
        count(&statistics::predicate_exact);
        if (precision(p).expanded && cd::supported_value(target, 150)
            && boundary(p, d, target, target, sign)) return sign;
        auto const native = ec::fixed_input_homogeneous(d == 0 ? target : 0.0, d == 1 ? target : 0.0);
        return ec::fixed_coordinate_compare(fixed(p), native, d);
    }

    int orientation_nni(coordinates const& a, coordinates const& b, point const& c)
    {
        int sign;
        if (uf::try_orientation_nni(a, b, c.approx, sign))
        {
            if constexpr (Instrument) ++counts.orientation_filters[0];
            return sign;
        }
        if constexpr (Instrument) ++counts.orientation_exact[0];
        return ec::fixed_orientation(ec::fixed_input_homogeneous(a[0], a[1]),
            ec::fixed_input_homogeneous(b[0], b[1]), fixed(c));
    }

    int orientation_nii(coordinates const& a, point const& b, point const& c)
    {
        int sign;
        if (uf::try_orientation_nii(a, b.approx, c.approx, sign))
        {
            if constexpr (Instrument) ++counts.orientation_filters[1];
            return sign;
        }
        if constexpr (Instrument) ++counts.orientation_exact[1];
        return ec::fixed_orientation(ec::fixed_input_homogeneous(a[0], a[1]), fixed(b), fixed(c));
    }

    int orientation_iii(point const& a, point const& b, point const& c)
    {
        int sign;
        if (uf::try_orientation_iii(a.approx, b.approx, c.approx, sign))
        {
            if constexpr (Instrument) ++counts.orientation_filters[2];
            return sign;
        }
        if constexpr (Instrument) ++counts.orientation_exact[2];
        // Existing fixed_orientation proves its six homogeneous monomial
        // capacity from the operand types. Here every implicit operand is a
        // native-segment intersection (degree three numerator, degree two
        // denominator), not an arbitrarily recursive construction.
        return ec::fixed_orientation(fixed(a), fixed(b), fixed(c));
    }

    std::array<double, 2> emit(point& p)
    {
        count(&statistics::emissions);
        if constexpr (Policy == policy::deferred)
        {
            for (unsigned d = 0; d < 2; ++d)
                if ((p.nearest_mask & (1u << d)) == 0)
                {
                    count(&statistics::emission_nearest);
                    refine(p, d, true);
                }
        }
        return p.approx;
    }

    // Approximate-policy output must do exactly this, with no refinement.
    static std::array<double, 2> emit_approximate(point const& p) { return p.approx; }
};

// Independent oracle: existing cpp_rational exact input conversion and an
// exact segment-parameter formula, not the tested homogeneous/expansion code.
// The conversion only seeds a candidate; exact rational midpoints and parity
// determine nearest. This code is used before/after timing, never within it.
inline std::array<rational, 2> reference(endpoints const& v)
{
    std::array<rational, 8> q;
    for (unsigned i = 0; i < 8; ++i) q[i] = ec::exact_value(v[i]);
    rational const ux = q[2] - q[0], uy = q[3] - q[1];
    rational const vx = q[6] - q[4], vy = q[7] - q[5];
    rational const denominator = ux * vy - uy * vx;
    if (denominator == 0) throw std::domain_error("parallel oracle input");
    rational const t = ((q[4] - q[0]) * vy - (q[5] - q[1]) * vx) / denominator;
    return {{q[0] + t * ux, q[1] + t * uy}};
}

inline double reference_nearest(rational const& q)
{
    if (q == 0) return 0;
    double a = q.convert_to<double>();
    for (unsigned i = 0; i < 8; ++i)
    {
        if (!std::isfinite(a) || a == 0 || !std::isnormal(a))
            throw std::domain_error("oracle corpus requires normal finite outputs or exact zero");
        rational const lo = (ec::exact_value(previous(a)) + ec::exact_value(a)) / 2;
        rational const hi = (ec::exact_value(a) + ec::exact_value(next(a))) / 2;
        if ((q > lo || (q == lo && even(a))) && (q < hi || (q == hi && even(a)))) return a;
        a = q <= lo ? previous(a) : next(a);
    }
    throw std::runtime_error("oracle conversion did not converge");
}

inline bool one_ulp(double a, rational const& q)
{
    if (!std::isfinite(a)) return false;
    if (a == 0) return q == 0;
    return ec::exact_value(previous(a)) <= q && q <= ec::exact_value(next(a));
}

struct query { unsigned dimension; double target; int expected; };
struct orientation_query
{
    coordinates first, second;
    std::size_t second_index, third_index;
    // NNI, NII, III, NNI on a producing line, NII with a source endpoint.
    std::array<int, 5> expected;
};
struct sample
{
    endpoints input;
    std::array<rational, 2> exact;
    std::array<double, 2> rounded;
    std::array<query, 6> queries;
    std::array<orientation_query, 3> orientation_queries{};
};

inline std::array<rational, 2> native_reference(coordinates const& p)
{ return {{ec::exact_value(p[0]), ec::exact_value(p[1])}}; }

inline int reference_orientation(std::array<rational, 2> const& a,
                                 std::array<rational, 2> const& b,
                                 std::array<rational, 2> const& c)
{
    rational const value = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
    return int(value > 0) - int(value < 0);
}

inline sample prepare(endpoints const& v)
{
    sample s;
    s.input = v;
    s.exact = reference(v);
    for (unsigned d = 0; d < 2; ++d)
    {
        s.rounded[d] = reference_nearest(s.exact[d]);
        for (unsigned j = 0; j < 3; ++j)
        {
            double target = s.rounded[d];
            if (target == 0 && j != 1) target = j == 0 ? -1 : 1;
            else for (unsigned k = 0; k < 8 && j != 1; ++k)
                target = j == 0 ? previous(target) : next(target);
            rational const exact_target = ec::exact_value(target);
            s.queries[3 * d + j] = {d, target,
                int(s.exact[d] > exact_target) - int(s.exact[d] < exact_target)};
        }
    }
    return s;
}

// Call after the whole corpus is assembled. Native query points are exact
// binary64 fixture values, irrespective of how they were chosen. Independent
// rational orientation determines every expectation outside the timed region.
inline void prepare_orientation_queries(std::vector<sample>& samples)
{
    for (std::size_t i = 0; i < samples.size(); ++i)
    {
        auto& s = samples[i];
        coordinates const source_a{{s.input[0], s.input[1]}}, source_b{{s.input[2], s.input[3]}};
        double const span = std::ldexp((std::max)(1.0,
            (std::max)(std::abs(s.rounded[0]), std::abs(s.rounded[1]))), -8);
        for (unsigned j = 0; j < s.orientation_queries.size(); ++j)
        {
            auto& q = s.orientation_queries[j];
            std::size_t const offset_b = samples.size() > 1 ? 1 + j % (samples.size() - 1) : 0;
            std::size_t const offset_c = samples.size() > 1 ? 1 + (j + 1) % (samples.size() - 1) : 0;
            q.second_index = (i + offset_b) % samples.size();
            q.third_index = (i + offset_c) % samples.size();
            // General-position native query lines, not the producing lines.
            // Include both sides; the rational oracle makes no assumption
            // about their position relative to the true implicit coordinate.
            if (j < 2)
            {
                double const y = s.rounded[1] + (j == 0 ? -span : span);
                q.first = {{s.rounded[0] - span, y}};
                q.second = {{s.rounded[0] + span, y}};
            }
            else
            {
                q.first = {{s.rounded[0] + span, s.rounded[1] - span}};
                q.second = {{s.rounded[0] + span, s.rounded[1] + span}};
            }
            auto const a = native_reference(q.first), b = native_reference(q.second);
            auto const& other = samples[q.second_index].exact;
            q.expected = {{reference_orientation(a, b, s.exact),
                reference_orientation(a, s.exact, other),
                reference_orientation(s.exact, other, samples[q.third_index].exact),
                reference_orientation(native_reference(source_a), native_reference(source_b), s.exact),
                reference_orientation(native_reference({{s.input[4], s.input[5]}}), s.exact, other)}};
        }
    }
}

inline std::vector<sample> corpus(std::string const& shape, std::size_t count, std::uint64_t seed)
{
    std::mt19937_64 random(seed);
    std::vector<sample> result;
    result.reserve(count);
    auto dyadic = [&]() { return double(std::int64_t(random() % 2000001) - 1000000) / 1024; };
    for (std::size_t i = 0; i < count; ++i)
    {
        endpoints v;
        if (shape == "ordinary")
        {
            // Mix short dyadics with general binary64 significands, rather
            // than timing only constructions of small integer coordinates.
            auto fraction = [&]() { return std::ldexp(double(random() >> 11), -53); };
            double const x = dyadic() + (i & 2 ? fraction() : 0);
            double const y = dyadic() + (i & 2 ? fraction() : 0);
            double const r = 1 + double(random() % 100000) / 1024 + (i & 2 ? fraction() : 0);
            v = {{x - r, y - r, x + r, y + r + 0.25,
                  x - r, y + r, x + r, y - r - 0.5}};
        }
        else if (shape == "shared_boundary")
        {
            double const x = dyadic(), y = dyadic(), r = 1 + double(random() % 1000);
            v = {{x, y, x + r, y + 2 * r, x, y, x - 2 * r, y + r}};
        }
        else if (shape == "near_parallel")
        {
            if (i % 2 == 0)
            {
                // Exact determinant -1; rounded products coincide (D=0).
                // Supporting lines intersect outside their finite segments.
                double const n = std::ldexp(1.0, 27) + 2 * double(random() % 512);
                v = {{0, 0, n + 1, n, 0, 1, n, n}};
            }
            else
            {
                // A real segment crossing, with a small nonzero determinant.
                double const e = std::ldexp(1.0, -48 - int(random() % 4));
                v = {{0, 0, 1, 1, 0, e, 1, 1 - e}};
            }
        }
        else throw std::invalid_argument("unknown corpus: " + shape);
        if (i & 1) { std::swap(v[0], v[2]); std::swap(v[1], v[3]); }
        result.push_back(prepare(v));
    }
    prepare_orientation_queries(result);
    return result;
}

template <policy Policy, bool Instrument>
inline int evaluate_query(operation<Policy, Instrument>& op, std::vector<point> const& points,
                           std::vector<sample> const& samples, std::size_t i, unsigned j, predicate kind)
{
    auto const& s = samples[i];
    if (kind == predicate::coordinate)
    {
        auto const& q = s.queries[j % s.queries.size()];
        return op.compare(points[i], q.dimension, q.target);
    }
    auto const& q = s.orientation_queries[j % s.orientation_queries.size()];
    switch (kind)
    {
    case predicate::nni: return op.orientation_nni(q.first, q.second, points[i]);
    case predicate::nii: return op.orientation_nii(q.first, points[i], points[q.second_index]);
    case predicate::iii: return op.orientation_iii(points[i], points[q.second_index], points[q.third_index]);
    case predicate::nni_support:
        return op.orientation_nni({{s.input[0], s.input[1]}}, {{s.input[2], s.input[3]}}, points[i]);
    case predicate::nii_source:
        return op.orientation_nii({{s.input[4], s.input[5]}}, points[i], points[q.second_index]);
    default: throw std::invalid_argument("unknown orientation predicate");
    }
}

inline int query_reference(sample const& s, unsigned j, predicate kind)
{
    if (kind == predicate::coordinate) return s.queries[j % s.queries.size()].expected;
    return s.orientation_queries[j % s.orientation_queries.size()].expected[unsigned(kind) - 1];
}

template <policy Policy>
inline statistics verify(std::vector<sample> const& samples, unsigned uses, unsigned emit_every,
                          predicate kind = predicate::coordinate)
{
    operation<Policy, true> op(samples.size());
    std::vector<point> points;
    points.reserve(samples.size());
    for (auto const& s : samples)
    {
        points.push_back(op.construct(s.input));
        for (unsigned d = 0; d < 2; ++d)
        {
            if (!one_ulp(points.back().approx[d], s.exact[d]))
                throw std::runtime_error("initial one-ULP certificate disagrees with oracle");
            if constexpr (Policy == policy::eager)
                if (bits(points.back().approx[d]) != bits(s.rounded[d]))
                    throw std::runtime_error("eager nearest disagrees with oracle");
        }
    }
    for (std::size_t i = 0; i < samples.size(); ++i)
    {
        auto const& s = samples[i];
        auto& p = points[i];
        for (unsigned j = 0; j < uses; ++j)
        {
            if (evaluate_query(op, points, samples, i, j, kind) != query_reference(s, j, kind))
                throw std::runtime_error(std::string("ULP-filtered ") + name(kind) + " disagrees with oracle");
        }
        if (emit_every != 0 && i % emit_every == 0)
        {
            auto const emitted = op.emit(p);
            for (unsigned d = 0; d < 2; ++d)
                if (bits(emitted[d]) != bits(s.rounded[d]))
                    throw std::runtime_error("final nearest disagrees with oracle");
            auto const before = op.counts.boundary_tests + op.counts.fixed_rounds;
            op.emit(p);
            if (before != op.counts.boundary_tests + op.counts.fixed_rounds)
                throw std::runtime_error("repeated output recomputed nearest");
        }
    }
    return op.counts;
}

} // namespace rounding_policies

#endif
