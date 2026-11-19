// Boost.Geometry
// Distributed under the Boost Software License, Version 1.0.

#ifndef BOOST_GEOMETRY_STRATEGIES_BUFFER_CARTESIAN_MP_HPP
#define BOOST_GEOMETRY_STRATEGIES_BUFFER_CARTESIAN_MP_HPP

#include <array>
#include <cmath>
#include <type_traits>
#include <boost/geometry/algorithms/buffer.hpp>
#include <boost/geometry/strategies/buffer/cartesian.hpp>
#include <boost/geometry/strategies/cartesian/buffer_join_round.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_int.hpp>

namespace boost { namespace geometry { namespace strategies { namespace buffer
{ namespace cartesian_mp_detail
{

template <typename Point>
using admitted_type = std::integral_constant<bool,
    std::is_same<coordinate_type_t<Point>, multiprecision::cpp_bin_float_100>::value
    && std::is_same<cs_tag_t<Point>, cartesian_tag>::value>;

using rational = multiprecision::cpp_rational;
struct ray { rational x, y; };

template <typename Point>
bool admissible(Point const& point)
{
    using C = coordinate_type_t<Point>;
    // Keep the dyadic-to-rational conversion bounded. No geometric
    // comparison within this band uses a tolerance.
    using std::ldexp;
    static C const high = ldexp(C(1), 512), low = ldexp(C(1), -512);
    auto in_band = [](C const& v)
    { return v == 0 || (v >= low && v <= high) || (v <= -low && v >= -high); };
    return in_band(geometry::get<0>(point)) && in_band(geometry::get<1>(point));
}

template <typename Point>
bool same(Point const& a, Point const& b)
{
    return geometry::get<0>(a) == geometry::get<0>(b)
        && geometry::get<1>(a) == geometry::get<1>(b);
}

template <typename Point>
ray direction(Point const& from, Point const& to)
{
    return {rational(geometry::get<0>(to)) - rational(geometry::get<0>(from)),
            rational(geometry::get<1>(to)) - rational(geometry::get<1>(from))};
}

inline bool nonzero(ray const& r) { return r.x != 0 || r.y != 0; }

inline int cross(ray const& a, ray const& b)
{
    rational const d = a.x * b.y - a.y * b.x;
    return (d > 0) - (d < 0);
}

inline int half(ray const& r)
{ return r.y < 0 || (r.y == 0 && r.x < 0) ? 1 : 0; }

inline int compare(ray const& a, ray const& b)
{ return half(a) == half(b) ? -cross(a, b) : half(a) - half(b); }

inline int octant(ray const& r)
{
    if (r.y == 0) return r.x > 0 ? 0 : r.x < 0 ? 4 : -1;
    if (r.x == 0) return r.y > 0 ? 2 : 6;
    if (r.x == r.y) return r.x > 0 ? 1 : 5;
    if (r.x == -r.y) return r.x > 0 ? 7 : 3;
    return -1;
}

// For distinct rays: strict membership in the CCW interval (first,last).
inline bool inside(ray const& first, ray const& last, ray const& r)
{
    bool const after = compare(first, r) < 0, before = compare(r, last) < 0;
    return compare(first, last) < 0 ? after && before : after || before;
}

} // namespace cartesian_mp_detail

// Opt-in repairs for represented binary multiprecision Cartesian buffers.
// All other strategy operations are those of the conventional umbrella.
struct cartesian_mp : cartesian<>
{
    // Return the outgoing side of a proved switch, or -1 if inconclusive.
    template <typename Range, typename Point>
    static typename std::enable_if<cartesian_mp_detail::admitted_type<Point>::value, int>::type
    exposed_vertex_exit(Range const& p, Range const& q, Point const& vertex)
    {
        namespace local = cartesian_mp_detail;
        if (!local::same(p.at(1), q.at(1)) || !local::same(p.at(1), vertex))
            return -1;

        auto const& origin = p.at(1);
        if (!local::admissible(origin) || !local::admissible(p.at(0)) || !local::admissible(p.at(2))
            || !local::admissible(q.at(0)) || !local::admissible(q.at(2)))
            return -1;
        std::array<local::ray, 4> const rays{{
            local::direction(origin, p.at(0)), local::direction(origin, p.at(2)),
            local::direction(origin, q.at(0)), local::direction(origin, q.at(2))}};
        for (std::size_t a = 0; a != rays.size(); ++a)
        {
            if (!local::nonzero(rays[a])) return -1;
            for (std::size_t b = 0; b != a; ++b)
                if (local::compare(rays[a], rays[b]) == 0) return -1;
        }

        // A clockwise ring has its interior CCW from the backward
        // incoming ray to the outgoing ray. Distinct rays make coverage
        // of each ray equivalent to coverage of its adjacent sector.
        bool const pin = !local::inside(rays[2], rays[3], rays[0]);
        bool const pout = !local::inside(rays[2], rays[3], rays[1]);
        bool const qin = !local::inside(rays[0], rays[1], rays[2]);
        bool const qout = !local::inside(rays[0], rays[1], rays[3]);
        if (pin == qin || pout == qout || pin == pout) return -1;
        return pout ? 0 : 1;
    }

    template <typename Join, typename Segment, typename Distance, typename Point,
              typename OutputPoint, typename Range>
    static typename std::enable_if<
        cartesian_mp_detail::admitted_type<Point>::value
        && cartesian_mp_detail::admitted_type<OutputPoint>::value
        && std::is_same<Join, strategy::buffer::join_round>::value
        && std::is_same<Segment, strategy::buffer::side_straight>::value
        && (std::is_same<Distance, strategy::buffer::distance_symmetric<coordinate_type_t<Point>>>::value
            || std::is_same<Distance, strategy::buffer::distance_asymmetric<coordinate_type_t<Point>>>::value),
        bool>::type
    apply_source_round_join(Join const& join, Segment const&,
        Distance const& distance, Point const& before, Point const& vertex, Point const& after,
        OutputPoint const& first, OutputPoint const& last,
        strategy::buffer::buffer_side_selector side, Range& out)
    {
        namespace local = cartesian_mp_detail;
        using namespace strategy::buffer;
        auto const radius = distance.apply(vertex, after, side);
        if (distance.negative() || !(radius > 0)
            || distance.apply(vertex, after, buffer_side_left) != distance.apply(vertex, after, buffer_side_right)
            || !local::admissible(before) || !local::admissible(vertex) || !local::admissible(after))
            return false;
        auto const a = local::direction(before, vertex), b = local::direction(vertex, after);
        // side_straight uses the left normal of the directed source edge.
        int const start = local::octant({-a.y, a.x}), finish = local::octant({-b.y, b.x});
        if (start < 0 || finish < 0) return false;
        int const span = (start + 8 - finish) % 8;
        if (span <= 0 || span >= 4 || join.m_points_per_circle > 4096
            || geometry::equal_to<OutputPoint>()(first, last))
            return false;

        // Derive the sample count without rounding an exact octant span
        // through atan2 before ceil; only replace fans with a cardinal sample.
        int const n = static_cast<int>((join.m_points_per_circle * span + 7) / 8);
        int const distance_to_cardinal = start % 2 == 0 ? 2 : 1;
        if (distance_to_cardinal >= span || (distance_to_cardinal * n) % span != 0)
            return false;
        int const cardinal_sample = distance_to_cardinal * n / span;

        using calculation = coordinate_type_t<OutputPoint>;
        out.push_back(first);
        for (int i = 1; i < n; ++i)
        {
            OutputPoint p(vertex);
            if (i == cardinal_sample)
            {
                switch ((start + 8 - distance_to_cardinal) % 8)
                {
                case 0: set<0>(p, get<0>(vertex) + radius); break;
                case 2: set<1>(p, get<1>(vertex) + radius); break;
                case 4: set<0>(p, get<0>(vertex) - radius); break;
                case 6: set<1>(p, get<1>(vertex) - radius); break;
                }
            }
            else
            {
                calculation const angle = geometry::math::pi<calculation>() * calculation(start * n - span * i)
                                        / calculation(4 * n);
                set<0>(p, get<0>(vertex) + radius * cos(angle));
                set<1>(p, get<1>(vertex) + radius * sin(angle));
            }
            out.push_back(p);
        }
        out.push_back(last);
        return true;
    }
};

}} // namespace strategies::buffer

// Concrete-geometry entry point. Including this header does not specialize
// default_strategy or change the meaning of the seven-argument overload.
template <typename Input, typename Output, typename Distance, typename Side,
          typename Join, typename End, typename Point>
inline void buffer(Input const& input, Output& output, Distance const& distance,
                   Side const& side, Join const& join, End const& end, Point const& point,
                   strategies::buffer::cartesian_mp const& strategy)
{
    static_assert(std::is_same<cs_tag_t<Input>, cartesian_tag>::value
                  && std::is_same<cs_tag_t<Output>, cartesian_tag>::value,
                  "cartesian_mp requires concrete Cartesian buffer geometries");
    concepts::check<Input const>();
    concepts::check<Output>();
    geometry::clear(output);
    dispatch::buffer_all<Input, Output>::apply(input, output, distance, side, join, end, point, strategy);
}

}} // namespace boost::geometry

#endif
