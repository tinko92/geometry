// Boost.Geometry
// Copyright (c) 2026
// Distributed under the Boost Software License, Version 1.0.

#ifndef BOOST_GEOMETRY_ALGORITHMS_BUFFER_LINESTRING_ARRANGEMENT_CERTIFIED_HPP
#define BOOST_GEOMETRY_ALGORITHMS_BUFFER_LINESTRING_ARRANGEMENT_CERTIFIED_HPP

#include <algorithm>
#include <map>
#include <set>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

#include <boost/geometry/algorithms/buffer_linestring_star_certified.hpp>
#include <boost/geometry/algorithms/buffer_polygon_certified.hpp>

namespace boost { namespace geometry
{

#ifndef DOXYGEN_NO_DETAIL
namespace detail { namespace certified_star_buffer
{

namespace rational_arr = rational_chord_arrangement;
namespace rational_line = certified_linestring_buffer;

// This arrangement uses the same exact producers and local predicate field
// as the star sweep, but has no radial ordering or common-witness assumption.
class algebraic_chord_union
{
    static constexpr id max_pieces = 16;
    static constexpr id max_edges = 384;
    static constexpr id max_pair_checks = 8192;
    static constexpr id max_atoms = 4096;
    static constexpr id max_cycles = 16;

    struct rectangle
    {
        double lx, ly, ux, uy;

        explicit rectangle(box_point const& p)
            : lx(p.x.lower()), ly(p.y.lower()), ux(p.x.upper()), uy(p.y.upper()) {}

        void extend(rectangle const& r)
        {
            lx = (std::min)(lx, r.lx); ly = (std::min)(ly, r.ly);
            ux = (std::max)(ux, r.ux); uy = (std::max)(uy, r.uy);
        }

        bool overlaps(rectangle const& r) const
        { return lx <= r.ux && r.lx <= ux && ly <= r.uy && r.ly <= uy; }
    };

    struct piece_range { id first, end; rectangle box; };
    struct edge_cuts
    {
        id owner;
        rectangle box;
        bool ascending;
        std::vector<knot> cuts;
    };

public:
    void append(capsule::rational_segment const& source, exact_number const& radius, exact_number const& sagitta)
    {
        if (m_pieces.size() == max_pieces || m_edges.size() == max_edges) { throw uncertified{}; }
        auto const span = m_geometry.append_capsule(source, radius, sagitta, max_edges - m_edges.size());
        rectangle box(vertex(m_geometry.m_edges[span.first].from).bounds);
        for (id i = span.first; i < span.second; ++i)
        {
            auto const& e = m_geometry.m_edges[i];
            auto const a = vertex(e.from), b = vertex(e.to);
            rectangle bounds(a.bounds);
            bounds.extend(rectangle(b.bounds));
            box.extend(bounds);
            m_edges.push_back({m_pieces.size(), bounds, compare(a, b) < 0, {a, b}});
        }
        m_pieces.push_back({span.first, span.second, box});
    }

    template <typename Polygon>
    bool build(rational_xy const& anchor, exact_number const& displacement, Polygon& output)
    {
        node();
        auto rings = cycles();
        return certify_and_export(rings, anchor, displacement, output);
    }

private:
    knot vertex(id i) const { return {knot_kind::vertex, i, 0, m_geometry.m_vertices[i].bounds, 0, false}; }

    int compare_axis(knot const& a, knot const& b, bool y)
    {
        if (a.kind == b.kind && a.a == b.a && a.b == b.b) { return 0; }
        auto const s = ec::sign(y ? a.bounds.y - b.bounds.y : a.bounds.x - b.bounds.x);
        if (s.certain) { return s.value; }
        m_geometry.fallback();
        envelope::local f(m_geometry);
        auto const p = f.evaluate(a), q = f.evaluate(b);
        return rn::compare(f.f, y ? p.y : p.x, y ? q.y : q.x);
    }

    int compare(knot const& a, knot const& b)
    {
        int const x = compare_axis(a, b, false);
        return x != 0 ? x : compare_axis(a, b, true);
    }

    void node()
    {
        id checks = 0;
        for (id i = 0; i < m_edges.size(); ++i)
        {
            for (id j = i + 1; j < m_edges.size(); ++j)
            {
                auto& a = m_edges[i];
                auto& b = m_edges[j];
                if (a.owner == b.owner || !a.box.overlaps(b.box)) { continue; }
                if (++checks > max_pair_checks) { throw uncertified{}; }
                envelope::local f(m_geometry);
                auto const s = f.segment(i), t = f.segment(j);
                auto const hit = rn::intersect(f.f, s, t);
                if (hit.kind == rn::intersection_kind::disjoint) { continue; }
                if (hit.kind == rn::intersection_kind::degenerate) { throw uncertified{}; }
                auto const insert = [&](rn::parameter const& u, rn::parameter const& v)
                {
                    // An overlap end is an original endpoint. Retain it as
                    // such; only an interior-interior hit needs two producers.
                    auto const& ab = m_geometry.m_edges[i];
                    auto const& cd = m_geometry.m_edges[j];
                    knot p = f.f.sign(u.n) == 0 ? vertex(ab.from)
                        : f.f.sign(f.f.subtract(u.d, u.n)) == 0 ? vertex(ab.to)
                        : f.f.sign(v.n) == 0 ? vertex(cd.from)
                        : f.f.sign(f.f.subtract(v.d, v.n)) == 0 ? vertex(cd.to)
                        : knot{knot_kind::crossing, i, j, m_geometry.crossing_bounds(i, j), 0, false};
                    a.cuts.push_back(p);
                    b.cuts.push_back(p);
                };
                insert(hit.t0, hit.u0);
                if (hit.kind == rn::intersection_kind::overlap) { insert(hit.t1, hit.u1); }
            }
        }
    }

    int midpoint_side(id e, knot const& a, knot const& b, box_point const& mid)
    {
        auto const& edge = m_geometry.m_edges[e];
        auto const& p = m_geometry.m_vertices[edge.from].bounds;
        auto const& q = m_geometry.m_vertices[edge.to].bounds;
        auto const s = ec::sign(det(q - p, mid - p));
        if (s.certain) { return s.value; }
        m_geometry.fallback();

        // Both nodes are on the same producer: they use at most three
        // roots together, plus the tested facet's root. Keep the midpoint
        // homogeneous; all denominators are positive.
        envelope::local f(m_geometry);
        auto const u = f.evaluate(a), v = f.evaluate(b);
        auto const edge_exact = f.segment(e);
        auto const offset = [&](rn::parameter const& x, rn::parameter const& y, rn::value const& base)
        {
            auto const n = f.f.add(f.f.multiply(x.n, y.d), f.f.multiply(y.n, x.d));
            auto const d = f.f.scale(f.f.multiply(x.d, y.d), 2);
            return f.f.subtract(n, f.f.multiply(d, base));
        };
        return f.f.sign(rn::cross(f.f, rn::subtract(f.f, edge_exact.to, edge_exact.from),
            {offset(u.x, v.x, edge_exact.from.x), offset(u.y, v.y, edge_exact.from.y)}));
    }

    bool exposed(id e, knot const& a, knot const& b)
    {
        ec::interval const half{0.5, 0};
        box_point const mid{(a.bounds.x + b.bounds.x) * half, (a.bounds.y + b.bounds.y) * half};
        rectangle const bounds(mid);
        for (id k = 0; k < m_pieces.size(); ++k)
        {
            auto const& p = m_pieces[k];
            if (k == m_edges[e].owner || !p.box.overlaps(bounds)) { continue; }
            bool covers_right = true;
            for (id i = p.first; i < p.end; ++i)
            {
                int const s = midpoint_side(i, a, b, mid);
                // Exact noding excludes a transverse boundary through the
                // midpoint. A zero incident facet only contains the right
                // neighborhood if it runs opposite to this atom.
                if (s < 0 || (s == 0 && m_edges[i].ascending == m_edges[e].ascending))
                {
                    covers_right = false;
                    break;
                }
            }
            if (covers_right) { return false; }
        }
        return true;
    }

    struct link { knot to; id edge; };

    std::vector<contour> cycles()
    {
        auto const less = [&](knot const& a, knot const& b) { return compare(a, b) < 0; };
        std::map<knot, link, decltype(less)> next(less);
        std::set<knot, decltype(less)> incoming(less);
        id atoms = 0;
        for (id i = 0; i < m_edges.size(); ++i)
        {
            auto& c = m_edges[i].cuts;
            std::sort(c.begin(), c.end(), less);
            c.erase(std::unique(c.begin(), c.end(), [&](knot const& a, knot const& b) { return compare(a, b) == 0; }),
                    c.end());
            for (id j = 1; j < c.size(); ++j)
            {
                if (++atoms > max_atoms) { throw uncertified{}; }
                auto const& a = m_edges[i].ascending ? c[j - 1] : c[j];
                auto const& b = m_edges[i].ascending ? c[j] : c[j - 1];
                if (!exposed(i, a, b)) { continue; }
                auto const added = next.emplace(a, link{b, i});
                if (added.second ? !incoming.insert(b).second : compare(added.first->second.to, b) != 0)
                {
                    throw uncertified{}; // Touching or branching boundary.
                }
            }
        }

        std::vector<contour> out;
        while (!next.empty())
        {
            if (out.size() == max_cycles) { throw uncertified{}; }
            contour ring;
            auto const first = next.begin()->first;
            auto cur = first;
            do
            {
                auto const it = next.find(cur);
                if (it == next.end()) { throw uncertified{}; }
                ring.push_back({it->first, it->second.edge});
                cur = it->second.to;
                next.erase(it);
            }
            while (compare(cur, first) != 0);
            if (ring.size() < 3) { throw uncertified{}; }
            out.push_back(std::move(ring));
        }
        return out;
    }

    int winding_order(contour const& ring)
    {
        // The lexicographically least vertex is a strict extreme turn of
        // a simple noded ring. The original directions suffice here; no
        // area sum involving all of the linestring's roots is required.
        envelope::local f(m_geometry);
        auto const a = f.segment(ring.back().edge), b = f.segment(ring.front().edge);
        return f.f.sign(rn::cross(f.f, rn::subtract(f.f, a.to, a.from), rn::subtract(f.f, b.to, b.from)));
    }

    bool inside(contour const& ring, knot const& p)
    {
        bool result = false;
        for (id i = 0; i < ring.size(); ++i)
        {
            auto const& a = ring[i];
            auto const& b = ring[(i + 1) % ring.size()];
            bool const high_a = compare_axis(a.at, p, true) > 0, high_b = compare_axis(b.at, p, true) > 0;
            if (high_a == high_b) { continue; }
            int const side = m_geometry.side(a.edge, p);
            if (side == 0) { throw uncertified{}; }
            if ((side > 0) == high_b) { result = !result; }
        }
        return result;
    }

    template <typename Polygon>
    bool certify_and_export(std::vector<contour>& ideal, rational_xy const& anchor,
                            exact_number const& displacement, Polygon& output)
    {
        using coordinate = coordinate_type_t<point_type_t<Polygon>>;
        static_assert(std::is_integral<coordinate>::value || std::is_same<coordinate, float>::value
                   || std::is_same<coordinate, double>::value, "Chord arrangement export requires native coordinates");
        id shell = ideal.size();
        for (id i = 0; i < ideal.size(); ++i)
        {
            int const turn = winding_order(ideal[i]);
            if (turn == 0) { return false; }
            if (turn < 0) { continue; }
            if (shell != ideal.size()) { return false; }
            shell = i;
        }
        if (shell == ideal.size()) { return false; }
        std::swap(ideal.front(), ideal[shell]);
        for (id i = 1; i < ideal.size(); ++i)
        {
            if (!inside(ideal.front(), ideal[i].front().at)) { return false; }
            for (id j = 1; j < i; ++j)
            {
                if (inside(ideal[j], ideal[i].front().at) || inside(ideal[i], ideal[j].front().at)) { return false; }
            }
        }

        std::vector<rational_arr::ring> stored;
        for (auto const& ring : ideal)
        {
            rational_line::piece vertices;
            for (auto const& m : ring)
            {
                coordinate x, y;
                if (!m_geometry.native_euclidean(m.at, anchor, displacement, x, y)) { return false; }
                rational_xy const p{common::rational_value(x), common::rational_value(y)};
                if (!vertices.empty() && capsule::same_point(vertices.back(), p)) { return false; }
                vertices.push_back(p);
            }
            if (capsule::same_point(vertices.back(), vertices.front())) { return false; }
            stored.emplace_back(std::move(vertices));
            auto const& area = stored.back().area2;
            if (stored.size() == 1 ? area <= 0 : area >= 0) { return false; }
        }
        // Independent exact rational read-back topology: simplicity, all
        // rings disjoint, all holes inside the shell and no nested holes.
        if (!certified_polygon_buffer::valid_source(stored)) { return false; }

        Polygon candidate;
        bool constexpr cw = geometry::point_order<Polygon>::value == clockwise;
        bool constexpr closed = closure<Polygon>::value == geometry::closed;
        if (!rational_line::write_ring_lossless(stored.front().points, exterior_ring(candidate), cw, closed))
        { return false; }
        auto& holes = interior_rings(candidate);
        range::resize(holes, stored.size() - 1);
        auto h = boost::begin(holes);
        for (id i = 1; i < stored.size(); ++i, ++h)
        {
            if (h == boost::end(holes) || !rational_line::write_ring_lossless(stored[i].points, *h, cw, closed))
            { return false; }
        }
        if (h != boost::end(holes)) { return false; }
        output = std::move(candidate);
        return true;
    }

    envelope m_geometry;
    std::vector<piece_range> m_pieces;
    std::vector<edge_cuts> m_edges;
};

}} // namespace detail::certified_star_buffer
#endif // DOXYGEN_NO_DETAIL

/*!
\brief Bounded opt-in native round buffer of a 2D Cartesian Linestring.
\details No common interior witness is required. Nodes the exact inscribed
    chord capsules before conversion. Requires one simple shell with disjoint
    (non-touching) holes, and certifies the same ring topology and corresponding
    boundary displacement at most max_displacement on the native grid.
    The vertex bound is Euclidean, checked against the exact algebraic contour.
    A zero displacement requires lossless export; nonzero displacement does not
    promise nearest-rounded vertices.
    Requires 0 < sagitta < radius and 0 <= max_displacement < radius-sagitta.
    The stored region contains the radius-(r-sagitta-max_displacement)
    neighborhood and is contained in the radius-(r+max_displacement)
    neighborhood. Limits apply. This does not certify smooth-buffer topology.
\return success, invalid_input, or cannot_certify. Failure leaves output unchanged.
*/
template <typename Linestring, typename Polygon, typename Radius, typename Error, typename Displacement>
[[nodiscard]] inline certified_buffer_status buffer_linestring_arrangement_certified(
    Linestring const& input, Polygon& output, Radius const& radius,
    Error const& sagitta, Displacement const& max_displacement)
{
    namespace impl = detail::certified_star_buffer;
    using status = certified_buffer_status;
    impl::capsule::check_types<linestring_tag, Linestring, Polygon>();
    auto it = boost::begin(input);
    auto const end = boost::end(input);
    if (it == end) { return status::invalid_input; }

    impl::rational_xy anchor, previous;
    impl::exact_number r, e, displacement;
    try
    {
        anchor = impl::common::read_point(*it++);
        r = impl::common::rational_value(radius);
        e = impl::common::rational_value(sagitta);
        displacement = impl::common::rational_value(max_displacement);
        if (r <= 0 || e < 0 || e >= r || displacement < 0 || displacement >= r - e)
        { return status::invalid_input; }
        if (e == 0) { return status::cannot_certify; }
    }
    catch (std::domain_error const&) { return status::invalid_input; }

    impl::algebraic_chord_union arrangement;
    for (; it != end; ++it)
    {
        impl::rational_xy current;
        try { current = impl::difference(impl::common::read_point(*it), anchor); }
        catch (std::domain_error const&) { return status::invalid_input; }
        if (impl::capsule::same_point(previous, current)) { continue; }
        try { arrangement.append(impl::capsule::rational_segment(previous, current), r, e); }
        catch (impl::uncertified const&) { return status::cannot_certify; }
        previous = std::move(current);
    }
    try { return arrangement.build(anchor, displacement, output) ? status::success : status::cannot_certify; }
    catch (impl::uncertified const&) { return status::cannot_certify; }
    catch (std::overflow_error const&) { return status::cannot_certify; }
    catch (std::domain_error const&) { return status::cannot_certify; }
}

}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_ALGORITHMS_BUFFER_LINESTRING_ARRANGEMENT_CERTIFIED_HPP
