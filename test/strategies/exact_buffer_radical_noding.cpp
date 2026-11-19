// Boost.Geometry
// Unit Test
// Copyright (c) 2026
// Distributed under the Boost Software License, Version 1.0.

#include <boost/geometry/strategies/cartesian/exact/radical_buffer_noding.hpp>
#include <boost/core/lightweight_test.hpp>

namespace ec = boost::geometry::strategies::detail::exact_cartesian;
namespace rb = ec::round_buffer;
namespace rn = rb::radical_noding;
using number = ec::exact_number;
using kind = rn::intersection_kind;

void test_independent_supports()
{
    rn::field f;
    rn::segment a, b;
    // Turn to the right. The finite inward offsets actually cross; their
    // exact node uses independent square classes 2 and 5.
    BOOST_TEST(rn::lift_offset(f, rb::straight_offset({0, 0}, {12, 12}, 4, rb::offset_side::right), a));
    BOOST_TEST(rn::lift_offset(f, rb::straight_offset({12, 12}, {24, 18}, 4, rb::offset_side::right), b));
    BOOST_TEST_EQ(f.rank(), 2u);
    auto const i = rn::intersect(f, a, b);
    BOOST_TEST(i.kind == kind::contact);
    if (i.kind != kind::contact) { return; }
    BOOST_TEST(f.sign(i.t0.n) > 0 && f.sign(f.subtract(i.t0.d, i.t0.n)) > 0);
    BOOST_TEST(f.sign(i.u0.n) > 0 && f.sign(f.subtract(i.u0.d, i.u0.n)) > 0);
    auto const p = rn::at(f, a, i.t0), q = rn::at(f, b, i.u0);
    BOOST_TEST_EQ(rn::compare(f, p.x, q.x), 0);
    BOOST_TEST_EQ(rn::compare(f, p.y, q.y), 0);

    rn::value sqrt2, sqrt5;
    BOOST_TEST(f.try_sqrt(2, sqrt2));
    BOOST_TEST(f.try_sqrt(5, sqrt5));
    // Solve y=x-4sqrt(2) and y=x/2+6-2sqrt(5) independently.
    auto const x = f.add(f.rational(12), f.subtract(f.scale(sqrt2, 8), f.scale(sqrt5, 4)));
    auto const y = f.add(f.rational(12), f.scale(f.subtract(sqrt2, sqrt5), 4));
    BOOST_TEST_EQ(rn::compare(f, p.x, {x, f.rational(1)}), 0);
    BOOST_TEST_EQ(rn::compare(f, p.y, {y, f.rational(1)}), 0);
    number rational = 101;
    BOOST_TEST_NOT(f.try_rational_quotient(p.x.n, p.x.d, rational));
    BOOST_TEST_EQ(rational, 101);

    // Appending a third square class does not change existing nodes. A
    // second pair's cut can be ordered on the shared first edge.
    rn::value sqrt3;
    BOOST_TEST(f.try_sqrt(3, sqrt3));
    auto const xx = f.add(f.rational(13), f.scale(sqrt3, number(1) / 100));
    rn::segment const vertical{{xx, f.rational(0)}, {xx, f.rational(30)}};
    auto const next = rn::intersect(f, a, vertical);
    BOOST_TEST(next.kind == kind::contact);
    if (next.kind == kind::contact) { BOOST_TEST(rn::compare(f, next.t0, i.t0) < 0); }
}

void test_contacts_and_overlaps()
{
    rn::field f;
    rn::segment a, b;
    BOOST_TEST(rn::lift_offset(f, rb::straight_offset({0, 0}, {4, 4}, 1, rb::offset_side::left), a));
    BOOST_TEST(rn::lift_offset(f, rb::straight_offset({3, 3}, {1, 1}, 1, rb::offset_side::right), b));
    BOOST_TEST_EQ(f.rank(), 1u); // rationally related lengths, exactly coincident
    auto const overlap = rn::intersect(f, a, b);
    BOOST_TEST(overlap.kind == kind::overlap);
    auto const rat = [&](number const& x) { return rn::parameter{f.rational(x), f.rational(1)}; };
    if (overlap.kind == kind::overlap)
    {
        BOOST_TEST_EQ(rn::compare(f, overlap.t0, rat(number(1) / 4)), 0);
        BOOST_TEST_EQ(rn::compare(f, overlap.t1, rat(number(3) / 4)), 0);
        BOOST_TEST_EQ(rn::compare(f, overlap.u0, rat(1)), 0);
        BOOST_TEST_EQ(rn::compare(f, overlap.u1, rat(0)), 0);
    }

    auto const end = a.to;
    auto const step = [&](rn::point const& p, number const& dx, number const& dy)
    { return rn::point{f.add(p.x, f.rational(dx)), f.add(p.y, f.rational(dy))}; };
    rn::segment later{end, step(end, 1, 1)};
    auto const touch = rn::intersect(f, a, later);
    BOOST_TEST(touch.kind == kind::contact);
    if (touch.kind == kind::contact) { BOOST_TEST_EQ(rn::compare(f, touch.t0, rat(1)), 0); }
    // Move just the supporting line of the second edge; no rounding or
    // tolerance is allowed to convert the parallel gap back into a contact.
    number const tiny = number(1) / (boost::multiprecision::cpp_int(1) << 300);
    later = {step(later.from, 0, tiny), step(later.to, 0, tiny)};
    BOOST_TEST(rn::intersect(f, a, later).kind == kind::disjoint);
    BOOST_TEST(rn::intersect(f, a, {end, end}).kind == kind::degenerate);
}

int main()
{
    test_independent_supports();
    test_contacts_and_overlaps();
    return boost::report_errors();
}
