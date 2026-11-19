// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.

#include <boost/core/lightweight_test.hpp>
#include <boost/geometry/strategies/cartesian/exact/borrowed_output.hpp>
#include <random>
#include <vector>

namespace ec = boost::geometry::strategies::detail::exact_cartesian;
using mode = ec::supported_append_mode;
using index_type = boost::geometry::signed_size_type;

struct point
{
    int x = 0, y = 0, producer = 0;
    point const* native = nullptr;
};

struct equal_points
{
    bool operator()(point const& a, point const& b) const
    { return a.x == b.x && a.y == b.y; }
};

struct side_predicate
{
    int apply(point const& a, point const& b, point const& p) const
    {
        int const cross = (b.x - a.x) * (p.y - a.y) - (b.y - a.y) * (p.x - a.x);
        return (cross > 0) - (cross < 0);
    }
    int direction_code(point const& a, point const& b, point const& p) const
    {
        int const dot = (b.x - a.x) * (p.x - b.x) + (b.y - a.y) * (p.y - b.y);
        return (dot > 0) - (dot < 0);
    }
};

// The production supported_vector contract, without Root or numerical headers.
struct range : std::vector<point>
{
    int pending = -1;
    std::vector<int> incoming;
    void push_back(point const& p)
    { std::vector<point>::push_back(p); incoming.push_back(pending); }
    void resize(std::size_t n)
    { std::vector<point>::resize(n); incoming.resize(n); }
};

void append(range& out, point const& p, int support, mode m = mode::spikes)
{
    out.pending = support;
    ec::append_supported(out, p, side_predicate{}, equal_points{}, m);
    BOOST_TEST_EQ(out.size(), out.incoming.size());
    BOOST_TEST_EQ(out.pending, support);
}

void check_labels(range const& out, std::vector<int> const& expected)
{
    BOOST_TEST(out.incoming == expected);
    BOOST_TEST_EQ(out.size(), out.incoming.size());
}

void legacy_append(std::vector<point>& out, point p, mode m)
{
    equal_points equal;
    side_predicate side;
    if (m != mode::none && !out.empty()
        && (m == mode::duplicates || out.size() == 1) && equal(out.back(), p)) return;
    out.push_back(p);
    if (m == mode::none || m == mode::duplicates) return;
    while (out.size() >= 3
           && side.apply(out[out.size() - 3], out[out.size() - 2], p) == 0
           && (m == mode::collinear
               || side.direction_code(out[out.size() - 3], out[out.size() - 2], p) < 1))
    {
        out.resize(out.size() - 2);
        out.push_back(p);
    }
}

void test_vertex_sequences()
{
    std::mt19937 random(7391);
    for (mode m : {mode::none, mode::duplicates, mode::spikes, mode::collinear})
    {
        for (int trial = 0; trial < 100; ++trial)
        {
            range out;
            std::vector<point> expected;
            for (int i = 0; i < 40; ++i)
            {
                point p{int(random() % 5) - 2, int(random() % 5) - 2, i};
                append(out, p, i, m);
                legacy_append(expected, p, m);
                BOOST_TEST_EQ(out.size(), expected.size());
                for (std::size_t j = 0; j < out.size() && j < expected.size(); ++j)
                {
                    BOOST_TEST(equal_points{}(out[j], expected[j]));
                    BOOST_TEST_EQ(out[j].producer, expected[j].producer);
                }
            }
        }
    }
}

void test_cancellation_and_concurrency()
{
    for (mode m : {mode::spikes, mode::collinear})
    {
        range out;
        append(out, {-1, 1, 1}, -1, m);
        append(out, {0, 0, 2}, 10, m);
        append(out, {2, 0, 3}, 20, m);
        append(out, {0, 0, 4}, 21, m);
        BOOST_TEST_EQ(out.size(), 2u);
        BOOST_TEST_EQ(out.back().producer, 4);
        check_labels(out, {-1, 10});
        append(out, {0, 0, 5}, 30, m);
        check_labels(out, {-1, 10});
        append(out, {0, 1, 6}, 40, m);
        check_labels(out, {-1, 10, 40});
    }

    // Three different line-pair producers at one concurrent intersection.
    range concurrent;
    append(concurrent, {0, 0, 12}, -1, mode::duplicates);
    append(concurrent, {0, 0, 13}, 13, mode::duplicates);
    append(concurrent, {0, 0, 23}, 23, mode::duplicates);
    BOOST_TEST_EQ(concurrent.size(), 1u);
    BOOST_TEST_EQ(concurrent.front().producer, 12);
    check_labels(concurrent, {-1});

    range closed;
    append(closed, {0, 0, 1}, -1);
    append(closed, {0, 2, 2}, 10);
    append(closed, {2, 2, 3}, 20);
    append(closed, {2, 0, 4}, 30);
    append(closed, closed.front(), 40);
    BOOST_TEST(equal_points{}(closed.front(), closed.back()));
    check_labels(closed, {-1, 10, 20, 30, 40});
}

void test_reversal()
{
    range out;
    for (int i = 0; i < 4; ++i) append(out, {i, i, i}, i * 10, mode::none);
    out.pending = 99;
    ec::reverse_supported<false>(out);
    check_labels(out, {0, 30, 20, 10});
    BOOST_TEST_EQ(out.front().producer, 3);
    BOOST_TEST_EQ(out.back().producer, 0);
    BOOST_TEST_EQ(out.pending, 99);
    ec::reverse_supported<false>(out);
    check_labels(out, {0, 10, 20, 30});

    range closed;
    append(closed, {0, 0, 1}, 40, mode::none);
    append(closed, {0, 2, 2}, 10, mode::none);
    append(closed, {2, 2, 3}, 20, mode::none);
    append(closed, {2, 0, 4}, 30, mode::none);
    append(closed, {0, 0, 5}, 40, mode::none);
    ec::reverse_supported<true>(closed);
    check_labels(closed, {10, 40, 30, 20, 10});
    BOOST_TEST_EQ(closed.front().producer, 5);
    BOOST_TEST_EQ(closed.back().producer, 1);
    ec::reverse_supported<true>(closed);
    check_labels(closed, {40, 10, 20, 30, 40});
    range empty;
    ec::reverse_supported<true>(empty);
    append(empty, {0, 0, 1}, -1);
    ec::reverse_supported<false>(empty);
    check_labels(empty, {-1});
}

template <bool Reverse, bool RemoveSpikes>
void copy_native(std::vector<point> const& input, bool circular,
                 index_type segment, index_type to_index, range& out)
{
    struct identifier { index_type segment_index; } id{segment};
    index_type const size = static_cast<index_type>(input.size());
    auto point_at = [&](index_type i, bool reverse)
    {
        auto const index = circular ? (reverse ? (size - i % size) % size : i % size)
                                    : (reverse ? size - i - 1 : i);
        auto const it = input.begin() + index;
        point result = *it;
        result.native = &*it;
        return result;
    };
    auto support_at = [&](index_type i, bool reverse)
    {
        BOOST_TEST(i >= 0 && i < size - (circular ? 0 : 1));
        auto const index = circular ? (reverse ? (size - (i + 1) % size) % size : i)
                                    : (reverse ? size - i - 2 : i);
        return 100 + int(index);
    };
    ec::copy_supported<Reverse, RemoveSpikes>(size + (circular ? 1 : 0), circular,
        id, to_index, out, side_predicate{}, equal_points{}, point_at, support_at);
}

void test_copy_native()
{
    std::vector<point> const square{{0, 0, 1}, {0, 4, 2}, {4, 4, 3}, {4, 0, 4}};
    range out;
    append(out, square.front(), -1);
    copy_native<false, true>(square, true, 0, 3, out);
    check_labels(out, {-1, 100, 101, 102});
    BOOST_TEST_EQ(out.pending, 103);
    for (std::size_t i = 1; i < out.size(); ++i) BOOST_TEST(out[i].native == &square[i]);
    ec::append_supported(out, square.front(), side_predicate{}, equal_points{}, mode::spikes);
    check_labels(out, {-1, 100, 101, 102, 103});

    range reversed;
    append(reversed, square.front(), -1);
    copy_native<true, true>(square, true, 0, 4, reversed);
    check_labels(reversed, {-1, 103, 102, 101, 100});
    BOOST_TEST_EQ(reversed.pending, 103);
    BOOST_TEST(reversed[1].native == &square[3]);
    BOOST_TEST(reversed.back().native == &square[0]);

    range wrapped;
    append(wrapped, square.back(), -1);
    copy_native<false, true>(square, true, 3, 1, wrapped);
    check_labels(wrapped, {-1, 103, 100});
    BOOST_TEST_EQ(wrapped.pending, 101);
    range round_trip;
    append(round_trip, square[1], -1);
    copy_native<false, true>(square, true, 1, 1, round_trip);
    check_labels(round_trip, {-1, 101, 102, 103, 100});
    BOOST_TEST_EQ(round_trip.pending, 101);

    std::vector<point> const line{{0, 0, 1}, {2, 0, 2}, {1, 0, 3}, {1, 1, 4}};
    range spikes, no_spikes;
    copy_native<false, false>(line, false, -1, 3, spikes);
    check_labels(spikes, {-1, 100, 101, 102});
    copy_native<false, true>(line, false, -1, 3, no_spikes);
    check_labels(no_spikes, {-1, 101, 102});
    BOOST_TEST_EQ(no_spikes.pending, 102);
    copy_native<false, true>(line, false, 3, 2, no_spikes);
    check_labels(no_spikes, {-1, 101, 102});
    BOOST_TEST_EQ(no_spikes.pending, 102);
}

void test_copy_working()
{
    range working;
    append(working, {0, 0, 1}, -1, mode::none);
    append(working, {0, 4, 2}, 101, mode::none);
    append(working, {4, 4, 3}, 202, mode::none);
    append(working, {4, 0, 4}, 303, mode::none);
    append(working, {0, 0, 5}, 404, mode::none);
    struct identifier { index_type segment_index; } id{0};
    auto point_at = [&](index_type i, bool reverse)
    { return working[reverse ? (4 - i % 4) % 4 : i % 4]; };
    auto support_at = [&](index_type i, bool reverse)
    {
        auto const edge = reverse ? (4 - (i + 1) % 4) % 4 : i;
        // Closing support lives at the appended endpoint, not the unused slot zero.
        return working.incoming[edge + 1];
    };
    range out;
    append(out, working.front(), -1);
    ec::copy_supported<true, true>(5, true, id, 4, out, side_predicate{}, equal_points{},
                                   point_at, support_at);
    check_labels(out, {-1, 404, 303, 202, 101});
    BOOST_TEST_EQ(out.pending, 404);
    BOOST_TEST_EQ(out[1].producer, 4);
    ec::reverse_supported<true>(out);
    check_labels(out, {404, 101, 202, 303, 404});
}

int main()
{
    test_vertex_sequences();
    test_cancellation_and_concurrency();
    test_reversal();
    test_copy_native();
    test_copy_working();
    return boost::report_errors();
}
