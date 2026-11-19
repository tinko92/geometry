// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.

#include <boost/core/lightweight_test.hpp>
#include <boost/geometry/algorithms/difference.hpp>
#include <boost/geometry/algorithms/intersection.hpp>
#include <boost/geometry/algorithms/is_valid.hpp>
#include <boost/geometry/algorithms/relation.hpp>
#include <boost/geometry/algorithms/sym_difference.hpp>
#include <boost/geometry/algorithms/union.hpp>
#include <boost/geometry/geometries/multi_point.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/strategies/relate/cartesian.hpp>
#include <boost/geometry/strategies/detail/point_conversion.hpp>
#include <boost/geometry/util/sorted_point_view.hpp>
#include <boost/iterator/iterator_adaptor.hpp>
#include <boost/iterator/iterator_facade.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <initializer_list>
#include <iterator>
#include <type_traits>

namespace bg = boost::geometry;

std::size_t source_copies = 0, source_moves = 0;

struct observed_copy
{
    observed_copy() = default;
    observed_copy(observed_copy const&) { ++source_copies; }
    observed_copy(observed_copy&&) { ++source_moves; }
    observed_copy& operator=(observed_copy const&) { ++source_copies; return *this; }
    observed_copy& operator=(observed_copy&&) { ++source_moves; return *this; }
};

struct no_copy
{
    no_copy() = default;
    no_copy(no_copy const&) = delete;
    no_copy& operator=(no_copy const&) = delete;
};

template <typename Copy>
struct input_point : Copy
{
    double x = 0, y = 0;
    input_point() = default;
    input_point(double x, double y) : x(x), y(y) {}
};

using observed_point = input_point<observed_copy>;
using xy = std::array<double, 2>;

template <typename Point, std::size_t N>
struct fixed_multipoint
{
    using value_type = Point;
    using iterator = Point*;
    using const_iterator = Point const*;
    std::array<Point, N> points;

    explicit fixed_multipoint(std::initializer_list<xy> coordinates)
    {
        BOOST_TEST_EQ(coordinates.size(), N);
        std::size_t i = 0;
        for (auto const& p : coordinates)
        {
            points[i].x = p[0];
            points[i++].y = p[1];
        }
    }
    iterator begin() { return points.data(); }
    iterator end() { return points.data() + N; }
    const_iterator begin() const { return points.data(); }
    const_iterator end() const { return points.data() + N; }
    std::size_t size() const { return N; }
};

struct proxy_point
{
    observed_point const* source;
    operator observed_point() const { return observed_point(source->x, source->y); }
};

// Deliberately not adapted as a Boost.Geometry point.
struct convertible_proxy
{
    observed_point const* source;
    operator observed_point() const { return observed_point(source->x, source->y); }
};

template <typename Iterator>
struct multipoint_view
{
    using value_type = typename std::iterator_traits<Iterator>::value_type;
    using iterator = Iterator;
    using const_iterator = Iterator;
    Iterator first, last;
    Iterator begin() const { return first; }
    Iterator end() const { return last; }
};

namespace boost { namespace geometry { namespace traits
{
template <typename C> struct tag<input_point<C>> { using type = point_tag; };
template <typename C> struct coordinate_type<input_point<C>> { using type = double; };
template <typename C> struct coordinate_system<input_point<C>> { using type = cs::cartesian; };
template <typename C> struct dimension<input_point<C>> : std::integral_constant<int, 2> {};
template <typename C, std::size_t D> struct access<input_point<C>, D>
{
    static double get(input_point<C> const& p) { return D == 0 ? p.x : p.y; }
    static void set(input_point<C>& p, double value) { (D == 0 ? p.x : p.y) = value; }
};
template <> struct tag<proxy_point> { using type = point_tag; };
template <> struct coordinate_type<proxy_point> { using type = double; };
template <> struct coordinate_system<proxy_point> { using type = cs::cartesian; };
template <> struct dimension<proxy_point> : std::integral_constant<int, 2> {};
template <std::size_t D> struct access<proxy_point, D>
{
    static double get(proxy_point const& p) { return geometry::get<D>(*p.source); }
};
template <typename P, std::size_t N> struct tag<fixed_multipoint<P, N>> { using type = multi_point_tag; };
template <typename I> struct tag<multipoint_view<I>> { using type = multi_point_tag; };
}}}

template <bool Proxy, typename ProxyPoint = proxy_point>
class reading_iterator
    : public boost::iterator_facade<reading_iterator<Proxy, ProxyPoint>, observed_point,
        boost::random_access_traversal_tag, std::conditional_t<Proxy, ProxyPoint, observed_point>>
{
    friend class boost::iterator_core_access;
    using reference_type = std::conditional_t<Proxy, ProxyPoint, observed_point>;
    observed_point const* m_point = nullptr;

    reference_type read(std::true_type) const { return ProxyPoint{m_point}; }
    reference_type read(std::false_type) const { return observed_point(m_point->x, m_point->y); }
    reference_type dereference() const { return read(std::integral_constant<bool, Proxy>()); }
    bool equal(reading_iterator const& other) const { return m_point == other.m_point; }
    void increment() { ++m_point; }
    void decrement() { --m_point; }
    void advance(std::ptrdiff_t count) { m_point += count; }
    std::ptrdiff_t distance_to(reading_iterator const& other) const { return other.m_point - m_point; }

public:
    reading_iterator() = default;
    explicit reading_iterator(observed_point const* p) : m_point(p) {}
};

template <typename Point>
class mutable_reading_iterator
    : public boost::iterator_adaptor<mutable_reading_iterator<Point>, Point const*,
        Point, boost::random_access_traversal_tag, Point const&>
{
    using base_type = boost::iterator_adaptor<mutable_reading_iterator<Point>, Point const*,
        Point, boost::random_access_traversal_tag, Point const&>;
public:
    mutable_reading_iterator() = default;
    explicit mutable_reading_iterator(Point const* p) : base_type(p) {}
    Point const& operator*() { return *this->base(); }
};

using point = bg::model::d2::point_xy<double>;
using less_type = bg::less<void, -1, bg::strategies::relate::cartesian<>>;

class cached_reading_iterator
    : public boost::iterator_adaptor<cached_reading_iterator, point const*,
        point, boost::random_access_traversal_tag, point const&>
{
    using base_type = boost::iterator_adaptor<cached_reading_iterator, point const*,
        point, boost::random_access_traversal_tag, point const&>;
    point m_cache;
public:
    cached_reading_iterator() = default;
    explicit cached_reading_iterator(point const* p) : base_type(p) {}
    point const& operator*() { m_cache = *this->base(); return m_cache; }
};

void reset_counters() { source_copies = source_moves = 0; }
void check_counters()
{
    BOOST_TEST_EQ(source_copies, 0u);
    BOOST_TEST_EQ(source_moves, 0u);
}

template <typename Range>
void check_input(Range const& range, std::initializer_list<xy> expected)
{
    BOOST_TEST_EQ(boost::size(range), expected.size());
    auto it = boost::begin(range);
    for (auto const& p : expected)
    {
        BOOST_TEST_EQ(bg::get<0>(*it), p[0]);
        BOOST_TEST_EQ(bg::get<1>(*it), p[1]);
        ++it;
    }
}

template <typename Range>
void check_sorted_view(Range const& source)
{
    reset_counters();
    less_type const less;
    bg::util::sorted_point_view<Range> const sorted(source, less);
    BOOST_TEST_EQ(sorted.size(), 4u);
    BOOST_TEST(!sorted.empty());
    double x = 0;
    for (auto it = sorted.begin(); it != sorted.end(); ++it, ++x)
    {
        BOOST_TEST_EQ(bg::get<0>(*it), x);
        BOOST_TEST_EQ(bg::get<1>(*it), 0.0);
    }
    BOOST_TEST(std::binary_search(sorted.begin(), sorted.end(), point(2, 0), less));
    BOOST_TEST(!std::binary_search(sorted.begin(), sorted.end(), point(2, 1), less));
    check_input(source, {{{3, 0}}, {{0, 0}}, {{2, 0}}, {{1, 0}}});
    check_counters();
}

template <typename First, typename Second>
void check_operations(First const& first, Second const& second)
{
    bg::model::multi_point<point> output;
    reset_counters();
    BOOST_TEST_EQ(bg::relation(first, second).str(), "0F0FFF0F2");
    check_counters();

    reset_counters();
    bg::intersection(first, second, output);
    check_input(output, {{{2, 0}}, {{0, 0}}});
    BOOST_TEST(bg::is_valid(output));
    check_counters();

    reset_counters();
    output.clear();
    bg::difference(first, second, output);
    check_input(output, {{{3, 0}}, {{1, 0}}});
    BOOST_TEST(bg::is_valid(output));
    check_counters();

    reset_counters();
    output.clear();
    bg::union_(first, second, output);
    check_input(output, {{{3, 0}}, {{0, 0}}, {{2, 0}}, {{1, 0}}, {{4, 0}}});
    BOOST_TEST(bg::is_valid(output));
    check_counters();

    reset_counters();
    output.clear();
    bg::sym_difference(first, second, output);
    check_input(output, {{{3, 0}}, {{1, 0}}, {{4, 0}}});
    BOOST_TEST(bg::is_valid(output));
    check_counters();

    check_input(first, {{{3, 0}}, {{0, 0}}, {{2, 0}}, {{1, 0}}});
    check_input(second, {{{4, 0}}, {{2, 0}}, {{0, 0}}});
}

template <typename CopyPolicy>
void test_native()
{
    using input = input_point<CopyPolicy>;
    using source = fixed_multipoint<input, 4>;
    source const first({{{3, 0}}, {{0, 0}}, {{2, 0}}, {{1, 0}}});
    fixed_multipoint<input, 3> const second({{{4, 0}}, {{2, 0}}, {{0, 0}}});
    using view = bg::util::sorted_point_view<source>;
    static_assert(std::is_same<decltype(*std::declval<view const&>().begin()), input const&>::value,
                  "Lvalue source access stays read-only and non-owning");
    static_assert(!std::is_constructible<view, source&&, less_type const&>::value,
                  "Do not borrow a temporary source");
    check_sorted_view(first);
    check_operations(first, second);
}

template <bool Proxy>
void test_proxy()
{
    fixed_multipoint<observed_point, 4> const first({{{3, 0}}, {{0, 0}}, {{2, 0}}, {{1, 0}}});
    fixed_multipoint<observed_point, 3> const second({{{4, 0}}, {{2, 0}}, {{0, 0}}});
    using iterator = reading_iterator<Proxy>;
    multipoint_view<iterator> const a{iterator(first.begin()), iterator(first.end())};
    multipoint_view<iterator> const b{iterator(second.begin()), iterator(second.end())};
    using sorted = bg::util::sorted_point_view<multipoint_view<iterator>>;
    static_assert(!std::is_reference<decltype(*std::declval<sorted const&>().begin())>::value,
                  "Value/proxy dereferences must not become dangling references");
    check_sorted_view(a);
    check_operations(a, b);
}

template <typename CopyPolicy>
void test_mutable_dereference()
{
    using input = input_point<CopyPolicy>;
    using iterator = mutable_reading_iterator<input>;
    fixed_multipoint<input, 4> const first({{{3, 0}}, {{0, 0}}, {{2, 0}}, {{1, 0}}});
    fixed_multipoint<input, 3> const second({{{4, 0}}, {{2, 0}}, {{0, 0}}});
    multipoint_view<iterator> const a{iterator(first.begin()), iterator(first.end())};
    multipoint_view<iterator> const b{iterator(second.begin()), iterator(second.end())};
    using sorted = bg::util::sorted_point_view<decltype(a)>;
    static_assert(std::is_same<decltype(*std::declval<sorted const&>().begin()), input const&>::value,
                  "Reading a mutable cursor still borrows the native point");
    check_sorted_view(a);
    check_operations(a, b);
}

void test_cached_dereference()
{
    std::array<point, 4> const source{{point(3, 0), point(0, 0), point(2, 0), point(1, 0)}};
    using iterator = cached_reading_iterator;
    multipoint_view<iterator> const range{iterator(source.data()), iterator(source.data() + source.size())};
    bg::util::sorted_point_view<decltype(range)> const sorted(range, less_type());
    auto const& first = *sorted.begin();
    auto const& last = *(sorted.end() - 1);
    BOOST_TEST(&first != &last);
    BOOST_TEST_EQ(bg::get<0>(first), 0.0);
    BOOST_TEST_EQ(bg::get<0>(last), 3.0);
    BOOST_TEST(std::binary_search(sorted.begin(), sorted.end(), point(2, 0), less_type()));
    BOOST_TEST_EQ(bg::get<0>(first), 0.0);
    BOOST_TEST_EQ(bg::get<0>(last), 3.0);
    check_input(range, {{{3, 0}}, {{0, 0}}, {{2, 0}}, {{1, 0}}});
}

void test_convertible_proxy()
{
    fixed_multipoint<observed_point, 4> const first({{{3, 0}}, {{0, 0}}, {{2, 0}}, {{1, 0}}});
    fixed_multipoint<observed_point, 3> const second({{{4, 0}}, {{2, 0}}, {{0, 0}}});
    using iterator = reading_iterator<true, convertible_proxy>;
    multipoint_view<iterator> const range{iterator(second.begin()), iterator(second.end())};
    using view = bg::util::sorted_point_view<decltype(range)>;
    static_assert(std::is_same<bg::tag_t<convertible_proxy>, void>::value, "Unadapted proxy");
    static_assert(std::is_same<decltype(*std::declval<view const&>().begin()), observed_point>::value,
                  "Unadapted proxies are converted per read, never retained by reference");
    reset_counters();
    view const sorted(range, less_type());
    check_input(sorted, {{{0, 0}}, {{2, 0}}, {{4, 0}}});
    BOOST_TEST(std::binary_search(sorted.begin(), sorted.end(), point(2, 0), less_type()));
    BOOST_TEST(!std::binary_search(sorted.begin(), sorted.end(), point(1, 0), less_type()));
    bg::detail::relate::matrix_handler<bg::de9im::matrix> result;
    BOOST_TEST((!bg::detail::relate::multipoint_multipoint<decltype(range), decltype(first)>::search<
        false, bg::strategies::relate::cartesian<>>(range, first, result)));
    BOOST_TEST_EQ(result.result().str(), "0FFFFF0FF");
#ifndef BOOST_NO_CXX17_IF_CONSTEXPR
    // C++14 also instantiates the reversed-input branch, which has never
    // supported unadapted proxies as the iterated (rather than sorted) input.
    bg::model::multi_point<point> output;
    bg::difference(first, range, output);
    check_input(output, {{{3, 0}}, {{1, 0}}});
    BOOST_TEST(bg::is_valid(output));
#endif
    check_input(first, {{{3, 0}}, {{0, 0}}, {{2, 0}}, {{1, 0}}});
    check_input(second, {{{4, 0}}, {{2, 0}}, {{0, 0}}});
}

namespace adl_conversion
{
struct source { int value; };
struct destination { int value = 0; };
struct strategy
{
    unsigned* member_calls;
    unsigned* adl_calls;
    void convert_geometry(source const& from, destination& to) const
    {
        ++*member_calls;
        to.value = from.value;
    }
    friend void convert_geometry(strategy const& s, source const&, destination& to)
    {
        ++*s.adl_calls;
        to.value = -1;
    }
};
}

void test_qualified_conversion()
{
    unsigned member_calls = 0, adl_calls = 0;
    adl_conversion::source const input{42};
    adl_conversion::destination output;
    bg::strategies::detail::assign_point_from_iterator(
        adl_conversion::strategy{&member_calls, &adl_calls},
        mutable_reading_iterator<adl_conversion::source>(&input), output);
    BOOST_TEST_EQ(output.value, 42);
    BOOST_TEST_EQ(member_calls, 1u);
    BOOST_TEST_EQ(adl_calls, 0u);

    input_point<no_copy> const native(3, 4);
    point converted;
    bg::strategies::detail::assign_point_from_iterator(bg::default_strategy(),
        mutable_reading_iterator<input_point<no_copy>>(&native), converted);
    BOOST_TEST_EQ(bg::get<0>(converted), 3.0);
    BOOST_TEST_EQ(bg::get<1>(converted), 4.0);
}

struct bucket_less
{
    std::size_t* calls;
    template <typename A, typename B>
    bool operator()(A const& a, B const& b) const
    {
        ++*calls;
        return int(bg::get<0>(a)) / 2 < int(bg::get<0>(b)) / 2;
    }
};

void test_empty_and_equivalence()
{
    fixed_multipoint<input_point<no_copy>, 0> const empty({});
    less_type const less;
    bg::util::sorted_point_view<decltype(empty)> const sorted(empty, less);
    BOOST_TEST(sorted.empty());
    BOOST_TEST(!std::binary_search(sorted.begin(), sorted.end(), point(0, 0), less));
    fixed_multipoint<input_point<no_copy>, 4> const source({{{3, 0}}, {{0, 0}}, {{2, 0}}, {{1, 0}}});
    std::size_t calls = 0;
    bucket_less const buckets{&calls};
    bg::util::sorted_point_view<decltype(source)> const grouped(source, buckets);
    BOOST_TEST(calls != 0);
    BOOST_TEST(std::binary_search(grouped.begin(), grouped.end(), point(1.75, 99), buckets));
    BOOST_TEST(!std::binary_search(grouped.begin(), grouped.end(), point(4, 0), buckets));
}

int main()
{
    test_native<no_copy>();
    test_native<observed_copy>();
    test_proxy<false>();
    test_proxy<true>();
    test_mutable_dereference<no_copy>();
    test_mutable_dereference<observed_copy>();
    test_cached_dereference();
    test_convertible_proxy();
    test_qualified_conversion();
    test_empty_and_equivalence();
    return boost::report_errors();
}
