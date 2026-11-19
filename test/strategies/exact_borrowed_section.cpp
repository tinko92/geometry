// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.

#define BOOST_GEOMETRY_EXACT_FILTER_STATISTICS
#include <boost/core/lightweight_test.hpp>
#include <boost/iterator/iterator_adaptor.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/geometry/algorithms/detail/sections/sectionalize.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/ring.hpp>
#include <boost/geometry/strategies/relate/cartesian_exact.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <limits>
#include <type_traits>

namespace bg = boost::geometry;
namespace ec = bg::strategies::detail::exact_cartesian;
namespace sd = bg::strategies::detail;

struct copy_counter
{
    inline static std::size_t copies = 0, moves = 0;
    copy_counter() = default;
    copy_counter(copy_counter const&) { ++copies; }
    copy_counter(copy_counter&&) { ++moves; }
    copy_counter& operator=(copy_counter const&) { ++copies; return *this; }
    copy_counter& operator=(copy_counter&&) { ++moves; return *this; }
};

struct no_copy
{
    no_copy() = default;
    no_copy(no_copy const&) = delete;
    no_copy& operator=(no_copy const&) = delete;
};

template <typename CopyPolicy>
struct source_point : CopyPolicy
{
    double coordinates[2];
    mutable unsigned reads[2]{};
    source_point(double x, double y) : coordinates{x, y} {}
};

using observed_point = source_point<copy_counter>;

struct point_proxy
{
    observed_point const* source;
    operator observed_point() const
    {
        return observed_point(source->coordinates[0], source->coordinates[1]);
    }
};

struct convertible_proxy
{
    observed_point const* source;
    operator observed_point() const
    {
        return observed_point(source->coordinates[0], source->coordinates[1]);
    }
};

template <typename Point, bool Closed>
struct source_ring
{
    using value_type = Point;
    using iterator = Point*;
    using const_iterator = Point const*;
    std::array<Point, 7> points{{Point(0, 0), Point(2, 0), Point(2, 0),
                               Point(3, 1), Point(1, 3), Point(-1, 1), Point(0, 0)}};
    iterator begin() { return points.data(); }
    iterator end() { return points.data() + (Closed ? 7 : 6); }
    const_iterator begin() const { return points.data(); }
    const_iterator end() const { return points.data() + (Closed ? 7 : 6); }
};

namespace boost { namespace geometry { namespace traits
{
template <typename P> struct tag<source_point<P>> { using type = point_tag; };
template <typename P> struct coordinate_type<source_point<P>> { using type = double; };
template <typename P> struct coordinate_system<source_point<P>> { using type = cs::cartesian; };
template <typename P> struct dimension<source_point<P>> : std::integral_constant<int, 2> {};
template <typename P, std::size_t D>
struct access<source_point<P>, D>
{
    static double get(source_point<P> const& point)
    {
        ++point.reads[D];
        return point.coordinates[D];
    }
};
template <> struct tag<point_proxy> { using type = point_tag; };
template <> struct coordinate_type<point_proxy> { using type = double; };
template <> struct coordinate_system<point_proxy> { using type = cs::cartesian; };
template <> struct dimension<point_proxy> : std::integral_constant<int, 2> {};
template <std::size_t D> struct access<point_proxy, D>
{
    static double get(point_proxy const& point) { return geometry::get<D>(*point.source); }
};
template <typename P, bool Closed> struct tag<source_ring<P, Closed>> { using type = ring_tag; };
template <typename P, bool Closed>
struct closure<source_ring<P, Closed>>
    : std::integral_constant<closure_selector, Closed ? closed : open> {};
template <typename P, bool Closed>
struct point_order<source_ring<P, Closed>> : std::integral_constant<order_selector, clockwise> {};
}}}

template <bool Proxy, typename ProxyPoint = point_proxy>
class reading_iterator
    : public boost::iterator_facade<reading_iterator<Proxy, ProxyPoint>, observed_point,
        boost::random_access_traversal_tag, std::conditional_t<Proxy, ProxyPoint, observed_point>>
{
    friend class boost::iterator_core_access;
    using reference_type = std::conditional_t<Proxy, ProxyPoint, observed_point>;
    observed_point const* m_point = nullptr;

    reference_type dereference() const
    {
        if constexpr (Proxy) return ProxyPoint{m_point};
        else return observed_point(m_point->coordinates[0], m_point->coordinates[1]);
    }
    bool equal(reading_iterator const& other) const { return m_point == other.m_point; }
    void increment() { ++m_point; }
    void decrement() { --m_point; }
    void advance(std::ptrdiff_t count) { m_point += count; }
    std::ptrdiff_t distance_to(reading_iterator const& other) const { return other.m_point - m_point; }

public:
    reading_iterator() = default;
    explicit reading_iterator(observed_point const* point) : m_point(point) {}
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
    explicit mutable_reading_iterator(Point const* point) : base_type(point) {}
    Point const& operator*() { return *this->base(); }
};

using kernel = bg::strategies::relate::cartesian_exact<double>;
using dimensions = std::integer_sequence<std::size_t, 0, 1>;
using sections = bg::sections<ec::search_box, 2>;

// Preserve the previous owning section-point behavior as the reference.
struct owning_reference : kernel
{
    template <typename Iterator>
    using section_iterator_point_type = ec::exact_point<double>;
    struct section_type : ec::static_section<double>
    {
        template <typename Iterator>
        static void assign_point(Iterator const& iterator, ec::exact_point<double>& point)
        {
            point = ec::make_point<double>(*iterator);
        }
    };
    static section_type section() { return {}; }
};

void compare_sections(sections const& actual, sections const& expected)
{
    BOOST_TEST_EQ(actual.size(), expected.size());
    for (std::size_t i = 0; i < actual.size() && i < expected.size(); ++i)
    {
        auto const& a = actual[i];
        auto const& b = expected[i];
        BOOST_TEST_EQ(a.begin_index, b.begin_index);
        BOOST_TEST_EQ(a.end_index, b.end_index);
        BOOST_TEST_EQ(a.count, b.count);
        BOOST_TEST_EQ(a.range_count, b.range_count);
        BOOST_TEST_EQ(a.non_duplicate_index, b.non_duplicate_index);
        BOOST_TEST_EQ(a.duplicate, b.duplicate);
        BOOST_TEST_EQ(a.is_non_duplicate_first, b.is_non_duplicate_first);
        BOOST_TEST_EQ(a.is_non_duplicate_last, b.is_non_duplicate_last);
        BOOST_TEST(a.ring_id == b.ring_id);
        for (std::size_t d = 0; d != 2; ++d) BOOST_TEST_EQ(a.directions[d], b.directions[d]);
        BOOST_TEST_EQ((bg::get<0, 0>(a.bounding_box)), (bg::get<0, 0>(b.bounding_box)));
        BOOST_TEST_EQ((bg::get<0, 1>(a.bounding_box)), (bg::get<0, 1>(b.bounding_box)));
        BOOST_TEST_EQ((bg::get<1, 0>(a.bounding_box)), (bg::get<1, 0>(b.bounding_box)));
        BOOST_TEST_EQ((bg::get<1, 1>(a.bounding_box)), (bg::get<1, 1>(b.bounding_box)));
    }
}

template <typename Apply>
sections check_borrowed(Apply apply, bool native = true)
{
    sections expected, actual;
    ec::filtering_statistics() = {};
    apply(expected, owning_reference());
    if (native) BOOST_TEST(ec::filtering_statistics().native_point_constructions > 0);
    ec::filtering_statistics() = {};
    copy_counter::copies = copy_counter::moves = 0;
    apply(actual, kernel());
    BOOST_TEST_EQ(copy_counter::copies, 0u);
    BOOST_TEST_EQ(copy_counter::moves, 0u);
    BOOST_TEST_EQ(ec::filtering_statistics().native_point_constructions, 0u);
    compare_sections(actual, expected);
    return actual;
}

template <typename CopyPolicy, bool Closed, bool Reverse>
void test_source_view()
{
    source_ring<source_point<CopyPolicy>, Closed> const ring;
    for (std::size_t max_count : {2u, 10u})
    {
        check_borrowed([&](auto& result, auto const& strategy)
        {
            bg::sectionalize<Reverse, dimensions>(ring, result, strategy, 0, max_count);
        });
    }
}

template <bool Proxy>
void test_value_dereference()
{
    source_ring<observed_point, true> const ring;
    using iterator = reading_iterator<Proxy>;
    using view = ec::iterator_point_view<iterator>;
    using selected = typename sd::section_iterator_point_type<kernel, iterator, ec::search_box>::type;
    static_assert(std::is_same<selected, view>::value, "Native input keeps only its iterator");
    static_assert(sizeof(view) == sizeof(iterator), "No coordinate snapshot");
    iterator begin(ring.begin()), end(ring.end());
    view point;
    sd::assign_point_from_iterator(kernel(), begin, point);
    ++begin;
    BOOST_TEST_EQ(bg::get<0>(point), 0.0);
    BOOST_TEST_EQ(bg::get<1>(point), 0.0);

    using part = bg::detail::sectionalize::sectionalize_part<dimensions>;
    begin = iterator(ring.begin());
    check_borrowed([&](auto& result, auto const& strategy)
    {
        part::apply(result, begin, end, strategy, bg::ring_identifier(), 2);
    });
    auto const reverse_begin = std::make_reverse_iterator(end);
    auto const reverse_end = std::make_reverse_iterator(begin);
    check_borrowed([&](auto& result, auto const& strategy)
    {
        part::apply(result, reverse_begin, reverse_end, strategy, bg::ring_identifier(), 2);
    });
}

template <typename Ring, typename Iterator>
void check_iterator_sections(Ring const& source, Iterator begin, Iterator end)
{
    using view = ec::iterator_point_view<Iterator>;
    static_assert(sizeof(view) == sizeof(Iterator), "No coordinate snapshot");
    view const sample(begin);
    ++begin;
    BOOST_TEST_EQ(bg::get<0>(sample), 0.0);
    BOOST_TEST_EQ(bg::get<1>(sample), 0.0);
    --begin;
    using part = bg::detail::sectionalize::sectionalize_part<dimensions>;
    auto check = [&](auto first, auto last, auto reference_first, auto reference_last)
    {
        sections expected, actual;
        part::apply(expected, reference_first, reference_last, owning_reference(), bg::ring_identifier(), 2);
        ec::filtering_statistics() = {};
        copy_counter::copies = copy_counter::moves = 0;
        part::apply(actual, first, last, kernel(), bg::ring_identifier(), 2);
        BOOST_TEST_EQ(ec::filtering_statistics().native_point_constructions, 0u);
        BOOST_TEST_EQ(copy_counter::copies, 0u);
        BOOST_TEST_EQ(copy_counter::moves, 0u);
        compare_sections(actual, expected);
    };
    check(begin, end, source.begin(), source.end());
    check(std::make_reverse_iterator(end), std::make_reverse_iterator(begin),
          std::make_reverse_iterator(source.end()), std::make_reverse_iterator(source.begin()));
}

template <typename Policy, bool Closed>
void test_mutable_dereference()
{
    source_ring<source_point<Policy>, Closed> const source;
    using iterator = mutable_reading_iterator<source_point<Policy>>;
    check_iterator_sections(source, iterator(source.begin()), iterator(source.end()));
}

void test_convertible_proxy()
{
    source_ring<observed_point, true> const source;
    using iterator = reading_iterator<true, convertible_proxy>;
    static_assert(std::is_same<bg::tag_t<convertible_proxy>, void>::value, "Unadapted proxy");
    using view = ec::iterator_point_view<iterator>;
    static_assert(std::is_same<typename view::reference, observed_point>::value, "Convert only per read");
    check_iterator_sections(source, iterator(source.begin()), iterator(source.end()));
}

void test_on_demand_access()
{
    std::array<source_point<no_copy>, 2> const source{{source_point<no_copy>(3, 4), source_point<no_copy>(5, 6)}};
    using iterator = decltype(source.data());
    ec::iterator_point_view<iterator> point;
    auto it = source.data();
    sd::assign_point_from_iterator(kernel(), it, point);
    auto const copy = point;
    ++it;
    BOOST_TEST_EQ(source[0].reads[0], 0u);
    BOOST_TEST_EQ(source[0].reads[1], 0u);
    BOOST_TEST_EQ(bg::get<0>(copy), 3.0);
    BOOST_TEST_EQ(bg::get<0>(copy), 3.0);
    BOOST_TEST_EQ(source[0].reads[0], 2u);
    BOOST_TEST_EQ(source[1].reads[0], 0u);
}

void test_short_ranges()
{
    source_ring<observed_point, true> const ring;
    using part = bg::detail::sectionalize::sectionalize_part<dimensions>;
    for (std::size_t count : {0u, 1u})
    {
        auto const result = check_borrowed([&](auto& sections, auto const& strategy)
        {
            part::apply(sections, ring.begin(), ring.begin() + count,
                        strategy, bg::ring_identifier(), 2);
        }, count != 0);
        BOOST_TEST(result.empty());
    }
}

void test_implicit_preservation()
{
    ec::point_cache_scope<double> scope;
    using point = ec::exact_point<double>;
    using iterator = bg::model::ring<point>::const_iterator;
    using selected = typename sd::section_iterator_point_type<kernel, iterator, ec::search_box>::type;
    static_assert(std::is_same<selected, point>::value, "Keep implicit point provenance");
    double const delta = std::ldexp(1.0, -54);
    std::array<point, 6> const source{{point(0, -1), point(2, 1), point(0, delta),
                                     point(2, delta), point(0, 0), point(0, 1)}};
    auto const shifted = ec::construct_intersection(source[0], source[1], source[2], source[3]);
    auto const low = point::splice_coordinates(shifted, 0, source[4], 1);
    auto const high = point::splice_coordinates(shifted, 0, source[5], 1);
    bg::model::ring<point> ring{{point(1, 0), point(1, 1), high, low, point(1, 0)}};
    point copied;
    sd::assign_point_from_iterator(kernel(), ring.cbegin() + 2, copied);
    BOOST_TEST(copied.shares_node(high));
    sd::assign_point_from_iterator(kernel(), mutable_reading_iterator<point>(ring.data() + 2), copied);
    BOOST_TEST(copied.shares_node(high));
    auto const result = check_borrowed([&](auto& sections, auto const& strategy)
    {
        bg::sectionalize<false, dimensions>(ring, sections, strategy);
    }, false);
    BOOST_TEST_EQ(result.size(), 4u);
    if (result.size() == 4)
    {
        int const expected[4][2] = {{0, 1}, {1, 0}, {0, -1}, {-1, 0}};
        for (std::size_t i = 0; i != 4; ++i)
        {
            BOOST_TEST(!result[i].duplicate);
            BOOST_TEST_EQ(result[i].directions[0], expected[i][0]);
            BOOST_TEST_EQ(result[i].directions[1], expected[i][1]);
        }
    }
}

struct legacy_selection
{
    template <typename> using section_point_type = bg::model::d2::point_xy<float>;
};

void test_legacy_selection()
{
    using point = bg::model::d2::point_xy<double>;
    using box = bg::model::box<point>;
    using iterator = point const*;
    using selected = typename sd::section_iterator_point_type<legacy_selection, iterator, box>::type;
    using previous = typename sd::section_point_type<legacy_selection, point, box>::type;
    static_assert(std::is_same<selected, previous>::value, "Old strategy point selection is unchanged");
    point const source(3, 4);
    point destination;
    sd::assign_point_from_iterator(bg::default_strategy(), &source, destination);
    BOOST_TEST_EQ(bg::get<0>(destination), 3.0);
    BOOST_TEST_EQ(bg::get<1>(destination), 4.0);
}

template <typename C>
void native_search_limits()
{
    using p = bg::model::d2::point_xy<C>;
    using ring = bg::model::ring<p>;
    using root = ec::borrowed_root<C, bg::strategies::relate::approximate_output, ring, ring>;
    using section = ec::borrowed_section<root>;
    double const x = 0x1p53;
    ec::search_box const limit(ec::search_point(x, -std::numeric_limits<double>::infinity()),
                              ec::search_point(x, std::numeric_limits<double>::infinity()));
    BOOST_TEST(!section::template preceding<0>(1, p(C(x), 0), limit, limit));
    BOOST_TEST(!section::template preceding<1>(1, p(C(x), 0), limit, limit));
    BOOST_TEST(!section::template preceding<1>(-1, p(C(x), 0), limit, limit));
    BOOST_TEST(section::template preceding<0>(1, p(C(x) - 1, 0), limit, limit));
    if constexpr (std::is_integral<C>::value)
        // Rounding x+1 to binary64 here would erase strict order.
        BOOST_TEST(section::template preceding<0>(-1, p(C(x) + 1, 0), limit, limit));
    else
        BOOST_TEST(section::template preceding<0>(-1, p(std::nextafter(x, x * 2), 0), limit, limit));
}

int main()
{
    native_search_limits<double>();
    native_search_limits<std::int64_t>();
    test_on_demand_access();
    test_short_ranges();
    test_legacy_selection();
    test_source_view<no_copy, false, false>();
    test_source_view<no_copy, false, true>();
    test_source_view<no_copy, true, false>();
    test_source_view<no_copy, true, true>();
    test_source_view<copy_counter, false, false>();
    test_source_view<copy_counter, false, true>();
    test_source_view<copy_counter, true, false>();
    test_source_view<copy_counter, true, true>();
    test_value_dereference<false>();
    test_value_dereference<true>();
    test_mutable_dereference<no_copy, false>();
    test_mutable_dereference<no_copy, true>();
    test_mutable_dereference<copy_counter, false>();
    test_mutable_dereference<copy_counter, true>();
    test_convertible_proxy();
    test_implicit_preservation();
    return boost::report_errors();
}
