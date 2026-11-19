// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.

#define BOOST_GEOMETRY_EXACT_FILTER_STATISTICS
#include <boost/core/lightweight_test.hpp>
#include <boost/iterator/iterator_adaptor.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/geometry/algorithms/detail/overlay/assign_parents.hpp>
#include <boost/geometry/algorithms/detail/overlay/ring_properties.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/ring.hpp>
#include <boost/geometry/strategies/cartesian/exact/static_ring_topology.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <type_traits>
#include <vector>

namespace bg = boost::geometry;
namespace ec = bg::strategies::detail::exact_cartesian;
namespace rt = ec::ring_topology_detail;
namespace ov = bg::detail::overlay;
namespace sd = bg::strategies::detail;

struct copied_ring_tag {};
struct copied_ring_source
{
    bg::model::ring<bg::model::d2::point_xy<double>> ring;
};

namespace boost { namespace geometry { namespace detail { namespace overlay
{
template <> struct get_ring<copied_ring_tag>
{
    static auto apply(ring_identifier const&, copied_ring_source const& source)
    {
        return source.ring;
    }
};
}}}}

struct ring_copy_counter
{
    inline static std::size_t copies = 0, moves = 0;
    ring_copy_counter() = default;
    ring_copy_counter(ring_copy_counter const&) { ++copies; }
    ring_copy_counter(ring_copy_counter&&) { ++moves; }
    ring_copy_counter& operator=(ring_copy_counter const&) { ++copies; return *this; }
    ring_copy_counter& operator=(ring_copy_counter&&) { ++moves; return *this; }
};

struct ring_no_copy
{
    ring_no_copy() = default;
    ring_no_copy(ring_no_copy const&) = delete;
    ring_no_copy& operator=(ring_no_copy const&) = delete;
};

template <typename Policy>
struct ring_source_point : Policy
{
    double x, y;
    ring_source_point(double x, double y) : x(x), y(y) {}
};

using observed_point = ring_source_point<ring_copy_counter>;

struct ring_proxy
{
    observed_point const* source;
    operator observed_point() const { return observed_point(source->x, source->y); }
};

struct convertible_proxy
{
    observed_point const* source;
    operator observed_point() const { return observed_point(source->x, source->y); }
};

template <typename Point, bool Closed>
struct spiked_ring
{
    using value_type = Point;
    using iterator = Point*;
    using const_iterator = Point const*;
    // The lexicographic minimum lies on a canceled, unequally subdivided spike.
    std::array<Point, 8> points{{Point(0, 0), Point(-4, 0), Point(-2, 0), Point(0, 0),
                               Point(0, 4), Point(4, 4), Point(4, 0), Point(0, 0)}};
    iterator begin() { return points.data(); }
    iterator end() { return points.data() + (Closed ? 8 : 7); }
    const_iterator begin() const { return points.data(); }
    const_iterator end() const { return points.data() + (Closed ? 8 : 7); }
};

namespace boost { namespace geometry { namespace traits
{
template <typename P> struct tag<ring_source_point<P>> { using type = point_tag; };
template <typename P> struct coordinate_type<ring_source_point<P>> { using type = double; };
template <typename P> struct coordinate_system<ring_source_point<P>> { using type = cs::cartesian; };
template <typename P> struct dimension<ring_source_point<P>> : std::integral_constant<int, 2> {};
template <typename P, std::size_t D> struct access<ring_source_point<P>, D>
{
    static double get(ring_source_point<P> const& point) { return D == 0 ? point.x : point.y; }
    static void set(ring_source_point<P>& point, double value)
    {
        (D == 0 ? point.x : point.y) = value;
    }
};
template <> struct tag<ring_proxy> { using type = point_tag; };
template <> struct coordinate_type<ring_proxy> { using type = double; };
template <> struct coordinate_system<ring_proxy> { using type = cs::cartesian; };
template <> struct dimension<ring_proxy> : std::integral_constant<int, 2> {};
template <std::size_t D> struct access<ring_proxy, D>
{
    static double get(ring_proxy const& point) { return geometry::get<D>(*point.source); }
};
template <typename P, bool Closed> struct tag<spiked_ring<P, Closed>> { using type = ring_tag; };
template <typename P, bool Closed> struct closure<spiked_ring<P, Closed>>
    : std::integral_constant<closure_selector, Closed ? closed : open> {};
template <typename P, bool Closed> struct point_order<spiked_ring<P, Closed>>
    : std::integral_constant<order_selector, clockwise> {};
}}}

template <bool Proxy, typename ProxyPoint = ring_proxy>
class ring_reading_iterator
    : public boost::iterator_facade<ring_reading_iterator<Proxy, ProxyPoint>, observed_point,
        boost::random_access_traversal_tag, std::conditional_t<Proxy, ProxyPoint, observed_point>>
{
    friend class boost::iterator_core_access;
    using reference_type = std::conditional_t<Proxy, ProxyPoint, observed_point>;
    observed_point const* m_point = nullptr;
    reference_type dereference() const
    {
        if constexpr (Proxy) return ProxyPoint{m_point};
        else return observed_point(m_point->x, m_point->y);
    }
    bool equal(ring_reading_iterator const& other) const { return m_point == other.m_point; }
    void increment() { ++m_point; }
    void decrement() { --m_point; }
    void advance(std::ptrdiff_t count) { m_point += count; }
    std::ptrdiff_t distance_to(ring_reading_iterator const& other) const { return other.m_point - m_point; }
public:
    ring_reading_iterator() = default;
    explicit ring_reading_iterator(observed_point const* point) : m_point(point) {}
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

template <typename Iterator, bool Closed = true>
struct ring_range
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
template <typename I, bool Closed> struct tag<ring_range<I, Closed>> { using type = ring_tag; };
template <typename I, bool Closed> struct closure<ring_range<I, Closed>>
    : std::integral_constant<closure_selector, Closed ? closed : open> {};
template <typename I, bool Closed> struct point_order<ring_range<I, Closed>>
    : std::integral_constant<order_selector, clockwise> {};
}}}

template <typename C>
struct topology_strategy
{
    static auto ring_topology() { return ec::ring_topology<C>(); }
};

void reset_counters()
{
    ec::filtering_statistics() = {};
    ring_copy_counter::copies = ring_copy_counter::moves = 0;
}

void check_counters()
{
    BOOST_TEST_EQ(ec::filtering_statistics().native_point_constructions, 0u);
    BOOST_TEST_EQ(ring_copy_counter::copies, 0u);
    BOOST_TEST_EQ(ring_copy_counter::moves, 0u);
}

template <typename Ring>
void check_borrowed_metadata(Ring const& ring, bool clockwise = true)
{
    using source_point = typename boost::range_value<Ring>::type;
    using sample = ec::iterator_point_view<typename boost::range_iterator<Ring const>::type>;
    using properties = ov::ring_properties<source_point, sd::topological_ring_tag<sample>>;
    reset_counters();
    properties const metadata(ring, topology_strategy<double>());
    BOOST_TEST(metadata.valid);
    BOOST_TEST_EQ(metadata.clockwise, clockwise);
    sample const first(boost::begin(ring));
    BOOST_TEST_EQ(bg::get<0>(metadata.point), bg::get<0>(first));
    BOOST_TEST_EQ(bg::get<1>(metadata.point), bg::get<1>(first));
    // analyze's boundary and residual have died. These samples must still reach
    // the original input, not either temporary vector (also checked under ASan).
    std::vector<std::array<double, 32>> churn(128);
    BOOST_TEST_EQ(bg::get<0>(metadata.boundary_point), 0.0);
    BOOST_TEST_EQ(bg::get<1>(metadata.boundary_point), 0.0);
    BOOST_TEST(ec::ring_topology<double>().analyze(ring).valid);
    check_counters();
}

template <typename Policy, bool Closed>
void test_native_ring()
{
    spiked_ring<ring_source_point<Policy>, Closed> const ring;
    check_borrowed_metadata(ring);
    using iterator = typename decltype(ring)::const_iterator;
    ec::iterator_point_view<iterator> representative;
    reset_counters();
    ec::ring_topology<double> topology;
    BOOST_TEST(topology.analyze(ring, representative).valid);
    BOOST_TEST(representative.base() == ring.begin() || representative.base() == ring.begin() + 3);
    BOOST_TEST(topology.representative(ring, representative));
    BOOST_TEST(representative.base() == ring.begin());
    check_counters();
}

template <bool Proxy, typename ProxyPoint = ring_proxy>
void test_proxy_ring()
{
    spiked_ring<observed_point, true> const source;
    using iterator = ring_reading_iterator<Proxy, ProxyPoint>;
    ring_range<iterator> const forward{iterator(source.begin()), iterator(source.end())};
    check_borrowed_metadata(forward);
    using reverse = std::reverse_iterator<iterator>;
    ring_range<reverse> const reversed{reverse(forward.end()), reverse(forward.begin())};
    check_borrowed_metadata(reversed, false);
}

template <typename Ring>
void check_borrowed_boundary(Ring const& ring, bool clockwise = true)
{
    using sample = ec::iterator_point_view<typename boost::range_iterator<Ring const>::type>;
    sample representative;
    reset_counters();
    {
        auto const boundary = rt::make_boundary<double>(ring);
        auto const result = rt::analyze(boundary, representative);
        BOOST_TEST(result.valid);
        BOOST_TEST_EQ(result.clockwise, clockwise);
        bg::model::d2::point_xy<double> const inside(2, 2), outside(-2, 2);
        BOOST_TEST_EQ(rt::locate(inside, boundary), 1);
        BOOST_TEST_EQ(rt::locate(outside, boundary), -1);
    }
    std::vector<std::array<double, 32>> churn(128);
    BOOST_TEST_EQ(bg::get<0>(representative), 0.0);
    BOOST_TEST_EQ(bg::get<1>(representative), 0.0);
    BOOST_TEST(ec::ring_topology<double>().representative(ring, representative));
    check_counters();
    check_borrowed_metadata(ring, clockwise);
}

template <typename Ring>
void check_borrowed_fast_analysis(Ring const& ring, bool clockwise = true)
{
    using sample = ec::iterator_point_view<typename boost::range_iterator<Ring const>::type>;
    sample representative;
    sd::ring_topology_result result;
    reset_counters();
    BOOST_TEST(rt::fast_analyze<double>(ring, representative, result, bg::ring_tag()));
    BOOST_TEST(result.valid);
    BOOST_TEST_EQ(result.clockwise, clockwise);
    BOOST_TEST_EQ(bg::get<0>(representative), 0.0);
    BOOST_TEST_EQ(bg::get<1>(representative), 0.0);
    check_counters();
}

void test_convertible_proxy()
{
    spiked_ring<observed_point, true> const source;
    using iterator = ring_reading_iterator<true, convertible_proxy>;
    ring_range<iterator> const forward{iterator(source.begin()), iterator(source.end())};
    check_borrowed_boundary(forward);
    using reverse = std::reverse_iterator<iterator>;
    ring_range<reverse> const reversed{reverse(forward.end()), reverse(forward.begin())};
    check_borrowed_boundary(reversed, false);
    ring_range<iterator> const simple{iterator(source.begin() + 3), iterator(source.end())};
    check_borrowed_fast_analysis(simple);
    ring_range<reverse> const simple_reversed{reverse(simple.end()), reverse(simple.begin())};
    check_borrowed_fast_analysis(simple_reversed, false);
}

template <typename Policy, bool Closed>
void test_mutable_ring()
{
    using point = ring_source_point<Policy>;
    spiked_ring<point, Closed> const source;
    using iterator = mutable_reading_iterator<point>;
    ring_range<iterator, Closed> const forward{iterator(source.begin()), iterator(source.end())};
    check_borrowed_boundary(forward);
    using reverse = std::reverse_iterator<iterator>;
    ring_range<reverse, Closed> const reversed{reverse(forward.end()), reverse(forward.begin())};
    check_borrowed_boundary(reversed, false);
    ring_range<iterator, Closed> const simple{iterator(source.begin() + 3), iterator(source.end())};
    check_borrowed_fast_analysis(simple);
    ring_range<reverse, Closed> const simple_reversed{reverse(simple.end()), reverse(simple.begin())};
    check_borrowed_fast_analysis(simple_reversed, false);
}

template <typename C>
void test_parenting_and_midpoints()
{
    ec::point_cache_scope<C> scope;
    using point = bg::model::d2::point_xy<C>;
    using ring = bg::model::ring<point>;
    using sample = ec::iterator_point_view<typename ring::const_iterator>;
    using properties = ov::ring_properties<point, sd::topological_ring_tag<sample>>;
    std::vector<ring> const roots{
        ring{{0, 0}, {0, 4}, {4, 4}, {4, 0}, {0, 0}},
        ring{{0, 0}, {4, 0}, {0, 4}, {0, 0}}
    };
    bg::model::polygon<point> const empty;
    std::map<bg::ring_identifier, properties> metadata;
    topology_strategy<C> strategy;
    reset_counters();
    for (std::size_t i = 0; i != roots.size(); ++i)
        metadata.try_emplace(bg::ring_identifier(2, i, -1), roots[i], strategy);
    // Every triangle vertex lies on the square boundary. Its diagonal midpoint
    // is strictly inside, so parent classification must test a true midpoint.
    ov::assign_parents<bg::overlay_union>(empty, empty, roots, metadata, strategy);
    BOOST_TEST(metadata.at(bg::ring_identifier(2, 1, -1)).parent == bg::ring_identifier(2, 0, -1));
    BOOST_TEST(metadata.at(bg::ring_identifier(2, 0, -1)).parent.source_index == -1);
    BOOST_TEST(!metadata.at(bg::ring_identifier(2, 1, -1)).discarded);
    auto const* owner = ec::operation_detail::current<C>();
    BOOST_TEST(owner != nullptr && !owner->storage().auxiliary.empty());
    auto const inner = rt::make_boundary<C>(roots[1]);
    auto const middle = rt::midpoint<C>(inner.points[1], inner.points[2]);
    BOOST_TEST_EQ((ec::static_compare_coordinate<0, C>(middle, point(2, 2))), 0);
    BOOST_TEST_EQ((ec::static_compare_coordinate<1, C>(middle, point(2, 2))), 0);
    check_counters();
}

void test_polygon_and_cancellation()
{
    ec::point_cache_scope<double> scope;
    using point = bg::model::d2::point_xy<double>;
    bg::model::polygon<point> polygon;
    polygon.outer() = {{0, 0}, {0, 8}, {8, 8}, {8, 0}, {0, 0}};
    polygon.inners().push_back({{2, 2}, {6, 2}, {6, 6}, {2, 6}, {2, 2}});
    using iterator = decltype(polygon.outer().cbegin());
    ec::iterator_point_view<iterator> sample;
    reset_counters();
    BOOST_TEST(ec::ring_topology<double>().analyze(polygon, sample).valid);
    BOOST_TEST(sample.base() == polygon.outer().cbegin());
    bg::model::ring<point> const canceled{{0, 0}, {0, 4}, {4, 4}, {4, 0}, {0, 0},
        {2, 0}, {4, 0}, {4, 2}, {4, 4}, {2, 4}, {0, 4}, {0, 2}, {0, 0}};
    BOOST_TEST(!ec::ring_topology<double>().analyze(canceled).valid);
    check_counters();
}

void test_implicit_points()
{
    ec::point_cache_scope<double> scope;
    using point = ec::exact_point<double>;
    double const delta = std::ldexp(1.0, -54);
    std::array<point, 6> const source{{point(0, -1), point(2, 1), point(0, delta),
                                     point(2, delta), point(0, 0), point(0, 1)}};
    auto const shifted = ec::construct_intersection(source[0], source[1], source[2], source[3]);
    auto const low = point::splice_coordinates(shifted, 0, source[4], 1);
    auto const high = point::splice_coordinates(shifted, 0, source[5], 1);
    bg::model::ring<point> const ring{{point(1, 0), point(1, 1), high, low, point(1, 0)}};
    auto const boundary = rt::make_boundary<double>(ring);
    static_assert(std::is_same<typename decltype(boundary)::point_type, point>::value,
                  "Implicit point representation is preserved");
    BOOST_TEST(boundary.points[2].shares_node(high));
    reset_counters();
    using properties = ov::ring_properties<point, ec::ring_metadata<double>>;
    properties const metadata(ring, topology_strategy<double>());
    BOOST_TEST(metadata.valid && metadata.clockwise);
    BOOST_TEST(metadata.point.shares_node(ring.front()));
    check_counters();
}

void test_internal_rational_points()
{
    ec::point_cache_scope<double> scope;
    using rational = ec::exact_number;
    using point = bg::model::d2::point_xy<rational>;
    using ring = bg::model::ring<point>;
    using properties = ov::ring_properties<point, ec::ring_metadata<double>>;
    using working = typename rt::boundary_point_type<double, ring>::type;
    static_assert(std::is_same<working, ec::exact_point<double>>::value,
                  "Internal rational adapters use the selected native kernel's bounded literals");
    rational const third = rational(1) / 3;
    rational const sixth = rational(1) / 6;
    // The canceled spike forces residual sorting, which must never infer a
    // kernel with cpp_rational as its native coordinate type.
    std::vector<ring> const roots{
        ring{{0, 0}, {-third, 0}, {-sixth, 0}, {0, 0},
             {0, third}, {third, third}, {third, 0}, {0, 0}},
        ring{{0, 0}, {third, 0}, {0, third}, {0, 0}}
    };
    bg::model::polygon<point> const empty;
    topology_strategy<double> strategy;
    std::map<bg::ring_identifier, properties> metadata;
    reset_counters();
    for (std::size_t i = 0; i != roots.size(); ++i)
        metadata.try_emplace(bg::ring_identifier(2, i, -1), roots[i], strategy);
    BOOST_TEST(metadata.at(bg::ring_identifier(2, 0, -1)).valid);
    BOOST_TEST(metadata.at(bg::ring_identifier(2, 0, -1)).clockwise);
    ov::assign_parents<bg::overlay_union>(empty, empty, roots, metadata, strategy);
    BOOST_TEST(metadata.at(bg::ring_identifier(2, 1, -1)).parent == bg::ring_identifier(2, 0, -1));
    auto const boundary = rt::make_boundary<double>(roots[0]);
    BOOST_TEST_EQ(rt::locate(point(sixth, sixth), boundary), 1);
    check_counters();
}

void test_default_metadata_barrier()
{
    using point = bg::model::d2::point_xy<double>;
    using properties = ov::ring_properties<point, ec::ring_metadata<double>>;
    bg::model::ring<point> const ring{{0, 0}, {0, 4}, {4, 4}, {4, 0}, {0, 0}};
    static_assert(std::is_same<typename properties::point_type, ec::exact_point<double>>::value,
                  "The stock metadata still requires owning exact points");
    reset_counters();
    properties const metadata(ring, topology_strategy<double>());
    BOOST_TEST(metadata.valid);
    // Deliberately expose, rather than hide via homogeneous snapshots, the two
    // conversions mandated by the existing default metadata selection.
    BOOST_TEST_EQ(ec::filtering_statistics().native_point_constructions, 2u);
}

void test_mixed_parenting_without_rewrapping()
{
    ec::point_cache_scope<double> scope;
    using point = bg::model::d2::point_xy<double>;
    using implicit = ec::exact_point<double>;
    using properties = ov::ring_properties<point, ec::ring_metadata<double>>;
    bg::model::polygon<point> outer;
    outer.outer() = {{0, 0}, {0, 10}, {10, 10}, {10, 0}, {0, 0}};
    bg::model::polygon<observed_point> subtractor;
    for (auto const& xy : std::array<std::array<double, 2>, 5>{{{{2, 2}}, {{2, 8}}, {{8, 8}}, {{8, 2}}, {{2, 2}}}})
        subtractor.outer().emplace_back(xy[0], xy[1]);
    std::vector<bg::model::ring<implicit>> const generated{
        bg::model::ring<implicit>{{3, 3}, {3, 7}, {7, 7}, {7, 3}, {3, 3}}
    };
    topology_strategy<double> strategy;
    std::map<bg::ring_identifier, properties> metadata;
    metadata.try_emplace(bg::ring_identifier(0, -1, -1), outer.outer(), strategy);
    metadata.try_emplace(bg::ring_identifier(1, -1, -1), subtractor.outer(), strategy);
    metadata.try_emplace(bg::ring_identifier(2, 0, -1), generated[0], strategy);
    metadata.at(bg::ring_identifier(1, -1, -1)).reversed = true;
    reset_counters();
    ov::assign_parents<bg::overlay_union>(outer, subtractor, generated, metadata, strategy);
    BOOST_TEST(metadata.at(bg::ring_identifier(1, -1, -1)).parent == bg::ring_identifier(0, -1, -1));
    BOOST_TEST(metadata.at(bg::ring_identifier(0, -1, -1)).parent.source_index == -1);
    BOOST_TEST(metadata.at(bg::ring_identifier(2, 0, -1)).parent.source_index == -1);
    check_counters();
}

void test_temporary_ring_accessor()
{
    using point = bg::model::d2::point_xy<double>;
    copied_ring_source const source{{{0, 0}, {0, 4}, {4, 4}, {4, 0}, {0, 0}}};
    auto const selected = rt::selected_source<double, copied_ring_tag>(bg::ring_identifier(), source);
    static_assert(std::is_same<typename decltype(selected)::point_type, ec::exact_point<double>>::value,
                  "Do not borrow an accessor's temporary owning ring");
    std::vector<std::array<double, 32>> churn(128);
    BOOST_TEST_EQ(rt::locate(point(2, 2), selected), 1);
    BOOST_TEST_EQ(rt::locate(point(0, 2), selected), 0);
    BOOST_TEST_EQ(rt::locate(point(-1, 2), selected), -1);
}

int main()
{
    test_native_ring<ring_no_copy, false>();
    test_native_ring<ring_no_copy, true>();
    test_native_ring<ring_copy_counter, false>();
    test_native_ring<ring_copy_counter, true>();
    test_proxy_ring<false>();
    test_proxy_ring<true>();
    test_convertible_proxy();
    test_mutable_ring<ring_no_copy, false>();
    test_mutable_ring<ring_no_copy, true>();
    test_mutable_ring<ring_copy_counter, false>();
    test_mutable_ring<ring_copy_counter, true>();
    test_parenting_and_midpoints<double>();
    test_parenting_and_midpoints<std::int_least64_t>();
    test_polygon_and_cancellation();
    test_implicit_points();
    test_internal_rational_points();
    test_default_metadata_barrier();
    test_mixed_parenting_without_rewrapping();
    test_temporary_ring_accessor();
    return boost::report_errors();
}
