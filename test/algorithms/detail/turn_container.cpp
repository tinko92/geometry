// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.

#include <algorithm>
#include <deque>
#include <functional>
#include <list>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include <geometry_test_common.hpp>
#include <boost/iterator/indirect_iterator.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/algorithms/detail/overlay/get_intersection_points.hpp>
#include <boost/geometry/algorithms/detail/turns/filter_continue_turns.hpp>
#include <boost/geometry/algorithms/detail/turns/remove_duplicate_turns.hpp>
#include <boost/geometry/strategies/detail/turn_container.hpp>

namespace sd = boost::geometry::strategies::detail;
namespace ov = boost::geometry::detail::overlay;

struct hook_state
{
    int factories = 0, sorts = 0, range_sorts = 0, erases = 0, uniques = 0;
    int depth = 0, scopes = 0, emissions = 0;
    bool reversals[2][2] = {};
    bg::segment_identifier ids[2];
    void const* geometries[2] = {};
};

// Only the active pointer sequence may move; records retain their identity.
template <typename Turn>
class mock_turn_container
{
    using pointers = std::vector<Turn*>;
    std::deque<Turn> m_records;
    pointers m_active;
    hook_state* m_state;

public:
    using value_type = Turn;
    using size_type = std::size_t;
    using iterator = boost::indirect_iterator<typename pointers::iterator>;
    using const_iterator = boost::indirect_iterator<typename pointers::const_iterator, Turn const>;

    explicit mock_turn_container(hook_state& state) : m_state(&state) {}
    mock_turn_container(mock_turn_container&&) = default;
    mock_turn_container(mock_turn_container const&) = delete;

    iterator begin() { return iterator(m_active.begin()); }
    iterator end() { return iterator(m_active.end()); }
    const_iterator begin() const { return const_iterator(m_active.begin()); }
    const_iterator end() const { return const_iterator(m_active.end()); }
    std::size_t size() const { return m_active.size(); }
    bool empty() const { return m_active.empty(); }
    Turn& operator[](std::size_t i) { return *m_active[i]; }
    Turn const& operator[](std::size_t i) const { return *m_active[i]; }
    void clear() { m_active.clear(); }

    template <typename... Args>
    Turn& emplace(Args&&... args)
    {
        m_records.emplace_back(std::forward<Args>(args)...);
        m_active.push_back(std::addressof(m_records.back()));
        return m_records.back();
    }

    void push_back(Turn const& turn)
    {
        BOOST_CHECK_EQUAL(m_state->depth, 1);
        BOOST_CHECK(turn.operations[0].seg_id == m_state->ids[0]);
        BOOST_CHECK(turn.operations[1].seg_id == m_state->ids[1]);
        ++m_state->emissions;
        emplace(turn);
    }

    template <typename Less>
    void sort(Less less)
    {
        ++m_state->sorts;
        sort(begin(), end(), less);
    }

    template <typename Less>
    void sort(iterator first, iterator last, Less less)
    {
        ++m_state->range_sorts;
        std::sort(first.base(), last.base(), [&](Turn const* a, Turn const* b)
        {
            return less(*a, *b);
        });
    }

    template <typename Predicate>
    void erase_if(Predicate predicate)
    {
        ++m_state->erases;
        m_active.erase(std::remove_if(m_active.begin(), m_active.end(),
            [&](Turn const* turn) { return predicate(*turn); }), m_active.end());
    }

    template <typename Predicate>
    void unique(Predicate predicate)
    {
        ++m_state->uniques;
        m_active.erase(std::unique(m_active.begin(), m_active.end(),
            [&](Turn const* a, Turn const* b) { return predicate(*a, *b); }), m_active.end());
    }
};

class observed_scope
{
    hook_state* m_state;
    bg::segment_identifier m_previous[2];

public:
    template <bool Reverse1, bool Reverse2, typename Turn>
    observed_scope(hook_state& state, Turn const& turn,
                   std::integral_constant<bool, Reverse1>, std::integral_constant<bool, Reverse2>)
        : m_state(&state), m_previous{state.ids[0], state.ids[1]}
    {
        BOOST_CHECK_EQUAL(state.depth, 0);
        ++state.depth;
        ++state.scopes;
        state.reversals[Reverse1][Reverse2] = true;
        for (int i = 0; i < 2; ++i)
        {
            BOOST_CHECK(turn.operations[i].seg_id.segment_index >= 0);
            state.ids[i] = turn.operations[i].seg_id;
        }
    }

    observed_scope(observed_scope const&) = delete;
    observed_scope(observed_scope&& other) : m_state(other.m_state),
        m_previous{other.m_previous[0], other.m_previous[1]}
    {
        other.m_state = nullptr;
    }

    ~observed_scope()
    {
        if (m_state != nullptr)
        {
            --m_state->depth;
            m_state->ids[0] = m_previous[0];
            m_state->ids[1] = m_previous[1];
        }
    }
};

struct observed_strategy : bg::strategies::relate::cartesian<>
{
    hook_state& state;
    explicit observed_strategy(hook_state& state) : state(state) {}

    template <typename DefaultContainer, typename Geometry1, typename Geometry2>
    auto make_turn_container(Geometry1 const& g1, Geometry2 const& g2) const
    {
        ++state.factories;
        if (state.geometries[0] != nullptr)
        {
            BOOST_CHECK(state.geometries[0] == std::addressof(g1));
            BOOST_CHECK(state.geometries[1] == std::addressof(g2));
        }
        return mock_turn_container<typename DefaultContainer::value_type>(state);
    }

    template <bool Reverse1, bool Reverse2, typename Turn>
    observed_scope model_turn_scope(Turn const& turn) const
    {
        return observed_scope(state, turn, std::integral_constant<bool, Reverse1>(),
                                           std::integral_constant<bool, Reverse2>());
    }
};

struct immovable_turn
{
    int key;
    explicit immovable_turn(int key) : key(key) {}
    immovable_turn(immovable_turn const&) = delete;
    immovable_turn(immovable_turn&&) = delete;
    immovable_turn& operator=(immovable_turn const&) = delete;
    immovable_turn& operator=(immovable_turn&&) = delete;
};

void test_services()
{
    struct ordinary_strategy {};
    ordinary_strategy strategy;
    static_assert(std::is_same<decltype(sd::make_turn_container<std::vector<int>>(strategy, 1, 2)),
                               std::vector<int>>::value, "Preserve vector defaults");
    static_assert(std::is_same<decltype(sd::make_turn_container<std::deque<int>>(strategy, 1, 2)),
                               std::deque<int>>::value, "Preserve deque defaults");
    auto values = sd::make_turn_container<std::vector<int>>(strategy, 1, 2);
    BOOST_CHECK(values.empty());
    values = {4, 2, 3, 2, 1};
    sd::sort_turns(values, std::less<int>());
    sd::unique_turns(values, std::equal_to<int>());
    sd::erase_turns_if(values, [](int value) { return value % 2 == 0; });
    BOOST_CHECK(values == std::vector<int>({1, 3}));
    sd::sort_turns(values, values.begin(), values.end(), std::greater<int>());
    BOOST_CHECK(values == std::vector<int>({3, 1}));
    auto empty_scope = sd::model_turn_scope(strategy, values);
    boost::ignore_unused(empty_scope);

    std::list<int> linked{3, 2, 2, 1};
    sd::sort_turns(linked, std::less<int>());
    sd::unique_turns(linked, std::equal_to<int>());
    BOOST_CHECK(linked == std::list<int>({1, 2, 3}));

    hook_state state;
    mock_turn_container<immovable_turn> turns(state);
    auto const* first = std::addressof(turns.emplace(3));
    auto const* second = std::addressof(turns.emplace(1));
    auto const* duplicate = std::addressof(turns.emplace(1));
    turns.emplace(2);
    auto less = [](immovable_turn const& a, immovable_turn const& b) { return a.key < b.key; };
    sd::sort_turns(turns, less);
    sd::sort_turns(turns, turns.begin() + 2, turns.end(),
        [&](immovable_turn const& a, immovable_turn const& b) { return less(b, a); });
    BOOST_CHECK_EQUAL(turns[2].key, 3);
    sd::unique_turns(turns,
        [](immovable_turn const& a, immovable_turn const& b) { return a.key == b.key; });
    sd::erase_turns_if(turns, [](immovable_turn const& turn) { return turn.key == 3; });
    BOOST_CHECK_EQUAL(turns.size(), 2u);
    BOOST_CHECK_EQUAL(first->key, 3);
    BOOST_CHECK_EQUAL(second->key, 1);
    BOOST_CHECK_EQUAL(duplicate->key, 1);
    BOOST_CHECK_EQUAL(state.sorts, 1);
    BOOST_CHECK_EQUAL(state.range_sorts, 2);
    BOOST_CHECK_EQUAL(state.erases, 1);
    BOOST_CHECK_EQUAL(state.uniques, 1);
}

using point = bg::model::d2::point_xy<double>;
using line = bg::model::linestring<point>;
using polygon = bg::model::polygon<point>;
using turn = ov::turn_info<point>;

template <typename Geometry>
Geometry from_wkt(char const* wkt)
{
    Geometry result;
    bg::read_wkt(wkt, result);
    return result;
}

template <typename Geometry>
std::string as_wkt(Geometry const& geometry)
{
    std::ostringstream output;
    output << bg::wkt(geometry);
    return output.str();
}

void test_filters()
{
    hook_state state;
    mock_turn_container<turn> turns(state);
    auto& continued = turns.emplace();
    continued.point = point(0, 0);
    continued.method = ov::method_collinear;
    continued.operations[0].operation = ov::operation_continue;
    continued.operations[1].operation = ov::operation_opposite;
    auto& retained = turns.emplace();
    retained.point = point(1, 1);
    turns.emplace(retained);
    bg::detail::turns::filter_continue_turns<decltype(turns), true>::apply(turns);
    bg::detail::turns::remove_duplicate_turns<decltype(turns), true>
        ::apply(turns, bg::strategies::relate::cartesian<>());
    BOOST_CHECK_EQUAL(turns.size(), 1u);
    BOOST_CHECK(std::addressof(turns[0]) == std::addressof(retained));
    BOOST_CHECK_EQUAL(continued.method, ov::method_collinear);
    BOOST_CHECK_EQUAL(state.erases, 1);
    BOOST_CHECK_EQUAL(state.uniques, 1);
}

template <typename Geometry1, typename Geometry2>
void test_relation(Geometry1 const& g1, Geometry2 const& g2, bool grouped_sort = false)
{
    hook_state state;
    state.geometries[0] = std::addressof(g1);
    state.geometries[1] = std::addressof(g2);
    observed_strategy strategy(state);
    BOOST_CHECK_EQUAL(bg::relation(g1, g2, strategy).str(), bg::relation(g1, g2).str());
    BOOST_CHECK_EQUAL(state.factories, 1);
    BOOST_CHECK(state.sorts > 0);
    BOOST_CHECK(state.scopes > 0);
    BOOST_CHECK(state.emissions > 0);
    if (grouped_sort)
    {
        BOOST_CHECK(state.range_sorts > state.sorts);
    }
    BOOST_CHECK_EQUAL(state.depth, 0);
}

void test_algorithms()
{
    auto const l1 = from_wkt<line>("LINESTRING(-1 2,5 2,2 2,2 5)");
    auto const l2 = from_wkt<line>("LINESTRING(0 2,4 2,2 0,2 4)");
    auto const p1 = from_wkt<polygon>("POLYGON((0 0,0 4,4 4,4 0,0 0))");
    auto const p2 = from_wkt<polygon>("POLYGON((2 -1,2 3,6 3,6 -1,2 -1))");
    auto const multi = from_wkt<bg::model::multi_polygon<polygon>>
        ("MULTIPOLYGON(((0 0,0 4,2 4,2 0,0 0)),((2 0,2 4,4 4,4 0,2 0)))");
    test_relation(l1, l2);
    test_relation(l1, p1);
    test_relation(l1, multi, true);
    test_relation(from_wkt<line>("LINESTRING(0 0,0 4,4 4,4 0,0 0)"), p1, true);
    test_relation(p1, p2);

    hook_state state;
    observed_strategy strategy(state);
    bg::model::multi_linestring<line> output, expected;
    bg::intersection(l1, l2, output, strategy);
    bg::intersection(l1, l2, expected);
    BOOST_CHECK_EQUAL(as_wkt(output), as_wkt(expected));
    output.clear(); expected.clear();
    bg::intersection(l1, p1, output, strategy);
    bg::intersection(l1, p1, expected);
    BOOST_CHECK_EQUAL(as_wkt(output), as_wkt(expected));

    std::vector<polygon> areas;
    bg::difference(p1, p2, areas, strategy);
    BOOST_CHECK(state.reversals[0][1]);
    BOOST_CHECK_EQUAL(areas.size(), 1u);
    BOOST_CHECK_CLOSE(bg::area(areas.front()), 10.0, 0.0001);
    BOOST_CHECK_EQUAL(bg::disjoint(l1, l2, strategy), bg::disjoint(l1, l2));
    BOOST_CHECK_EQUAL(bg::touches(p1, p2, strategy), bg::touches(p1, p2));
    BOOST_CHECK_EQUAL(bg::is_simple(l1, strategy), bg::is_simple(l1));
    BOOST_CHECK_EQUAL(bg::is_valid(p1, strategy), bg::is_valid(p1));
    BOOST_CHECK_EQUAL(state.depth, 0);
}

struct throwing_turn_policy
{
    template <typename Range1, typename Range2, typename Turn, typename Strategy, typename Output>
    static Output apply(Range1 const&, Range2 const&, Turn const&, Strategy const& strategy, Output)
    {
        BOOST_CHECK_EQUAL(strategy.state.depth, 1);
        throw std::runtime_error("turn policy failure");
    }
};

struct checking_interrupt_policy
{
    static bool const enabled = true;
    hook_state& state;
    int calls = 0;

    template <typename Range>
    bool apply(Range const&)
    {
        BOOST_CHECK_EQUAL(state.depth, 0);
        ++calls;
        return true;
    }
};

void test_model_scopes()
{
    hook_state state;
    observed_strategy strategy(state);
    auto const l1 = from_wkt<line>("LINESTRING(0 0,2 2)");
    auto const l2 = from_wkt<line>("LINESTRING(0 2,2 0)");
    auto turns = sd::make_turn_container<std::deque<turn>>(strategy, l1, l2);
    bg::get_intersection_points(l1, l2, turns, strategy);
    BOOST_CHECK_EQUAL(turns.size(), 1u);
    BOOST_CHECK_EQUAL(turns[0].operations[0].seg_id.source_index, 0);
    BOOST_CHECK_EQUAL(turns[0].operations[1].seg_id.source_index, 1);
    BOOST_CHECK_EQUAL(bg::get<0>(turns[0].point), 1.0);
    BOOST_CHECK_EQUAL(state.depth, 0);

    checking_interrupt_policy interrupt{state};
    bg::dispatch::get_turns<bg::linestring_tag, bg::linestring_tag,
        line, line, false, false, ov::get_turn_info<ov::assign_null_policy>>
        ::apply(0, l1, 1, l2, strategy, turns, interrupt);
    BOOST_CHECK_EQUAL(interrupt.calls, 1);
    BOOST_CHECK_EQUAL(state.depth, 0);

    auto const self = from_wkt<bg::model::ring<point>>("POLYGON((0 0,2 2,0 2,2 0,0 0))");
    bg::detail::get_turns::no_interrupt_policy policy;
    turns.clear();
    bg::detail::self_get_turn_points::get_turns<true, ov::get_turn_info<ov::assign_null_policy>>
        ::apply(self, strategy, turns, policy, 0, true);
    BOOST_CHECK(!turns.empty());
    BOOST_CHECK(state.reversals[1][1]);
    BOOST_CHECK_EQUAL(state.depth, 0);

    auto const ring = from_wkt<bg::model::ring<point>>("POLYGON((0 0,0 4,4 4,4 0,0 0))");
    bg::model::box<point> const box(point(2, -1), point(5, 2));
    turns.clear();
    bg::dispatch::get_turns<bg::ring_tag, bg::box_tag, decltype(ring), decltype(box),
        true, false, ov::get_turn_info<ov::assign_null_policy>>
        ::apply(0, ring, 1, box, strategy, turns, policy);
    BOOST_CHECK(!turns.empty());
    BOOST_CHECK(state.reversals[1][0]);
    BOOST_CHECK_EQUAL(state.depth, 0);

    auto const other = from_wkt<bg::model::ring<point>>("POLYGON((2 -1,2 3,6 3,6 -1,2 -1))");
    BOOST_CHECK_THROW((bg::dispatch::get_turns<bg::ring_tag, bg::ring_tag,
        bg::model::ring<point>, bg::model::ring<point>, false, true, throwing_turn_policy>
        ::apply(0, ring, 1, other, strategy, turns, policy)), std::runtime_error);
    BOOST_CHECK_EQUAL(state.depth, 0);
    BOOST_CHECK_EQUAL(state.ids[0].source_index, -1);
    BOOST_CHECK_EQUAL(state.ids[1].source_index, -1);
}

int test_main(int, char*[])
{
    test_services();
    test_filters();
    test_algorithms();
    test_model_scopes();
    return 0;
}
