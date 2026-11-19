// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.

// clang++ -std=c++17 -O0 -g -ffp-contract=off -Iinclude -I../..
//   test/strategies/exact_borrowed_turns.cpp -o /tmp/bg-borrowed-turns
#include <boost/core/lightweight_test.hpp>
#include <boost/geometry/strategies/cartesian/exact/borrowed_turns.hpp>
#include <boost/geometry/strategies/detail/turn_container.hpp>

#include <array>
#include <cstddef>
#include <cstdlib>
#include <memory>
#include <new>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <variant>

namespace bg = boost::geometry;
namespace ec = bg::strategies::detail::exact_cartesian;
namespace sd = bg::strategies::detail;
namespace ov = bg::detail::overlay;

namespace
{
bool count_allocations = false;
std::size_t allocations = 0;
}

void* operator new(std::size_t size)
{
    if (count_allocations) ++allocations;
    if (void* p = std::malloc(size ? size : 1)) return p;
    throw std::bad_alloc();
}
void operator delete(void* p) noexcept { std::free(p); }
void* operator new[](std::size_t size) { return ::operator new(size); }
void operator delete[](void* p) noexcept { std::free(p); }
using native_point = bg::model::d2::point_xy<double>;
using native_polygon = bg::model::polygon<native_point>;
using native_ring = bg::model::ring<native_point>;

struct mock_root
{
    using coordinate_type = double;
    using output_policy = bg::strategies::relate::approximate_output;
    using native_geometries = boost::mp11::mp_list<native_polygon, native_ring>;
    using native_ranges = boost::mp11::mp_list<native_ring>;
    using native_boundaries = boost::mp11::mp_list<native_ring>;
    using native_access = boost::mp11::mp_list<native_point const*,
        ec::iterator_point_view<native_ring::const_iterator>>;
    using point_type = ec::borrowed_point<mock_root>;
    using operation_type = ov::turn_operation<point_type, ec::borrowed_segment_position_tag>;
    using numeric_state = ec::construction_rounding::point_state<double, output_policy>;
};

using point = mock_root::point_type;
using record = ec::borrowed_record<mock_root>;
using stage = ec::borrowed_stage<mock_root>;
using context = ec::borrowed_context<mock_root>;
using working_polygon = ec::borrowed_geometry_types<mock_root>::polygon<true, true>;

struct observed_operation : mock_root::operation_type { int metadata = 7; };
using turn_base = ov::turn_info<point, mock_root::operation_type::segment_ratio_type, observed_operation>;

struct observed_turn : turn_base
{
    inline static int copies = 0, moves = 0, assignments = 0;
    observed_turn() = default;
    observed_turn(observed_turn const& other) : turn_base(other) { ++copies; }
    observed_turn(observed_turn&& other) : turn_base(std::move(other)) { ++moves; }
    observed_turn& operator=(observed_turn const& other)
    {
        turn_base::operator=(other);
        ++assignments;
        return *this;
    }
    observed_turn& operator=(observed_turn&& other)
    {
        turn_base::operator=(std::move(other));
        ++assignments;
        return *this;
    }
};

using turns = ec::borrowed_turn_container<mock_root, observed_turn>;
using ordinary_turn = ov::turn_info<native_point>;

struct original_layout
{
    native_point point;
    ov::method_type method;
    bool touch_only;
    bg::signed_size_type cluster_id;
    bool discarded;
    bool is_traversable;
    std::array<ordinary_turn::turn_operation_type, 2> operations;
};

static_assert(std::is_empty<ov::turn_numeric_storage<native_point>>::value, "Default storage is empty");
static_assert(sizeof(ordinary_turn) == sizeof(original_layout), "Preserve default turn size");
static_assert(alignof(ordinary_turn) == alignof(original_layout), "Preserve default turn alignment");
static_assert(std::is_standard_layout<ordinary_turn>::value, "Preserve default layout category");
static_assert(offsetof(ordinary_turn, point) == offsetof(original_layout, point), "Preserve point offset");
static_assert(offsetof(ordinary_turn, operations) == offsetof(original_layout, operations), "Preserve operations offset");
static_assert(std::is_base_of<record, observed_turn>::value, "Embed the Root record in the turn");
static_assert(!std::is_polymorphic<record>::value && !std::is_polymorphic<observed_turn>::value,
              "No virtual numerical dispatch");
static_assert(std::is_same<decltype(*std::declval<turns const&>().begin()), observed_turn const&>::value,
              "Const traversal must not mutate records");

native_polygon first_polygon()
{
    native_polygon result;
    result.outer() = native_ring{{0, 0}, {0, 6}, {6, 6}, {6, 0}, {0, 0}};
    return result;
}

native_polygon second_polygon()
{
    native_polygon result;
    result.outer() = native_ring{{-1, 2}, {-1, 4}, {7, 4}, {7, 2}, {-1, 2}};
    return result;
}

bg::segment_identifier id(int source, int segment)
{
    return bg::segment_identifier(source, -1, -1, segment);
}

void seed_candidate(observed_turn& candidate, record& seed, stage& source,
                    int key, int first = 0, int second = 3)
{
    candidate.cluster_id = key;
    candidate.operations[0].seg_id = id(0, first);
    candidate.operations[1].seg_id = id(1, second);
    seed.stage = &source;
    seed.is_seed = true;
    auto const a = source.context.support_coordinates(source.support(candidate.operations[0].seg_id));
    auto const b = source.context.support_coordinates(source.support(candidate.operations[1].seg_id));
    seed.numeric.init({{a[0], a[1], a[2], a[3], b[0], b[1], b[2], b[3]}});
    candidate.point = point(seed);
}

struct saved_point
{
    point value;
    observed_turn* owner;
    record* storage;
    std::array<bg::segment_identifier, 2> original_ids;
    std::array<ec::original_edge, 2> edges;
    double x, y;

    explicit saved_point(observed_turn& turn)
        : value(turn.point), owner(&turn), storage(&turn.point.record())
        , original_ids{{turn.operations[0].seg_id, turn.operations[1].seg_id}}
        , edges{{storage->stage->support(original_ids[0]), storage->stage->support(original_ids[1])}}
        , x(value.emit<double>(0)), y(value.emit<double>(1))
    {}

    void check(bool finalized) const
    {
        BOOST_TEST(&value.record() == storage);
        BOOST_TEST(static_cast<record*>(owner) == storage);
        BOOST_TEST(!storage->is_seed);
        BOOST_TEST(storage->operations[0] == &owner->operations[0]);
        BOOST_TEST(storage->operations[1] == &owner->operations[1]);
        BOOST_TEST_EQ(storage->stage->finalized, finalized);
        BOOST_TEST_EQ(value.emit<double>(0), x);
        BOOST_TEST_EQ(value.emit<double>(1), y);
        BOOST_TEST_EQ(value.exact_coordinate(0), ec::exact_number(x));
        BOOST_TEST_EQ(value.exact_coordinate(1), ec::exact_number(y));
        for (int i = 0; i != 2; ++i)
        {
            auto const expected = finalized
                ? bg::segment_identifier(static_cast<bg::signed_size_type>(edges[i].range), -1, -1,
                                         static_cast<bg::signed_size_type>(edges[i].segment))
                : original_ids[i];
            BOOST_TEST(owner->operations[i].seg_id == expected);
            BOOST_TEST_EQ(owner->operations[i].metadata, 7);
        }
    }
};

void check_finished(stage const& source)
{
    BOOST_TEST(source.finalized);
    BOOST_TEST(std::holds_alternative<std::monostate>(source.sources[0]));
    BOOST_TEST(std::holds_alternative<std::monostate>(source.sources[1]));
}

void test_stable_storage()
{
    auto const a = first_polygon(), b = second_polygon();
    context root;
    root.register_range(b.outer());
    root.register_range(a.outer());
    turns container(root, a, b);
    stage& source = *root.active_stage;
    record seed;
    observed_turn candidate;
    seed_candidate(candidate, seed, source, 42);
    container.push_back(candidate);
    saved_point const first(container.back());
    BOOST_TEST(first.storage != &seed);
    seed_candidate(candidate, seed, source, 13, 2, 1);
    container.push_back(candidate);
    saved_point const second(container.back());
    first.check(false);

    // Existing, non-seed records stay shared even within the current stage.
    candidate.point = first.value;
    candidate.cluster_id = 99;
    container.push_back(candidate);
    BOOST_TEST(&container.back().point.record() == first.storage);
    BOOST_TEST(static_cast<record*>(&container.back()) != first.storage);

    observed_turn native;
    native.point = point::borrow(a.outer()[2]);
    container.push_back(native);
    observed_turn* const native_owner = &container.back();
    BOOST_TEST(native_owner->point.is_native());

    for (int i = 0; i != 2048; ++i)
    {
        seed_candidate(candidate, seed, source, i % 23);
        container.push_back(candidate);
    }
    int const copies = observed_turn::copies, moves = observed_turn::moves;
    int const assignments = observed_turn::assignments;
    BOOST_TEST_EQ(copies, 2052);
    BOOST_TEST_EQ(moves, 0);
    BOOST_TEST_EQ(assignments, 0);
    auto less = [](observed_turn const& x, observed_turn const& y) { return x.cluster_id < y.cluster_id; };
    sd::sort_turns(container, less);
    auto const* front = &container.front();
    auto const* back = &container.back();
    sd::sort_turns(container, container.begin() + 1, container.begin() + 4,
        [&](observed_turn const& x, observed_turn const& y) { return less(y, x); });
    BOOST_TEST(&container.front() == front);
    BOOST_TEST(&container.back() == back);
    first.check(false); second.check(false);
    sd::erase_turns_if(container, [](observed_turn const& turn) { return turn.cluster_id == 42; });
    sd::sort_turns(container, less);
    sd::unique_turns(container,
        [](observed_turn const& x, observed_turn const& y) { return x.cluster_id == y.cluster_id; });
    BOOST_TEST_EQ(observed_turn::copies, copies);
    BOOST_TEST_EQ(observed_turn::moves, moves);
    BOOST_TEST_EQ(observed_turn::assignments, assignments);
    first.check(false); second.check(false);
    container.clear();
    BOOST_TEST(container.empty());
    first.check(false); second.check(false);
    container.finish();
    BOOST_TEST(root.active_stage == nullptr);
    first.check(true); second.check(true);
    check_finished(*first.storage->stage);
    BOOST_TEST(native_owner->operations[0].seg_id == bg::segment_identifier());
    BOOST_TEST(native_owner->operations[1].seg_id == bg::segment_identifier());
    container.finish();
    first.check(true);
}

void test_stage_lifetimes()
{
    auto const a = first_polygon(), b = second_polygon();
    context root;
    root.register_range(b.outer());
    root.register_range(a.outer());
    std::unique_ptr<saved_point> first, second;
    record foreign_seed;
    stage* initial_stage = nullptr;
    {
        turns original(root, a, b);
        initial_stage = root.active_stage;
        record seed;
        observed_turn candidate;
        seed_candidate(candidate, seed, *initial_stage, 1);
        original.push_back(candidate);
        initial_stage = root.active_stage;
        first = std::make_unique<saved_point>(original.back());
    }
    first->check(true);
    check_finished(*initial_stage);
    {
        working_polygon intermediate;
        auto& ring = intermediate.outer();
        std::size_t const range = root.register_range(a.outer());
        auto append = [&](point const& p, std::size_t segment)
        {
            ring.pending = ec::original_edge{range, segment};
            ring.push_back(p);
        };
        append(point::borrow(a.outer()[0]), 3);
        append(first->value, 0);
        append(point::borrow(a.outer()[1]), 0);
        append(point::borrow(a.outer()[2]), 1);
        append(point::borrow(a.outer()[3]), 2);
        append(point::borrow(a.outer()[0]), 3);

        turns next(root, intermediate, b);
        stage& source = *root.active_stage;
        record seed;
        observed_turn candidate;
        seed_candidate(candidate, seed, source, 2, 1, 1);
        next.push_back(candidate);
        second = std::make_unique<saved_point>(next.back());
        BOOST_TEST_EQ(second->x, 0.0);
        BOOST_TEST_EQ(second->y, 4.0);

        candidate.point = first->value;
        next.push_back(candidate);
        BOOST_TEST(&next.back().point.record() == first->storage);
        BOOST_TEST(next.back().point.record().stage == initial_stage);

        // A seed belonging to another stage is not cloned either.
        foreign_seed = *first->storage;
        foreign_seed.is_seed = true;
        candidate.point = point(foreign_seed);
        next.push_back(candidate);
        BOOST_TEST(&next.back().point.record() == &foreign_seed);
        next.erase_if([&](observed_turn const& turn) { return &turn.point.record() == &foreign_seed; });
        second->check(false);
    }
    BOOST_TEST(root.active_stage == nullptr);
    first->check(true);
    second->check(true);
    check_finished(*second->storage->stage);
}

void test_reverse_and_move()
{
    auto const a = first_polygon(), b = second_polygon();
    context root;
    turns original(root, a, b);
    stage* source = root.active_stage;
    source->reverse = {{true, true}};
    record seed;
    observed_turn candidate;
    seed_candidate(candidate, seed, *source, 1, 3, 0);
    original.push_back(candidate);
    source = root.active_stage;
    saved_point const saved(original.back());
    turns moved(std::move(original));
    BOOST_TEST(root.active_stage == source);
    saved.check(false);
    moved.finish();
    saved.check(true);
    BOOST_TEST_EQ(saved.owner->operations[0].seg_id.segment_index, 0);
    BOOST_TEST_EQ(saved.owner->operations[1].seg_id.segment_index, 3);
    BOOST_TEST(root.active_stage == nullptr);
}

void test_unwinding()
{
    auto const a = first_polygon(), b = second_polygon();
    context root;
    turns outer(root, a, b);
    stage* const previous = root.active_stage;
    stage* abandoned = nullptr;
    observed_turn* retained = nullptr;
    try
    {
        turns inner(root, a, b);
        abandoned = root.active_stage;
        record seed;
        observed_turn candidate;
        seed_candidate(candidate, seed, *abandoned, 1);
        inner.push_back(candidate);
        abandoned = root.active_stage;
        retained = &inner.back();
        throw std::runtime_error("abandon the Root operation");
    }
    catch (std::runtime_error const&) {}
    BOOST_TEST(root.active_stage == previous);
    BOOST_TEST(!abandoned->finalized);
    BOOST_TEST(!std::holds_alternative<std::monostate>(abandoned->sources[0]));
    BOOST_TEST(retained->operations[0].seg_id == id(0, 0));
    BOOST_TEST(retained->operations[1].seg_id == id(1, 3));
    outer.finish();
    BOOST_TEST(root.active_stage == nullptr);
}

void test_empty_and_moved_promotion()
{
    auto const a = first_polygon(), b = second_polygon();
    context root;
    stage previous(root, a, b);
    root.active_stage = &previous;
    allocations = 0;
    count_allocations = true;
    {
        turns original(root, a, b);
        auto const* old_stage = root.active_stage;
        turns moved(std::move(original));
        BOOST_TEST(root.active_stage != old_stage);
        auto* moved_stage = root.active_stage;
        original.finish();
        BOOST_TEST(root.active_stage == moved_stage);
        moved.finish();
        BOOST_TEST(root.active_stage == &previous);
        check_finished(*moved_stage);
    }
    count_allocations = false;
    BOOST_TEST_EQ(allocations, 0u);
    root.active_stage = nullptr;

    std::unique_ptr<saved_point> saved;
    {
        turns original(root, a, b);
        turns moved(std::move(original));
        stage* const stack_stage = root.active_stage;
        record seed;
        observed_turn candidate;
        seed_candidate(candidate, seed, *stack_stage, 1);
        stack_stage->model = {{&candidate.operations[0], &candidate.operations[1]}};
        seed.operations = stack_stage->model;
        moved.push_back(candidate);
        BOOST_TEST(root.active_stage != stack_stage);
        BOOST_TEST(moved.back().point.record().stage == root.active_stage);
        saved = std::make_unique<saved_point>(moved.back());
        stack_stage->model = {{nullptr, nullptr}};
        // The same policy call may still insert another seed from its old frame.
        moved.push_back(candidate);
        BOOST_TEST(&moved.back().point.record() != &seed);
        BOOST_TEST(moved.back().point.record().stage == saved->storage->stage);
        BOOST_TEST(moved.back().point.record().operations[0] == &moved.back().operations[0]);
        saved->check(false);
    }
    BOOST_TEST(root.active_stage == nullptr);
    saved->check(true);
}

int main()
{
    test_stable_storage();
    test_stage_lifetimes();
    test_reverse_and_move();
    test_unwinding();
    test_empty_and_moved_promotion();
    return boost::report_errors();
}
