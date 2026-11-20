// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.

#include <boost/core/lightweight_test.hpp>
#include <boost/geometry/geometries/pointing_segment.hpp>
#include <boost/geometry/strategies/cartesian/exact/borrowed_bounds.hpp>

#include <array>
#include <cstdint>
#include <limits>

namespace bg = boost::geometry;
namespace ec = bg::strategies::detail::exact_cartesian;

template <typename Policy>
void test_segment_box()
{
    using C = std::int64_t;
    using native = bg::model::d2::point_xy<C>;
    using native_segment = bg::model::segment<native>;
    using native_box = bg::model::box<native>;
    using metadata_point = bg::model::point<double, 2, bg::cs::cartesian>;
    using metadata_box = bg::model::box<metadata_point>;
    using root = ec::borrowed_root<C, Policy, native_segment, native_box>;
    using point = ec::borrowed_point<root>;
    using predicate = ec::root_segment_box<root>;
    using root_capacity = ec::construction_capacity<C, ec::construction_family::original_support_intersection>;
    using double_grid = ec::construction_capacity<double, ec::construction_family::native_vertex>;
    static_assert(ec::construction_predicate_capacity<root_capacity, root_capacity, double_grid>::orientation_bits
        == 2427, "Metadata orientation must retain the actual double grid, not the integral Root format");
    ec::borrowed_context<root> context;
    metadata_box const bounds(metadata_point(-0.5, -0.5), metadata_point(0.5, 0.5));
    ec::search_box const search(ec::search_point(-0.5, -0.5), ec::search_point(0.5, 0.5));
    std::array<native, 9> const values{{{-2, -2}, {-2, 0}, {-2, 2}, {0, -2}, {0, 0},
                                      {0, 2}, {2, -2}, {2, 0}, {2, 2}}};
    for (auto const& first : values)
    for (auto const& second : values)
    {
        point const a = point::borrow(first), b = point::borrow(second);
        bg::model::referring_segment<point const> const referred(a, b);
        bg::model::pointing_segment<point const> const pointed(a, b);
        bg::model::referring_segment<native const> const original(first, second);
        bool const expected = ec::segment_box<C>::apply(original, bounds);
        BOOST_TEST_EQ(predicate::apply(referred, bounds), expected);
        BOOST_TEST_EQ(predicate::apply(pointed, search), expected);
    }
    native const x0(-2, 0), x1(0, 2);
    point const a = point::borrow(x0), b = point::borrow(x1);
    bg::model::referring_segment<point const> const diagonal(a, b);
    BOOST_TEST(predicate::apply(diagonal, bounds));
    native_box const original_box(native(-1, -1), native(1, 1));
    ec::root_bounds<root> origins;
    ec::root_envelope<root>::apply(original_box, origins);
    BOOST_TEST(!predicate::apply(diagonal, origins)); // Closed contact at (-1,1).
    BOOST_TEST(!predicate::apply(diagonal, original_box));
    BOOST_TEST(predicate::get_disjoint_point_box_strategy().apply(point::borrow(original_box.min_corner()), origins));

    double const inf = (std::numeric_limits<double>::infinity)();
    for (double lx : {-inf, 0.0}) for (double ux : {0.0, inf})
    for (double ly : {-inf, 0.0}) for (double uy : {0.0, inf})
    {
        metadata_box const infinite(metadata_point(lx, ly), metadata_point(ux, uy));
        native_box const clipped(native(std::isinf(lx) ? -3 : 0, std::isinf(ly) ? -3 : 0),
                                 native(std::isinf(ux) ? 3 : 0, std::isinf(uy) ? 3 : 0));
        for (auto const& first : values) for (auto const& second : values)
        {
            point const a = point::borrow(first), b = point::borrow(second);
            bg::model::pointing_segment<point const> const segment(a, b);
            bg::model::referring_segment<native const> const original(first, second);
            BOOST_TEST_EQ(predicate::apply(segment, infinite), ec::segment_box<C>::apply(original, clipped));
        }
    }
    ec::root_bounds<root> const empty;
    BOOST_TEST(predicate::apply(diagonal, empty));
    native_segment const left(native(0, 0), native(1, 1));
    native_segment const right(native(0, 1), native(1, 0));
    ec::borrowed_stage<root> stage(context, left, right);
    typename root::operation_type operations[2];
    operations[0].seg_id = bg::segment_identifier(0, -1, -1, 0);
    operations[1].seg_id = bg::segment_identifier(1, -1, -1, 0);
    ec::borrowed_record<root> record;
    record.stage = &stage;
    record.operations = {{&operations[0], &operations[1]}};
    point const crossing(record), origin = point::borrow(left.first);
    record.numeric.init(crossing.producers());
    bg::model::pointing_segment<point const> const implicit(crossing, origin);
    metadata_box const below(metadata_point(0.3, 0.1), metadata_point(0.4, 0.2));
    metadata_box const contact(metadata_point(0.25, 0.25), metadata_point(0.5, 0.5));
    BOOST_TEST(predicate::apply(implicit, below));
    BOOST_TEST(!predicate::apply(implicit, contact));
    BOOST_TEST(ec::operation_detail::require<C>().storage().auxiliary.empty());
}

int main()
{
    test_segment_box<bg::strategies::relate::approximate_output>();
    test_segment_box<bg::strategies::relate::nearest_output>();
    return boost::report_errors();
}
