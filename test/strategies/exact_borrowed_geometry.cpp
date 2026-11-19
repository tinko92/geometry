// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.

#include <boost/core/lightweight_test.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/strategies/cartesian/exact/borrowed_predicates.hpp>
#include <boost/geometry/strategies/cartesian/exact/borrowed_geometry.hpp>
#include <boost/geometry/views/detail/boundary_view.hpp>
#include <boost/tuple/tuple.hpp>
#include <algorithm>
#include <cstdint>
#include <type_traits>

namespace bg = boost::geometry;
namespace ec = bg::strategies::detail::exact_cartesian;
namespace policy = bg::strategies::relate;

template <typename Ring>
void fill_ring(Ring& ring, int offset, int extent)
{
    using point = typename Ring::value_type;
    ring = {point(offset, offset), point(offset, offset + extent),
            point(offset + extent, offset + extent), point(offset + extent, offset)};
    if (bg::point_order<Ring>::value == bg::counterclockwise)
        std::reverse(ring.begin(), ring.end());
    if (bg::closure<Ring>::value == bg::closed) ring.push_back(ring.front());
}

template <typename Ring, typename Line>
void check_vertices(Ring const& original, Line const& output)
{
    auto const size = ec::borrowed_detail::range_size(original);
    BOOST_TEST_EQ(output.size(), size == 0 ? 0u : size + 1);
    if (size == 0) return;
    for (std::size_t i = 0; i < output.size(); ++i)
    {
        auto const index = i % size;
        BOOST_TEST_EQ(bg::get<0>(output[i]), bg::get<0>(original[index]));
        BOOST_TEST_EQ(bg::get<1>(output[i]), bg::get<1>(original[index]));
    }
}

template <typename Root, typename Ring, typename Line>
void check_boundary(ec::borrowed_context<Root>& context, Ring const& original, Line const& output)
{
    check_vertices(original, output);
    BOOST_TEST_EQ(output.incoming.size(), output.size());
    if (output.empty()) return;
    auto const id = context.register_range(original);
    BOOST_TEST(!output.incoming.front().valid());
    for (std::size_t i = 1; i < output.size(); ++i)
        BOOST_TEST((output.incoming[i] == ec::original_edge{id, i - 1}));
    BOOST_TEST(output.pending == output.incoming.back());
    using accessor = ec::iterator_point_view<typename Ring::const_iterator>;
    auto const size = ec::borrowed_detail::range_size(original);
    for (std::size_t i = 0; i < output.size(); ++i)
    {
        // Accessors must refer to original const iterators after the view and
        // any intermediate supported ring have both been destroyed.
        BOOST_TEST(output[i].is_native());
        output[i].visit([&](auto const& source)
        {
            if constexpr (std::is_same<std::decay_t<decltype(source)>, accessor>::value)
                BOOST_TEST(source.base() == original.cbegin() + i % size);
            else BOOST_TEST(false);
        });
    }
    ec::borrowed_detail::validate_supports<Root>(context, output);
}

template <typename C, bool Clockwise, bool Closed, typename Policy>
void test_boundary_conversion()
{
    using native = bg::model::d2::point_xy<C>;
    using polygon = bg::model::polygon<native, Clockwise, Closed>;
    using ring = typename polygon::ring_type;
    using polygons = bg::model::multi_polygon<polygon>;
    using root = ec::borrowed_root<C, Policy, polygons, polygons>;
    using working_polygon = typename ec::borrowed_working<root, polygon>::type;
    using line = typename ec::borrowed_geometry_types<root>::line;
    using lines = bg::model::multi_linestring<line>;
    polygons input;
    input.resize(2);
    for (std::size_t i = 0; i < input.size(); ++i)
    {
        int const offset = int(i) * 20;
        fill_ring(input[i].outer(), offset, 8);
        input[i].inners().resize(1);
        fill_ring(input[i].inners().front(), offset + 2, 2);
        std::reverse(input[i].inners().front().begin(), input[i].inners().front().end());
    }
    ec::borrowed_context<root> context;
    ec::borrowed_geometry_converter<root> converter(context);
    line boundary;
    {
        bg::detail::boundary_view<ring const> view(input.front().outer());
        converter.apply(view, boundary);
    }
    check_boundary(context, input.front().outer(), boundary);
    {
        bg::detail::boundary_view<ring> view(input.front().outer());
        converter.apply(view, boundary);
    }
    check_boundary(context, input.front().outer(), boundary);

    // A polygon/multipolygon view owns its nested boundary-view objects.
    // Destroying those must not invalidate the copied point or support metadata.
    lines native_boundaries;
    {
        bg::detail::boundary_view<polygons const> view(input);
        converter.apply(view, native_boundaries);
    }
    BOOST_TEST_EQ(native_boundaries.size(), 4u);
    if (native_boundaries.size() == 4)
        for (std::size_t i = 0; i < input.size(); ++i)
        {
            check_boundary(context, input[i].outer(), native_boundaries[2 * i]);
            check_boundary(context, input[i].inners().front(), native_boundaries[2 * i + 1]);
        }
    lines working_boundaries;
    {
        working_polygon working;
        converter.apply(input.front(), working);
        if (Closed)
        {
            working.outer().incoming.front() = {};
            working.inners().front().incoming.front() = {};
        }
        bg::detail::boundary_view<working_polygon const> view(working);
        converter.apply(view, working_boundaries);
    }
    BOOST_TEST_EQ(working_boundaries.size(), 2u);
    if (working_boundaries.size() == 2)
    {
        check_boundary(context, input.front().outer(), working_boundaries[0]);
        check_boundary(context, input.front().inners().front(), working_boundaries[1]);
    }

    using exported_point = bg::model::d2::point_xy<int>;
    using exported_line = bg::model::linestring<exported_point>;
    using exported_polygon = bg::model::polygon<exported_point>;
    // The target alternatives differ from the working alternatives and their order.
    boost::variant<line, ec::borrowed_point<root>> dynamic(boundary);
    boost::variant<exported_point, exported_polygon, exported_line> exported;
    converter.apply(dynamic, exported);
    auto const* exported_boundary = boost::get<exported_line>(&exported);
    BOOST_TEST(exported_boundary != nullptr);
    if (exported_boundary) check_vertices(input.front().outer(), *exported_boundary);

    using source_tuple = boost::tuple<polygon, native>;
    using working_tuple = typename ec::borrowed_working<root, source_tuple>::type;
    static_assert(std::is_same<typename boost::tuples::element<0, working_tuple>::type,
                              working_polygon>::value, "Retain Boost tuple working specializations");
    source_tuple const tuple(input.front(), native(3, 4));
    working_tuple working;
    converter.apply(tuple, working);
    boost::tuple<exported_polygon, exported_point> tuple_output;
    converter.apply(working, tuple_output);
    BOOST_TEST_EQ(bg::get<0>(boost::get<1>(tuple_output)), 3);
    BOOST_TEST_EQ(bg::get<1>(boost::get<1>(tuple_output)), 4);

    ring const empty;
    bg::detail::boundary_view<ring const> empty_view(empty);
    converter.apply(empty_view, boundary);
    BOOST_TEST(boundary.empty() && boundary.incoming.empty() && !boundary.pending.valid());
}

template <typename C, typename Policy>
void test_policy()
{
    test_boundary_conversion<C, true, true, Policy>();
    test_boundary_conversion<C, true, false, Policy>();
    test_boundary_conversion<C, false, true, Policy>();
    test_boundary_conversion<C, false, false, Policy>();
}

int main()
{
    test_policy<double, policy::approximate_output>();
    test_policy<double, policy::eager_nearest_output>();
    test_policy<double, policy::deferred_nearest_output>();
    test_policy<std::int64_t, policy::approximate_output>();
    test_policy<std::int64_t, policy::eager_nearest_output>();
    test_policy<std::int64_t, policy::deferred_nearest_output>();
    return boost::report_errors();
}
