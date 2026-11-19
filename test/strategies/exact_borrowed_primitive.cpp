// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.

#include <boost/core/lightweight_test.hpp>
#include <boost/geometry/algorithms/disjoint.hpp>
#include <boost/geometry/algorithms/intersection.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/segment.hpp>
#include <boost/geometry/strategies/relate/cartesian.hpp>

#include <array>
#include <iterator>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

#if __cplusplus >= 201703L
#include <boost/geometry/strategies/cartesian/exact/borrowed_primitive.hpp>
#endif

namespace bg = boost::geometry;
using point = bg::model::d2::point_xy<double>;
using segment = bg::model::segment<point>;
using primitive = bg::detail::intersection::intersection_segment_segment_point<point>;

void check(point const& actual, point const& expected)
{
    BOOST_TEST_EQ(bg::get<0>(actual), bg::get<0>(expected));
    BOOST_TEST_EQ(bg::get<1>(actual), bg::get<1>(expected));
}

struct unavailable_hook : bg::strategies::relate::cartesian<>
{
    template <typename S1, typename S2, typename Output, typename Convert,
              std::enable_if_t<sizeof(S1) == 0, int> = 0>
    Output intersect_segments(S1 const&, S2 const&, Output, Convert&&) const;
    template <typename S1, typename S2, std::enable_if_t<sizeof(S1) == 0, int> = 0>
    bool disjoint_segments(S1 const&, S2 const&) const;
};

template <typename Strategy>
void test_defaults(Strategy const& strategy)
{
    auto test = [&](segment const& a, segment const& b, std::vector<point> const& expected)
    {
        std::vector<point> result;
        bg::intersection(a, b, result, strategy);
        BOOST_TEST_EQ(bg::disjoint(a, b, strategy), expected.empty());
        BOOST_TEST_EQ(result.size(), expected.size());
        for (std::size_t i = 0; i < result.size() && i < expected.size(); ++i)
        {
            check(result[i], expected[i]);
        }
    };
    test(segment(point(0, 0), point(3, 3)), segment(point(0, 3), point(3, 0)), {{1.5, 1.5}});
    test(segment(point(0, 0), point(4, 0)), segment(point(2, 0), point(6, 0)), {{2, 0}, {4, 0}});
    test(segment(point(4, 0), point(0, 0)), segment(point(2, 0), point(6, 0)), {{4, 0}, {2, 0}});
    test(segment(point(0, 0), point(1, 1)), segment(point(1, 1), point(2, 2)), {{1, 1}});
    test(segment(point(0, 0), point(1, 0)), segment(point(2, 0), point(3, 0)), {});
    test(segment(point(1, 0), point(1, 0)), segment(point(0, 0), point(2, 0)), {{1, 0}});
}

struct hooked_strategy : bg::strategies::relate::cartesian<>
{
    mutable int calls = 0;
    mutable int disjoint_calls = 0;
    void relate() const = delete;

    template <typename S1, typename S2>
    bool disjoint_segments(S1 const&, S2 const&) const
    {
        ++disjoint_calls;
        return true;
    }

    template <typename S1, typename S2, typename Output, typename Convert>
    Output intersect_segments(S1 const&, S2 const&, Output out, Convert&& convert) const
    {
        ++calls;
        *out++ = convert(point(9, 7));
        return out;
    }
};

void test_hook_dispatch()
{
    segment const a(point(0, 0), point(3, 3)), b(point(0, 3), point(3, 0));
    hooked_strategy strategy;
    std::array<point, 2> result;
    auto* end = primitive::apply(a, b, result.data(), strategy);
    BOOST_TEST(end == result.data() + 1);
    BOOST_TEST_EQ(strategy.calls, 1);
    check(result[0], point(9, 7));
    BOOST_TEST(bg::disjoint(a, b, strategy));
    BOOST_TEST_EQ(strategy.disjoint_calls, 1);
}

#if __cplusplus >= 201703L
namespace ec = bg::strategies::detail::exact_cartesian;

template <typename Root>
struct bound_strategy : bg::strategies::relate::cartesian<>
{
    ec::borrowed_context<Root>& context;
    ec::borrowed_stage<Root>* expected_active = nullptr;
    mutable int calls = 0, conversions = 0;

    explicit bound_strategy(ec::borrowed_context<Root>& context) : context(context) {}
    void relate() const = delete;

    template <typename S1, typename S2, typename Output, typename Convert>
    Output intersect_segments(S1 const& first, S2 const& second, Output out, Convert&& convert) const
    {
        ++calls;
        return ec::borrowed_intersect_segments(context, first, second, out, std::forward<Convert>(convert));
    }

    template <typename Source, typename Destination>
    void convert_geometry(Source const& source, Destination& destination) const
    {
        ++conversions;
        BOOST_TEST(context.active_stage == expected_active);
        ec::borrowed_geometry_converter<Root>(context).apply(source, destination);
    }
};

template <typename Root>
void check_finalized(ec::borrowed_point<Root> const& p,
                     ec::original_edge const& first, ec::original_edge const& second)
{
    BOOST_TEST(!p.is_native());
    auto const& record = p.record();
    BOOST_TEST(!record.is_seed);
    BOOST_TEST(record.stage->finalized);
    BOOST_TEST(std::holds_alternative<std::monostate>(record.stage->sources[0]));
    BOOST_TEST(std::holds_alternative<std::monostate>(record.stage->sources[1]));
    BOOST_TEST(record.operations == record.stage->model);
    BOOST_TEST(record.operations[0]->seg_id == bg::segment_identifier(first.range, -1, -1, first.segment));
    BOOST_TEST(record.operations[1]->seg_id == bg::segment_identifier(second.range, -1, -1, second.segment));
    BOOST_TEST_EQ(p.exact_coordinate(0), ec::exact_number(3) / 2);
    BOOST_TEST_EQ(p.exact_coordinate(1), ec::exact_number(3) / 2);
    BOOST_TEST_EQ(p.template emit<double>(0), 1.5);
    BOOST_TEST_EQ(p.refine(1).approximate, 1.5);
}

template <typename Policy>
void test_root_primitive()
{
    using root = ec::borrowed_root<double, Policy, segment, segment>;
    using borrowed_point = ec::borrowed_point<root>;
    segment const a(point(0, 0), point(3, 3)), b(point(0, 3), point(3, 0));
    ec::borrowed_context<root> context;
    auto const second_range = context.register_range(b);
    auto const first_range = context.register_range(a);
    ec::original_edge const first{first_range, 0}, second{second_range, 0};
    bound_strategy<root> strategy(context);

    std::vector<point> crossing;
    bg::intersection(a, b, crossing, strategy);
    BOOST_TEST_EQ(crossing.size(), 1u);
    if (!crossing.empty()) check(crossing[0], point(1.5, 1.5));
    BOOST_TEST_EQ(strategy.calls, 1);
    BOOST_TEST_EQ(strategy.conversions, 1);
    BOOST_TEST(context.active_stage == nullptr);

    std::vector<borrowed_point> saved;
    bg::intersection(a, b, saved, strategy);
    BOOST_TEST_EQ(saved.size(), 1u);
    if (saved.empty()) return;
    check_finalized(saved[0], first, second);
    auto const* const retained_record = &saved[0].record();

    BOOST_TEST(!ec::borrowed_segments_disjoint<root>(a, b));
    segment const far(point(10, 10), point(11, 11));
    BOOST_TEST(ec::borrowed_segments_disjoint<root>(a, far));
    segment const touching(point(3, 3), point(4, 4));
    BOOST_TEST(!ec::borrowed_segments_disjoint<root>(a, touching));
    segment const degenerate(point(1, 1), point(1, 1));
    BOOST_TEST(!ec::borrowed_segments_disjoint<root>(a, degenerate));
    bg::model::segment<borrowed_point> const constructed
        (saved[0], borrowed_point::borrow(a.first));
    BOOST_TEST(!ec::borrowed_segments_disjoint<root>(constructed, b));
    BOOST_TEST(ec::borrowed_segments_disjoint<root>(constructed, far));

    // No turn container owns the result; context retention alone preserves it.
    for (int i = 0; i != 128; ++i)
    {
        std::array<point, 2> output;
        auto* end = primitive::apply(a, b, output.data(), strategy);
        BOOST_TEST(end == output.data() + 1);
        check(output[0], point(1.5, 1.5));
    }
    BOOST_TEST(&saved[0].record() == retained_record);
    check_finalized(saved[0], first, second);

    ec::borrowed_stage<root> outer(context, a, b);
    context.active_stage = &outer;
    strategy.expected_active = &outer;
    std::vector<point> nested;
    bg::intersection(a, b, nested, strategy);
    BOOST_TEST(context.active_stage == &outer);
    BOOST_TEST(!outer.finalized);
    check(nested.front(), point(1.5, 1.5));

    borrowed_point converted_before_throw;
    try
    {
        ec::borrowed_intersect_segments(context, a, b, std::back_inserter(nested),
            [&](borrowed_point const& p) -> point
            {
                BOOST_TEST(context.active_stage == &outer);
                converted_before_throw = p;
                throw std::runtime_error("converter failure");
            });
        BOOST_ERROR("Converter should throw");
    }
    catch (std::runtime_error const&) {}
    BOOST_TEST(context.active_stage == &outer);
    check_finalized(converted_before_throw, first, second);
    context.active_stage = nullptr;

    segment const c(point(0, 0), point(4, 0)), d(point(2, 0), point(6, 0));
    ec::borrowed_context<root> overlap_context;
    bound_strategy<root> overlap_strategy(overlap_context);
    std::vector<point> overlap;
    bg::intersection(c, d, overlap, overlap_strategy);
    BOOST_TEST_EQ(overlap.size(), 2u);
    if (overlap.size() == 2)
    {
        check(overlap[0], point(2, 0));
        check(overlap[1], point(4, 0));
    }
    std::vector<borrowed_point> selected;
    bg::intersection(c, d, selected, overlap_strategy);
    BOOST_TEST_EQ(selected.size(), 2u);
    if (selected.size() == 2)
    {
        BOOST_TEST(selected[0].is_native());
        BOOST_TEST(selected[1].is_native());
        BOOST_TEST_EQ(selected[0].template emit<double>(0), 2.0);
        BOOST_TEST_EQ(selected[1].template emit<double>(0), 4.0);
    }
    BOOST_TEST(overlap_context.active_stage == nullptr);
}
#endif

int main()
{
    test_defaults(bg::strategies::relate::cartesian<>());
    test_defaults(unavailable_hook());
    test_hook_dispatch();
#if __cplusplus >= 201703L
    test_root_primitive<bg::strategies::relate::approximate_output>();
    test_root_primitive<bg::strategies::relate::eager_nearest_output>();
    test_root_primitive<bg::strategies::relate::deferred_nearest_output>();
#endif
    return boost::report_errors();
}
