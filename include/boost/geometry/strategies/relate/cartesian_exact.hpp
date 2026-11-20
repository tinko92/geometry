// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Licensed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_GEOMETRY_STRATEGIES_RELATE_CARTESIAN_EXACT_HPP
#define BOOST_GEOMETRY_STRATEGIES_RELATE_CARTESIAN_EXACT_HPP

#include <type_traits>
#include <tuple>
#include <vector>
#include <boost/range/size.hpp>
#include <boost/geometry/algorithms/convert.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/ring.hpp>
#include <boost/geometry/geometries/linestring.hpp>
#include <boost/geometry/geometries/multi_point.hpp>
#include <boost/geometry/geometries/multi_linestring.hpp>
#include <boost/geometry/geometries/multi_polygon.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/geometry_collection.hpp>
#include <boost/geometry/geometries/adapted/boost_variant.hpp>
#include <boost/geometry/core/geometry_types.hpp>
#include <boost/geometry/core/visit.hpp>
#include <boost/geometry/util/sequence.hpp>
#include <boost/geometry/util/tuples.hpp>
#include <boost/geometry/strategies/cartesian/exact/iterator_point_view.hpp>
#include <boost/geometry/strategies/cartesian/exact/native_predicates.hpp>
#include <boost/geometry/strategies/cartesian/exact/construction_rounding.hpp>
#include <boost/geometry/strategies/cartesian/exact/static_box.hpp>
#include <boost/geometry/strategies/cartesian/exact/static_predicates.hpp>
#include <boost/geometry/strategies/cartesian/exact/static_ring_topology.hpp>
#include <boost/geometry/strategies/cartesian/exact/segment.hpp>
#include <boost/geometry/strategies/relate/cartesian.hpp>
#include <boost/geometry/strategies/detail/distance_measure.hpp>
#include <boost/geometry/strategies/detail/ring_properties.hpp>
#include <boost/geometry/strategies/detail/section.hpp>
#include <boost/geometry/strategies/detail/box_intersection.hpp>
#include <boost/geometry/strategies/detail/box_topology.hpp>

namespace boost { namespace geometry { namespace strategies
{
namespace detail { namespace exact_cartesian
{

using bounding_box = search_box;

template <typename Coordinate>
struct side
{
    using cs_tag = cartesian_tag;

    template <typename A, typename B, typename P>
    static int apply(A const& a, B const& b, P const& p)
    {
        check_coordinate_type<Coordinate, A>();
        check_coordinate_type<Coordinate, B>();
        check_coordinate_type<Coordinate, P>();
        return static_orientation<Coordinate>(a, b, p);
    }

    template <typename A, typename B, typename P>
    static int direction_code(A const& first, B const& second, P const& point)
    {
        return static_direction<Coordinate>(first, second, point);
    }
};

template <typename Coordinate>
struct point_equal
{
    using cs_tag = cartesian_tag;

    template <typename A, typename B>
    static bool apply(A const& a, B const& b)
    {
        check_coordinate_type<Coordinate, A>();
        check_coordinate_type<Coordinate, B>();
        return static_equal_points<Coordinate>(a, b);
    }
};

template <typename Coordinate>
struct winding
{
    using cs_tag = cartesian_tag;
    struct state_type { int count = 0; bool boundary = false; };

    template <typename P, typename FirstPoint, typename SecondPoint = FirstPoint>
    static bool apply(P const& point, FirstPoint const& first, SecondPoint const& second,
                      state_type& state)
    {
        check_coordinate_type<Coordinate, P>();
        check_coordinate_type<Coordinate, FirstPoint>();
        check_coordinate_type<Coordinate, SecondPoint>();
        int const ay = static_compare_coordinate<1, Coordinate>(first, point);
        int const by = static_compare_coordinate<1, Coordinate>(second, point);
        // A segment strictly above or below the query cannot cross its ray or
        // contain it. Equal ordinates need only the exact abscissa comparisons.
        if (ay * by > 0) { return true; }
        if (ay == 0 && by == 0)
        {
            if (static_compare_coordinate<0, Coordinate>(first, point)
                * static_compare_coordinate<0, Coordinate>(second, point) <= 0)
            {
                state.boundary = true;
                return false;
            }
            return true;
        }
        int const side = static_orientation<Coordinate>(first, second, point);
        if (side == 0 && static_compare_coordinate<0, Coordinate>(first, point)
                * static_compare_coordinate<0, Coordinate>(second, point) <= 0
            && ay * by <= 0)
        {
            state.boundary = true;
            return false;
        }
        if (ay <= 0 && by > 0 && side > 0)
        {
            ++state.count;
        }
        else if (by <= 0 && ay > 0 && side < 0)
        {
            --state.count;
        }
        return true;
    }

    static int result(state_type const& state)
    {
        return state.boundary ? 0 : state.count == 0 ? -1 : 1;
    }
};

template <typename Coordinate, typename ComparePolicy>
struct point_compare
{
    template <typename A, typename B>
    static bool apply(A const& first, B const& second)
    {
        check_coordinate_type<Coordinate, A>();
        check_coordinate_type<Coordinate, B>();
        int const x = static_compare_coordinate<0, Coordinate>(first, second);
        int const order = x == 0 ? static_compare_coordinate<1, Coordinate>(first, second) : x;
        return std::is_same<ComparePolicy, strategy::compare::equal_to>::value
            ? order == 0 : ComparePolicy::apply(order, 0);
    }
};

template <typename Coordinate>
struct topology
{
    template <typename A, typename B>
    bool cluster_equal(A const& a, B const& b) const
    {
        return point_equal<Coordinate>::apply(a, b);
    }

    template <typename A, typename B>
    int compare_y(A const& a, B const& b) const
    {
        return static_compare_coordinate<1, Coordinate>(a, b);
    }

    template <typename A, typename B>
    bool cluster_window_exceeded(A const& upper, B const& lower) const
    {
        return compare_y(upper, lower) > 0;
    }

    template <typename Indices, typename Turns>
    void colocate_cluster(Indices const& indices, Turns& turns) const
    {
        auto const point = turns[*indices.begin()].point;
        for (auto const index : indices)
        {
            turns[index].point = point;
        }
    }

    template <typename A, typename B>
    bool same_direction_point(A const& a, B const& b) const
    {
        return point_equal<Coordinate>::apply(a, b);
    }

    template <typename A, typename B>
    auto squared_distance(A const& first, B const& second) const
    {
        return static_squared_distance<Coordinate>(first, second);
    }

    template <typename Distance>
    bool distance_is_zero(Distance const& distance) const { return distance.is_zero(); }

    template <typename A, typename B, typename P>
    auto distance_measure(A const& first, B const& second, P const& point) const
    {
        return static_distance_measure<Coordinate>(first, second, point);
    }
};

template <typename Coordinate>
struct comparable_squared_distance
{
    template <typename A, typename B>
    static auto apply(A const& a, B const& b)
    { return static_squared_distance<Coordinate>(a, b); }
};

template <typename Coordinate>
struct static_section : strategies::detail::strict_cartesian_section
{
private:
    template <typename DimensionVector, std::size_t Index, typename Point>
    static int direction(Point const& previous, Point const& current)
    {
        using dimension = typename util::sequence_element<Index, DimensionVector>::type;
        return static_box_detail::compare_coordinate<dimension::value, Coordinate>(current, previous);
    }
    template <typename DimensionVector, typename Point, std::size_t... Index>
    static void directions_impl(Point const& previous, Point const& current, int* result,
                                std::index_sequence<Index...>)
    {
        int const unused[] = {0, (result[Index] = direction<DimensionVector, Index>(previous, current), 0)...};
        (void)unused;
    }
public:
    template <typename Iterator>
    static void assign_point(Iterator const& iterator, iterator_point_view<Iterator>& point)
    {
        point = iterator_point_view<Iterator>(iterator);
    }

    template <typename Iterator, typename Point>
    static auto assign_point(Iterator iterator, Point& point)
        -> std::enable_if_t<std::is_same<typename std::iterator_traits<Iterator>::value_type, Point>::value>
    {
        point = *iterator;
    }

    template <typename DimensionVector, typename Point>
    static void directions(Point const& previous, Point const& current, int* result)
    {
        directions_impl<DimensionVector>(previous, current, result,
            std::make_index_sequence<util::sequence_size<DimensionVector>::value>());
    }

    template <typename Point>
    static bool duplicate(Point const& previous, Point const& current)
    {
        return point_equal<Coordinate>::apply(previous, current);
    }

    template <std::size_t Dimension, typename Point, typename Box>
    static bool preceding(int direction, Point const& point, Box const&, Box const& other)
    {
        return (direction == 1 && static_box_detail::compare_coordinate<Dimension, Coordinate>
                    (point, static_box_detail::corner<min_corner>(other)) < 0)
            || (direction == -1 && static_box_detail::compare_coordinate<Dimension, Coordinate>
                    (point, static_box_detail::corner<max_corner>(other)) > 0);
    }
};

template <typename Coordinate, typename Geometry, typename Tag = geometry::tag_t<Geometry>>
struct working_geometry { using type = Geometry; };

template <typename Coordinate, typename Geometry>
struct working_geometry<Coordinate, Geometry, linestring_tag>
{
    using type = model::linestring<exact_point<Coordinate>>;
};

template <typename Coordinate, typename Geometry>
struct working_geometry<Coordinate, Geometry, point_tag>
{
    using type = exact_point<Coordinate>;
};

template <typename Coordinate, typename Geometry>
struct working_geometry<Coordinate, Geometry, box_tag>
{
    using type = model::box<exact_point<Coordinate>>;
};

template <typename Coordinate, typename Geometry>
struct working_geometry<Coordinate, Geometry, ring_tag>
{
    using type = model::ring<exact_point<Coordinate>,
        geometry::point_order<Geometry>::value == clockwise,
        geometry::closure<Geometry>::value == closed>;
};

template <typename Coordinate, typename Geometry>
struct working_geometry<Coordinate, Geometry, polygon_tag>
{
    using type = model::polygon<exact_point<Coordinate>,
        geometry::point_order<Geometry>::value == clockwise,
        geometry::closure<Geometry>::value == closed>;
};

template <typename Coordinate, typename Geometry>
struct working_geometry<Coordinate, Geometry, multi_point_tag>
{
    using type = model::multi_point<exact_point<Coordinate>>;
};

template <typename Coordinate, typename Geometry>
struct working_geometry<Coordinate, Geometry, multi_linestring_tag>
{
    using type = model::multi_linestring<model::linestring<exact_point<Coordinate>>>;
};

template <typename Coordinate, typename Geometry>
struct working_geometry<Coordinate, Geometry, multi_polygon_tag>
{
    using type = model::multi_polygon<typename working_geometry<Coordinate,
        typename boost::range_value<Geometry>::type>::type>;
};

template <typename Coordinate, typename Sequence>
struct working_dynamic;

template <typename Coordinate, typename... Geometry>
struct working_dynamic<Coordinate, util::type_sequence<Geometry...>>
{
    using type = boost::variant<typename working_geometry<Coordinate, Geometry>::type...>;
};

template <typename Coordinate, typename Geometry>
struct working_geometry<Coordinate, Geometry, dynamic_geometry_tag>
    : working_dynamic<Coordinate, typename traits::geometry_types<Geometry>::type>
{};

template <typename Coordinate, typename Geometry>
struct working_geometry<Coordinate, Geometry, geometry_collection_tag>
{
    using type = model::geometry_collection<typename working_dynamic<Coordinate,
        typename traits::geometry_types<Geometry>::type>::type>;
};

template <typename Coordinate, typename Geometry>
struct working_geometry<Coordinate, Geometry, void>
{
    using type = std::vector<typename working_geometry<Coordinate,
        typename boost::range_value<Geometry>::type>::type>;
};

template <typename Coordinate, typename... Geometry>
struct working_geometry<Coordinate, std::tuple<Geometry...>, void>
{
    using type = std::tuple<typename working_geometry<Coordinate, Geometry>::type...>;
};

template <typename Coordinate, typename First, typename Second>
struct working_geometry<Coordinate, std::pair<First, Second>, void>
{
    using type = std::pair<typename working_geometry<Coordinate, First>::type,
                           typename working_geometry<Coordinate, Second>::type>;
};

template <typename Coordinate>
struct working_geometry<Coordinate, boost::tuples::null_type, void>
{
    using type = boost::tuples::null_type;
};

template <typename Coordinate, typename Head, typename Tail>
struct working_geometry<Coordinate, boost::tuples::cons<Head, Tail>, void>
{
    using type = boost::tuples::cons<typename working_geometry<Coordinate, Head>::type,
                                    typename working_geometry<Coordinate, Tail>::type>;
};

template <typename Coordinate, typename... Geometry>
struct working_geometry<Coordinate, boost::tuples::tuple<Geometry...>, void>
{
    using type = boost::tuples::tuple<typename working_geometry<Coordinate, Geometry>::type...>;
};

template <typename Coordinate, typename Collection, bool IsTuple = tuples::is_tuple<Collection>::value>
struct working_output : working_geometry<Coordinate, Collection> {};

template <typename Coordinate, typename Collection>
struct working_output<Coordinate, Collection, true>
{
    template <std::size_t... I>
    static auto type_of(std::index_sequence<I...>)
        -> std::tuple<typename working_output<Coordinate,
            typename tuples::element<I, Collection>::type>::type...>;
    using type = decltype(type_of(std::make_index_sequence<tuples::size<Collection>::value>()));
};

template <typename Coordinate, typename OutputPolicy = strategies::relate::approximate_output>
struct point_converter
{
    template <typename Source>
    static void apply(Source const& source, exact_point<Coordinate>& destination)
    {
        destination = make_point<Coordinate>(source);
    }

    template <typename Source, typename Destination>
    static void apply(Source const& source, Destination& destination)
    {
        if constexpr (static_predicates_detail::exact_input<std::remove_cv_t<Source>>::value)
        {
            assign_output<OutputPolicy>(source, destination);
        }
        else
        {
            using target = geometry::coordinate_type_t<Destination>;
            auto const x = geometry::util::numeric_cast<target>(geometry::get<0>(source));
            auto const y = geometry::util::numeric_cast<target>(geometry::get<1>(source));
            geometry::set<0>(destination, x);
            geometry::set<1>(destination, y);
        }
    }

};

// Most overlay calls request implicit points. Direct primitive callers may
// request native output instead, so keep that conversion at the same boundary.
template <typename Coordinate, typename OutputPolicy>
struct output_segment_intersection : segment_intersection<Coordinate>
{
    using base_type = segment_intersection<Coordinate>;
    using base_type::apply;
    using exact_result = segment_intersection_result<exact_point<Coordinate>>;

    template <typename Result>
    static Result convert(exact_result const& exact)
    {
        Result result;
        result.count = exact.count;
        for (std::size_t i = 0; i != exact.count; ++i)
            point_converter<Coordinate, OutputPolicy>::apply(exact.intersections[i], result.intersections[i]);
        return result;
    }

    template <typename First, typename Second, typename Result>
    static Result apply(First const& first, Second const& second,
                        policies::relate::segments_intersection_points<Result> const& policy)
    {
        if constexpr (std::is_same<typename Result::point_type, exact_point<Coordinate>>::value)
        {
            return base_type::apply(first, second, policy);
        }
        else
        {
            auto const exact = base_type::apply(first, second,
                policies::relate::segments_intersection_points<exact_result>());
            return convert<Result>(exact);
        }
    }

    template <typename First, typename Second, typename Result>
    static auto apply(First const& first, Second const& second,
                      policies::relate::segments_intersection_policy<Result> const& policy)
        -> typename policies::relate::segments_intersection_policy<Result>::return_type
    {
        if constexpr (std::is_same<typename Result::point_type, exact_point<Coordinate>>::value)
        {
            return base_type::apply(first, second, policy);
        }
        else
        {
            auto const exact = base_type::apply(first, second,
                policies::relate::segments_intersection_policy<exact_result>());
            return {convert<Result>(exact.intersection_points), exact.direction};
        }
    }
};

template <typename Coordinate, typename Tag, typename OutputPolicy = strategies::relate::approximate_output>
struct geometry_converter
{
    template <typename Source, typename Destination>
    static void apply(Source const& source, Destination& destination,
        std::enable_if_t<!util::is_polygon<Destination>::value>* = nullptr)
    {
        geometry::detail::conversion::range_to_range::apply(source, destination,
                                                           point_converter<Coordinate, OutputPolicy>());
    }

    template <typename Source, typename Destination>
    static void apply(Source const& source, Destination& destination,
        std::enable_if_t<util::is_polygon<Destination>::value>* = nullptr)
    {
        geometry::clear(destination);
        geometry::detail::conversion::range_to_range::apply(source, exterior_ring(destination),
                                                           point_converter<Coordinate, OutputPolicy>());
    }
};

template <typename Coordinate, typename OutputPolicy>
struct geometry_converter<Coordinate, point_tag, OutputPolicy> : point_converter<Coordinate, OutputPolicy> {};

template <typename Coordinate, typename OutputPolicy>
struct geometry_converter<Coordinate, ring_tag, OutputPolicy>
{
    template <typename FirstIterator, typename Destination>
    struct ring_point_converter
    {
        Destination* ring;
        FirstIterator first;
        std::size_t closure_index;
        std::size_t index = 0;

        template <typename SourcePoint, typename OutputPoint>
        void apply(SourcePoint const& source, OutputPoint& destination)
        {
            // Closing the ring reuses its first emitted vertex, even if the
            // source's closing occurrence has been refined independently.
            if (index == 0)
            {
                point_converter<Coordinate, OutputPolicy>::apply(*first, destination);
            }
            else if (index == closure_index)
            {
                destination = *boost::begin(*ring);
            }
            else
            {
                point_converter<Coordinate, OutputPolicy>::apply(source, destination);
            }
            ++index;
        }
    };

    template <typename Source, typename Destination>
    static void apply(Source const& source, Destination& destination,
        std::enable_if_t<!util::is_polygon<Destination>::value>* = nullptr)
    {
        if constexpr (static_predicates_detail::exact_input<geometry::point_type_t<Destination>>::value)
        {
            geometry::detail::conversion::range_to_range::apply(source, destination,
                point_converter<Coordinate, OutputPolicy>());
            return;
        }
        auto const source_first = boost::begin(source);
        if (source_first == boost::end(source))
        {
            geometry::clear(destination);
            return;
        }
        using policy = ring_point_converter<decltype(source_first), Destination>;
        geometry::detail::conversion::range_to_range::apply(source, destination,
            policy{&destination, source_first,
                boost::size(source) - (geometry::closure<Source>::value == closed ? 1 : 0)});
    }

    template <typename Source, typename Destination>
    static void apply(Source const& source, Destination& destination,
        std::enable_if_t<util::is_polygon<Destination>::value>* = nullptr)
    {
        geometry::clear(destination);
        apply(source, exterior_ring(destination));
    }
};

template <typename Coordinate, typename OutputPolicy>
struct geometry_converter<Coordinate, polygon_tag, OutputPolicy>
{
    template <typename Source, typename Destination>
    static void apply(Source const& source, Destination& destination,
        std::enable_if_t<util::is_polygon<Destination>::value>* = nullptr)
    {
        geometry_converter<Coordinate, ring_tag, OutputPolicy>::apply(exterior_ring(source), exterior_ring(destination));
        auto const& src = interior_rings(source);
        auto&& dst = interior_rings(destination);
        using rings_type = std::remove_reference_t<decltype(dst)>;
        traits::resize<rings_type>::apply(dst, boost::size(src));
        auto output = boost::begin(dst);
        for (auto const& ring : src)
        {
            geometry_converter<Coordinate, ring_tag, OutputPolicy>::apply(ring, *output++);
        }
    }

    template <typename Source, typename Destination>
    static void apply(Source const& source, Destination& destination,
        std::enable_if_t<util::is_ring<Destination>::value>* = nullptr)
    {
        geometry_converter<Coordinate, ring_tag, OutputPolicy>::apply(exterior_ring(source), destination);
    }
};

template <typename Coordinate, typename OutputPolicy>
struct multi_geometry_converter
{
    template <typename Source, typename Destination>
    static void apply(Source const& source, Destination& destination)
    {
        traits::resize<Destination>::apply(destination, boost::size(source));
        auto output = boost::begin(destination);
        for (auto const& item : source)
        {
            geometry_converter<Coordinate, geometry::tag_t<decltype(item)>, OutputPolicy>::apply(item, *output++);
        }
    }
};

template <typename Coordinate, typename OutputPolicy>
struct geometry_converter<Coordinate, multi_point_tag, OutputPolicy> : multi_geometry_converter<Coordinate, OutputPolicy> {};
template <typename Coordinate, typename OutputPolicy>
struct geometry_converter<Coordinate, multi_linestring_tag, OutputPolicy> : multi_geometry_converter<Coordinate, OutputPolicy> {};
template <typename Coordinate, typename OutputPolicy>
struct geometry_converter<Coordinate, multi_polygon_tag, OutputPolicy> : multi_geometry_converter<Coordinate, OutputPolicy> {};

template <typename Coordinate, typename OutputPolicy>
struct geometry_converter<Coordinate, geometry_collection_tag, OutputPolicy> : multi_geometry_converter<Coordinate, OutputPolicy> {};

template <typename Geometry>
struct matching_geometry_tag
{
    template <typename Candidate>
    struct apply : std::is_same<geometry::tag_t<Geometry>, geometry::tag_t<Candidate>> {};
};

template <typename Coordinate, typename OutputPolicy>
struct geometry_converter<Coordinate, dynamic_geometry_tag, OutputPolicy>
{
    template <typename Source, typename Destination>
    static void apply(Source const& source, Destination& destination)
    {
        traits::visit<Source>::apply([&](auto const& geometry)
        {
            using source_type = util::remove_cref_t<decltype(geometry)>;
            using destination_type = typename util::sequence_find_if
                <typename traits::geometry_types<Destination>::type,
                 matching_geometry_tag<source_type>::template apply>::type;
            static_assert(!std::is_void<destination_type>::value,
                          "Exact output collection must support the result geometry tag");
            destination_type result;
            geometry_converter<Coordinate, geometry::tag_t<source_type>, OutputPolicy>::apply(geometry, result);
            destination = std::move(result);
        }, source);
    }
};

template <typename Coordinate, typename OutputPolicy>
struct geometry_converter<Coordinate, box_tag, OutputPolicy>
{
    template <typename Source, typename Destination>
    static void apply(Source const& source, Destination& destination,
        std::enable_if_t<util::is_box<Destination>::value>* = nullptr)
    {
        using target = typename geometry::point_type<Destination>::type;
        if constexpr (std::is_same<target, exact_point<Coordinate>>::value)
        {
            static_box_detail::convert_box<Coordinate>(source, destination);
        }
        else
        {
            target first, second;
            point_converter<Coordinate, OutputPolicy>::apply(static_box_detail::borrowed_corner<min_corner>(source), first);
            point_converter<Coordinate, OutputPolicy>::apply(static_box_detail::borrowed_corner<max_corner>(source), second);
            geometry::set<min_corner, 0>(destination, geometry::get<0>(first));
            geometry::set<min_corner, 1>(destination, geometry::get<1>(first));
            geometry::set<max_corner, 0>(destination, geometry::get<0>(second));
            geometry::set<max_corner, 1>(destination, geometry::get<1>(second));
        }
    }

    template <typename Source, typename Destination>
    static void apply(Source const& source, Destination& destination,
        std::enable_if_t<!util::is_box<Destination>::value>* = nullptr)
    {
        convert_boundary(source, destination, geometry::tag_t<Destination>());
    }

private:
    template <typename Source, typename Ring>
    static void convert_boundary(Source const& source, Ring& ring, ring_tag)
    {
        auto const lower = make_point<Coordinate>(static_box_detail::corner<min_corner>(source));
        auto const upper = make_point<Coordinate>(static_box_detail::corner<max_corner>(source));
        using point = exact_point<Coordinate>;
        auto left = point::splice_coordinates(lower, 0, upper, 1);
        auto right = point::splice_coordinates(upper, 0, lower, 1);
        if (geometry::point_order<Ring>::value != clockwise) std::swap(left, right);
        std::array<point, 5> const corners{{lower, left, upper, right, lower}};
        geometry::clear(ring);
        auto const count = geometry::closure<Ring>::value == closed ? 5u : 4u;
        for (unsigned i = 0; i != count; ++i)
        {
            typename boost::range_value<Ring>::type output;
            point_converter<Coordinate, OutputPolicy>::apply(corners[i], output);
            traits::push_back<Ring>::apply(ring, std::move(output));
        }
    }

    template <typename Source, typename Polygon>
    static void convert_boundary(Source const& source, Polygon& polygon, polygon_tag)
    {
        geometry::clear(polygon);
        convert_boundary(source, exterior_ring(polygon), ring_tag());
    }
};

enum class box_relation { within, covered_by, disjoint, touches, overlaps };

template <typename Coordinate, box_relation Relation>
struct box_predicate
{
    using cs_tag = cartesian_tag;
private:
    template <typename Point, typename Box>
    static bool apply_impl(Point const& point, Box const& box, point_tag)
    {
        using predicates = static_box_predicates<Coordinate>;
        if (Relation == box_relation::within) return predicates::point_covered_by(point, box, true);
        if (Relation == box_relation::covered_by) return predicates::point_covered_by(point, box);
        if (Relation == box_relation::disjoint) return !predicates::point_covered_by(point, box);
        if (Relation == box_relation::touches) return predicates::point_covered_by(point, box)
            && !predicates::point_covered_by(point, box, true);
        return false;
    }
    template <typename First, typename Second>
    static bool apply_impl(First const& first, Second const& second, box_tag)
    {
        using predicates = static_box_predicates<Coordinate>;
        if (Relation == box_relation::within) return predicates::within(first, second);
        if (Relation == box_relation::covered_by) return predicates::covered_by(first, second);
        if (Relation == box_relation::disjoint) return predicates::disjoint(first, second);
        if (Relation == box_relation::touches) return predicates::touches(first, second);
        return predicates::overlaps(first, second);
    }
public:
    template <typename Geometry1, typename Geometry2>
    static bool apply(Geometry1 const& first, Geometry2 const& second)
    {
        return apply_impl(first, second, geometry::tag_t<Geometry1>());
    }
};

template <typename Coordinate>
struct box_topology
{
    static auto touches()
    {
        return box_predicate<Coordinate, box_relation::touches>();
    }

    static auto overlaps()
    {
        return box_predicate<Coordinate, box_relation::overlaps>();
    }
};

template <typename Coordinate, typename OutputPolicy = strategies::relate::approximate_output>
struct box_intersection
{
private:
    template <std::size_t Index, std::size_t Dimension, typename Point, typename Output>
    static void assign_bound(Point const& source, Output& output)
    {
        using target_point = typename geometry::point_type<Output>::type;
        if constexpr (std::is_same<target_point, exact_point<Coordinate>>::value)
        {
            static_box_detail::assign_bound<Index, Dimension, Coordinate>(source, output);
        }
        else
        {
            using target = geometry::coordinate_type_t<Output>;
            if constexpr (static_predicates_detail::exact_input<std::remove_cv_t<Point>>::value)
                geometry::set<Index, Dimension>(output, output_coordinate<target, OutputPolicy>(source, Dimension));
            else
                geometry::set<Index, Dimension>(output, geometry::util::numeric_cast<target>(geometry::get<Dimension>(source)));
        }
    }

    template <std::size_t Dimension, typename Box1, typename Box2, typename Output>
    static bool apply_dimension(Box1 const& first, Box2 const& second, Output& output)
    {
        auto const& max1 = static_box_detail::borrowed_corner<max_corner>(first);
        auto const& min2 = static_box_detail::borrowed_corner<min_corner>(second);
        if (static_box_detail::compare_coordinate<Dimension, Coordinate>(max1, min2) < 0)
        {
            return false;
        }
        auto const& max2 = static_box_detail::borrowed_corner<max_corner>(second);
        auto const& min1 = static_box_detail::borrowed_corner<min_corner>(first);
        if (static_box_detail::compare_coordinate<Dimension, Coordinate>(max2, min1) < 0)
        {
            return false;
        }

        if (static_box_detail::compare_coordinate<Dimension, Coordinate>(min1, min2) < 0)
            assign_bound<min_corner, Dimension>(min2, output);
        else assign_bound<min_corner, Dimension>(min1, output);
        if (static_box_detail::compare_coordinate<Dimension, Coordinate>(max1, max2) > 0)
            assign_bound<max_corner, Dimension>(max2, output);
        else assign_bound<max_corner, Dimension>(max1, output);
        return true;
    }

public:
    template <typename Box1, typename Box2, typename Output>
    static bool apply(Box1 const& first, Box2 const& second, Output& output)
    {
        // Commit each overlapping dimension without reading untouched output.
        return apply_dimension<0>(first, second, output)
            && apply_dimension<1>(first, second, output);
    }
};

template <typename Coordinate, typename OutputPolicy = strategies::relate::approximate_output,
          typename Source, typename Destination>
inline void append_output(Source const& source, Destination& destination,
    std::enable_if_t<!tuples::is_tuple<Destination>::value && !util::is_box<Destination>::value>* = nullptr)
{
    for (auto const& item : source)
    {
        typename boost::range_value<Destination>::type result;
        geometry_converter<Coordinate, geometry::tag_t<decltype(item)>, OutputPolicy>::apply(item, result);
        traits::push_back<Destination>::apply(destination, std::move(result));
    }
}

template <typename Coordinate, typename OutputPolicy = strategies::relate::approximate_output,
          typename Source, typename Destination>
inline void append_output(Source const& source, Destination& destination,
    std::enable_if_t<util::is_box<Destination>::value>* = nullptr)
{
    geometry_converter<Coordinate, box_tag, OutputPolicy>::apply(source, destination);
}

template <typename Coordinate, typename OutputPolicy, typename Source, typename Destination, std::size_t... I>
inline void append_tuple_output(Source const& source, Destination& destination, std::index_sequence<I...>)
{
    using expand = int[];
    (void)expand{0, (append_output<Coordinate, OutputPolicy>(tuples::get<I>(source), tuples::get<I>(destination)), 0)...};
}

template <typename Coordinate, typename OutputPolicy = strategies::relate::approximate_output,
          typename Source, typename Destination>
inline void append_output(Source const& source, Destination& destination,
    std::enable_if_t<tuples::is_tuple<Destination>::value>* = nullptr)
{
    append_tuple_output<Coordinate, OutputPolicy>(source, destination,
        std::make_index_sequence<tuples::size<Destination>::value>());
}

template <typename Coordinate, typename OutputPolicy = strategies::relate::approximate_output>
struct clipping
{
    using point_type = exact_point<Coordinate>;

    static auto get_equals_point_point_strategy() { return point_equal<Coordinate>(); }

    template <typename Source, typename Destination>
    static void convert_geometry(Source const& source, Destination& destination)
    {
        point_converter<Coordinate, OutputPolicy>::apply(source, destination);
    }

    template <bool Construct = true, typename Box>
    bool clip_segment(Box const& box, model::referring_segment<point_type>& segment,
                      bool& first_clipped, bool& second_clipped) const
    {
        auto a = segment.first;
        auto b = segment.second;
        auto const lower = make_point<Coordinate>(static_box_detail::corner<min_corner>(box));
        auto const upper = make_point<Coordinate>(static_box_detail::corner<max_corner>(box));
        bool clipped_a = false, clipped_b = false;
        if (!clip_halfplane<0, 1>(lower, segment.first, segment.second, a, b, clipped_a, clipped_b)
            || !clip_halfplane<0, -1>(upper, segment.first, segment.second, a, b, clipped_a, clipped_b)
            || !clip_halfplane<1, 1>(lower, segment.first, segment.second, a, b, clipped_a, clipped_b)
            || !clip_halfplane<1, -1>(upper, segment.first, segment.second, a, b, clipped_a, clipped_b)) return false;
        if (Construct)
        {
            first_clipped = clipped_a;
            second_clipped = clipped_b;
            if (clipped_a) segment.first = std::move(a);
            if (clipped_b) segment.second = std::move(b);
        }
        return true;
    }

private:
    template <std::size_t Dimension, int Direction>
    static bool clip_halfplane(point_type const& bound, point_type const& support_a,
                               point_type const& support_b, point_type& a, point_type& b,
                               bool& clipped_a, bool& clipped_b)
    {
        int const ca = Direction * static_compare_coordinate<Dimension>(a, bound);
        int const cb = Direction * static_compare_coordinate<Dimension>(b, bound);
        if (ca < 0 && cb < 0) return false;
        if (ca >= 0 && cb >= 0) return true;
        if (ca == 0) { b = a; clipped_b = true; return true; }
        if (cb == 0) { a = b; clipped_a = true; return true; }
        point_type first, second;
        Coordinate native;
        if (bound.native_coordinate(Dimension, native))
        {
            first = Dimension == 0 ? point_type(native, 0) : point_type(0, native);
            second = Dimension == 0 ? point_type(native, 1) : point_type(1, native);
        }
        else
        {
            point_type const zero(0, 0), one(1, 1);
            first = Dimension == 0 ? point_type::splice_coordinates(bound, 0, zero, 1)
                                   : point_type::splice_coordinates(zero, 0, bound, 1);
            second = Dimension == 0 ? point_type::splice_coordinates(bound, 0, one, 1)
                                    : point_type::splice_coordinates(one, 0, bound, 1);
        }
        auto intersection = construct_intersection(support_a, support_b, first, second);
        if (ca < 0) { a = std::move(intersection); clipped_a = true; }
        else { b = std::move(intersection); clipped_b = true; }
        return true;
    }

public:

    template <typename Linestring, typename Iterator>
    static void apply(Linestring& line, Iterator& output)
    {
        if (boost::size(line) > 1)
        {
            *output++ = line;
        }
        geometry::clear(line);
    }
};

template <typename Coordinate>
struct segment_box
{
    static auto get_disjoint_point_box_strategy()
    {
        return box_predicate<Coordinate, box_relation::covered_by>();
    }

    template <typename Segment, typename Box>
    static bool apply(Segment const& segment, Box const& box)
    {
        auto const& a = endpoint<0>(segment, 0);
        auto const& b = endpoint<1>(segment, 0);
        auto const& lower = static_box_detail::borrowed_corner<min_corner>(box);
        auto const& upper = static_box_detail::borrowed_corner<max_corner>(box);
        auto const compare_x = [](auto const& first, auto const& second)
        { return static_box_detail::compare_coordinate<0, Coordinate>(first, second); };
        auto const compare_y = [](auto const& first, auto const& second)
        { return static_box_detail::compare_coordinate<1, Coordinate>(first, second); };
        if (compare_x(lower, upper) > 0 || compare_y(lower, upper) > 0
            || (compare_x(a, lower) < 0 && compare_x(b, lower) < 0)
            || (compare_x(a, upper) > 0 && compare_x(b, upper) > 0)
            || (compare_y(a, lower) < 0 && compare_y(b, lower) < 0)
            || (compare_y(a, upper) > 0 && compare_y(b, upper) > 0))
            return true;

        if constexpr (std::is_floating_point<geometry::coordinate_type_t<Box>>::value)
        {
            if (static_box_detail::infinity<0>(lower, std::true_type{})
                || static_box_detail::infinity<1>(lower, std::true_type{})
                || static_box_detail::infinity<0>(upper, std::true_type{})
                || static_box_detail::infinity<1>(upper, std::true_type{}))
            {
                int const dx = compare_x(b, a), dy = compare_y(b, a);
                if (dx == 0 || dy == 0) return false;
                // The linear orientation functional reaches its extrema at
                // opposite corners. An unbounded extremum has a known sign.
                if (dy > 0)
                    return dx > 0 ? unbounded_disjoint<max_corner, min_corner>(a, b, box)
                                  : unbounded_disjoint<max_corner, max_corner>(a, b, box);
                return dx > 0 ? unbounded_disjoint<min_corner, min_corner>(a, b, box)
                              : unbounded_disjoint<min_corner, max_corner>(a, b, box);
            }
        }

        // The remaining separating axis is normal to the supporting line.
        // Strict signs retain closed-boundary contacts, including degeneracies.
        int const side = static_box_detail::corner_orientation<Coordinate>(a, b,
            static_box_detail::indexed_point<Box, min_corner, min_corner>{box});
        return side != 0
            && static_box_detail::corner_orientation<Coordinate>(a, b,
                static_box_detail::indexed_point<Box, max_corner, min_corner>{box}) == side
            && static_box_detail::corner_orientation<Coordinate>(a, b,
                static_box_detail::indexed_point<Box, max_corner, max_corner>{box}) == side
            && static_box_detail::corner_orientation<Coordinate>(a, b,
                static_box_detail::indexed_point<Box, min_corner, max_corner>{box}) == side;
    }

private:
    template <std::size_t Index, typename Segment>
    static auto endpoint(Segment const& segment, int)
        -> std::enable_if_t<std::is_pointer<decltype(segment.first)>::value, decltype(*segment.first)>
    {
        // segment_iterator yields pointing_segment: retain the actual implicit
        // point instead of extracting its coordinates into a rational snapshot.
        BOOST_ASSERT(segment.first != nullptr && segment.second != nullptr);
        return Index == 0 ? *segment.first : *segment.second;
    }

    template <std::size_t Index, typename Segment>
    static decltype(auto) endpoint(Segment const& segment, long)
    { return static_box_detail::borrowed_corner<Index>(segment); }

    template <std::size_t X, std::size_t Y, typename A, typename B, typename Box>
    static int extreme_side(A const& a, B const& b, Box const& box, int infinite_sign)
    {
        static_box_detail::indexed_point<Box, X, Y> const corner{box};
        return static_box_detail::infinity<0>(corner, std::true_type{})
            || static_box_detail::infinity<1>(corner, std::true_type{})
            ? infinite_sign : static_box_detail::corner_orientation<Coordinate>(a, b, corner);
    }

    template <std::size_t MinX, std::size_t MinY, typename A, typename B, typename Box>
    static bool unbounded_disjoint(A const& a, B const& b, Box const& box)
    {
        return extreme_side<MinX, MinY>(a, b, box, -1) > 0
            || extreme_side<1 - MinX, 1 - MinY>(a, b, box, 1) < 0;
    }
};

}} // namespace detail::exact_cartesian

namespace relate
{

/*!
\brief Exact Cartesian predicates and implicit straight-segment intersections.
\tparam Coordinate Native input coordinate type: integral, float or double.
\tparam OutputPolicy approximate_output (default), a nearest policy, bounded_ulp_output<N>,
        or lean_bounded_ulp_output<N>.
\details Finite binary floating-point inputs denote their exact binary values.
         Input promotion must preserve the complete source coordinate type.
         Constructed points use forward error filters
         with exact fallback. Overlay keeps exact
         working geometries until the complete result is converted to output.
\note Requires two-dimensional Cartesian input, IEEE binary arithmetic, and no
      fast-math/reassociation. Floating overflow and underflow are out of scope.
      Round-to-nearest and controlled contraction are required.
      Long double input/output is not supported. approximate_output emits the
      currently stored approximation without output-only refinement.
      nearest_output and deferred_nearest_output keep a one-ULP approximation
      and round floating output to nearest, ties to even, at emission.
      eager_nearest_output rounds at construction. These policies omit stored
      coordinate error fields. approximate_output retains forward error bounds.
      bounded_ulp_output<N, Retry> certifies binary64 intersections within N
      representable steps of nearest, ties to even, and retains predicate
      error bounds. lean_bounded_ulp_output has the same emission contract;
      its root-bound constructions omit both error fields and its filters
      specialize on N. Both bounded policies accept a third template
      argument to enable or disable selective determinant tightening after
      a miss (enabled by default at every budget). If disabled, Retry
      (default false) enables only an FMA
      quotient-residual retry. The coordinate guarantee is unchanged.
      Uncertified axes, non-double Coordinate types, float output and legacy
      exact_point conversion use nearest rounding for floating constructions.
      Integral output uses BG numeric_cast conversion, which truncates floating
      values toward zero, not the displacement rounding of native integer segment
      construction. No policy guarantees validity after output conversion;
      no topology decision is made on the converted result.
*/
template <typename Coordinate, typename OutputPolicy = approximate_output>
class cartesian_exact : public cartesian<double>
{
    static_assert(std::is_same<OutputPolicy, approximate_output>::value
        || std::is_same<OutputPolicy, nearest_output>::value
        || std::is_same<OutputPolicy, eager_nearest_output>::value
        || std::is_same<OutputPolicy, deferred_nearest_output>::value
        || strategies::detail::exact_cartesian::is_bounded_ulp_output<OutputPolicy>::value,
        "Unknown exact Cartesian output policy");
    using primitives = strategies::detail::exact_cartesian::point_equal<Coordinate>;

public:
    using base_type = cartesian<double>;
    using output_policy = OutputPolicy;
    using operation_scope = strategies::detail::exact_cartesian::point_cache_scope<Coordinate>;
    using base_type::area;
    using base_type::disjoint;

    template <typename A, typename B, typename Apply>
    decltype(auto) with_sources(A const& a, B const& b, Apply&& apply) const;

    template <typename Geometry, typename Box>
    static auto envelope(Geometry const&, Box const&)
    {
        using namespace strategies::detail::exact_cartesian;
        using tag = geometry::tag_t<Geometry>;
        using selected = std::conditional_t<std::is_same<tag, point_tag>::value, static_envelope_point<Coordinate>,
            std::conditional_t<std::is_same<tag, box_tag>::value, static_envelope_box<Coordinate>,
            std::conditional_t<std::is_same<tag, segment_tag>::value, static_envelope_segment<Coordinate>,
            std::conditional_t<std::is_same<tag, multi_point_tag>::value
                || std::is_same<tag, linestring_tag>::value || std::is_same<tag, ring_tag>::value
                || std::is_same<tag, polygon_tag>::value, static_envelope_range<Coordinate>, static_envelope_boxes<Coordinate>>>>>;
        return selected();
    }

    template <typename Box, typename Geometry>
    static auto expand(Box const&, Geometry const&)
    {
        using namespace strategies::detail::exact_cartesian;
        using tag = geometry::tag_t<Geometry>;
        using selected = std::conditional_t<std::is_same<tag, point_tag>::value, static_expand_point<Coordinate>,
            std::conditional_t<std::is_same<tag, segment_tag>::value, static_expand_segment<Coordinate>, static_expand_box<Coordinate>>>;
        return selected();
    }

    template <typename Point, typename Box>
    static auto within(Point const&, Box const&,
        std::enable_if_t<util::is_pointlike<Point>::value && util::is_box<Box>::value>* = nullptr)
    {
        return strategies::detail::exact_cartesian::box_predicate<Coordinate,
            strategies::detail::exact_cartesian::box_relation::within>();
    }

    template <typename Box1, typename Box2>
    static auto within(Box1 const&, Box2 const&,
        std::enable_if_t<util::is_box<Box1>::value && util::is_box<Box2>::value>* = nullptr)
    {
        return strategies::detail::exact_cartesian::box_predicate<Coordinate,
            strategies::detail::exact_cartesian::box_relation::within>();
    }

    template <typename Point, typename Box>
    static auto covered_by(Point const&, Box const&,
        std::enable_if_t<util::is_pointlike<Point>::value && util::is_box<Box>::value>* = nullptr)
    {
        return strategies::detail::exact_cartesian::box_predicate<Coordinate,
            strategies::detail::exact_cartesian::box_relation::covered_by>();
    }

    template <typename Box1, typename Box2>
    static auto covered_by(Box1 const&, Box2 const&,
        std::enable_if_t<util::is_box<Box1>::value && util::is_box<Box2>::value>* = nullptr)
    {
        return strategies::detail::exact_cartesian::box_predicate<Coordinate,
            strategies::detail::exact_cartesian::box_relation::covered_by>();
    }

    template <typename Box1, typename Box2>
    static auto disjoint(Box1 const&, Box2 const&,
        std::enable_if_t<util::is_box<Box1>::value && util::is_box<Box2>::value>* = nullptr)
    {
        return strategies::detail::exact_cartesian::box_predicate<Coordinate,
            strategies::detail::exact_cartesian::box_relation::disjoint>();
    }

    template <typename Segment, typename Box>
    static auto disjoint(Segment const&, Box const&,
        std::enable_if_t<util::is_segment<Segment>::value && util::is_box<Box>::value>* = nullptr)
    {
        return strategies::detail::exact_cartesian::segment_box<Coordinate>();
    }
    template <typename Point>
    using intersection_point_type = strategies::detail::exact_cartesian::exact_point<Coordinate>;

    using segment_position_type = strategies::detail::exact_cartesian::segment_position_tag;

    template <typename Point>
    using comparable_distance_type = strategies::detail::exact_cartesian::squared_distance_value<Coordinate>;

    template <typename Point>
    using box_type = strategies::detail::exact_cartesian::bounding_box;

    template <typename Point>
    using geometry_box_type = model::box<strategies::detail::exact_cartesian::exact_point<Coordinate>>;

    template <typename Point>
    using section_point_type = Point;

    template <typename Iterator>
    using section_iterator_point_type = std::conditional_t<
        strategies::detail::exact_cartesian::fixed_native_supported<
            geometry::coordinate_type_t<typename std::iterator_traits<Iterator>::value_type>>::value,
        strategies::detail::exact_cartesian::iterator_point_view<Iterator>,
        typename std::iterator_traits<Iterator>::value_type>;

#if defined(BOOST_GEOMETRY_EXACT_DISABLE_ALGORITHM_OPTIMIZATIONS)
    static constexpr bool exact_predicates = false;
#else
    static constexpr bool exact_predicates = true;
#endif
    using exact_cluster_compare = strategies::detail::exact_cartesian::point_compare
        <Coordinate, strategy::compare::less>;

    template <typename Point>
    using intersection_result_type = strategies::detail::exact_cartesian::segment_intersection_result
        <intersection_point_type<Point>>;

    template <typename Geometry>
    using overlay_geometry_type = typename strategies::detail::exact_cartesian::working_geometry
        <Coordinate, Geometry>::type;

    // Box failure preserves untouched coordinates; its primitive writes in place.
    template <typename Collection>
    using overlay_output_type = std::conditional_t<util::is_box<Collection>::value,
        Collection, typename strategies::detail::exact_cartesian::working_output
            <Coordinate, Collection>::type>;

    static auto side() { return strategies::detail::exact_cartesian::side<Coordinate>(); }
    static auto topology() { return strategies::detail::exact_cartesian::topology<Coordinate>(); }
    template <typename A, typename B>
    static auto comparable_distance(A const&, B const&)
    { return strategies::detail::exact_cartesian::comparable_squared_distance<Coordinate>(); }
    static auto segment_position() { return strategies::detail::exact_cartesian::segment_position<Coordinate>(); }
    static auto ring_topology() { return strategies::detail::exact_cartesian::ring_topology<Coordinate>(); }
    static auto section() { return strategies::detail::exact_cartesian::static_section<Coordinate>(); }
    static auto box_intersection() { return strategies::detail::exact_cartesian::box_intersection<Coordinate, OutputPolicy>(); }
    static auto box_topology() { return strategies::detail::exact_cartesian::box_topology<Coordinate>(); }
    template <typename Point, typename Box>
    static auto clipping(Box const&) { return strategies::detail::exact_cartesian::clipping<Coordinate, OutputPolicy>(); }
    static auto relate()
    {
        return strategies::detail::exact_cartesian::output_segment_intersection<Coordinate, OutputPolicy>();
    }

    template <typename A, typename B>
    static auto relate(A const&, B const&,
        std::enable_if_t<util::is_pointlike<A>::value && util::is_pointlike<B>::value>* = nullptr)
    {
        return primitives();
    }

    template <typename A, typename B>
    static auto relate(A const&, B const&,
        std::enable_if_t<util::is_pointlike<A>::value
            && (util::is_linear<B>::value || util::is_polygonal<B>::value)>* = nullptr)
    {
        return strategies::detail::exact_cartesian::winding<Coordinate>();
    }

    template <typename ComparePolicy, typename EqualsPolicy>
    using compare_type = strategies::detail::exact_cartesian::point_compare<Coordinate, ComparePolicy>;

    template <typename Geometry>
    using ring_measure_type = strategies::detail::exact_cartesian::ring_metadata<Coordinate>;

    template <typename Source, typename Destination>
    static void convert_geometry(Source const& source, Destination& destination,
        std::enable_if_t<!std::is_same<Source, Destination>::value || std::is_array<Source>::value>* = nullptr)
    {
        strategies::detail::exact_cartesian::geometry_converter<Coordinate, geometry::tag_t<Source>, OutputPolicy>
            ::apply(source, destination);
    }

    template <typename Geometry>
    static void convert_geometry(Geometry const& source, Geometry& destination,
        std::enable_if_t<!std::is_array<Geometry>::value>* = nullptr)
    {
        destination = source;
    }

    template <typename Source, typename Destination>
    static void append_output(Source const& source, Destination& destination)
    {
        strategies::detail::exact_cartesian::append_output<Coordinate, OutputPolicy>(source, destination);
    }
};

} // namespace relate

}}} // namespace boost::geometry::strategies

#include <boost/geometry/strategies/cartesian/exact/borrowed_strategy.hpp>

namespace boost { namespace geometry { namespace strategies { namespace relate
{
template <typename Coordinate, typename OutputPolicy>
template <typename A, typename B, typename Apply>
decltype(auto) cartesian_exact<Coordinate, OutputPolicy>::with_sources(
    A const& a, B const& b, Apply&& apply) const
{
    namespace exact = strategies::detail::exact_cartesian;
    if constexpr (exact::fixed_native_supported<geometry::coordinate_type_t<A>>::value
                  && exact::fixed_native_supported<geometry::coordinate_type_t<B>>::value)
    {
        exact::check_coordinate_type<Coordinate, A>();
        exact::check_coordinate_type<Coordinate, B>();
        using root = exact::borrowed_root<Coordinate, OutputPolicy, A, B>;
        exact::borrowed_context<root> context(a, b);
        exact::borrowed_strategy<root> bound(context);
        return std::forward<Apply>(apply)(bound);
    }
    else return std::forward<Apply>(apply)(*this);
}
}}}}

#endif
