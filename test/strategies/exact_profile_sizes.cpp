// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.
#include <boost/geometry.hpp>
#include <boost/geometry/strategies/relate/cartesian_exact.hpp>
#include <boost/geometry/algorithms/detail/overlay/traversal_info.hpp>
#include <cstdint>
#include <iostream>

template <typename C>
void sizes(char const* name)
{
    namespace bg = boost::geometry;
    namespace ec = bg::strategies::detail::exact_cartesian;
    namespace services = bg::strategies::relate::services;
    using point = bg::model::d2::point_xy<C>;
    using exact = bg::strategies::relate::cartesian_exact<C>;
    using ordinary = bg::strategies::relate::cartesian<>;
    using implicit = typename services::intersection_point_type<exact,point>::type;
    using distance = typename services::comparable_distance_type<exact,point>::type;
    using exact_turn = bg::detail::overlay::traversal_turn_info<implicit,
        typename services::segment_ratio_type<exact,point>::type,distance>;
    using ordinary_turn = bg::detail::overlay::traversal_turn_info<point,
        typename services::segment_ratio_type<ordinary,point>::type,
        typename services::comparable_distance_type<ordinary,point>::type>;
    std::cout << name << ',' << sizeof(point) << ',' << sizeof(ec::native_point_view<C, point>) << ','
        << sizeof(implicit) << ',' << sizeof(distance) << ',' << sizeof(exact_turn) << ','
        << sizeof(ordinary_turn) << ',' << sizeof(typename ec::operation_detail::storage_for<C>::native_entry) << '\n';
}

int main()
{
    std::cout << "coordinate,native_point,native_view,implicit_point,squared_distance,exact_traversal_turn,ordinary_traversal_turn,native_precision_entry\n";
    sizes<double>("double"); sizes<std::int64_t>("int64");
}
