// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.

#ifndef BOOST_GEOMETRY_TEST_STRATEGIES_EXACT_TEST_STRATEGY_HPP
#define BOOST_GEOMETRY_TEST_STRATEGIES_EXACT_TEST_STRATEGY_HPP

#if defined(BOOST_GEOMETRY_TEST_EXACT_CARTESIAN)
#include <cmath>
#include <cstddef>
#include <iostream>
#include <type_traits>
#include <boost/geometry/core/coordinate_dimension.hpp>
#include <boost/geometry/core/coordinate_type.hpp>
#include <boost/geometry/core/cs.hpp>
#include <boost/geometry/algorithms/detail/visit.hpp>
#include <boost/geometry/algorithms/for_each.hpp>
#include <boost/geometry/strategies/relate/cartesian_exact.hpp>
#endif

namespace bg_test
{

struct default_strategy_policy
{
    static constexpr bool is_default = true;

    template <typename Function>
    decltype(auto) apply(Function const& function) const
    {
        return function();
    }
};

#if defined(BOOST_GEOMETRY_TEST_EXACT_CARTESIAN)

struct exact_strategy_statistics
{
    std::size_t cases = 0;
    std::size_t nonfinite = 0;

    ~exact_strategy_statistics()
    {
        std::cout << "\nExact Cartesian cases attempted: " << cases << '\n';
        std::cout << "Exact Cartesian cases skipped (nonfinite): " << nonfinite << '\n';
    }
};

inline exact_strategy_statistics& exact_statistics()
{
    static exact_strategy_statistics statistics;
    return statistics;
}

template <typename Strategy>
struct explicit_strategy_policy
{
    static constexpr bool is_default = false;

    template <typename Function>
    decltype(auto) apply(Function const& function) const
    {
        return function(Strategy());
    }
};

template <typename Geometry1, typename Geometry2>
struct supports_exact_strategy
{
    using coordinate_type = boost::geometry::coordinate_type_t<Geometry1>;
    // The kernel requires two-dimensional inputs with one native coordinate type.
    static constexpr bool value =
        std::is_same<boost::geometry::cs_tag_t<Geometry1>, boost::geometry::cartesian_tag>::value
        && std::is_same<boost::geometry::cs_tag_t<Geometry2>, boost::geometry::cartesian_tag>::value
        && boost::geometry::dimension<Geometry1>::value == 2
        && boost::geometry::dimension<Geometry2>::value == 2
        && std::is_same<coordinate_type, boost::geometry::coordinate_type_t<Geometry2>>::value
        && (std::is_integral<coordinate_type>::value
            || std::is_same<coordinate_type, float>::value
            || std::is_same<coordinate_type, double>::value);
};

template <typename Geometry>
bool finite_indexed_geometry(Geometry const& geometry)
{
    return std::isfinite(boost::geometry::get<0, 0>(geometry))
        && std::isfinite(boost::geometry::get<0, 1>(geometry))
        && std::isfinite(boost::geometry::get<1, 0>(geometry))
        && std::isfinite(boost::geometry::get<1, 1>(geometry));
}

template <typename Geometry, typename Tag>
bool finite_leaf_geometry(Geometry const& geometry, Tag)
{
    return boost::geometry::all_points_of(geometry, [](auto const& point)
    {
        return std::isfinite(boost::geometry::get<0>(point))
            && std::isfinite(boost::geometry::get<1>(point));
    });
}

template <typename Geometry>
bool finite_leaf_geometry(Geometry const& geometry, boost::geometry::box_tag)
{
    return finite_indexed_geometry(geometry);
}

template <typename Geometry>
bool finite_leaf_geometry(Geometry const& geometry, boost::geometry::segment_tag)
{
    return finite_indexed_geometry(geometry);
}

template <typename Geometry>
bool finite_geometry(Geometry const& geometry)
{
    bool finite = true;
    boost::geometry::detail::visit_breadth_first([&](auto const& leaf)
    {
        using leaf_type = typename std::decay<decltype(leaf)>::type;
        finite = finite_leaf_geometry(leaf, boost::geometry::tag_t<leaf_type>());
        return finite;
    }, geometry);
    return finite;
}

template <typename Geometry1, typename Geometry2, typename Function>
void for_exact_strategy(Geometry1 const& geometry1, Geometry2 const& geometry2,
                        Function const& function, std::true_type)
{
    using strategy = boost::geometry::strategies::relate::cartesian_exact
        <boost::geometry::coordinate_type_t<Geometry1>>;
    auto& statistics = exact_statistics();
    // Non-finite inputs are outside the exact kernel's domain; default checks still run.
    if (!finite_geometry(geometry1) || !finite_geometry(geometry2))
    {
        ++statistics.nonfinite;
        return;
    }
    ++statistics.cases;
    function(explicit_strategy_policy<strategy>());
}

template <typename Geometry1, typename Geometry2, typename Function>
void for_exact_strategy(Geometry1 const&, Geometry2 const&, Function const&, std::false_type)
{}

template <typename Geometry1, typename Geometry2, typename Function>
void for_exact_strategy(Geometry1 const& geometry1, Geometry2 const& geometry2,
                        Function const& function)
{
    for_exact_strategy(geometry1, geometry2, function,
        std::integral_constant<bool, supports_exact_strategy<Geometry1, Geometry2>::value>());
}

#else

template <typename Geometry1, typename Geometry2, typename Function>
void for_exact_strategy(Geometry1 const&, Geometry2 const&, Function const&)
{}

#endif

} // namespace bg_test

#endif // BOOST_GEOMETRY_TEST_STRATEGIES_EXACT_TEST_STRATEGY_HPP
