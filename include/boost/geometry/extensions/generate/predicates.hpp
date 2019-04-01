// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2019 Tinko Bartels, Berlin, Germany.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_EXTENSIONS_GENERATE_PREDICATES_HPP
#define BOOST_GEOMETRY_EXTENSIONS_GENERATE_PREDICATES_HPP

#include <boost/geometry/algorithms/covered_by.hpp>
#include <boost/geometry/algorithms/crosses.hpp>
#include <boost/geometry/algorithms/disjoint.hpp>
#include <boost/geometry/algorithms/equals.hpp>
#include <boost/geometry/algorithms/intersects.hpp>
#include <boost/geometry/algorithms/is_simple.hpp>
#include <boost/geometry/algorithms/overlaps.hpp>
#include <boost/geometry/algorithms/touches.hpp>
#include <boost/geometry/algorithms/within.hpp>

namespace boost { namespace geometry { namespace generate
{

namespace detail { namespace predicates
{

    template <typename UnaryPredicate, typename Geometry>
    class satisfies
    {
        public:
        typedef Geometry geo;
        satisfies(UnaryPredicate const& pred) : pred(pred) {}
        satisfies(satisfies<UnaryPredicate, geo> const& s) : pred(s.pred) {}
        bool operator()(geo const& g) const { return pred(g); }
        bool operator()(geo const* g) const { return pred(*g); }
        private:
        UnaryPredicate pred;
    };

}} // namespace detail::predicates

    template<typename Geometry, typename UnaryPredicate>
    detail::predicates::satisfies<UnaryPredicate, Geometry>
    satisfies(UnaryPredicate const& pred) { return detail::predicates::satisfies<UnaryPredicate, Geometry>(pred); }

    template<typename Geometry1, typename Geometry2>
    auto within(const Geometry2& g2) { return satisfies<Geometry1>([g2](Geometry1 const& g1) { return geometry::within(g1,g2); }); }
    template<typename Geometry1, typename Geometry2>
    auto touches(const Geometry2& g2) { return satisfies<Geometry1>([g2](Geometry1 const& g1) { return geometry::touches(g1,g2); }); }
    template<typename Geometry1>
    auto touches() { return satisfies<Geometry1>([](Geometry1 const& g1) { return geometry::touches(g1); }); }
    template<typename Geometry1, typename Geometry2>
    auto equals(const Geometry2& g2) { return satisfies<Geometry1>([g2](Geometry1 const& g1) { return geometry::equals(g1,g2); }); }
    template<typename Geometry1, typename Geometry2>
    auto disjoint(const Geometry2& g2) { return satisfies<Geometry1>([g2](Geometry1 const& g1) { return geometry::disjoint(g1,g2); }); }
    template<typename Geometry1, typename Geometry2>
    auto intersects(const Geometry2& g2) { return satisfies<Geometry1>([g2](Geometry1 const& g1) { return geometry::intersects(g1,g2); }); }
    template<typename Geometry1>
    auto intersects() { return satisfies<Geometry1>([](Geometry1 const& g1) { return geometry::intersects(g1); }); }
    template<typename Geometry1, typename Geometry2>
    auto covered_by(const Geometry2& g2) { return satisfies<Geometry1>([g2](Geometry1 const& g1) { return geometry::covered_by(g1,g2); }); }
    template<typename Geometry1>
    auto is_simple() { return satisfies<Geometry1>([](Geometry1 const& g1) { return geometry::is_simple(g1); }); }
    template<typename Geometry1, typename Geometry2>
    auto overlaps(const Geometry2& g2) { return satisfies<Geometry1>([g2](Geometry1 const& g1) { return geometry::overlaps(g1,g2); }); }
    template<typename Geometry1, typename Geometry2>
    auto crosses(const Geometry2& g2) { return satisfies<Geometry1>([g2](Geometry1 const& g1) { return geometry::crosses(g1,g2); }); }
    template<typename Predicate1, typename Predicate2>
    auto operator&&(const Predicate1& pred1, const Predicate2& pred2) { return satisfies<typename Predicate1::geo>([pred1,pred2](typename Predicate1::geo const& g){return pred1(g) && pred2(g);}); }
    template<typename Predicate1, typename Predicate2>
    auto operator||(const Predicate1& pred1, const Predicate2& pred2) { return satisfies<typename Predicate1::geo>([pred1,pred2](typename Predicate1::geo const& g) {return pred1(g) || pred2(g);}); }
    template<typename Predicate>
    auto operator!(const Predicate& pred) { return satisfies<typename Predicate::geo>([pred](typename Predicate::geo const& g) {return !pred(g);}); }

}}} // namespace boost::geometry::generate

#endif // BOOST_GEOMETRY_EXTENSIONS_GENERATE_PREDICATES_HPP