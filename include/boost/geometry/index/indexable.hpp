// Boost.Geometry Index
//
// Copyright (c) 2011-2019 Adam Wulkiewicz, Lodz, Poland.
//
// This file was modified by Oracle on 2020.
// Modifications copyright (c) 2020 Oracle and/or its affiliates.
// Contributed and/or modified by Adam Wulkiewicz, on behalf of Oracle
//
// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_INDEX_INDEXABLE_HPP
#define BOOST_GEOMETRY_INDEX_INDEXABLE_HPP

#include <boost/geometry/core/static_assert.hpp>

#include <boost/geometry/index/detail/is_indexable.hpp>

#include <boost/geometry/util/type_traits.hpp>

#include <tuple>

namespace boost { namespace geometry { namespace index { namespace detail
{

template <typename From, typename To>
struct is_referenceable
    : std::is_same
        <
            typename util::remove_cref<From>::type,
            typename util::remove_cref<To>::type
        >
{};

template <typename Indexable, typename V>
inline Indexable const& indexable_prevent_any_type(V const& )
{
    BOOST_GEOMETRY_STATIC_ASSERT_FALSE("Unexpected type.", V);
    return Indexable();
}

/*!
\brief The function object extracting Indexable from Value.

It translates Value object to Indexable object. The default version handles Values which are Indexables.
This template is also specialized for std::pair<Indexable, T2> and
std::tuple<Indexable, ...>.

\tparam Value       The Value type which may be translated directly to the Indexable.
\tparam IsIndexable If true, the const reference to Value is returned.
*/
template <typename Value, bool IsIndexable = is_indexable<Value>::value>
struct indexable
{
    BOOST_GEOMETRY_STATIC_ASSERT(
        (detail::is_indexable<Value>::value),
        "Value has to be an Indexable.",
        Value);

    /*! \brief The type of result returned by function object. */
    typedef Value const& result_type;

    /*!
    \brief Return indexable extracted from the value.

    \param v The value.
    \return The indexable.
    */
    inline result_type operator()(Value const& v) const
    {
        return v;
    }

    /*!
    \brief Prevent reference to temporary for types convertible to Value.
    */
    template <typename V>
    inline result_type operator()(V const& v) const
    {
        return indexable_prevent_any_type<Value>(v);
    }
};

/*!
\brief The function object extracting Indexable from Value.

This specialization translates from std::pair<Indexable, T2>.

\tparam Indexable       The Indexable type.
\tparam Second          The second type.
*/
template <typename Indexable, typename Second>
struct indexable<std::pair<Indexable, Second>, false>
{
    typedef std::pair<Indexable, Second> value_type;

    BOOST_GEOMETRY_STATIC_ASSERT(
        (detail::is_indexable<Indexable>::value),
        "The first type of std::pair has to be an Indexable.",
        Indexable);

    /*! \brief The type of result returned by function object. */
    typedef Indexable const& result_type;

    /*!
    \brief Return indexable extracted from the value.

    \param v The value.
    \return The indexable.
    */
    inline result_type operator()(value_type const& v) const
    {
        return v.first;
    }

    /*!
    \brief Return indexable extracted from compatible type different than value_type.

    \param v The value.
    \return The indexable.
    */
    template <typename I, typename S>
    inline result_type operator()(std::pair<I, S> const& v) const
    {
        BOOST_GEOMETRY_STATIC_ASSERT(
            (is_referenceable<I, result_type>::value),
            "Unexpected type.",
            std::pair<I, S>);
        return v.first;
    }

    /*!
    \brief Prevent reference to temporary for types convertible to Value.
    */
    template <typename V>
    inline result_type operator()(V const& v) const
    {
        return indexable_prevent_any_type<Indexable>(v);
    }
};

}}}} // namespace boost::geometry::index::detail

namespace boost { namespace geometry { namespace index { namespace detail {

/*!
\brief The function object extracting Indexable from Value.

This specialization translates from std::tuple<Indexable, Args...>.
It's defined if the compiler supports tuples and variadic templates.

\tparam Indexable   The Indexable type.
*/
template <typename Indexable, typename ...Args>
struct indexable<std::tuple<Indexable, Args...>, false>
{
    typedef std::tuple<Indexable, Args...> value_type;

    BOOST_GEOMETRY_STATIC_ASSERT(
        (detail::is_indexable<Indexable>::value),
        "The first type of std::tuple has to be an Indexable.",
        Indexable);

    /*! \brief The type of result returned by function object. */
    typedef Indexable const& result_type;

    /*!
    \brief Return indexable extracted from the value.

    \param v The value.
    \return The indexable.
    */
    result_type operator()(value_type const& v) const
    {
        return std::get<0>(v);
    }

    /*!
    \brief Return indexable extracted from compatible type different than value_type.

    \param v The value.
    \return The indexable.
    */
    template <typename I, typename ...A>
    inline result_type operator()(std::tuple<I, A...> const& v) const
    {
        BOOST_GEOMETRY_STATIC_ASSERT(
            (is_referenceable<I, result_type>::value),
            "Unexpected type.",
            std::tuple<I, A...>);
        return std::get<0>(v);
    }

    /*!
    \brief Prevent reference to temporary for types convertible to Value.
    */
    template <typename V>
    inline result_type operator()(V const& v) const
    {
        return indexable_prevent_any_type<Indexable>(v);
    }
};

}}}} // namespace boost::geometry::index::detail


namespace boost { namespace geometry { namespace index {

/*!
\brief The function object extracting Indexable from Value.

It translates Value object to Indexable object. By default, it can handle Values which are Indexables,
std::pair<Indexable, T2> and std::tuple<Indexable, ...>.

\tparam Value       The Value type which may be translated directly to the Indexable.
*/
template <typename Value>
struct indexable
    : detail::indexable<Value>
{
    /*! \brief The type of result returned by function object. It should be const Indexable reference. */
    typedef typename detail::indexable<Value>::result_type result_type;

    /*!
    \brief Return indexable extracted from the value.

    \param v The value.
    \return The indexable.
    */
    inline result_type operator()(Value const& v) const
    {
        return detail::indexable<Value>::operator()(v);
    }

    /*!
    \brief Return indexable extracted from the value. Overload for types
           compatible with Value but different yet holding referenceable
           Indexable, e.g. tuple containing a reference.

    \param v The value.
    \return The indexable.
    */
    template <typename V>
    inline result_type operator()(V const& v) const
    {
        return detail::indexable<Value>::operator()(v);
    }
};

}}} // namespace boost::geometry::index

#endif // BOOST_GEOMETRY_INDEX_INDEXABLE_HPP
