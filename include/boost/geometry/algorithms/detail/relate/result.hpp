// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2007-2015 Barend Gehrels, Amsterdam, the Netherlands.
// Copyright (c) 2017-2023 Adam Wulkiewicz, Lodz, Poland.

// This file was modified by Oracle on 2013-2022.
// Modifications copyright (c) 2013-2022 Oracle and/or its affiliates.
// Contributed and/or modified by Adam Wulkiewicz, on behalf of Oracle

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_ALGORITHMS_DETAIL_RELATE_RESULT_HPP
#define BOOST_GEOMETRY_ALGORITHMS_DETAIL_RELATE_RESULT_HPP

#include <algorithm>
#include <cstddef>
#include <cstring>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>

#include <boost/throw_exception.hpp>

#include <boost/geometry/core/assert.hpp>
#include <boost/geometry/core/coordinate_dimension.hpp>
#include <boost/geometry/core/exception.hpp>
#include <boost/geometry/core/static_assert.hpp>
#include <boost/geometry/util/sequence.hpp>

namespace boost { namespace geometry {

#ifndef DOXYGEN_NO_DETAIL
namespace detail { namespace relate {

enum field { interior = 0, boundary = 1, exterior = 2 };

// TODO: IF THE RESULT IS UPDATED WITH THE MAX POSSIBLE VALUE FOR SOME PAIR OF GEOEMTRIES
// THE VALUE ALREADY STORED MUSN'T BE CHECKED
// update() calls chould be replaced with set() in those cases
// but for safety reasons (STATIC_ASSERT) we should check if parameter D is valid and set() doesn't do that
// so some additional function could be added, e.g. set_dim()


template <typename MatrixOrMask, field F1, field F2>
using fields_in_bounds = util::bool_constant
    <
        (F1 < MatrixOrMask::static_height && F2 < MatrixOrMask::static_width)
    >;

// --------------- MATRIX ----------------

// matrix

template <std::size_t Height, std::size_t Width = Height>
class matrix
{
public:
    typedef char value_type;
    typedef std::size_t size_type;
    typedef const char * const_iterator;
    typedef const_iterator iterator;

    static const std::size_t static_width = Width;
    static const std::size_t static_height = Height;
    static const std::size_t static_size = Width * Height;

    inline matrix()
    {
        std::fill_n(m_array, static_size, 'F');
    }

    template <field F1, field F2>
        requires fields_in_bounds<matrix, F1, F2>::value
    inline char get() const
    {
        static const std::size_t index = F1 * Width + F2;
        static_assert(index < static_size, "Invalid index pair (out of bounds).");
        return m_array[index];
    }

    template <field F1, field F2, char V>
        requires fields_in_bounds<matrix, F1, F2>::value
    inline void set()
    {
        static const std::size_t index = F1 * Width + F2;
        static_assert(index < static_size, "Invalid index pair (out of bounds).");
        m_array[index] = V;
    }

    inline char operator[](std::size_t index) const
    {
        BOOST_GEOMETRY_ASSERT(index < static_size);
        return m_array[index];
    }

    inline const_iterator begin() const
    {
        return m_array;
    }

    inline const_iterator end() const
    {
        return m_array + static_size;
    }

    inline static std::size_t size()
    {
        return static_size;
    }

    inline const char * data() const
    {
        return m_array;
    }

    inline std::string str() const
    {
        return std::string(m_array, static_size);
    }

private:
    char m_array[static_size];
};

// matrix_handler

template <typename Matrix>
class matrix_handler
{
public:
    typedef Matrix result_type;

    static const bool interrupt = false;

    matrix_handler()
    {}

    result_type const& result() const
    {
        return m_matrix;
    }

    result_type const& matrix() const
    {
        return m_matrix;
    }

    result_type & matrix()
    {
        return m_matrix;
    }

    template <field F1, field F2, char D>
    inline bool may_update() const
    {
        static_assert('0' <= D && D <= '9', "D must be digit.");
        char const c = m_matrix.template get<F1, F2>();
        return D > c || c > '9';
    }

    template <field F1, field F2, char V>
    inline void update()
    {
        static_assert(('0' <= V && V <= '9') || V == 'T', "D must be digit or T");
        char const c = m_matrix.template get<F1, F2>();
        // If c == T and V == T it will be set anyway but that's fine
        if (V > c || c > '9')
        {
            m_matrix.template set<F1, F2, V>();
        }
    }

    template <field F1, field F2, char V>
    inline void set()
    {
        static_assert(('0' <= V && V <= '9') || V == 'T', "V must be digit or T");
        m_matrix.template set<F1, F2, V>();
    }

    template <field F1, field F2>
    inline char get() const
    {
        return m_matrix.template get<F1, F2>();
    }

private:
    Matrix m_matrix;
};

// --------------- RUN-TIME MASK ----------------

// run-time mask

template <std::size_t Height, std::size_t Width = Height>
class mask
{
public:
    static const std::size_t static_width = Width;
    static const std::size_t static_height = Height;
    static const std::size_t static_size = Width * Height;

    inline mask(const char * s)
    {
        char * it = m_array;
        char * const last = m_array + static_size;
        for ( ; it != last && *s != '\0' ; ++it, ++s )
        {
            char c = *s;
            check_char(c);
            *it = c;
        }
        if ( it != last )
        {
            std::fill(it, last, '*');
        }
    }

    inline mask(const char * s, std::size_t count)
    {
        if ( count > static_size )
        {
            count = static_size;
        }
        if ( count > 0 )
        {
            std::for_each(s, s + count, check_char);
            std::copy_n(s, count, m_array);
        }
        if ( count < static_size )
        {
            std::fill_n(m_array + count, static_size - count, '*');
        }
    }

    template <field F1, field F2>
        requires fields_in_bounds<mask, F1, F2>::value
    inline char get() const
    {
        static const std::size_t index = F1 * Width + F2;
        static_assert(index < static_size, "Index pair invalid (out of bounds).");
        return m_array[index];
    }

private:
    static inline void check_char(char c)
    {
        bool const is_valid = c == '*' || c == 'T' || c == 'F'
                         || ( c >= '0' && c <= '9' );
        if ( !is_valid )
        {
            BOOST_THROW_EXCEPTION(geometry::invalid_input_exception());
        }
    }

    char m_array[static_size];
};

// interrupt()

template <typename Mask, bool InterruptEnabled>
struct interrupt_dispatch
{
    template <field F1, field F2, char V>
    static inline bool apply(Mask const&)
    {
        return false;
    }
};

template <typename Mask>
struct interrupt_dispatch<Mask, true>
{
    template <field F1, field F2, char V>
    static inline bool apply(Mask const& mask)
    {
        char m = mask.template get<F1, F2>();
        return check_element<V>(m);
    }

    template <char V>
    static inline bool check_element(char m)
    {
        if constexpr (V >= '0' && V <= '9')
        {
            return m == 'F' || ( m < V && m >= '0' && m <= '9' );
        }
        else if constexpr (V == 'T')
        {
            return m == 'F';
        }
        else
        {
            return false;
        }
    }
};

template <typename Masks, int I = 0, int N = std::tuple_size<Masks>::value>
struct interrupt_dispatch_tuple
{
    template <field F1, field F2, char V>
    static inline bool apply(Masks const& masks)
    {
        typedef typename std::tuple_element<I, Masks>::type mask_type;
        mask_type const& mask = std::get<I>(masks);
        return interrupt_dispatch<mask_type, true>::template apply<F1, F2, V>(mask)
            && interrupt_dispatch_tuple<Masks, I+1>::template apply<F1, F2, V>(masks);
    }
};

template <typename Masks, int N>
struct interrupt_dispatch_tuple<Masks, N, N>
{
    template <field F1, field F2, char V>
    static inline bool apply(Masks const& )
    {
        return true;
    }
};

template <typename ...Masks>
struct interrupt_dispatch<std::tuple<Masks...>, true>
{
    typedef std::tuple<Masks...> mask_type;

    template <field F1, field F2, char V>
    static inline bool apply(mask_type const& mask)
    {
        return interrupt_dispatch_tuple<mask_type>::template apply<F1, F2, V>(mask);
    }
};

template <field F1, field F2, char V, bool InterruptEnabled, typename Mask>
inline bool interrupt(Mask const& mask)
{
    return interrupt_dispatch<Mask, InterruptEnabled>
                ::template apply<F1, F2, V>(mask);
}

// may_update()

template <typename Mask>
struct may_update_dispatch
{
    template <field F1, field F2, char D, typename Matrix>
    static inline bool apply(Mask const& mask, Matrix const& matrix)
    {
        static_assert('0' <= D && D <= '9', "D must be digit.");

        char const m = mask.template get<F1, F2>();

        if ( m == 'F' )
        {
            return true;
        }
        else if ( m == 'T' )
        {
            char const c = matrix.template get<F1, F2>();
            return c == 'F'; // if it's T or between 0 and 9, the result will be the same
        }
        else if ( m >= '0' && m <= '9' )
        {
            char const c = matrix.template get<F1, F2>();
            return D > c || c > '9';
        }

        return false;
    }
};

template <typename Masks, int I = 0, int N = std::tuple_size<Masks>::value>
struct may_update_dispatch_tuple
{
    template <field F1, field F2, char D, typename Matrix>
    static inline bool apply(Masks const& masks, Matrix const& matrix)
    {
        typedef typename std::tuple_element<I, Masks>::type mask_type;
        mask_type const& mask = std::get<I>(masks);
        return may_update_dispatch<mask_type>::template apply<F1, F2, D>(mask, matrix)
            || may_update_dispatch_tuple<Masks, I+1>::template apply<F1, F2, D>(masks, matrix);
    }
};

template <typename Masks, int N>
struct may_update_dispatch_tuple<Masks, N, N>
{
    template <field F1, field F2, char D, typename Matrix>
    static inline bool apply(Masks const& , Matrix const& )
    {
        return false;
    }
};

template <typename ...Masks>
struct may_update_dispatch<std::tuple<Masks...>>
{
    typedef std::tuple<Masks...> mask_type;

    template <field F1, field F2, char D, typename Matrix>
    static inline bool apply(mask_type const& mask, Matrix const& matrix)
    {
        return may_update_dispatch_tuple<mask_type>::template apply<F1, F2, D>(mask, matrix);
    }
};

template <field F1, field F2, char D, typename Mask, typename Matrix>
inline bool may_update(Mask const& mask, Matrix const& matrix)
{
    return may_update_dispatch<Mask>
                ::template apply<F1, F2, D>(mask, matrix);
}

// check_matrix()

template <typename Mask>
struct check_dispatch
{
    template <typename Matrix>
    static inline bool apply(Mask const& mask, Matrix const& matrix)
    {
        return per_one<interior, interior>(mask, matrix)
            && per_one<interior, boundary>(mask, matrix)
            && per_one<interior, exterior>(mask, matrix)
            && per_one<boundary, interior>(mask, matrix)
            && per_one<boundary, boundary>(mask, matrix)
            && per_one<boundary, exterior>(mask, matrix)
            && per_one<exterior, interior>(mask, matrix)
            && per_one<exterior, boundary>(mask, matrix)
            && per_one<exterior, exterior>(mask, matrix);
    }

    template <field F1, field F2, typename Matrix>
    static inline bool per_one(Mask const& mask, Matrix const& matrix)
    {
        const char mask_el = mask.template get<F1, F2>();
        const char el = matrix.template get<F1, F2>();

        if ( mask_el == 'F' )
        {
            return el == 'F';
        }
        else if ( mask_el == 'T' )
        {
            return el == 'T' || ( el >= '0' && el <= '9' );
        }
        else if ( mask_el >= '0' && mask_el <= '9' )
        {
            return el == mask_el;
        }

        return true;
    }
};

template <typename Masks, int I = 0, int N = std::tuple_size<Masks>::value>
struct check_dispatch_tuple
{
    template <typename Matrix>
    static inline bool apply(Masks const& masks, Matrix const& matrix)
    {
        typedef typename std::tuple_element<I, Masks>::type mask_type;
        mask_type const& mask = std::get<I>(masks);
        return check_dispatch<mask_type>::apply(mask, matrix)
            || check_dispatch_tuple<Masks, I+1>::apply(masks, matrix);
    }
};

template <typename Masks, int N>
struct check_dispatch_tuple<Masks, N, N>
{
    template <typename Matrix>
    static inline bool apply(Masks const&, Matrix const&)
    {
        return false;
    }
};

template <typename ...Masks>
struct check_dispatch<std::tuple<Masks...>>
{
    typedef std::tuple<Masks...> mask_type;

    template <typename Matrix>
    static inline bool apply(mask_type const& mask, Matrix const& matrix)
    {
        return check_dispatch_tuple<mask_type>::apply(mask, matrix);
    }
};

template <typename Mask, typename Matrix>
inline bool check_matrix(Mask const& mask, Matrix const& matrix)
{
    return check_dispatch<Mask>::apply(mask, matrix);
}

// matrix_width

template <typename MatrixOrMask>
struct matrix_width
{
    static const std::size_t value = MatrixOrMask::static_width;
};

template <typename Tuple,
          int I = 0,
          int N = std::tuple_size<Tuple>::value>
struct matrix_width_tuple
{
    static const std::size_t
        current = matrix_width<typename std::tuple_element<I, Tuple>::type>::value;
    static const std::size_t
        next = matrix_width_tuple<Tuple, I+1>::value;

    static const std::size_t
        value = current > next ? current : next;
};

template <typename Tuple, int N>
struct matrix_width_tuple<Tuple, N, N>
{
    static const std::size_t value = 0;
};

template <typename ...Masks>
struct matrix_width<std::tuple<Masks...>>
{
    static const std::size_t
        value = matrix_width_tuple<std::tuple<Masks...>>::value;
};

// mask_handler

template <typename Mask, bool Interrupt>
class mask_handler
    : private matrix_handler
        <
            relate::matrix<matrix_width<Mask>::value>
        >
{
    typedef matrix_handler
        <
            relate::matrix<matrix_width<Mask>::value>
        > base_t;

public:
    typedef bool result_type;

    bool interrupt;

    inline explicit mask_handler(Mask const& m)
        : interrupt(false)
        , m_mask(m)
    {}

    result_type result() const
    {
        return !interrupt
            && check_matrix(m_mask, base_t::matrix());
    }

    template <field F1, field F2, char D>
    inline bool may_update() const
    {
        return detail::relate::may_update<F1, F2, D>(m_mask, base_t::matrix());
    }

    template <field F1, field F2, char V>
    inline void update()
    {
        if (relate::interrupt<F1, F2, V, Interrupt>(m_mask))
        {
            interrupt = true;
        }
        else
        {
            base_t::template update<F1, F2, V>();
        }
    }

    template <field F1, field F2, char V>
    inline void set()
    {
        if (relate::interrupt<F1, F2, V, Interrupt>(m_mask))
        {
            interrupt = true;
        }
        else
        {
            base_t::template set<F1, F2, V>();
        }
    }

    template <field F1, field F2>
    inline char get() const
    {
        return base_t::template get<F1, F2>();
    }

private:
    Mask const& m_mask;
};

// --------------- FALSE MASK ----------------

struct false_mask {};

// --------------- COMPILE-TIME MASK ----------------

// static_check_characters
template <typename Seq>
struct static_check_characters {};

template <char C, char ...Cs>
struct static_check_characters<std::integer_sequence<char, C, Cs...>>
    : static_check_characters<std::integer_sequence<char, Cs...>>
{
    typedef std::integer_sequence<char, C, Cs...> type;
    static const bool is_valid = (C >= '0' && C <= '9')
                               || C == 'T' || C == 'F' || C == '*';
    BOOST_GEOMETRY_STATIC_ASSERT((is_valid),
                                 "Invalid static mask character",
                                 type);
};

template <char ...Cs>
struct static_check_characters<std::integral_constant<char, Cs...>>
{};

// static_mask

template <typename Seq, std::size_t Height, std::size_t Width = Height>
struct static_mask
{
    static const std::size_t static_width = Width;
    static const std::size_t static_height = Height;
    static const std::size_t static_size = Width * Height;

    static_assert(std::size_t(util::sequence_size<Seq>::value) == static_size,
        "Sequence must match matrix size.");

    template <detail::relate::field F1, detail::relate::field F2>
    struct static_get
    {
        static_assert(std::size_t(F1) < static_height, "Index F1 out of bounds.");
        static_assert(std::size_t(F2) < static_width, "Index F2 out of bounds.");

        static const char value
            = util::sequence_element<F1 * static_width + F2, Seq>::value;
    };

private:
    // check static_mask characters
    enum { mask_check = sizeof(static_check_characters<Seq>) };
};

// static_should_handle_element

template <typename StaticMask, field F1, field F2>
consteval bool static_should_handle_element_value()
{
    if constexpr (util::is_sequence<StaticMask>::value)
    {
        return []<std::size_t... I>(std::index_sequence<I...>)
        {
            return (static_should_handle_element_value
                <typename util::sequence_element<I, StaticMask>::type,
                 F1, F2>() || ...);
        }(std::make_index_sequence<util::sequence_size<StaticMask>::value>{});
    }
    else
    {
        constexpr char mask_el = StaticMask::template static_get<F1, F2>::value;
        return mask_el == 'F' || mask_el == 'T'
            || (mask_el >= '0' && mask_el <= '9');
    }
}

template <typename StaticMask, field F1, field F2>
struct static_should_handle_element
    : std::bool_constant<static_should_handle_element_value<StaticMask, F1, F2>()>
{
};

// static_interrupt

template <typename StaticMask, char V, field F1, field F2,
          bool InterruptEnabled>
consteval bool static_interrupt_value()
{
    if constexpr (! InterruptEnabled)
    {
        return false;
    }
    else if constexpr (util::is_sequence<StaticMask>::value)
    {
        return []<std::size_t... I>(std::index_sequence<I...>)
        {
            return (static_interrupt_value
                <typename util::sequence_element<I, StaticMask>::type,
                 V, F1, F2, true>() && ...);
        }(std::make_index_sequence<util::sequence_size<StaticMask>::value>{});
    }
    else
    {
        constexpr char mask_el = StaticMask::template static_get<F1, F2>::value;
        if constexpr (V >= '0' && V <= '9')
        {
            return mask_el == 'F'
                || (mask_el < V && mask_el >= '0' && mask_el <= '9');
        }
        else
        {
            return V == 'T' && mask_el == 'F';
        }
    }
}

template <typename StaticMask, char V, field F1, field F2, bool EnableInterrupt>
struct static_interrupt
    : std::bool_constant
        <static_interrupt_value<StaticMask, V, F1, F2, EnableInterrupt>()>
{
};

// static_may_update

template <typename StaticMask, char D, field F1, field F2, typename Matrix>
inline bool static_may_update_apply(Matrix const& matrix)
{
    if constexpr (util::is_sequence<StaticMask>::value)
    {
        return [&]<std::size_t... I>(std::index_sequence<I...>)
        {
            return (static_may_update_apply
                <typename util::sequence_element<I, StaticMask>::type,
                 D, F1, F2>(matrix) || ...);
        }(std::make_index_sequence<util::sequence_size<StaticMask>::value>{});
    }
    else
    {
        constexpr char mask_el = StaticMask::template static_get<F1, F2>::value;
        if constexpr (mask_el == 'F')
        {
            return true;
        }
        else if constexpr (mask_el == 'T')
        {
            return matrix.template get<F1, F2>() == 'F';
        }
        else if constexpr (mask_el >= '0' && mask_el <= '9')
        {
            char const c = matrix.template get<F1, F2>();
            return D > c || c > '9';
        }
        else
        {
            return false;
        }
    }
}

template <typename StaticMask, char D, field F1, field F2>
struct static_may_update
{
    template <typename Matrix>
    static inline bool apply(Matrix const& matrix)
    {
        return static_may_update_apply<StaticMask, D, F1, F2>(matrix);
    }
};

// static_check_matrix

template <typename StaticMask, field F1, field F2, typename Matrix>
inline bool static_check_element(Matrix const& matrix)
{
    constexpr char mask_el = StaticMask::template static_get<F1, F2>::value;
    char const element = matrix.template get<F1, F2>();
    if constexpr (mask_el == 'F')
    {
        return element == 'F';
    }
    else if constexpr (mask_el == 'T')
    {
        return element == 'T' || (element >= '0' && element <= '9');
    }
    else if constexpr (mask_el >= '0' && mask_el <= '9')
    {
        return element == mask_el;
    }
    else
    {
        return true;
    }
}

template <typename StaticMask, typename Matrix>
inline bool static_check_apply(Matrix const& matrix)
{
    if constexpr (util::is_sequence<StaticMask>::value)
    {
        return [&]<std::size_t... I>(std::index_sequence<I...>)
        {
            return (static_check_apply
                <typename util::sequence_element<I, StaticMask>::type>(matrix)
                || ...);
        }(std::make_index_sequence<util::sequence_size<StaticMask>::value>{});
    }
    else
    {
        return static_check_element<StaticMask, interior, interior>(matrix)
            && static_check_element<StaticMask, interior, boundary>(matrix)
            && static_check_element<StaticMask, interior, exterior>(matrix)
            && static_check_element<StaticMask, boundary, interior>(matrix)
            && static_check_element<StaticMask, boundary, boundary>(matrix)
            && static_check_element<StaticMask, boundary, exterior>(matrix)
            && static_check_element<StaticMask, exterior, interior>(matrix)
            && static_check_element<StaticMask, exterior, boundary>(matrix)
            && static_check_element<StaticMask, exterior, exterior>(matrix);
    }
}

template <typename StaticMask>
struct static_check_matrix
{
    template <typename Matrix>
    static inline bool apply(Matrix const& matrix)
    {
        return static_check_apply<StaticMask>(matrix);
    }
};

// static_mask_handler

template <typename StaticMask, bool Interrupt>
class static_mask_handler
    : private matrix_handler< matrix<3> >
{
    typedef matrix_handler< relate::matrix<3> > base_type;

public:
    typedef bool result_type;

    bool interrupt;

    inline static_mask_handler()
        : interrupt(false)
    {}

    inline explicit static_mask_handler(StaticMask const& /*dummy*/)
        : interrupt(false)
    {}

    result_type result() const
    {
        return (!Interrupt || !interrupt)
            && static_check_matrix<StaticMask>::apply(base_type::matrix());
    }

    template <field F1, field F2, char D>
    inline bool may_update() const
    {
        return static_may_update<StaticMask, D, F1, F2>::
                    apply(base_type::matrix());
    }

    template <field F1, field F2, char V>
    inline void update()
    {
        if constexpr (static_interrupt<StaticMask, V, F1, F2, Interrupt>::value)
        {
            interrupt = true;
        }
        else if constexpr (static_should_handle_element<StaticMask, F1, F2>::value)
        {
            base_type::template update<F1, F2, V>();
        }
    }

    template <field F1, field F2, char V>
    inline void set()
    {
        if constexpr (static_interrupt<StaticMask, V, F1, F2, Interrupt>::value)
        {
            interrupt = true;
        }
        else
        {
            base_type::template set<F1, F2, V>();
        }
    }

    template <field F1, field F2>
    inline char get() const
    {
        return base_type::template get<F1, F2>();
    }

};

// --------------- UTIL FUNCTIONS ----------------

// update

template <field F1, field F2, char D, typename Result>
inline void update(Result & res)
{
    res.template update<F1, F2, D>();
}

template <field F1, field F2, char D, bool Transpose, typename Result>
inline void update(Result & res)
{
    if constexpr (Transpose)
    {
        res.template update<F2, F1, D>();
    }
    else
    {
        res.template update<F1, F2, D>();
    }
}

// may_update

template <field F1, field F2, char D, typename Result>
inline bool may_update(Result const& res)
{
    return res.template may_update<F1, F2, D>();
}

template <field F1, field F2, char D, bool Transpose, typename Result>
inline bool may_update(Result const& res)
{
    if constexpr (Transpose)
    {
        return res.template may_update<F2, F1, D>();
    }
    else
    {
        return res.template may_update<F1, F2, D>();
    }
}

// result_dimension

template <typename Geometry>
struct result_dimension
{
    static const std::size_t dim = geometry::dimension<Geometry>::value;
    static_assert(dim >= 0, "dim must be non-negative.");
    static const char value = (dim <= 9) ? ('0' + dim) : 'T';
};

}} // namespace detail::relate
#endif // DOXYGEN_NO_DETAIL

}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_ALGORITHMS_DETAIL_RELATE_RESULT_HPP
