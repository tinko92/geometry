// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2007-2012 Barend Gehrels, Amsterdam, the Netherlands.
// Copyright (c) 2008-2012 Bruno Lalande, Paris, France.
// Copyright (c) 2009-2012 Mateusz Loskot, London, UK.

// Parts of Boost.Geometry are redesigned from Geodan's Geographic Library
// (geolib/GGL), copyright (c) 1995-2010 Geodan, Amsterdam, the Netherlands.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_STRATEGIES_SIDE_INFO_HPP
#define BOOST_GEOMETRY_STRATEGIES_SIDE_INFO_HPP

#include <cmath>
#include <utility>

#include <boost/geometry/strategies/side.hpp>

#if defined(BOOST_GEOMETRY_DEBUG_INTERSECTION) || defined(BOOST_GEOMETRY_DEBUG_ROBUSTNESS)
#  include <iostream>
#endif

namespace boost { namespace geometry
{

// Silence warning C4127: conditional expression is constant
#if defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable : 4127)
#endif

/*!
\brief Class side_info: small class wrapping for sides (-1,0,1)
*/
class side_info
{
public :
    inline side_info(side_type side_a1 = side_type::collinear, side_type side_a2 = side_type::collinear,
                     side_type side_b1 = side_type::collinear, side_type side_b2 = side_type::collinear)
        : sides{{side_a1, side_a2}, {side_b1, side_b2}} {}

    template <int Which>
    inline void set(side_type first, side_type second)
    {
        sides[Which].first = first;
        sides[Which].second = second;
    }

    template <int Which, int Index>
    inline void correct_to_zero()
    {
        if (Index == 0)
        {
            sides[Which].first = side_type::collinear;
        }
        else
        {
            sides[Which].second = side_type::collinear;
        }
    }

    template <int Which, int Index>
    inline side_type get() const
    {
        return Index == 0 ? sides[Which].first : sides[Which].second;
    }


    // Returns true if both lying on the same side WRT the other
    // (so either 1,1 or -1-1)
    template <int Which>
    inline bool same() const
    {
        return static_cast<int>(sides[Which].first) * static_cast<int>(sides[Which].second) == 1;
    }

    inline bool collinear() const
    {
        return sides[0].first  == side_type::collinear
            && sides[0].second == side_type::collinear
            && sides[1].first  == side_type::collinear
            && sides[1].second == side_type::collinear;
    }

    inline bool crossing() const
    {
        return static_cast<int>(sides[0].first) * static_cast<int>(sides[0].second) == -1
            && static_cast<int>(sides[1].first) * static_cast<int>(sides[1].second) == -1;
    }

    inline bool touching() const
    {
        return   (static_cast<int>(sides[0].first) * static_cast<int>(sides[1].first) == -1
               && sides[0].second == side_type::collinear && sides[1].second == side_type::collinear)
            ||   (static_cast<int>(sides[1].first) * static_cast<int>(sides[0].first) == -1
               && sides[1].second == side_type::collinear && sides[0].second == side_type::collinear);
    }

    template <int Which>
    inline bool one_touching() const
    {
        // This is normally a situation which can't occur:
        // If one is completely left or right, the other cannot touch
        return one_zero<Which>()
            && static_cast<int>(sides[1 - Which].first) * static_cast<int>(sides[1 - Which].second) == 1;
    }

    inline bool meeting() const
    {
        // Two of them (in each segment) zero, two not
        return one_zero<0>() && one_zero<1>();
    }

    template <int Which>
    inline bool zero() const
    {
        return sides[Which].first  == side_type::collinear
            && sides[Which].second == side_type::collinear;
    }

    template <int Which>
    inline bool one_zero() const
    {
        return (sides[Which].first == side_type::collinear && sides[Which].second != side_type::collinear)
            || (sides[Which].first != side_type::collinear && sides[Which].second == side_type::collinear);
    }

    inline bool one_of_all_zero() const
    {
        int const sum = std::abs(static_cast<int>(sides[0].first))
                      + std::abs(static_cast<int>(sides[0].second))
                      + std::abs(static_cast<int>(sides[1].first))
                      + std::abs(static_cast<int>(sides[1].second));
        return sum == 3;
    }

#if defined(BOOST_GEOMETRY_DEBUG_INTERSECTION) || defined(BOOST_GEOMETRY_DEBUG_ROBUSTNESS)
    inline void debug() const
    {
        std::cout << static_cast<int>(sides[0].first) << " "
                  << static_cast<int>(sides[0].second) << " "
                  << static_cast<int>(sides[1].first) << " "
                  << static_cast<int>(sides[1].second)
                  << std::endl;
    }
#endif

    inline void reverse()
    {
        std::swap(sides[0], sides[1]);
    }

//private :
    std::pair<side_type, side_type> sides[2];

};

#if defined(_MSC_VER)
#pragma warning(pop)
#endif

}} // namespace boost::geometry


#endif // BOOST_GEOMETRY_STRATEGIES_SIDE_INFO_HPP
