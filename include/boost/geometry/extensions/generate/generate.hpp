// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2019 Tinko Bartels, Berlin, Germany.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_EXTENSIONS_GENERATE_GENERATE_HPP
#define BOOST_GEOMETRY_EXTENSIONS_GENERATE_GENERATE_HPP

namespace boost { namespace geometry { namespace generate
{

    template< typename OutputIt, typename Size, typename Generator, typename Pred>
    OutputIt generate_n(OutputIt const& first, Size const& count, Generator g, Pred const& pred)
    {
        OutputIt it = first;
        for( Size i = 0; i < count; ) {
            *it = g();
            if(pred(*it)) {
                ++it;
                ++i;
            }
        }
        return it;
    }

    template< typename OutputIt, typename Size, typename Generator, typename Pred>
    OutputIt generate_n(OutputIt const& first, Size const& count, Generator g, Pred& pred)
    {
        OutputIt it = first;
        for( Size i = 0; i < count; ) {
            *it = g();
            if(pred(*it)) {
                ++it;
                ++i;
            }
        }
        return it;
    }

    template<typename ForwardIt, typename Generator, typename Pred>
    void generate(ForwardIt const& first, ForwardIt const& last, Generator g, Pred& pred)
    {
        ForwardIt it = first;
        while (it != last) {
            *it = g();
            if(pred(*it)) {
                ++it;
            }
        }
    }

    template<typename ForwardIt, typename Generator, typename Pred>
    void generate(ForwardIt const& first, ForwardIt const& last, Generator g, Pred const& pred)
    {
        ForwardIt it = first;
        while (it != last) {
            *it = g();
            if(pred(*it)) {
                ++it;
            }
        }
    }

}}} // namespace boost::geometry::generate

#endif // BOOST_GEOMETRY_EXTENSIONS_GENERATE_GENERATE_HPP
