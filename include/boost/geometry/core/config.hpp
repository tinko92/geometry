// Boost.Geometry

// Copyright (c) 2018 Barend Gehrels, Amsterdam, the Netherlands.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_CORE_CONFIG_HPP
#define BOOST_GEOMETRY_CORE_CONFIG_HPP

#define BOOST_GEOMETRY_USE_RESCALING
#define BOOST_GEOMETRY_USE_KRAMER_RULE
#define BOOST_GEOMETRY_USE_COMPLEX_SEGMENT_RATIO

// For backward compatibility, deprecated
#if defined(BOOST_GEOMETRY_USE_RESCALING)
#define BOOST_GEOMETRY_NO_ROBUSTNESS
#endif

#endif // BOOST_GEOMETRY_CORE_CONFIG_HPP
