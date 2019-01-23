// Boost.Geometry

// Copyright (c) 2018 Barend Gehrels, Amsterdam, the Netherlands.
// Copyright (c) 2018 Oracle and/or its affiliates.
// Contributed and/or modified by Adam Wulkiewicz, on behalf of Oracle

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_CORE_CONFIG_HPP
#define BOOST_GEOMETRY_CORE_CONFIG_HPP

#include <boost/config.hpp>

// NOTE: workaround for VC++ 12 (aka 2013): cannot specify explicit initializer for arrays
#if !defined(BOOST_NO_CXX11_UNIFIED_INITIALIZATION_SYNTAX) && (!defined(_MSC_VER) || (_MSC_VER >= 1900))
#define BOOST_GEOMETRY_CXX11_ARRAY_UNIFIED_INITIALIZATION
#endif

//#define BOOST_GEOMETRY_USE_RESCALING
//#define BOOST_GEOMETRY_USE_KRAMER_RULE
//#define BOOST_GEOMETRY_USE_COMPLEX_SEGMENT_RATIO

#define BOOST_GEOMETRY_GENERAL_INT_SUPPORT_DOUBT1

// For backward compatibility, deprecated
#if defined(BOOST_GEOMETRY_USE_RESCALING)
#define BOOST_GEOMETRY_NO_ROBUSTNESS
#endif

#endif // BOOST_GEOMETRY_CORE_CONFIG_HPP
