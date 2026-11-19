# Boost.Geometry
# Copyright (c) 2024 Barend Gehrels, Amsterdam, the Netherlands.
# Use, modification and distribution is subject to the Boost Software License,
# Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
# http://www.boost.org/LICENSE_1_0.txt)

foreach(item IN ITEMS
    andoyer
    buffer_join
    buffer_join_geo
    buffer_side_straight
    cross_track
    crossings_multiply
    distance_default_result
    distance_cross_track
    douglas_peucker
    envelope_segment
    exact_arithmetic
    exact_construction_filter
    exact_forward_error
    exact_fixed_arithmetic
    exact_fixed_measures
    exact_fixed_rounding
    exact_point
    exact_borrowed_native
    exact_construction_capacity
    exact_construction_rounding
    exact_rounding_policies
    exact_ulp_filter
    exact_typed_cache
    franklin
    haversine
    point_in_box
    projected_point
    pythagoras
    pythagoras_point_box
    segment_intersection
    segment_intersection_collinear
    segment_intersection_geo
    segment_intersection_sph
    spherical_side
    side_rounded_input
    thomas
    transform_cs
    transformer
    matrix_transformer
    vincenty
    winding
  )
  boost_geometry_add_unit_test("strategies" ${item})
  if (item MATCHES "^exact_")
    target_compile_features(${BOOST_GEOMETRY_UNIT_TEST_NAME} PRIVATE cxx_std_17)
  endif()
  if (item MATCHES "^exact_")
    if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
      target_compile_options(${BOOST_GEOMETRY_UNIT_TEST_NAME} PRIVATE -ffp-contract=off)
    elseif (MSVC)
      target_compile_options(${BOOST_GEOMETRY_UNIT_TEST_NAME} PRIVATE /fp:strict)
    endif()
  endif()
endforeach()


find_package(Threads REQUIRED)
# Exact point-cache regressions exercise threads. Older Unix toolchains need
# the native thread link flags even when a particular test is single-threaded.
# Keep this dependency on these tests, not on the production geometry target.
foreach(item IN ITEMS exact_point exact_construction_filter exact_forward_error exact_fixed_arithmetic exact_fixed_measures exact_fixed_rounding exact_borrowed_native exact_construction_capacity exact_construction_rounding exact_rounding_policies exact_ulp_filter exact_typed_cache)
  target_link_libraries(boost_geometry_strategies_${item} PRIVATE Threads::Threads)
endforeach()
