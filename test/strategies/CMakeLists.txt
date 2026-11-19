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
    exact_sterbenz
    exact_algorithm_optimization
    exact_benchmark_check
    exact_additional_issue_api
    exact_additional_issues
    exact_gc
    exact_construction_filter
    exact_coordinate
    exact_default_pilot
    exact_filter_routing
    exact_static_box
    exact_static_clipping
    exact_filtered_segment
    exact_forward_error
    exact_fixed_arithmetic
    exact_fixed_measures
    exact_fixed_rounding
    exact_issues
    exact_native_predicates
    exact_point
    exact_segment_position
    exact_relate_overlay
    exact_segment
    exact_segment_box
    exact_shared_boundary
    exact_static_predicates
    exact_static_ring_topology
    exact_borrowed_alias
    exact_borrowed_bounds
    exact_borrowed_capacity
    exact_borrowed_clipping
    exact_borrowed_gc
    exact_borrowed_geometry
    exact_borrowed_mixed
    exact_borrowed_native
    exact_borrowed_operations
    exact_borrowed_output
    exact_borrowed_predicates
    exact_borrowed_primitive
    exact_borrowed_public_predicates
    exact_borrowed_query
    exact_borrowed_records
    exact_borrowed_reverse_dispatch
    exact_borrowed_ring
    exact_borrowed_root_ring
    exact_borrowed_section
    exact_borrowed_segment_box
    exact_borrowed_tuple
    exact_borrowed_turns
    exact_gc_boxes
    exact_gc_conversion_hooks
    exact_construction_capacity
    exact_construction_rounding
    exact_double_qualification
    exact_cluster
    exact_rounding_policies
    exact_output_policy
    exact_output_support_hooks
    exact_ulp_filter
    exact_typed_cache
    exact_tagged_lines
    exact_tagged_dispatch
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
    simple_cartesian
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
  if (item MATCHES "^exact_" OR item STREQUAL "simple_cartesian")
    if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
      target_compile_options(${BOOST_GEOMETRY_UNIT_TEST_NAME} PRIVATE -ffp-contract=off)
    elseif (MSVC)
      target_compile_options(${BOOST_GEOMETRY_UNIT_TEST_NAME} PRIVATE /fp:strict)
    endif()
  endif()
endforeach()

target_link_libraries(boost_geometry_strategies_exact_additional_issues PRIVATE Boost::property_tree)
target_compile_definitions(boost_geometry_strategies_exact_benchmark_check PRIVATE BG_BENCH_EXACT)
target_compile_definitions(boost_geometry_strategies_exact_tagged_dispatch PRIVATE
    BOOST_GEOMETRY_EXACT_TAGGED_ORIENTATION=1 BOOST_GEOMETRY_EXACT_FILTER_STATISTICS)
target_compile_definitions(boost_geometry_strategies_exact_default_pilot PRIVATE
    BOOST_GEOMETRY_ENABLE_CARTESIAN_EXACT_DEFAULT_PILOT)

find_package(Threads REQUIRED)
# Exact point-cache regressions exercise threads. Older Unix toolchains need
# the native thread link flags even when a particular test is single-threaded.
# Keep this dependency on these tests, not on the production geometry target.
foreach(item IN ITEMS exact_additional_issue_api exact_additional_issues exact_algorithm_optimization exact_benchmark_check
    exact_gc exact_issues exact_native_predicates exact_point exact_segment_position exact_relate_overlay exact_segment exact_segment_box exact_shared_boundary
    exact_construction_filter exact_coordinate exact_default_pilot exact_filter_routing exact_static_box exact_static_clipping
    exact_filtered_segment exact_forward_error exact_fixed_arithmetic exact_fixed_measures exact_fixed_rounding
    exact_static_predicates exact_static_ring_topology
    exact_borrowed_alias exact_borrowed_bounds exact_borrowed_capacity exact_borrowed_clipping exact_borrowed_gc
    exact_borrowed_geometry exact_borrowed_mixed exact_borrowed_native exact_borrowed_operations
    exact_borrowed_output exact_borrowed_predicates exact_borrowed_primitive exact_borrowed_public_predicates exact_borrowed_query
    exact_borrowed_records exact_borrowed_reverse_dispatch exact_borrowed_ring exact_borrowed_root_ring exact_borrowed_section
    exact_borrowed_segment_box exact_gc_boxes exact_gc_conversion_hooks
    exact_borrowed_tuple exact_borrowed_turns exact_construction_capacity exact_construction_rounding exact_double_qualification
    exact_cluster exact_rounding_policies exact_output_policy exact_output_support_hooks exact_ulp_filter exact_typed_cache
    exact_tagged_dispatch)
  target_link_libraries(boost_geometry_strategies_${item} PRIVATE Threads::Threads)
endforeach()
