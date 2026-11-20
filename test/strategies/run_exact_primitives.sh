#!/bin/sh
# Copyright (c) 2026 Tinko Bartels.
# Distributed under the Boost Software License, Version 1.0.
set -u
repo=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
out=${1:-/tmp/bg-exact-primitives}
boost=${BOOST_ROOT:-/Users/tinko/dev/boost}
cxx=${CXX:-clang++}
mkdir -p "$out"
snapshot_headers()
{
    find "$repo/include" -type f -exec shasum -a 256 {} + | LC_ALL=C sort
}
snapshot_headers > "$out/start-headers.sha256"
printf 'test\tcompile\trun\n' > "$out/results.tsv"
failed=0
selected=0
for name in arithmetic round_buffer round_buffer_api convex_erosion buffer_radical_field buffer_radical_noding construction_filter coordinate default_pilot decimal_overlay decimal_issues sterbenz forward_error fixed_arithmetic fixed_measures \
    construction_capacity construction_rounding borrowed_records borrowed_predicates \
    borrowed_query borrowed_alias borrowed_bounds borrowed_root_ring borrowed_turns borrowed_output \
    borrowed_primitive borrowed_clipping borrowed_operations output_support_hooks \
    borrowed_capacity borrowed_gc borrowed_geometry borrowed_tuple borrowed_reverse_dispatch \
    borrowed_segment_box borrowed_public_predicates gc_conversion_hooks gc_boxes \
    fixed_rounding typed_cache tagged_lines tagged_dispatch borrowed_native borrowed_mixed borrowed_section borrowed_ring segment_box \
    native_predicates point segment_position \
    segment filtered_segment static_predicates static_box static_clipping \
    static_ring_topology relate_overlay shared_boundary cluster rounding_policies output_policy ulp_filter filter_routing \
    lean_bounded_ulp algorithm_optimization benchmark_check additional_issue_api additional_issues gc issues sorted_point_view
do
    source=$repo/test/strategies/exact_$name.cpp
    standard=c++17
    if [ "$name" = sorted_point_view ]; then
        source=$repo/test/util/sorted_point_view.cpp
        standard=c++14
    fi
    case "$name" in ${CASE_FILTER:-*}) ;; *) continue ;; esac
    selected=$((selected + 1))
    compile=0
    define=
    if [ "$name" = benchmark_check ]; then define=-DBG_BENCH_EXACT; fi
    if [ "$name" = tagged_dispatch ]; then
        define='-DBOOST_GEOMETRY_EXACT_TAGGED_ORIENTATION=1 -DBOOST_GEOMETRY_EXACT_FILTER_STATISTICS'
    fi
    if [ "$name" = default_pilot ]; then
        define='-DBOOST_GEOMETRY_ENABLE_CARTESIAN_EXACT_DEFAULT_PILOT -DBOOST_GEOMETRY_TEST_EXPECT_EXACT_DEFAULT_PILOT'
    fi
    "$cxx" "-std=$standard" -O1 -g0 -DNDEBUG -ffp-contract=off -pthread \
        -Wno-constant-conversion $define -I"$repo/include" -I"$boost" -I"$repo/test" \
        "$source" -o "$out/$name" > "$out/$name.compile.log" 2>&1 || compile=$?
    run=SKIP
    if [ "$compile" -eq 0 ]; then
        run=0
        perl -e 'alarm shift; exec @ARGV or die $!' "${RUN_TIMEOUT:-180}" "$out/$name" \
            > "$out/$name.run.log" 2>&1 || run=$?
    fi
    printf '%s\t%s\t%s\n' "$name" "$compile" "$run" >> "$out/results.tsv"
    printf '%s compile=%s run=%s\n' "$name" "$compile" "$run"
    if [ "$compile" -ne 0 ] || [ "$run" != 0 ]; then failed=1; fi
done
if [ "$selected" -eq 0 ]; then
    printf 'No test matches shell pattern: %s\n' "${CASE_FILTER:-*}" >&2
    exit 2
fi
snapshot_headers > "$out/end-headers.sha256"
if ! cmp -s "$out/start-headers.sha256" "$out/end-headers.sha256"; then
    printf '%s\n' 'Headers changed during qualification; rerun against a frozen tree.' >&2
    failed=1
fi
exit "$failed"
