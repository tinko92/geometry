#!/bin/sh
# Copyright (c) 2026 Tinko Bartels.
# Distributed under the Boost Software License, Version 1.0.
# Compile the default-strategy Cartesian set-operation and relate regressions.
# Override CXX, CXXSTD, BOOST_ROOT, OPT, JOBS (1 or 2), or RUN_TIMEOUT as needed.
# CASE_FILTER selects paths; EXTRA_CXXFLAGS adds flags; BASELINE pins reference headers.
set -u

repo=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
out=${1:-/tmp/bg-exact-regression-default}
boost_root=${BOOST_ROOT:-$(CDPATH= cd -- "$repo/../.." && pwd)}
cxx=${CXX:-clang++}
standard=${CXXSTD:-c++17}
opt=${OPT:--O0}
jobs=${JOBS:-2}
run_timeout=${RUN_TIMEOUT:-180}
extra_flags=${EXTRA_CXXFLAGS:-}
baseline=${BASELINE:-72b2d832}

case "$jobs" in 1|2) ;; *) printf '%s\n' 'JOBS must be 1 or 2' >&2; exit 2 ;; esac
case "$opt" in -O0|-O1) ;; *) printf '%s\n' 'OPT must be -O0 or -O1' >&2; exit 2 ;; esac
mkdir -p "$out/cases" "$out/base"
if [ ! -f "$out/base/include/boost/geometry.hpp" ]; then
    git -C "$repo" archive "$baseline" include | tar -xf - -C "$out/base"
fi
if [ ! -f "$out/base/include/boost/geometry.hpp" ]; then
    printf '%s\n' 'Could not extract isolated baseline headers' >&2
    exit 2
fi

find "$repo/test/algorithms/set_operations" "$repo/test/algorithms/relate" \
    -name '*.cpp' ! -name '*_sph.cpp' ! -name '*_geo.cpp' | sort > "$out/all-sources.txt"
{
    # Exercise output-iterator adaptations before the larger areal suites.
    for file in intersection/intersection_tupled difference/difference_tupled \
                union/union_tupled sym_difference/sym_difference_tupled \
                intersection/intersection_pl_pl intersection/intersection_pl_a \
                intersection/intersection_pl_l intersection/intersection_areal_areal_linear; do
        printf '%s/test/algorithms/set_operations/%s.cpp\n' "$repo" "$file"
    done
    cat "$out/all-sources.txt"
} | awk '!seen[$0]++' | awk -v filter="${CASE_FILTER:-.}" '$0 ~ filter' > "$out/sources.txt"

{
    date -u
    printf 'repo=%s\nboost_root=%s\nbaseline=%s\njobs=%s\n' "$repo" "$boost_root" "$baseline" "$jobs"
    printf 'flags=-std=%s %s -g0 -ffp-contract=off -pthread\n' "$standard" "$opt"
    printf 'extra_flags=%s\n' "$extra_flags"
    printf 'include_order=<selected>/include %s %s/test\n' "$boost_root" "$repo"
    "$cxx" --version
} > "$out/config.txt"
git -C "$repo" diff "$baseline" -- include > "$out/start-headers.diff"

compile()
{
    "$cxx" "-std=$standard" "$opt" -g0 -ffp-contract=off -pthread $extra_flags \
        -I"$1/include" -I"$boost_root" -I"$repo/test" \
        -I"$repo/test/algorithms/set_operations/intersection" \
        -I"$repo/test/algorithms/set_operations/union" \
        -I"$repo/test/algorithms/set_operations/difference" \
        -I"$repo/test/algorithms/set_operations/sym_difference" \
        "$source" -o "$case_dir/$2" \
        > "$case_dir/$2.compile.log" 2>&1
}

run_binary()
{
    # Bound a hung test without requiring GNU timeout on macOS.
    perl -e 'alarm shift; exec @ARGV or die $!' "$run_timeout" "$case_dir/$1" \
        > "$case_dir/$1.run.log" 2>&1
}

run_one()
{
    source=$1
    relative=${source#"$repo/test/algorithms/"}
    name=$(printf '%s' "$relative" | tr / _ | sed 's/\.cpp$//')
    case_dir=$out/cases/$name
    mkdir -p "$case_dir"
    printf 'START %s\n' "$relative"
    current_compile=0
    compile "$repo" current || current_compile=$?
    printf '%s\n' "$current_compile" > "$case_dir/current.compile.rc"
    current_run=SKIP
    if [ "$current_compile" -eq 0 ]; then
        current_run=0
        run_binary current || current_run=$?
        printf '%s\n' "$current_run" > "$case_dir/current.run.rc"
    fi
    classification=PASS
    if [ "$current_compile" -ne 0 ] || [ "$current_run" -ne 0 ]; then
        printf 'CHECK_BASE compile=%s run=%s %s\n' "$current_compile" "$current_run" "$relative"
        base_compile=0
        compile "$out/base" baseline || base_compile=$?
        printf '%s\n' "$base_compile" > "$case_dir/baseline.compile.rc"
        if [ "$base_compile" -ne 0 ]; then
            classification=BASE_COMPILE_FAILURE
            if [ "$current_compile" -eq 0 ]; then classification=BASE_COMPILE_CURRENT_RUNTIME_FAILURE; fi
        else
            base_run=0
            run_binary baseline || base_run=$?
            printf '%s\n' "$base_run" > "$case_dir/baseline.run.rc"
            if [ "$current_compile" -ne 0 ]; then
                classification=NEW_COMPILE_FAILURE
            elif [ "$base_run" -eq 0 ]; then
                classification=NEW_RUNTIME_FAILURE
            elif diff -u "$case_dir/baseline.run.log" "$case_dir/current.run.log" > "$case_dir/runtime.diff"; then
                classification=KNOWN_RUNTIME_FAILURE
            else
                classification=RUNTIME_DIFFERENCE
            fi
        fi
    fi
    printf '%s\t%s\t%s\t%s\n' "$relative" "$classification" "$current_compile" "$current_run" > "$case_dir/result.tsv"
    printf 'RESULT %s %s\n' "$classification" "$relative"
}

pending=
while IFS= read -r source; do
    run_one "$source" &
    current_pid=$!
    if [ "$jobs" -eq 1 ]; then
        wait "$current_pid"
    elif [ -n "$pending" ]; then
        wait "$pending"
        wait "$current_pid"
        pending=
    else
        pending=$current_pid
    fi
done < "$out/sources.txt"
if [ -n "$pending" ]; then wait "$pending"; fi

printf 'source\tclassification\tcompile_rc\trun_rc\n' > "$out/results.tsv"
while IFS= read -r source; do
    relative=${source#"$repo/test/algorithms/"}
    name=$(printf '%s' "$relative" | tr / _ | sed 's/\.cpp$//')
    cat "$out/cases/$name/result.tsv" >> "$out/results.tsv"
done < "$out/sources.txt"
awk -F '\t' 'NR > 1 { counts[$2]++; total++ } END { for (c in counts) print c, counts[c]; print "TOTAL", total }' \
    "$out/results.tsv" | sort
awk -F '\t' 'NR > 1 && ($2 ~ /^NEW_/ || $2 == "RUNTIME_DIFFERENCE" || $2 == "BASE_COMPILE_CURRENT_RUNTIME_FAILURE") { bad = 1 } END { exit bad }' \
    "$out/results.tsv"
