#!/usr/bin/env bash
# Run in the experiment directory with the two frozen include trees.
set -eu -o pipefail
cc=$1; out=timing/$2; core=$3
mkdir -p timing
mkdir "$out"
src=before/test/strategies/simple_cartesian_benchmark.cpp
boost=${BOOST_ROOT:?set BOOST_ROOT to the assembled Boost headers}
flags=(-std=c++17 -O3 -DNDEBUG -ffp-contract=off -fno-fast-math -fno-omit-frame-pointer -march=x86-64-v2 -pthread
       -DBOOST_GEOMETRY_DISABLE_DEPRECATED_14_WARNING -DBG_BENCH_COMPARE_BOUNDED -DBG_BENCH_EXACT -DBG_BENCH_LEAN)
budgets=(1 2 4 8 100)
{
    "$cc" --version
    echo "flags=${flags[*]}; cpu=$core; seconds=.15; repetitions=5; baseline=811718814"
    sha256sum "$src" before/test/strategies/projected_benchmark_fixture.hpp {before,after}/include/boost/geometry/strategies/cartesian/exact/{ulp_filter,borrowed_predicates}.hpp
} >"$out/metadata.txt"
for n in "${budgets[@]}"; do
    for version in before after; do
        b="$out/$version-$n"
        "$cc" "${flags[@]}" "-DBG_BENCH_ULPS=$n" -I"$version/include" -I"$boost" "$src" -o "$b" >"$b.build" 2>&1
    done
done
for n in "${budgets[@]}"; do
    for version in before after; do
        b="$out/$version-$n"
        taskset -c "$core" "$b" --label "$version-$n" --check 1 >"$b.check"
    done
done
test ! -e "$out/samples.csv"
k=0
for rep in 0 1 2 3 4; do
    echo "repeat $rep"
    for w in degenerate general decimal projected projected-local; do
        for op in union intersection difference sym_difference relate_aa relate_la relate_ll covered_by clip; do
            for j in 0 1 2 3 4; do
                n=${budgets[(k+j+rep)%5]}
                for i in 0 1; do
                    if (((k+j+i)%2)); then version=before; else version=after; fi
                    taskset -c "$core" "$out/$version-$n" --label "$version-$n" \
                        --workload "$w" --operation "$op" --seconds .15 --repeat 1 --repeat-offset "$rep" \
                        >"$out/window.csv" 2>>"$out/runner.log"
                    if [[ ! -f "$out/samples.csv" ]]; then head -1 "$out/window.csv" >"$out/samples.csv"; fi
                    tail -n +2 "$out/window.csv" >>"$out/samples.csv"
                done
            done
            ((++k))
        done
    done
done
