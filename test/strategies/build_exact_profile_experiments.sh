#!/bin/sh
# Copyright (c) 2026 Tinko Bartels.
# Distributed under the Boost Software License, Version 1.0.
set -eu
repo=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
out=${1:-/tmp/bg-exact-profile-20261113}
boost=${BOOST_ROOT:-/Users/tinko/dev/boost}
for coordinate in double int64; do
    flag=
    if [ "$coordinate" = int64 ]; then flag=-DBG_BENCH_INT64; fi
    for variant in ${VARIANTS:-compare interval combined validation provenance endpoints rounding all}; do
        clang++ -std=c++14 -O3 -DNDEBUG -gline-tables-only -fno-omit-frame-pointer \
            -ffp-contract=fast -pthread -DBOOST_GEOMETRY_DISABLE_DEPRECATED_14_WARNING \
            -DBG_BENCH_EXACT $flag -I"$out/experiments/$variant" -I"$repo/include" -I"$boost" \
            "$repo/test/strategies/exact_cartesian_benchmark.cpp" -o "$out/$variant-$coordinate" \
            > "$out/$variant-$coordinate.build.log" 2>&1
    done
done
