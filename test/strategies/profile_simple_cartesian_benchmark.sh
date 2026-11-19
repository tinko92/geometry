#!/bin/sh
# CPU instruction profiling with an isolated Callgrind, no system changes.
set -eu
here=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
action=${1:-build}
: "${PROFILE_ROOT:?Set PROFILE_ROOT to a fresh output directory}"
: "${BOOST_ROOT:?Set BOOST_ROOT to the existing Boost umbrella}"
: "${EXACT_INCLUDE:?Set EXACT_INCLUDE to the frozen geometry include directory}"
: "${VALGRIND_ROOT:?Set VALGRIND_ROOT to the extracted package directory}"
CXX=${CXX:-g++}
cpu=${TASKSET_CPU:-2}
rounds=${PROFILE_ROUNDS:-500}
case "$rounds" in ''|*[!0-9]*|0) echo 'Invalid PROFILE_ROUNDS' >&2; exit 1;; esac
export VALGRIND_LIB=$VALGRIND_ROOT/usr/libexec/valgrind
mkdir -p "$PROFILE_ROOT"

fingerprint()
{
    sha256sum "$here/simple_cartesian_benchmark.cpp" "$here/profile_simple_cartesian_calls.cpp"
    find "$EXACT_INCLUDE" -type f -exec sha256sum {} + | LC_ALL=C sort
}

build()
{
    fingerprint > "$PROFILE_ROOT/sources.sha256"
    {
        date -u
        "$CXX" --version
        "$VALGRIND_ROOT/usr/bin/valgrind" --version
        uname -a
        printf 'cpu=%s rounds=%s\n' "$cpu" "$rounds"
        printf 'flags=%s\n' '-std=c++17 -O3 -g -DNDEBUG -ffp-contract=off -fno-omit-frame-pointer -pthread'
    } > "$PROFILE_ROOT/build.txt"
    for variant in exact simple; do
        case "$variant" in exact) define=-DBG_BENCH_EXACT;; simple) define=-DBG_BENCH_SIMPLE;; esac
        "$CXX" -std=c++17 -O3 -g -DNDEBUG -ffp-contract=off -fno-omit-frame-pointer \
            -pthread -DBOOST_GEOMETRY_DISABLE_DEPRECATED_14_WARNING "$define" \
            -I"$EXACT_INCLUDE" -I"$BOOST_ROOT" -I"$VALGRIND_ROOT/usr/include" \
            "$here/profile_simple_cartesian_calls.cpp" -o "$PROFILE_ROOT/$variant"
        sha256sum "$PROFILE_ROOT/$variant" >> "$PROFILE_ROOT/build.txt"
    done
}

verify()
{
    fingerprint | cmp - "$PROFILE_ROOT/sources.sha256"
}

collect()
{
    verify
    date -u > "$PROFILE_ROOT/collection.txt"
    n=0
    for category in degenerate general decimal; do
        for operation in union intersection difference sym_difference relate_aa relate_la relate_ll covered_by clip; do
            if [ $((n % 2)) = 0 ]; then order='exact simple'; else order='simple exact'; fi
            for variant in $order; do
                name=$variant-$category-$operation
                [ ! -e "$PROFILE_ROOT/$name.callgrind" ] || { echo "Profile already exists: $name" >&2; exit 1; }
                echo "$name" >&2
                taskset -c "$cpu" "$VALGRIND_ROOT/usr/bin/valgrind" --tool=callgrind \
                    --instr-atstart=no --collect-jumps=yes --read-inline-info=yes \
                    --callgrind-out-file="$PROFILE_ROOT/$name.callgrind" \
                    "$PROFILE_ROOT/$variant" "$category" "$operation" "$rounds" \
                    > "$PROFILE_ROOT/$name.csv" 2> "$PROFILE_ROOT/$name.log"
            done
            n=$((n + 1))
        done
    done
    verify
    date -u >> "$PROFILE_ROOT/collection.txt"
}

reports()
{
    verify
    for data in "$PROFILE_ROOT/"*.callgrind; do
        "$VALGRIND_ROOT/usr/bin/callgrind_annotate" --inclusive=no --show-percs=no \
            --threshold=100 --auto=no "$data" > "$data-functions.txt"
    done
}

case "$action" in
    build) build;;
    collect) collect;;
    reports) reports;;
    *) echo "Usage: $0 build|collect|reports" >&2; exit 1;;
esac
