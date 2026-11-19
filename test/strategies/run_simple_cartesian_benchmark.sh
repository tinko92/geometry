#!/bin/sh
# Serial, case-paired benchmark. Only BENCH_COMPARE_BOUNDED=1 needs Python 3.
# EXACT_INCLUDE=... DEVELOP_INCLUDE=... BOOST_ROOT=... CXX=clang++ TASKSET_CPU=2 \
#   sh run_simple_cartesian_benchmark.sh build|check|run|all [output-directory]
# BEFORE_INCLUDE optionally adds three old-header policies; BENCH_ARCH=x86-64-v2
# applies the same -march target to every variant, including develop.
# For exported headers, set DEVELOP_REVISION to the exact archived commit.
# BENCH_COMPARE_BOUNDED=1 selects only develop and five stored/lean pairs.
set -eu

here=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
action=${1:-build}
out=${2:-"$here/simple-cartesian-results"}
source=${BENCH_SOURCE:-"$here/simple_cartesian_benchmark.cpp"}
CXX=${CXX:-clang++}
seconds=${BENCH_SECONDS:-0.1}
repetitions=${BENCH_REPETITIONS:-3}
arch=${BENCH_ARCH:-}
variants="develop simple exact eager deferred"
if [ "${BENCH_BOUNDED_ULP:-0}" = 1 ]; then variants="$variants ulp1 ulp4 ulp100"; fi
categories=${BENCH_CATEGORIES:-"degenerate general decimal"}
if [ -n "${BEFORE_INCLUDE:-}" ]; then variants="$variants before before-eager before-deferred"; fi
compare=${BENCH_COMPARE_BOUNDED:-0}
if [ "$compare" = 1 ]; then
    variants=develop
    for n in 1 2 4 8 100; do variants="$variants ulp$n lean$n"; done
    categories=${BENCH_CATEGORIES:-"degenerate general decimal projected projected-local"}
fi
: "${EXACT_INCLUDE:?Set EXACT_INCLUDE to the current geometry include directory}"
: "${DEVELOP_INCLUDE:?Set DEVELOP_INCLUDE to the pinned develop include directory}"
: "${BOOST_ROOT:?Set BOOST_ROOT to the umbrella Boost include directory}"
case "$repetitions" in ''|*[!0-9]*|0) echo "Invalid BENCH_REPETITIONS" >&2; exit 1;; esac
if [ -n "${TASKSET_CPU:-}" ]; then
    command -v taskset >/dev/null 2>&1 || { echo "taskset unavailable" >&2; exit 1; }
fi
mkdir -p "$out"
out=$(CDPATH= cd -- "$out" && pwd)

variant_config()
{
    policy=
    case "$1" in
        develop) include=$DEVELOP_INCLUDE; define=;;
        simple) include=${SIMPLE_INCLUDE:-$EXACT_INCLUDE}; define=-DBG_BENCH_SIMPLE;;
        exact|eager|deferred|ulp*|lean*) include=$EXACT_INCLUDE; define=-DBG_BENCH_EXACT;;
        before|before-eager|before-deferred) include=$BEFORE_INCLUDE; define=-DBG_BENCH_EXACT;;
        *) echo "Unknown variant: $1" >&2; exit 1;;
    esac
    case "$1" in
        eager|before-eager) policy=-DBG_BENCH_EAGER;;
        deferred|before-deferred) policy=-DBG_BENCH_DEFERRED;;
        ulp*) policy="-DBG_BENCH_ULPS=${1#ulp}";;
        lean*) policy="-DBG_BENCH_ULPS=${1#lean}";;
    esac
}

build_config()
{
    printf 'CXX=%s\nBENCH_ARCH=%s\nvariants=%s\n' "$CXX" "$arch" "$variants"
    printf 'develop_revision=%s\n' "${DEVELOP_REVISION:-unspecified}"
    printf 'EXACT_INCLUDE=%s\nDEVELOP_INCLUDE=%s\nSIMPLE_INCLUDE=%s\nBEFORE_INCLUDE=%s\nBOOST_ROOT=%s\n' \
        "$EXACT_INCLUDE" "$DEVELOP_INCLUDE" "${SIMPLE_INCLUDE:-$EXACT_INCLUDE}" "${BEFORE_INCLUDE:-}" "$BOOST_ROOT"
}

invoke()
{
    if [ -n "${TASKSET_CPU:-}" ]; then taskset -c "$TASKSET_CPU" "$@";
    else "$@"; fi
}

build()
{
    for old in "$out/results.csv" "$out/run-metadata.txt" "$out/raw/"*.csv; do
        if [ -f "$old" ]; then echo "Results already exist; build in a fresh output directory" >&2; exit 1; fi
    done
    build_config > "$out/build-config.txt"
    {
        date -u '+%Y-%m-%dT%H:%M:%SZ'
        uname -a
        "$CXX" --version
        build_config
        printf 'flags=%s' '-std=c++17 -O3 -DNDEBUG -fno-fast-math -ffp-contract=off -fno-omit-frame-pointer -pthread'
        if [ -n "$arch" ]; then printf ' -march=%s' "$arch"; fi
        printf '\n'
        cksum "$source"
    } > "$out/build-metadata.txt"
    for variant in $variants; do
        variant_config "$variant"
        echo "Building $variant" >&2
        set -- -std=c++17 -O3 -DNDEBUG -fno-fast-math -ffp-contract=off \
            -fno-omit-frame-pointer -pthread -DBOOST_GEOMETRY_DISABLE_DEPRECATED_14_WARNING
        if [ -n "$arch" ]; then set -- "$@" "-march=$arch"; fi
        if [ -n "$define" ]; then set -- "$@" "$define"; fi
        if [ -n "$policy" ]; then set -- "$@" "$policy"; fi
        if [ "$compare" = 1 ]; then set -- "$@" -DBG_BENCH_COMPARE_BOUNDED; fi
        case "$variant" in lean*) set -- "$@" -DBG_BENCH_LEAN;; esac
        printf 'variant=%s define=%s policy=%s include=%s\n' "$variant" "$define" "$policy" "$include" \
            >> "$out/build-metadata.txt"
        if [ "$compare" = 1 ]; then printf 'compile_args=%s\n' "$*" >> "$out/build-metadata.txt"; fi
        "$CXX" "$@" -I"$include" -I"$BOOST_ROOT" "$source" -o "$out/$variant"
        cksum "$out/$variant" >> "$out/build-metadata.txt"
        find "$include" -type f -exec cksum {} + | LC_ALL=C sort > "$out/$variant-headers.cksum"
    done
    source_checksums > "$out/source.cksum"
}

source_checksums()
{
    cksum "$source"
    if [ "$compare" = 1 ]; then cksum "$(dirname -- "$source")/projected_benchmark_fixture.hpp"; fi
}

verify_source()
{
    build_config | cmp - "$out/build-config.txt" || {
        echo "Compiler, target, variants, or include configuration changed since build" >&2; exit 1;
    }
    source_checksums | cmp - "$out/source.cksum" || {
        echo "Source changed since build" >&2; exit 1;
    }
    for variant in $variants; do
        variant_config "$variant"
        find "$include" -type f -exec cksum {} + | LC_ALL=C sort | \
            cmp - "$out/$variant-headers.cksum" || {
                echo "$variant headers changed since build" >&2; exit 1;
            }
    done
}

check()
{
    verify_source
    for variant in $variants; do
        invoke "$out/$variant" --label "$variant" --digest-self-test 1
        invoke "$out/$variant" --label "$variant" --check 1 > "$out/$variant-check.csv"
    done
    if [ "$compare" = 1 ]; then python3 "$here/analyze_bounded_comparison.py" "$out" check; fi
}

run()
{
    if [ "$compare" != 1 ]; then verify_source; fi
    mkdir -p "$out/raw"
    for old in "$out/raw/"*.csv; do
        if [ -f "$old" ]; then echo "Raw CSVs already exist; use a fresh output directory" >&2; exit 1; fi
    done
    if [ "$compare" = 1 ]; then check; fi
    {
        date -u '+started=%Y-%m-%dT%H:%M:%SZ'
        printf 'seconds=%s\nrepetitions=%s\nTASKSET_CPU=%s\nvariants=%s\ncorpus_per_group=3\n' \
            "$seconds" "$repetitions" "${TASKSET_CPU:-unbound}" "$variants"
        printf 'BENCH_ARCH=%s\n' "$arch"
        if [ "$compare" = 1 ]; then printf 'categories=%s\n' "$categories"; fi
    } > "$out/run-metadata.txt"
    order=$variants
    repetition=0
    while [ "$repetition" -lt "$repetitions" ]; do
        for category in $categories; do
            for operation in union intersection difference sym_difference relate_aa relate_la relate_ll covered_by clip; do
                for variant in $order; do
                    echo "$repetition $category $operation $variant" >&2
                    set --
                    if [ "$compare" = 1 ]; then set -- --repeat-offset "$repetition"; fi
                    invoke "$out/$variant" --label "$variant" --workload "$category" \
                        --operation "$operation" --seconds "$seconds" --repeat 1 "$@" \
                        > "$out/raw/$variant-$category-$operation-$repetition.csv"
                done
                if [ "$compare" = 1 ]; then
                    set -- $order; first=$1; shift; order="$* $first"
                fi
            done
        done
        # Rotate starting variant each repetition, without concurrent timings.
        if [ "$compare" != 1 ]; then
            set -- $order; first=$1; shift; order="$* $first"
        fi
        repetition=$((repetition + 1))
    done
    verify_source
    first=1
    for file in "$out/raw/"*.csv; do
        if [ "$first" = 1 ]; then cat "$file"; first=0;
        else sed '1d' "$file"; fi
    done > "$out/results.csv"
    if [ "$compare" = 1 ]; then
        python3 "$here/analyze_bounded_comparison.py" "$out" timed "$repetitions" $categories
    fi
    date -u '+finished=%Y-%m-%dT%H:%M:%SZ' >> "$out/run-metadata.txt"
}

case "$action" in
    build) build;;
    check) check;;
    run) run;;
    all) build; if [ "$compare" != 1 ]; then check; fi; run;;
    *) echo "Usage: $0 build|check|run|all [output-directory]" >&2; exit 1;;
esac
