# Operation-scoped Sterbenz orientation

The exact Cartesian strategy scans the original pair once for a combined
binary64 box. Each axis must lie in `[m, 2m]` or its negative reflection,
with all magnitudes in `[2^-150, 2^150]`. The two axes may have different
signs and scales. The scan includes all components and interior rings.

The certificate applies to all-native orientation triples in borrowed segment
intersection stages sourced from those same immutable input objects (including
boundary views and reverse dispatch). Distinct rounded products decide the
sign; tied nonzero products use exact FMA tails. Empty combined boxes,
zeros, mixed signs, nonfinite coordinates, failed exponent tests and boxes
wider than a factor of two do not qualify. Constructed/query points and stages
with different source objects keep the existing exact predicates. This does not
change coordinate order, construction, winding, or traversal filters.

## Reproduce the A/B run

Use Bash (or Zsh) from the Sterbenz branch root. The benchmark driver and
projected fixture are the same for both builds. Only the headers differ:
`off` uses the committed series HEAD, not any working tree changes.
Set `BOOST_ROOT` to a Boost umbrella include directory and use the same
compiler and target flags for both; GCC 13 is preferred for comparisons
with the existing performance notes.

```sh
export BOOST_ROOT=/Users/tinko/dev/boost
CXX=${CXX:-clang++}
base=423afeb62536192e7d05388b95b16017126dd717
scratch=$(mktemp -d /tmp/bg-sterbenz-ab.XXXXXX)
mkdir "$scratch/base"
git archive "$base" include | tar -x -C "$scratch/base"

"$CXX" -std=c++20 -O3 -DNDEBUG -ffp-contract=off -DBG_BENCH_EXACT \
  -I"$scratch/base/include" -isystem "$BOOST_ROOT" \
  test/strategies/exact_cartesian_benchmark.cpp -o "$scratch/off"
"$CXX" -std=c++20 -O3 -DNDEBUG -ffp-contract=off -DBG_BENCH_EXACT \
  -I"$PWD/include" -isystem "$BOOST_ROOT" \
  test/strategies/exact_cartesian_benchmark.cpp -o "$scratch/on"

for variant in off on; do
  "$scratch/$variant" --projected-self-test 1 --vertices 96 --corpus 8
  "$scratch/$variant" --label "$variant" --workload projected-pair \
    --vertices 96 --corpus 8 --seconds 0.2 --repeat 5 > "$scratch/$variant.csv"
done
diff -u \
  <(awk -F, 'NR > 1 { print $3, $5, $14, $15 }' "$scratch/off.csv" | sort -u) \
  <(awk -F, 'NR > 1 { print $3, $5, $14, $15 }' "$scratch/on.csv" | sort -u)
```

Compare median `ns_per_call` for each operation and frame. `projected`
qualifies; `projected-local` straddles zero and measures the admission
overhead without the fast path. Keep other jobs stopped, repeat in reverse
binary order, and vary vertices/corpus before drawing conclusions. The CSV
timing checksum depends on the number of calls; compare the corpus and
geometry digests instead. Leave filter counters disabled in timed builds.
The cross-frame self-test can fail with other fixture sizes even on the
baseline (for example, 128 vertices and eight cases with stored approximation
output). Check for agreement within each frame before attributing a mismatch
to the specialization.
The `strategies_exact_sterbenz` test separately verifies that public
intersection hits the filter and that constructed operands do not.
