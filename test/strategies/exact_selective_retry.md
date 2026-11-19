# Selective determinant retry for one-step output

Candidate based on `12bb0ff9a4d1dd71a36818d8d846806bc46aa55c`. The initial
certificate is unchanged. The new retry is used only by
`bounded_ulp_output<1>` (also when its existing `Retry` argument is true).
There is no new production selector. Larger budgets, uncertified denominators
and inputs outside the existing binary64 exponent guard keep their old
paths. Unlike the compensated-head experiment, this retry never replaces a
head, quotient, candidate coordinate or denominator lower bound.

## Bounds and admission

Use the existing notation: `n,d` are the rounded determinant heads,
`t=RN(n/d)`, `L>0` is the certified lower bound for the exact `|D|`,
`a=|t|`, and `z=fma(t,h,origin)` on an axis of rounded direction `h`.
Write `E_T,E_D` for determinant error bounds **in units of** `u=2^-53`.
The original bounds are `8*MT,8*MD`. Use the existing exact FMA residual
scaled by `1/u`, `R=|n-t*d|/u`.

- **Preserve the original determinant head.** If all four differences for
  that determinant have zero TwoDiff tails, let `p,q` be the same rounded
  products that formed `k=RN(p-q)`. With exact product tails `ep,eq` and
  exact subtraction tail `ek`, the true error is `ek+ep-eq`. Put
  `s=RN(ep-eq)` and `r=RN(ek+s)`. Then
  `|det-k| <= (|r|+u*|s|)/(1-u)`: the two posterior rounding errors sum
  to at most `u*(|r|+|s|)/(1-u)`. The implemented scaled bound
  `2*RN(|r|/u+|s|)` suffices because `2*(1-u)^2 > 1`.
  If any endpoint subtraction is inexact, retain the old bound. In the
  admitted band, differences and subtraction tails are on the `2^-116`
  grid; products and determinant residuals are on the `2^-232` grid.
  These terms and their power-of-two scalings remain finite and normal
  when nonzero.
- **Reuse the existing rounding-cell certificate.** The quotient relation
  gives `|T-t*D| <= u*(R+E_T+a*E_D)`. With `U` the exact direction,
  `|U| <= |h|/(1-u)` and `|U-h| <= u*|h|/(1-u)`. Hence the
  unrounded coordinate error is bounded by
  `u*|h|*(a*L+R+E_T+a*E_D)/(L*(1-u))`.
  Evaluating its nonnegative height in two FMAs and one addition loses at
  most `(1-u)^3`, just as in the initial certificate. The existing
  coefficient and strict rounded-product comparison therefore still apply.
- **Prove when work cannot help.** Keep `L,t,z` fixed: the rounded height
  and its product with `|h|` are monotone in all three error terms. First
  screen with direction alone, then evaluate the residual-only retry with
  both old determinant bounds. Recover tails only if setting both determinant
  errors to zero would pass on an axis still failing. Try the larger
  contribution first; after that, recover the second set only if erasing
  its remaining contribution could pass another failed axis. This is an
  optimistic admission test, not a claim that a passed gate guarantees a
  rescue. The new bounds are clamped to the old ones. Passing axes never
  enter the retry and retain their stored values.

## Build the pair

Use the same benchmark source, compiler and architecture on both builds.
Only the include tree changes. For example, from the candidate clone
(Bash or Zsh):

```sh
export BOOST_ROOT=/Users/tinko/dev/boost
CXX=${CXX:-clang++}
out=$(mktemp -d /tmp/bg-selective-retry-ab.XXXXXX)
mkdir "$out/base"
git archive 12bb0ff9a include | tar -x -C "$out/base"
build() {
  "$CXX" -std=c++17 -O3 -DNDEBUG -fno-fast-math -ffp-contract=off -pthread \
    -DBOOST_GEOMETRY_DISABLE_DEPRECATED_14_WARNING -DBG_BENCH_EXACT -DBG_BENCH_ULPS=1 \
    -I"$1" -isystem "$BOOST_ROOT" test/strategies/simple_cartesian_benchmark.cpp -o "$2"
}
build "$out/base/include" "$out/before"
build "$PWD/include" "$out/after"
for variant in before after; do
  "$out/$variant" --label "$variant" --check 1 > "$out/$variant-check.csv"
done
diff -u \
  <(awk -F, 'NR>1 {print $3,$4,$11,$12,$13,$14,$15,$16}' "$out/before-check.csv" | sort) \
  <(awk -F, 'NR>1 {print $3,$4,$11,$12,$13,$14,$15,$16}' "$out/after-check.csv" | sort)
```

For each individual operation and workload, alternate before/after and
after/before, at least five times, with `--seconds 0.045` or longer. Example:
`"$out/after" --workload projected-local --operation intersection --seconds 0.1`.
Take medians per cell, then geometric means of the matched ratios. Ratios
below one are faster. Do not compare time-dependent checksums.

The focused tests in `exact_construction_rounding.cpp` check the new budget
against rational arithmetic, including a rounded-product tie, and verify
numerator-only, denominator-only, both, residual-only and no-retry cases.
They check rejected axes, unchanged passing candidates, and the exact
one-step and predicate-enclosure guarantees against independent rational
rounding.

## Local results

Apple ARM64, Apple Clang 21 and GCC 16, using the flags above. These are
candidate / `12bb0ff9a` geometric means over the nine workload cells
(five alternating median runs per cell, 45 ms minimum per run):

| Workload | Clang | GCC |
| --- | ---: | ---: |
| degenerate | 0.9771 | 0.9880 |
| general | 0.8438 | 0.8643 |
| decimal | 0.9626 | 0.9564 |
| projected | 1.0049 | 0.9996 |
| projected-local | 0.8401 | 0.8415 |

In a seven-repeat, 140 ms recheck of slower cells, Clang projected
linear/areal relate was 0.9933x; GCC projected difference was 1.0030x and
GCC degenerate difference and symmetric difference were 1.0097x and 1.0164x.
The Clang degenerate `covered_by` control was 1.0561x despite not exercising
intersection construction. Treat small individual changes as local
measurements, not guarantees for other targets.

Both compilers matched all 45 corpus/output digests and shape counts against
the baseline. Construction-rounding tests passed under both compilers:
N=1 certified 381/1054 axes versus 372/1054 at baseline. N=4 and N=100
certificate counts did not change. The Clang public output-policy test also
passed. No server or full-suite run.
