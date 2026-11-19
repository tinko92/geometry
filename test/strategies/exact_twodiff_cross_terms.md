# Nonzero TwoDiff tails in the bounded-output retry

Experiment on `dfb6f0449`: `0b30b3206` retains the original heads and
tightens their budgets if all relevant subtraction tails are zero.
`3615f8d11` instead tried compensated replacement heads with a coarse
rounded-difference bound or a zero-tail gate. This experiment tightens the
**original** heads with nonzero tails. Define
`BOOST_GEOMETRY_EXACT_TWODIFF_CROSS_TERMS=1` consistently across translation
units to try it. The default is 0. Only the existing N<=4 weighted retry
changes; the first attempt, positive lower-bound gate, fixed quotient,
candidate, crossing construction, and exact fallback are unchanged. On
`cc7de2684` a failed lower-bound gate still enters denominator rescue,
independently of the TwoDiff switch. A positive lower bound allows the
same-head retry without rebuilding either determinant. The benchmark
numbers below are from the earlier `dfb6f0449` experiment, not this composition.

## Certificate

For rounded factors `x,y,vx,vy`, write the exact differences as
`x+xt, y+yt, vx+vxt, vy+vyt`. Recover the tails by TwoDiffTail. If
`k=RN(p-q)`, `p=RN(x*vy)`, and `q=RN(y*vx)`, the exact error of `k` for
the rounded factors is `ek+ep-eq`, from TwoDiffTail and the two exact FMA
product tails. To bound the original determinant, evaluate the signed
linear change `x*vyt + xt*vy - y*vxt - yt*vx` with four FMAs appended to
that residual. Bound `xt*vyt - yt*vxt` without evaluating it.

Let `u=2^-53` and `M=RN(|p|+|q|)`. All operations are within the existing
binary64 admission band: endpoint differences and tails are on the
`2^-116` grid, and the determinant terms and accumulated correction are
on the `2^-232` grid. Each nonzero value is normal. With round-to-nearest
and gradual underflow, the sums of the two exact product magnitudes and
the two rounded magnitudes are each less than `2*M` when `M>0`. Thus:

- `|ek|+|ep|+|eq| <= 4*u*M`.
- `|tail| <= u*|head|/(1-u)` gives a linear-term absolute sum `< 5*u*M`
  and a quadratic absolute sum `< 3*u*u*M`.
- For each of the two initial sums and four FMAs, the exact intermediate
  is bounded by `9*u*M` plus the previous rounding error, so induction
  gives `< 16*u*M` at every step and a total error `< 96*u*u*M`.

For `M=0` the weighted terms and residual vanish. Otherwise the final
corrected residual `r` gives the budget, in units of `u`,
`(1+4*u)*fma(128*u, M, abs(r)/u)`. The factor covers the two final rounds
since `(1+4*u)*(1-u)^2 > 1`. The existing zero-tail bound is unchanged,
and the new budget is capped by the already-certified `8*M`. Shared
`vx,vy` tails are recovered once within the retry. Numerator budget has
weight 1 and denominator budget weight `|t|`, as before; the second
determinant is skipped if it cannot recover another axis. No new positive
denominator lower bound is attempted.

## Verification and timing

`exact_construction_rounding` checks the determinant budgets against
original endpoint rationals, and each N=1,2,3,4 coordinate against the
independent exact line-equation oracle. It checks the stored enclosure,
unchanged already-certified values and untouched refused axes. Fixtures
include a proper crossing with a lost subtraction tail and a cancellation
with only a quadratic term left and no positive denominator certificate.
Across 2,561 mixed decimal/translated and exponent-stress constructions,
fallback counts at N=1,2,3,4 fell respectively from
`2387/2382/2375/2370` to `2360/2159/1775/1590`. This corpus includes line
intersections outside the segments and is not a whole-overlay workload.

Local Darwin arm64, Apple Clang 21, `-O3 -DNDEBUG -ffp-contract=off`.
Ratios are experiment / baseline wall time. Construction uses the existing
`--benchmark-bounded` mixed-decimal corpus, three interleaved runs per mode.
Polygon intersection uses six 0.20s repetitions per mode in ABBA batches:

| Workload | N=1 | N=4 |
| --- | ---: | ---: |
| Construction median | .836 | .488 |
| Decimal polygon intersection | .908 | .885 |
| Projected local polygon intersection | .985 | .981 |
| Matched translated polygon intersection | .980 | 1.031 |

Construction cache entries (512 inputs per sweep) dropped from 338 to 317
at N=1 and from 323 to 208 at N=4. Unchanged N=100 and ordinary policies
were within 1%. The projected pair both admit exact differences (common
lattice locally, Sterbenz in the translated frame) and are controls, not
evidence of useful nonzero-tail work. All 54 paired untimed operation
checks over the three polygon workloads matched in quantized geometry
fingerprints and output counts; this is not a bit-identical coordinate
promise. Keep the switch off by default pending broader measurements.

To reproduce, build `exact_construction_rounding.cpp` at switch values 0
and 1 and run normally and with `--benchmark-bounded`. Build
`simple_cartesian_benchmark.cpp` with `-DBG_BENCH_EXACT -DBG_BENCH_ULPS=1`
or `4`, each at both switch values. Run `--operation intersection
--workload decimal --seconds 0.20 --repeat 3`, then repeat for
`projected-local` and `projected`. Add `--check 1 --operation all` for
untimed fingerprint comparison.
