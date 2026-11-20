# Producer certificates for first-stage binary64 predicates

The construction policy publishes two independent facts: an arithmetic
certificate for the error relative to the true coordinate, and an ordering
radius measured in representable steps from its correctly rounded value.
The determinant, orientation, dot product and squared-distance kernels take
the arithmetic certificate as a type. They do not test the policy at runtime.

| Source | Arithmetic factor rho | Bilinear coefficient (units of u) | Order radius |
| --- | ---: | ---: | ---: |
| Lossless native | 0 | 12 | 0, exact |
| Eager `nearest_output` | 1 | 12 | 0, nearest |
| Lean bounded producer, K > 0 | K+2 | 2K+16 | K |
| Only a K-step promise, K > 0 | 2K+2 | 4K+16 | K |

Here is the ordering contract, independently on each axis. A native input is the exact real number
represented by its double. For a constructed input with head `z`, let `r`
be the finite binary64 nearest-even rounding of the true coordinate `Z`.
The number of adjacent representable steps from `z` to `r` is at most
`K = min(N, 2^24)`, with the two signed zeros counted once. For policies
larger than the cap, the **producer** must certify the capped distance.
Silently capping a merely N-step external guarantee would be unsound.

Kernels instantiated with `step_certificate<N>` do not require the tighter absolute error certificate available
on successful borrowed-record constructions. In particular, an arbitrary
K-step input need not satisfy a bound of `(K+2)u |z|`: near the bottom of a
normal binade an outward step can be almost twice this unit of relative error.
Eager correctly rounded means exactly **K=0 steps from RN(Z)**. It does not
mean half a step from RN(Z); the half-cell between RN(Z) and Z is already
accounted for below. A rounded zero need not mean Z=0.

The typed visitor uses this path only for eager-nearest or lean-bounded
Root policies. Each operand must be a losslessly binary64-representable
native type, the selected crossing record of that Root, or its fixed query
type (exact-native shortcut or direct nearest rounding of the fixed quotient).
The latter two are treated as uncertain even if a particular head is exact.
An ineligible combination uses the ordinary forward-error predicate with
the operands' own enclosures and the usual exact fallback. In particular
legacy points and wide integer conversions do not inherit eager nearest.
The stored-error policies keep their own enclosure path. Native-only
orientation and scalar comparison retain typed dispatch.

All operations below use binary64 round-to-nearest/ties-even, gradual
underflow, no excess precision, reassociation, contraction or fast-math.
Use `-ffp-contract=off`. These are environment preconditions, not runtime
tests. Write `u=2^-53`, `m=2^-1022`, `eta=2^-1074=2u m`.

## Scalar order

For finite IEEE binary64, the signless encoding is the ascending step index
of the magnitude, starting at zero. Adding that index to unsigned `2^63`
for positive numbers and subtracting it for negative numbers gives one
monotone lattice through both zeros. The difference of the larger and smaller
unsigned indices fits in 64 bits even for opposite largest finite inputs.
Reject nonfinite heads before decoding. The usual matching IEEE binary64
and uint64 byte order is required.

Give each constructed operand radius K in index units and each exact native
operand radius zero. If the index separation strictly exceeds the sum of
radii, the two RN values are strictly ordered in the direction of the heads.
Monotonicity of RN then strictly orders the true values. This proof needs no
extra half-step: disjoint RN cells cannot contain a common true value, even at
a midpoint tie. For overlapping index ranges the filter declines; only two
exact natives can certify equality. No neighbour beyond max-double is needed.
For K=0 the implementation just compares finite doubles: any two strictly
unequal heads suffice, with no bit decoding and no subtraction.

## Coordinate and difference bounds

Let `H(t)=max(|t|,m)`. Either adjacent spacing at a finite interior head t
is at most `2u H(t)` (including a binade boundary, the subnormal range and
zero). After i steps, H is therefore at most `(1+2u)^i H(z)`.
Include the last half-spacing from r to Z. Summing the bounds for K steps
and this half-cell gives

```text
|Z-z| <= (2K+1)u (1+2u)^K H(z) <= C H(z),   C=(2K+2)u.
```

Indeed, `(1+2u)^K <= (1-2u)^(-K) <= 1/(1-2Ku)` by Bernoulli,
and `(2K+2)2Ku <= 1` at the cap. The measure range guard
`|z| <= 2^500` keeps the entire K-step path away from overflow by this
same growth bound (less than a factor of two). Tests of this form reject
NaN and infinity. Neither zero nor any subnormal head is assumed exact.

For K=0 use the sharper coordinate bound `u H(z)` directly from the last
half-cell. For a factor `d=RN(b-a)` define a real scale S:

* If either endpoint is uncertain, `S=H(a)+H(b)`. Then `|b-a| <= S`.
  Subtraction rounding costs at most `u S`: a tiny difference is exact
  because binary64 operands are integer multiples of eta. The combined
  displacement of both endpoints costs at most `C S` (only `u S` for K=0).
  No extra `|d|` or individual coordinate error is added to the computed scale.
* For two exact endpoints use `S=|d|`, including zero. Subtraction need not
  be exact; its posterior rounding error is at most `u/(1-u) |d| <= 2u S`.

In both cases the true difference D and its head satisfy

```text
|d| <= (1+u) S,    |D-d| <= c u S,
c = 2 (K=0),   c = 2K+3 (K>0).
```

## Bilinear measures and signs

Use two products with the following pairing and final operation:

| Kernel | Real expression |
| --- | --- |
| Signed determinant / orientation | `(b-a)x(p-a)` (2D cross product) |
| Dot sign | `(b-a).(p-b)` |
| Squared distance | `(b-a).(b-a)` |

Let M be the sum of the two products of the corresponding real S factors.
The two linear perturbations and quadratic term cost at most
`[2(1+u)c u + (c u)^2] M`, including for repeated operands. The products
of the rounded differences have absolute sum at most `P=(1+u)^2 M` before
product rounding. The two multiplications cost at most `u P + eta`, and
the final addition/subtraction costs at most `u [(1+u)P + eta] + eta/2`.
Thus the total error is at most `R u M + 2 eta`, where

```text
R = 2(1+u)c + c^2 u + (2+u)(1+u)^2
  < 2c + 2 + 1/2.
```

For the strict inequality the excess above `2c+2` is increasing in c and
is less than 1/2 already at the worst case `c=2^25+3`. Hence R is below
6.5 for nearest, and below `4K+9` otherwise.

Each computed S is at least `(1-u)S` (zero or one addition). Each magnitude
product can lose a factor `1-u` and an absolute `eta/2`; their nonnegative
sum, if subnormal, is exact. The computed sum floored at m is T:

```text
T >= (1-u)^4 M - eta,    eta <= 2u T,
M <= ((1+2u)/(1-u)^4) T <= (1+8u) T.
```

Use `F=L u`, with L=12 for nearest (also for all-exact operands), and
`L=4K+16` otherwise. These static dyadic coefficients are exactly
representable. Even when subnormal, the returned error `RN(F T)` is at
least `(L(1-u)-1)u T`. This dominates the true error: for K=0,
`6.5(1+8u)+4 < 12(1-u)-1`; for positive K up to the cap,
`(4K+9)(1+8u)+4 < (4K+16)(1-u)-1` (the slack decreases with K).

The coordinate guard keeps every difference and computed S at most
`2^501`, and the sum of magnitude products at most `2^1003`.
All heads and errors are finite, including at the admitted
endpoints. Underflow of a value product, magnitude product, or final error
is covered above, not excluded by a relative-error assumption.

A measure success certifies `|true - value| <= error`, not necessarily a
sign. Orientation and dot certify only for `value > error` or
`value < -error`; a zero head, touching bound or other ambiguous case
requests fallback. On every false return the output is unchanged.

## Stronger certificate from the bounded producer

The generic `step_certificate<N>` accepts *any* K-step certificate. The
same kernels instantiated with `producer_certificate<N>` may use the
following additional fact for **every** uncertain ordinate:

```text
|Z-z| <= (K+2)u H(z),     1 <= K = min(N, 2^24).
```

This is a bound in real arithmetic; it does not require reading or rounding
a per-coordinate error field. Here is why the producer supplies it:

* For a successful `bounded_intersection` axis let h be the rounded
  direction, t the rounded quotient, y the real `origin+t*h`, and z=RN(y)
  its single FMA head. L is a positive lower bound for the true absolute
  denominator. The construction's determinant and division budgets give
  `|Z-y| <= u |h| height / (L (1-u)^4)`: one factor covers replacing the
  real direction by h, and at most three rounded nonnegative height operations
  cover the budget evaluation. See the [retry and rescue derivation](exact_selective_retry.md).
  The coordinate band keeps all nonzero comparison operands and products
  normal and finite. The strict comparison of rounded products implies
  `|h| height < |z| threshold` by monotonicity. The dyadic
  `1-(2K+16)u` is exact; rounding its product with K, then with L, gives
  `threshold <= KL [1-(2K+16)u](1+u)^2`. Since
  `(1-18u)(1+u)^2 < (1-u)^4`, acceptance implies `|Z-y| < Ku |z|`.
  Lastly `|y-z| <= u |z|/(1-u)`, so the stated enclosure follows.
  A zero candidate cannot pass a strict test against zero.
* Replacing the coarse division budget by the FMA residual, or tightening
  original determinant budgets with certified tails, changes only the height.
  Optimistic retry screens do **not** accept an axis. Denominator rescue
  requires a newly separated determinant and recomputes L, t, the residual,
  both candidates and the threshold before applying the same strict test.
* The bounded policy only attempts this construction for binary64 source
  coordinates in its exponent band. For each rejected axis (including a
  zero candidate), and for other supported native formats, `point_state`
  goes through `refine_value`. A cached midpoint-certified expansion
  quotient or directly rounded fixed quotient gives z=RN(Z), hence the
  stronger `|Z-z| <= u H(z)` inside the predicate guard. This includes a
  nonzero real rounding to zero. Lossy native conversions cannot enter
  expansion refinement as if exact. Later refinement or float emission can
  replace the binary64 head only by nearest binary64, not by a widened
  float. A failed construction is not a certificate. Nonfinite or out-of-
  measure-range heads are never usable here.
* Matching fixed queries have an exact or independently nearest head. An
  arbitrary auxiliary point, foreign record, non-exact native
  conversion, or a mere lattice-distance promise is not a producer certificate.

Now redo the **difference** estimate with the unchanged scale S. Endpoint
displacements cost at most `(K+2)u S` and subtraction costs at most `u S`.
Thus `c=K+3` suffices; two exact endpoints still use `S=|d|` and cost at
most `2u S`. In the bilinear proof above the same monotone remainder bound
gives `R < 2c+2.5 < 2K+9`, with no assumption of independent errors for
repeated coordinates. The computed magnitude T and subnormal losses are
unchanged. In units of `u T` the true error is therefore below
`(2K+9)(1+8u)+4`, whereas the returned `RN((2K+16)u T)` is at least
`(2K+16)(1-u)-1`. The latter is larger even at the cap (where the slack
is smallest). All-exact operands and eager nearest keep the separate 12u
proof; coordinate comparison remains a lattice-order test.

Two witnesses at K=100 illustrate the distinction. For native a=(0,0),
b=(1,1) and a constructed p=(1,1+2^-44), the positive determinant is
512u; the scale is about 2, so 216u times the scale certifies and 416u
times the scale does not. Conversely take a native origin and a constructed
(1,0) head whose true x is exactly `1+100*2^-52`. This is 100 legal
outward steps. Its squared distance differs from 1 by **more than 400u**,
while the producer-only bound on those heads is 216u. The generic 416u
bound encloses it. Thus the tighter coefficient is not a replacement for
the public K-step theorem.
