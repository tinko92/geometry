# First-stage binary64 predicates for step-certified coordinates

The `ulp_filter::try_*_n` functions in `ulp_filter.hpp` use the following
contract, independently on each axis. A native input is the exact real number
represented by its double. For a constructed input with head `z`, let `r`
be the finite binary64 nearest-even rounding of the true coordinate `Z`.
The number of adjacent representable steps from `z` to `r` is at most
`K = min(N, 2^24)`, with the two signed zeros counted once. For policies
larger than the cap, the **producer** must certify the capped distance.
Silently capping a merely N-step external guarantee would be unsound.

These kernels do not require the tighter absolute error certificate available
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
Deferred and legacy nearest-output policies retain their adjacent-value
kernels; native-only orientation and scalar comparison retain typed dispatch.

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
