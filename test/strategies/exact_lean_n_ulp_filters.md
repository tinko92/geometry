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
Even K=0 leaves half a rounding cell of uncertainty for constructed inputs.
Lossy native conversions cannot use the exact flag; they need their own
certificate before using the constructed flag.

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

For a component `d=RN(b-a)` let S be the real sum of `|d|` and H of each
constructed endpoint; omit both H terms for two native endpoints. A normal
subtraction result has posterior rounding error at most `2u |d|`; a tiny
difference is exact because both operands are integer multiples of eta.
Thus `|D-d| <= C S` for the true difference D. This is an absolute bound:
near cancellation cannot hide the uncertainties of the endpoints.

## Bilinear measures and signs

Use two products with the following pairing and final operation:

| Kernel | Real expression |
| --- | --- |
| Signed determinant / orientation | `(b-a)x(p-a)` (2D cross product) |
| Dot sign | `(b-a).(p-b)` |
| Squared distance | `(b-a).(b-a)` |

Let M be the sum of the two products of the corresponding real S factors.
For each product the two linear perturbations and their quadratic term cost
at most `(2C+C^2) S1 S2`. This uses no independence assumption when a
point, or an entire vector in squared distance, occurs twice.

Each rounded operation on real t has error at most `u |t| + eta/2` in
range. The two rounded products have total error at most `u M + eta`, so
the absolute sum of their heads is at most `(1+u)M + eta`. The final
addition/subtraction therefore brings the arithmetic error to at most
`(2u+u^2)M + (1.5+u)eta`. Since `C < 1/2`, the total absolute error
of the returned head is at most

```text
(5K+8)u M + 2 eta.
```

It remains to bound M using the *computed* magnitude. Nonnegative sums of
doubles, if subnormal, are exact. Thus each computed S is at least
`(1-u)^2 S` (at most two additions). Each magnitude product can lose a
further factor `1-u` and an absolute `eta/2`; the final nonnegative sum
loses at most a factor `1-u`. Its computed result, floored at m, is T:

```text
T >= (1-u)^6 M - eta,    eta <= 2u T,
M <= ((1+2u)/(1-u)^6) T < (9/8) T.
```

The coefficient `F=(8K+16)u` is exactly representable and less than 1.
The returned error `RN(F T)` is at least `(1-u)F T - eta/2`, hence at
least `(8K+14)u T`. This dominates the true error because
`(5K+8)*9/8 + 4 < 8K+14`. With all operands native the same proof uses
K=0. The coordinate guard keeps every difference at most `2^501`, each
computed S at most `2^502`, and the sum of magnitude products at most
`2^1005`. All heads and errors are finite, including at the admitted
endpoints. Underflow of a value product, magnitude product, or final error
is covered above, not excluded by a relative-error assumption.

A measure success certifies `|true - value| <= error`, not necessarily a
sign. Orientation and dot certify only for `value > error` or
`value < -error`; a zero head, touching bound or other ambiguous case
requests fallback. On every false return the output is unchanged.
