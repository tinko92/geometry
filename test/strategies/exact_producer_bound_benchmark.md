# Using the constructor's direct bound in lean predicates

Baseline: `811718814`. Same construction algorithm, retry policy, scalar
ordering and output contract. Only the bound on a bilinear result changes.
Where every uncertain operand has the matching producer certificate, its
coefficient falls from `(4K+16)u` to `(2K+16)u`. See the
[bound, admission rules and counterexample for raw step certificates](exact_lean_n_ulp_filters.md).

## Public operations

Intel Xeon Gold 6530, CPU only, GCC 13.3.0 and Clang 22.1.8, C++17, -O3,
x86-64-v2, no fast-math or contraction. Baseline and candidate use the
same public driver and three binary64 input pairs per cell. The 45 cells
are nine operations by five workloads: degenerate, general, decimal,
projected and projected translated near the origin. Five interleaved
0.15-second samples per binary per cell; same physical core for the pairs,
alternating order. Each public ratio is the geometric mean over cells of
the ratio of median wall times, candidate/baseline. No counters in timed
binaries.

| N | GCC time ratio | Clang time ratio |
|---:|---:|---:|
| 1 | .9997 | .9996 |
| 2 | 1.0007 | .9997 |
| 4 | 1.0014 | .9988 |
| 8 | 1.0008 | 1.0002 |
| 100 | 1.0000 | 1.0006 |

There is no measurable public-operation speedup on this corpus. Every ratio
is within 0.14% of parity; individual repeat ratios vary more than this.
Both compiler runs had CPU/wall time .9999. All 2,250 samples per compiler
passed input and output fingerprint checks against the untimed public walk.
Recompiling with the final proxy admission change gave identical executable
text sections and matching check results for all ten candidate binaries.

In a **separate instrumented walk**, both versions and both compilers made
the same decisions in every cell at every N. These are raw counts for one
complete corpus walk at one budget, not totals over the timed repetitions:

| Decision | Initial filter calls | Initial misses | Refinements | Exact exits |
|---|---:|---:|---:|---:|
| Orientation (excluding native-only bypass) | 1,642 | 128 | 0 | 0 |
| Dot sign | 10 | 0 | 0 | 0 |
| Squared-distance pair comparison | 174 | 105 | 26 | 105 |
| Signed or mixed measure pair comparison | 0 | 0 | 0 | 0 |
| Coordinate order | 24,083 | 231 | 0 | 0 |
| Coordinate equality | 1,821 | 0 | 0 | 0 |

Of the 105 exact squared-distance decisions, 79 used the small-integer
shortcut and 26 needed exact pairs after refinement. Initial sign or
coordinate misses resolved by identity or another shortcut are not exact
arithmetic fallbacks. Zero signed comparisons here means unexercised.

## Near-collinear sweep

For each integer m from -768 to 768 put s=2^-52. Intersect the segment
from (0,1+ms) to (3,1+(m+1)s) with the vertical at x=1. The exact crossing
P is (1,1+(m+1/3)s); every source coordinate is binary64. Write
O=(0,0), D=(1,1), U=(1,0), V=(0,2). Test side(O,D,P), the negative of
(D-V) dot (P-D), det(O,U,P) against det(O,U,D), and |P-O|^2 against |D-O|^2.
Each answer has the sign of 3m+1. Give every probe a freshly constructed
point and operation scope, so a previous exact fallback cannot warm it.

| N | Exact fallback count: old to new | GCC time ratio | Clang time ratio |
|---:|---:|---:|---:|
| 1 | 188 to 172 | .978 | .985 |
| 2 | 220 to 188 | .923 | .938 |
| 4 | 284 to 220 | .869 | .882 |
| 8 | 412 to 284 | .789 | .803 |
| 100 | 3,228 to 1,756 | .559 | .566 |

Counts are for the **four requested probes**, 6,148 in total per N.
Do not count nested calls made while resolving their misses as additional
probes. At N=100 side falls from 833 to 433, dot from 1,537 to 865, and
each measure comparison from 429 to 229. Counts agree on GCC and Clang.
Timings are ratios of median time per probe across seven interleaved
repetitions, including cold construction and the exact-sign check, without
instrumentation. This corpus deliberately populates the decision threshold;
its speedup is not a prediction for typical geometry.

## Correctness

Eleven focused executables passed on each compiler using the final source:
ULP filters, lean bounded output, construction rounding, rounding and output
policies, filter routing, borrowed predicates, bounds, segment/box,
borrowed operations, and shared boundary.

Independent rational oracles check both the constructor's direct bound and
the resulting signs and measure enclosures, including capped budgets, rounded
zero, determinant retries and denominator rescue. The negative control is a
legal 100-step coordinate for which the producer-only squared-distance
bound is false; it must stay on the generic kernel. Output and scalar-order
rules did not change.

The operand audit also found a separate wrong-sign case for tagged iterator
proxies: the view declares float ordinates but the proxy exposes doubles.
Taking a proxy directly bypassed the view's conversion. For the regression
triangle, the correct determinant is -2^-28; bypassing conversion gives
+2^-28. The lean reader now respects the declared coordinate type, and the
stored lvalue shortcut runs only when it preserves that type.
