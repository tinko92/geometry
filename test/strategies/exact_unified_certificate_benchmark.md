# Eager nearest and producer-certified arithmetic filters

Compared `cbe774cd9` against its parent `f6b03306a`. The old
`nearest_output` used deferred rounding; old explicit eager is the control
for the filter refactoring. On root-bound constructions the new
`nearest_output` and its compatibility alias `eager_nearest_output` are
the same type. See [filter certificates and scalar order](exact_lean_n_ulp_filters.md).

## Public operations

CPU-only Intel Xeon Gold 6530; GCC 13.3.0 and Clang 22.1.8. C++17,
-O3, -march=x86-64-v2, -fno-fast-math, -ffp-contract=off. Each compiler
used one pinned logical CPU for all timed runs. Three double input pairs
per cell; five workloads (degenerate, general, three-decimal, projected,
and the projected fixture translated near the origin); nine operations.
Each of 13 binaries ran five interleaved 0.14-second windows per cell,
2,925 windows per compiler. Input and output checks preceded timing.
Summaries are equal-cell geometric means of median paired wall-time ratios;
below one means faster. These results are against the previous experiment
branch, not develop.

| Comparison | GCC | Clang |
| --- | ---: | ---: |
| New nearest / old nearest (deferred) | 0.9183 | 0.9278 |
| New nearest / old explicit eager | 0.9968 | 0.9996 |
| Lean 1 / previous lean 1 | 0.9997 | 0.9998 |
| Lean 2 / previous lean 2 | 1.0027 | 0.9992 |
| Lean 4 / previous lean 4 | 1.0016 | 1.0005 |
| Lean 8 / previous lean 8 | 1.0010 | 1.0013 |
| Lean 100 / previous lean 100 | 0.9992 | 0.9999 |

| New nearest / old deferred, by operation | GCC | Clang |
| --- | ---: | ---: |
| Union | 0.882 | 0.892 |
| Intersection | 0.854 | 0.868 |
| Difference | 0.878 | 0.887 |
| Symmetric difference | 0.839 | 0.861 |
| Relate, area/area | 0.982 | 1.002 |
| Relate, line/area | 1.016 | 1.028 |
| Relate, line/line | 0.985 | 1.010 |
| Covered by | 1.000 | 1.001 |
| Clip | 0.852 | 0.828 |

The gain is from choosing eager instead of deferred, not from making an
already-eager policy appreciably faster. The largest workload difference
was on projected coordinates (0.768 GCC, 0.763 Clang); translated near the
origin, the same fixture was close to parity (0.992 and 1.003). Eager is
not universally faster: line/area relate took longer on both compilers.

## Correctness

The contract is round-to-nearest, ties-to-even. The two compilers passed
25 focused test executables plus the fixed-rounding boundary cases under
that mode. The suite includes 7,998 construction checks, 1,064 rational
predicate triples, and 192 measure comparisons. It tests binade boundaries,
subnormals, midpoint ties, and rounding directly to float to avoid double
rounding. It also exercises the scalar shortcut for unequal nearest heads,
exact fallback on equal heads, and rejection of the producer certificate
for foreign records and lossy integer conversions.

Separately, all 135 public workload/input/operation combinations produced
byte-identical lossless hexfloat transcripts across old eager, old deferred,
and new nearest on *both* compilers, including output structure and DE-9IM
strings. The timing driver's lean comparisons check same-budget structure
and spatial digests (quantized), not bitwise equality across budgets.
