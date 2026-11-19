# CPU qualification on develop plus PR 1517

## Final results

The measured implementation is `49a77482d` (locally retained as
`refs/experiments/20261119/xy-combined-final`). Its headers and test
sources match the PR-ordered experiment branch. The two CPU-only x86-64-v2
runs used strict FP, identical C++17 flags and native double geometries.
Ratios are time / same-compiler develop (`082b99181`), lower is better.
Each of 45 operation/workload cells has equal weight; each time is the
median of three samples. N counts representable steps from the correctly
rounded result, not decimal digits.

| Cartesian, binary64 | GCC 13 | Clang 22.1.8 |
|---|---:|---:|
| develop | 1.000 | 1.000 |
| Simple (uncertified) | .599 | .685 |
| Exact topology, emit stored approximation | .856 | 1.128 |
| Exact topology, 100 ULP steps | .890 | 1.179 |
| Exact topology, 4 ULP steps | .950 | 1.255 |
| Exact topology, 1 ULP step | 1.054 | 1.408 |
| Exact topology, correct rounding, eager | 1.237 | 1.654 |
| Exact topology, correct rounding, deferred | 1.321 | 1.779 |

All exact policy fingerprints agreed. Simple disagreed on three XOR cells
and develop on one; timings alone do not establish numerical correctness.
The 100-step policy costs about 4% over stored output on either compiler.
On decimal inputs its ratio to develop was 1.165/1.386 (GCC/Clang);
tightening to four steps gave 1.368/1.638 and to one step 1.691/2.071.
On the projected corpus all three budgets took about .63/1.00.

With exact stored output, absolute Clang time was .973 times GCC time;
on develop it was .738. Union is the largest remaining same-compiler
penalty (1.100/1.441). In profiles, search and section construction, turn
generation, supported coordinate ordering and vertex copying all matter;
the samples do not implicate a single cold exact fallback. Forced inlining
and a direct min/max specialization did not shift the box trial by even .4%.
The source-box Sterbenz dispatch did pay off: .930/.930 time vs disabled
on the projected parcel from #1281; the road-window example took .844/.924.
The translated parcel was a negative control, at 1.017/1.015.

The three-hour spherical areal run passed 390,775 pairs (3,126,200 calls);
the separate mixed run passed 52,456 pairs and the final-header mixed replay
added 6,586 fresh pairs. The tests check exact rays and the directed region
result, not a lossless angular encoding of every generated vertex. The
certified Cartesian buffer reproduced nesting for #578 and removed the two
spurious components in the #1281 dilation/erosion pipeline. It may decline
unrepresentable native exports. The explicit multiprecision buffer fixes
the full `rt_w6` and `nores_wt_1` regressions. Ordinary and opt-in
paths are identified in the separate issue replays; these statements are not
unconditional guarantees about the ordinary buffer.

Registered GCC 13 release B2 relate/set/buffer tests passed in C++14 and
C++17 (71 and 75 rebuilt executables on the respective qualification runs).
The combined pre-X/Y snapshot passed 78 explicit GCC strategy tests.
Four changed-path executables on the final X/Y snapshot passed on each
compiler, including the 872 rational predicate triples. The separate
46-source non-Jamfile sweep retained its baseline classification: 39 pass,
three baseline compile failures, one identical runtime failure and three
changed runtime logs. Native float difference on `case_precision_m2_s`
already changed on PR 1517; the explicit exact result recovers the expected
area. No unconditional library default has been changed.

After building the benchmark, use `run_simple_cartesian_benchmark.sh` to
repeat the policy matrix (categories `degenerate general decimal projected
projected-local`, all bounded policies enabled). The 45 cells comprise
union, intersection, difference, symmetric difference, line clipping,
covered-by and three relation calls. See `exact_optimization_audit.md`
for the per-optimization control experiments and proof boundaries and
`exact_integration_order.md` for proposed review units.
