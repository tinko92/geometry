# Borrowed side and dot: shared X/Y reads

Base: `5a3405c9a`. This is a local experiment for the stored-error binary64
first filter, not a replacement for numerical variants throughout Cartesian
predicates. No server measurements were taken.

## Route and numerical scope

The public `cartesian_exact<double>` umbrella binds native geometry sources
to a borrowed Root. Its `side().apply` and `side().direction_code` call
`static_orientation` and `static_direction`, respectively, and then the
typed Root visitor. The experiment is **inside** that visitor. Its supported
triples have at least one selected construction record; the other operands
are selected records or native values losslessly representable as binary64.
It supports `approximate_output` and the stored-error bounded-ULP policies.
Lean policies, native-only triples, wide integer conversion, value-producing
iterator proxies and auxiliary queries still use their old first filter.

The old filter already short-circuited to a native bilinear bound when all
six errors were zero. It also already coalesced X/Y loads for the ordinary
native pointer on Apple Clang. It did **not** eliminate the construction's
extra representation test, its generated native-coordinate visitor path,
or the three by-value 32-byte approximation packets passed to the general
sign filter. Check the construction alternative of `side_nni`, not just
the wrapper's outer variant dispatch, in optimized assembly.

The new filter retains that outer typed dispatch. It obtains one reference
to the selected record's numeric state (with the existing checked record
access), and dereferences a lvalue source iterator once, before fetching the
two axes. It passes no owning point or approximation packet to the kernel.
It uses the same rounded differences, `8u` exact-input bilinear bound, and
inflated uncertain-input error polynomial; both routes share the numerical
helpers. A published stored error of zero only skips error propagation when
all relevant input errors are zero. It does not change classification of
the point for original-support, intersection identity, or topology tests.
The stored errors are nonnegative; any nonfinite ordinate propagates through
a difference and a product, making the sign filter decline. The existing FP
contract, including its arithmetic range restriction, still applies.

After a miss, the old support/identity checks, refinement, retry and fixed
arithmetic remain in the same order. No state or descriptor layout changed.
The experiment does not address shared reads for lexicographic comparison
or equality, and does not change standalone legacy exact points.

## Local checks

Clang, GCC 16 and Clang ASan/UBSan passed. `exact_borrowed_predicates.cpp`
calls the public umbrella's side/dot in the 872 rational-oracle triples and
retains the 160 measure checks. Additional
tests assert one dereference per lvalue source iterator on a first-filter
hit, with a certified X and an uncertain Y, for both stored-error policy
families. The same cases exercise a rounded coordinate that masks a nonzero
exact result, a true zero and two certified zero coordinates. The 37
allocation checks recorded zero allocations. The test also reports the
unchanged borrowed point/record sizes. `exact_forward_error` and
`exact_borrowed_public_predicates` also passed.

## Local A/B

Apple Clang 21, arm64, `-O3 -DNDEBUG -ffp-contract=off`, no statistics.
The predicate driver calls the public strategy with a polygon Root. It uses
256 triples per loop, 450,000 calls per sample, seven serial repetitions
with alternating base/candidate order. Median ns per full sign:

| Operands | Side base / candidate | Dot base / candidate |
| --- | ---: | ---: |
| Native pointers and one uncertain construction | 9.56 / 3.76 | 9.20 / 4.26 |
| Lvalue iterator, borrowed native, uncertain construction | 13.45 / 5.63 | 10.78 / 5.80 |
| Three uncertain constructions | 12.47 / 7.69 | 11.94 / 8.06 |
| Native pointers and one zero-error construction | 8.67 / 3.12 | 11.18 / 6.66 |
| Lvalue iterator, borrowed native, zero-error construction | 13.85 / 4.86 | 12.97 / 8.15 |
| Three borrowed natives (unchanged code path) | 4.85 / 4.82 | 8.66 / 8.97 |

Every per-case checksum matched. The uncertain fixtures use constructed
points with stored positive error; the zero-error fixture refines the
origin through the normal point state. These timings include visitation,
not construction; the arrow does not represent an end-to-end speedup.

`exact_cartesian_benchmark.cpp` (`-std=c++20 -DBG_BENCH_EXACT`, 96 vertices,
8 corpus members) matched both corpus and geometry digests in all 27
projected, projected-local and shared operation cells. Both digest and
projected self-tests passed on the candidate. The second 55 ms window of
five serial, order-alternating pairs yielded these median us per call:

| Workload | Base | Candidate |
| --- | ---: | ---: |
| Projected union | 42.46 | 41.83 |
| Projected symmetric difference | 290.43 | 266.18 |
| Shared union | 54.72 | 54.28 |
| Shared intersection | 52.17 | 51.81 |
| Projected areal/areal relate | 13.12 | 13.10 |
| Projected linear/areal relate | 9.72 | 9.42 |

The small deltas are not a general speedup claim. In particular, the
unchanged native control is a measure of local variation. The full
Boost.Geometry suite and other CPU architectures were not tested.

Local artifacts in `/private/tmp`: `bg1517-xy-filter-probe.cpp`,
`bg1517-xy-filter-base.s`, `bg1517-xy-filter-new.s`, `bg1517-xy-filter-micro.cpp`,
`bg1517-xy-filter-run.sh` and the three `bg1517-xy-filter-*.tsv` result files.
