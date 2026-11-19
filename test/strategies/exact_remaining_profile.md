# Remaining Exact-Strategy Costs

Profiled `ad5c0f0cb` against local develop
`aaf0560e520041a833c11ad5d61017b684b46f9a`, Apple M5/Apple Clang 21.
The production headers are unchanged. All optimization experiments are isolated
header overlays, not applied fixes.

## Method

32 serial macOS `sample` profiles: eight workload/operation pairs, double/int64,
exact/develop, three seconds at 1 ms. The binaries are the unchanged final
operation-cache benchmark builds. The parser subtracts child counts to obtain
exclusive samples and checks conservation of the main-thread sample total.
It reports primitive groups, algorithm stages and inclusive subtrees separately;
these views must not be added together. Percentages below are approximate.

Separate instrumented headers count one eight-input corpus per operation. They
are never used for timings. The double control applies an exactly representable
dyadic affine transformation, so inputs have fractional coordinates without
changing their topology. Counts, sizes, raw-profile hashes and exploratory
timings are saved alongside this report. Raw stacks remain in
`/tmp/bg-exact-profile-20261113`.

## Common Cases

There is little final exact arithmetic to optimize in the ordinary convex calls:

| Per Call | Convex Union | Areal Relate | Linear Relate | Point/Polygon |
|---|---:|---:|---:|---:|
| Native side predicates | 266 | 264 | 256 | 2 |
| First-filter successes, double | 266 | 264 | 256 | 2 |
| Constructed intersections | 2 | 2 | 64 | 0 |
| Owning native-point constructions | 299 | 138 | 352 | 0 |
| Native floating-view validations | 798 | 792 | 800 | 192 |
| Initial native/approximate coordinate comparisons | 411 | 256 | 2711 | 0 |
| Coordinate refinements | 0 | 0 | 0 | 0 |
| Precision-cache entries | 0 | 0 | 0 | 0 |

The fractional control also has 100% native-side first-filter success and no
coordinate refinement in these cases. Int64 native sides use the bounded native
integer path, not a floating-point orientation filter. Output rounding is
separate: it still evaluates expansions even when all predicates finish early.

Important observations:

- **Linear relate:** point comparison accounts for about 20-23% of samples,
  scalar-bound comparison 12-13%, `fegetround` 15-20%, and stack probing 4-6%.
  The common comparisons are nevertheless resolved without refinement.
- **Point/polygon:** 192 floating-point views are validated although only two
  orientation tests are needed. The validator itself accounts for about 62% of
  sampled call time. Merely forcing it inline did not improve the timings.
- **Disjoint union:** zero intersections, refinements or exact materializations,
  but 650 owning point constructions and 1151 native coordinate comparisons.
  Exact is about nine times slower than develop here. Exact ring processing is
  about 35-39% and sectionalization about 21% of sampled time.
- **Output:** identifiable output-conversion subtrees account for about 18-19%
  of convex union and 39% of linear/areal intersection.

The representation also has a measurable cost:

| Object | Bytes, Both Coordinate Types |
|---|---:|
| Native point / native view | 16 / 16 |
| Implicit owning point | 128 |
| Exact squared-distance descriptor | 280 |
| Exact / ordinary overlay traversal turn | 984 / 376 |

The distance descriptor owns two complete points. Two such descriptors are part
of each enriched overlay turn. These sizes do not imply that all the extra bytes
are actively used in every operation.

## Small Proven Opportunities

Each isolated patch preserves the existing exact primitive or adds an exact
shortcut. The final experiment uses three serial, rotated 0.2-second samples
per case: 384 timing records. All exact variants match the original benchmark
digests. These are exploratory timings, not a full regression qualification.

| Change | Why Exactness Is Unchanged | Measured Benefit |
|---|---|---|
| Remove the unreachable final coordinate fallback | `compare_points` already returns an exact result or throws; every return is certain. Its result is retained unchanged. | Coordinate comparator stack frame drops from 8672 to 96 bytes; linear relate about 7% faster for double, 10% for int64. |
| Query the rounding mode only once on a successful scalar filter | The scalar filter already checks the mode. On failure, the wrapper still selects the same uncertain/interval path. | Removes a duplicated query on every successful implicit-coordinate filter; about 4%/2% in linear relate in this run. |
| Add the missing support-line first filter | A certified nonzero sign proves non-membership. Ambiguous cases retain the old expansion/fixed fallback. | Star symmetric difference about 1.25x/1.27x faster. |
| Recognize a supporting segment's own endpoints | After excluding a degenerate line, either endpoint belongs exactly to that line. No arithmetic is needed. | Star symmetric difference about 1.39x/1.45x faster. |
| Skip native orientation's intermediate interval evaluation when its expansion stage is supported | The same exact expansion supplies the answer, including zero. Unsupported ranges/modes retain the original path. | Shared-boundary double relate about 1.76x faster. |

The support-line helper currently sends even small integer inputs through
floating expansions. Reusing the existing checked native-integer side primitive
there is another exactness-preserving opportunity; it was not benchmarked here.

The experimental bundle gives:

| Case | Double Speedup | Int64 Speedup |
|---|---:|---:|
| Convex union | 1.01x | 1.02x |
| Disjoint union | 1.05x | 1.07x |
| Convex linear relate | 1.14x | 1.12x |
| Star symmetric difference | 2.15x | 2.28x |
| Shared-boundary areal relate | 1.80x | 1.10x |

The bundle includes the unsuccessful validation-inlining experiment; it is not
a recommended final patch set. Point queries are about 3-6% slower in that bundle.
Do not add individual speedups. Small differences are within desktop timing
variation, and code layout changes affect cases outside the edited path.
Point and native-predicate unit suites pass on the bundle, including their
rounding-mode, nonfinite-input and full-range controls. The full algorithm suite
and frozen issue replay have **not** been rerun on these prototypes.

## Higher-Value Follow-Ups

1. **Keep native inputs native in sectioning and ring scans.** The strategy
   currently chooses a 128-byte exact point even for native section temporaries;
   `static_compare_coordinate` also wraps native inputs in owning points.
   Select small native temporaries/views when coordinates are native, while
   preserving implicit points and mixed-type exact comparisons. This changes
   representation, not any predicate. Adapted/non-default-constructible points,
   closure traits, and temporary ranges need coverage before applying it.
   Unmodified source rings also retain their valid-input orientation and hole
   hierarchy; avoid re-proving those properties by scanning every coordinate.
   A strict conservative-box-disjoint certificate can select direct copying,
   concatenation or an empty result, rather than building and analyzing exact
   output rings, under the appropriate valid-input/output-type contracts.
2. **Do not look up an unbound cache handle.** `cached_bounds` currently asks for
   thread-local operation state before checking the point's index. A sentinel
   index cannot have cached bounds. Checking it first avoids TLS lookups in the
   ordinary no-refinement case; `_tlv_get_addr` alone is about 4% in int64 linear
   relate. This was identified but not separately benchmarked.
3. **Validate once at a suitable boundary, then use checked native views.**
   Preserve standalone predicate validation and invalid-input behavior. Reuse
   validated query/endpoints rather than repeating finite checks for every edge.
   Do not remove validation or assume a rounding mode globally.
4. **Avoid repeated output work.** Compute both numerators with one shared
   denominator, reuse certified rounded cache results, and copy a closed ring's
   already converted closing coordinate according to its closure traits.
   Integer-output rounding cells can often be certified from existing bounds.
   Final nearest/ties-to-even conversion must remain correct; merely copying an
   approximate intersection is not an optimization with the same contract.
5. **Reduce unused turn payload.** Avoid default-initializing/copying unused
   distance descriptors, or use operation-lifetime references/indices where
   lifetime is proved. Native-support outputs must remain self-contained.
   Do not replace general distance order with coordinate order without proving
   the compared displacements are on the same line.

### Symmetric Difference

`sym_difference_areal_areal` currently computes two differences and unions their
results. Exactness gives a useful conditional optimization. For regular closed
sets, write `DA = closure(interior(A \\ B))` and analogously `DB`. Then
`DA intersection DB` is a subset of `boundary(A) intersection boundary(B)`:
a point in the interior of either original set cannot be approached by the
other difference. If the original boundaries intersect only at finitely many
points, the two valid difference outputs have disjoint interiors and no shared
edge. Their multipolygon concatenation is already the regularized union.

Thus a certificate of **no positive-length original boundary overlap** could
avoid the third overlay for exact strategies. It must be carried from actual
turn classification, not guessed from approximate bounds. Simply concatenating
in the shared-edge case is wrong, and rounded intermediate geometries do not
have the same guarantee. This is an algorithm-side proposal, not an implemented
or timed optimization in this report.

### Boost.Graph

Graph detection/search occupies about 5-6% of exact convex union and 2-3% of exact
star symmetric difference in these samples. Its absolute cost is similar to
develop's, whose graph percentage is higher because the rest is faster. It is
not the main source of the remaining multi-fold gap; replacing it is not the
first performance change to make.

## Evidence And Reproduction

- [Profiles and hashes](exact_remaining_profiles.json)
- [Integer/fractional double counts](exact_remaining_counts_double.csv), [int64 counts](exact_remaining_counts_int64.csv)
- [Object sizes](exact_remaining_sizes.csv)
- [Comparator assembly](exact_remaining_comparator_asm.txt)
- [Experimental medians and checks](exact_remaining_experiments.json), [all timing records](exact_remaining_experiment_records.csv)
- [Unapplied overlays](exact_profile_patches)

Recreate overlays from this commit's unchanged production headers, e.g. for
`compare`, and repeat for the other patch names:

```sh
out=/tmp/bg-profile-new
mkdir -p "$out/experiments/compare"
cp -R include/boost "$out/experiments/compare/"
patch -p1 -d "$out/experiments/compare" < test/strategies/exact_profile_patches/compare.patch
```

Use the same compiler flags as the benchmark manifest. The instrumented overlay
is separate (`$out/instrumented`) and is used only by `exact_profile_counts.cpp`.
`build_exact_profile_experiments.sh` builds the exploratory binaries;
`profile_exact_cartesian.js` collects samples, summarizes them, runs serial
experiments and exports the overlays. The original benchmark directory must
contain the matching exact/develop binaries. Do not compile during timed runs.
