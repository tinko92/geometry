# Exact Cartesian Redesign Checkpoint

This is an implementation checkpoint, not completion of the agreed redesign.
All work is local on `codex/exact-cartesian-static-kernel`, based on `2886f2acd`.
The independent develop baseline is `aaf0560e5`; the PR1491 baseline is `72b2d832f`.

## Contract

Native integral coordinates up to 64 bits, float and double; immutable, finite,
valid Cartesian inputs. Floating overflow/underflow and long double are outside
scope. IEEE round-to-nearest, controlled contraction, no reassociation/fast-math.
Lossy coordinate-domain promotion is rejected. Internal auxiliaries belong to
the root operation; there is no persistent implicit-result contract.

The default exports the approximation currently stored in each constructed point.
Native vertex selections convert the original coordinates directly. Export does
not refine for output. Integral destinations use BG `numeric_cast` of
that approximation, including under `nearest_output`; endpoint-relative integer
displacement rounding is not reproduced. `nearest_output` correctly rounds
floating destinations. Neither policy promises validity of the exported native
geometry when distinct exact coordinates cannot be represented distinctly.
Ring closure reuses the first emitted native point. Internal exact closing
occurrences retain their own supporting-line provenance.

## Predicates And Constructions

Normalize homogeneous points as \(P=(X_P,Y_P,W_P)\), \(W_P>0\).
An input double denotes its exact binary rational value.

| Operation | Exact Expression / Rule |
|---|---|
| Coordinate order | \(\operatorname{sign}(X_PW_Q-X_QW_P)\), or the analogous Y expression |
| Point equality | Both coordinate comparisons are zero |
| Orientation | Sign of the determinant with rows \((X_P,Y_P,W_P)\), \((X_Q,Y_Q,W_Q)\), \((X_R,Y_R,W_R)\) |
| Order along a segment | Select a nonconstant coordinate and reverse its comparison if the segment runs in the decreasing direction |
| Direction | \(\operatorname{sign}((B-A)\mathbin{\cdot}(P-B))\) |
| Squared-distance order | Compare \(\lVert B-A\rVert^2\) and \(\lVert D-C\rVert^2\) by positive-denominator cross multiplication; no square roots |
| Signed distance measure | Compare determinants \(\det(B-A,P-A)\), including requested sign reversal |
| Segment classification | Exact orientation, coordinate order and endpoint equality; overlaps select existing endpoints |
| Proper intersection | \((A\times B)\times(C\times D)\) in homogeneous coordinates |
| Midpoint query | \((X_PW_Q+X_QW_P,Y_PW_Q+Y_QW_P,2W_PW_Q)\) |
| Coordinate splice | \((X_PW_Q,Y_QW_P,W_PW_Q)\) |
| Box/envelope predicates | Coordinate extrema and orientation against original supports; conservative search boxes are not result geometry |
| Ring topology | Exact boundary cancellation, orientation and point location; no signed-area ordering for parents |

Common predicate schedule:

```text
native/native coordinate order: compare original coordinates directly
native floating orientation: existing side_robust, exact equality, robustness 3
native integral orientation: compile-time safe width, checked wider width, fixed fallback
mixed/implicit predicate:
    evaluate the approximation and specialized forward-error bound
    return a sign only when certified
    exploit exact identity/support incidence when applicable
    refine reusable numerator/denominator expansions through a cache index
    retry the filter; otherwise evaluate the fixed homogeneous expression
```

There is no general interval-arithmetic orientation stage or runtime expression
tree. Native/native/implicit orientation can certify zero on an original support
without owning native wrappers or numerical-cache growth. This does not treat
an arbitrary third concurrent line as a recorded producing support.

For construction, let \(u=B-A,v=D-C,w=C-A\),
\(N=\det(w,v),D=\det(u,v),t=\mathrm{fl}(N/D)\).
Forward bounds \(E_N,E_D\) give
\[
E_t\ \leq\ \frac{E_N+|t|E_D+u_{\rm fp}|N|}{L},
\qquad 0<L\leq |D|-E_D.
\]
The coordinate approximation is `fma(t,u[k],A[k])`; its forward bound includes
endpoint error, direction error, parameter error and the fused rounding error.
The code inflates nonnegative bound evaluation by a proved static coefficient.
A zero/uncertain rounded denominator refines at construction, not first output.

## Current Storage

Members only, omitting access control and aliases:

```cpp
// Borrowed native operands; no cached coordinate snapshot.
Point const& native_point_view::m_point;
Iterator iterator_point_view::m_iterator;

// Transitional production implicit point: STILL OWNS PRODUCER COORDINATES.
bool exact_point<C>::m_intersection;
std::array<C, 8> exact_point<C>::m_values;
mutable std::array<interval, 2> exact_point<C>::m_bounds;
point_handle<C> exact_point<C>::m_precision; // one size_t index
mutable native_cache_handle exact_point<C>::m_cache; // one size_t index

// One typed numerical workspace per root operation.
storage_for<C> context<C>::m_storage;
std::vector<native_entry> point_storage<N, D>::native;
std::vector<fixed_homogeneous<N + D + 1, 2 * D + 1>> point_storage<N, D>::auxiliary;

// native_entry
fixed_homogeneous<N, D> value;
bool exact_ready;
std::array<expanded_coordinate, 2> expanded;

// expanded_coordinate
expansion<32> numerator;
expansion<16> denominator;
bool attempted, valid, bounded;
interval bounds;
```

Expansion constructors/copies touch only active components. Fixed integer
temporaries use statically bounded storage and active-limb arithmetic. Cache
growth may allocate; numerical predicates do not allocate arbitrary big integers.
Independent callers use independent workspaces. Thread-local state holds only
the active-context pointer and optional diagnostic counters, not numeric storage.

For native grid bounds \(|xQ|<2^b,Q=2^s\), native intersections need
\(N=3b+3,D=2b+s+3\) magnitude bits. One midpoint/splice of native intersections
fits the declared auxiliary budget. That is not a proof for arbitrary recursively
constructed auxiliaries; narrowing is checked rather than silently truncated.

## Outstanding Requirements

1. Bind immutable root source types and retain actual turn/clip owners through
   final emission. Remove owning native leaves, endpoint tuples and distance-key
   copies. A pointer substitution alone is unsafe: turns move and intermediate
   geometries are destroyed during composed operations.
2. Preserve original outgoing edge supports, not only vertex producer pairs.
   Equal coordinates at triple concurrency can have different useful supports.
3. Integrate policy-specific early-ULP point layouts into the algorithm. Current
   `nearest_output` still uses the stored-error production point layout. The
   eager/deferred ULP experiment is separate; it is not a production comparison.
4. Finish construction-origin/capacity proofs across composed overlays and
   collections. No unrestricted rational-geometry API is claimed.
5. Qualify remaining shared legacy failures with spatial oracles. New diagnostic
   differences have been classified, but this is not a clean legacy exact suite.

BGL remains: exact predicates do not remove genuine articulation/sector work.
No local replacement is claimed without an independent equivalence proof.

## Qualified Checkpoint Results

Final headers match the frozen issue-replay snapshot, SHA256
`ef0814503177bf6d514e9a08eeaed61aa4b2b70f7ede9fb5441ee636e124976b`
(the runner hashes relative filenames and contents).

- All 34 primitive/test-driver sources pass; complete start/end header hashes
  match. Final mixed-support and borrowed-ring tests pass ASan/UBSan, including
  noncopyable inputs, open/reversed ranges, mutable-dereference cursors and
  unadapted convertible proxies.
- Ordinary C++14 suite: 37 passed, four runtime failures identical to baseline,
  three matching baseline compile failures. No new default-strategy regressions.
- Issue replay: 36 successful compilations, 656 records, 82 fixtures, 47 issues.
  No execution/evaluator failures or current internal-exact observations. Both
  policies' 82 retained records and 82 native records are unchanged by the last
  support/iterator fixes. Native representation/feedback observations remain
  separate: 27 approximate, 23 nearest. These are not counts of closed issues.
- The legacy exact-enabled report has 33 passing sources, three shared compile
  failures, two unchanged runtime failures, five runtime differences and one new
  failing source. Of 56 added diagnostic messages, 37 concern integral conversion,
  16 naive floating export and three ordinary-validator false positives. All 54
  representative exact outputs and four additional int64 outputs are spatially
  identical to the previous exact baseline. Invalid `ticket_9942` inputs remain
  flagged; shared baseline assertions have not all been independently requalified.
- An independent BigInt rational arrangement oracle checks all 86 open 2D cells
  for #876: both exact winding variants equal A-minus-B. Naive CCW export adds
  about 7.09168107 area units; nearest avoids this large discrepancy. That is an
  export-accuracy difference, not a changed exact result.

Apple M5, Apple Clang 21, C++17, O3, NDEBUG, contraction off; five serial rotated
rounds, 0.1 seconds per case. Median per-case timings, then equal-weight geometric
means. There are 45 workload/operation cases per coordinate/control family.
Allocation counting is a separate run. These are local ARM measurements.

| Default Output | Speedup / Previous Exact | Time / Develop |
|---|---:|---:|
| Double, integer-valued inputs | 2.83x | 1.66x |
| Int64 inputs | 1.59x | 2.30x |
| Double, dyadic fractional control | 2.76x | 1.75x |

For double, polygon set operations improve 2.60x and relation matrices 2.51x
against previous exact. Final nearest-output adds about 11% across the double
corpus; this does not measure either pending early-ULP production layout.
Allocation counts are unchanged in 44/45 double cases and decrease by one in
the remaining case; no case adds allocations. Star symmetric-difference traffic
falls from 6,627,329 to 4,699,889 allocated bytes per call (not peak live memory).

Profiling caught and fixed a borrowed-NNI support-incidence regression. Star
polygon relate falls from 246 to 43 microseconds; convex linear/areal relate from
46 to 14. Sampled memory-copy/clear self time in the former falls from 38.7% to
2.5%. In overlay, remaining costs include point comparisons, ring-parenting work,
11-13% memory-copy/clear and about 12% allocation self time. BGL accounts for
7-11% inclusive time in the two sampled exact overlay cases, so eliminating it
alone would not remove most of the remaining slowdown. Sampling has unclassified
and folded symbols; these percentages are diagnostic, not cycle-exact attribution.

The standalone eager/deferred ULP experiment has 480 timed records and independent
oracle checks. Eager construction is cheaper when most points are emitted;
deferred construction wins when few are. Both experiment variants eagerly create
precision entries and reserve cache capacity before timing. Those limitations
preclude choosing the production policy from this experiment alone.

Compact manifests/tables are retained outside the worktree under
`/Users/tinko/dev/boost/.geometry-plans/exact-redesign-evidence-20261115`.
Raw snapshots/logs are under `/tmp/bg-exact-redesign-20261115`; the classified
legacy replay and independent oracle sources are retained under
`/Users/tinko/dev/boost/.geometry-plans/exact-postclosure-qualification-20261115`.

## Reproduction

`run_exact_primitives.sh` records start/end header hashes. The issue replay
snapshots headers, fixtures and the evaluator before compiling. Benchmark
manifests pin source, flags, compiler, baselines and header hashes. Paired timings
run serially; allocation instrumentation is measured separately.

```sh
export BOOST_ROOT=/Users/tinko/dev/boost
sh test/strategies/run_exact_primitives.sh /tmp/bg-exact-primitives
node test/strategies/run_exact_redesign_issues.js --prepare /tmp/bg-exact-issues
node test/strategies/run_exact_redesign_issues.js --run /tmp/bg-exact-issues
export BENCH_VARIANTS=develop,previous,optimized,nearest
node test/strategies/run_exact_benchmark.js build /tmp/bg-exact-bench
node test/strategies/run_exact_benchmark.js paired /tmp/bg-exact-bench
node test/strategies/run_exact_benchmark.js allocations /tmp/bg-exact-bench
node test/strategies/run_exact_benchmark.js summary /tmp/bg-exact-bench
```

Benchmark digests check component/hole/point counts and relation matrices, not
coordinate or spatial equivalence. Native-output representation/feedback checks
are reported separately from internal exactness and independent exact oracles.
