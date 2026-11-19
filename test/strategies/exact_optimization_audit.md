# Exact optimization TODO audit

Scope: all lines of `exact_optimization_todo.txt` at
`b7daf7d8ff579b512ae97f7b7585a7ca72061332` (use
`git show b7daf7d8f:test/strategies/exact_optimization_todo.txt`
for stable line numbers). Read-only inputs: the local series and its
`refs/experiments/20261119/`, their scratch checkouts, and
`/Users/tinko/dev/boost/.geometry-results/exact-1517-cpu-20261119`
(abbreviated **CPU**). No remote access. A measured experiment is not
called a candidate merely because it is absent from the baseline.
"Open" means not found in these refs/evidence, not a claim about other work.

## Predicate and operation backlog

| TODO lines | Status and branch evidence | Measurement or actual gap |
| --- | --- | --- |
| 1-9 | Scope and priorities, not an implementation. RN/exact fallback and borrowed support lifetimes remain requirements. | Not a pending experiment. |
| 10-18 | Landed: `fb4dab4a9`, rounded difference vs rounded nonnegative error sum, direct exact order and inconclusive ties. Earlier scalar pre-visitor `6987622b4` is distinct. | CPU `scalar-v1` / REPORT: the landed monotone change saved ~2-5% typically on GCC; pre-visitor probe lost ~9-10% on degeneracies. |
| 20-23 | Implemented in `ff60d958d` and value-category fix `13fe33ba0`: stored sign, negation, exact miss and sign-only callers. Not landed. | CPU `alg-v1` compares side variants on that common base; local `exact_first_stage_filters.md` checks signs and two sign-only sites. No clean independent x86 signed-measure-only number identified. |
| 25-29 | Implemented: deferred/refined modes of `4b782a55b` retain decided signs and reverse-first rejection. | CPU `alg-v1`, prepared/deferred/reuse not consistently faster; targeted local avoided one adaptive call. No need to recreate the variant. |
| 31-37 | Landed: `7d8d1980f`, `9c24240ae`; native source admission, strict product order and product-tail ties. Non-native/uncertified paths unchanged. | CPU original projected ratios .926 GCC/.922 Clang vs unoptimized exact; combined stage/strict-side .924/.909 in compact projected. Translated control first A/B ~1. |
| 39-43 | Implemented: known common support `333cfeb22` includes off-line probe and direction/order checks, experiment only. | CPU `next-v1` common1, mostly within 1%. No broad support search warranted. |
| 45-60 | `6f820e54a`, `8abe94acc`: stored vs axis/absolute-products/displacement, independently certified N=100 layouts; native Shewchuk stage A and error-based paths. Only expression/displacement reached public side/dot before this audit. `a04487baf` tests native numerical classification. | Local first-stage report has proof, zero/exponent guards, box setup and construction-only timing; CPU `relative-v1`/ `tight-v1` cover public relative/displacement/compensated variations, **not** an operation box. |
| 61-65 | First-stage corpus and separate counters/timers exist. Orders and squared distances measured as primitives; no separate full equality, structural on-support or along-segment composition timings. | Partial coverage, not a mandate to retime every operand mix. Decimal work in CPU policy matrices uses native binary64 decimal-looking inputs, not the authoritative decimal wrapper. `d7577e8b1`/baseline `cbf0ce0ba` checks wrapper-native certificates vs rational preparation separately. |
| 67-73 | Certified factor signs (including useful zero) exercised on side/dot in filter/error-dispatch trials. | Computing signs independently was tried. Cross-call reuse of signs already computed for a different predicate is not established as a general dispatch. Needs a concrete producer/consumer pair. |
| 74-77 | Sterbenz and zero-error dispatch trials compare rounded products; nonexact operand construction is not admitted by an N-step guarantee alone. | Landed native orientation only, broader zero-error side/dot mixed. |
| 78-81 | Final-rounding-free strict monotone form covered by scalar/order and expression filters. | A new retry that manufactures separate product certificates is not justified. Reusing already-available separate budgets in other operations is a distinct, unmeasured minor variant. |
| 82-97 | Expression magnitude and axis alternative proved and timed; reference stage A vs dot/stored, misses and fallback observations in `exact_first_stage_filters.md`; directional boxes and axis scale can become pessimistic. | Local dispatch and x86 `relative-v1` are not universal wins. In the local primitive study, translating to 2^30 took expression/axis/box side from passes to zero while displacement and stored still passed. |
| 99-108 | **Actual public gap before this audit:** global box and cascade were only in `exact_bounded_predicate_filters.cpp`; the report explicitly withheld public producer provenance. | This branch tests certified box-then-local vs local alone and stored bounds, on the same eager N=100 layout. Only side and direction are wired; other candidate predicates remain primitive-only or unchanged. |
| 110-115 | Implemented independently and combined in `ee4a211d2`: terminal relate and source-edge lookup memoization. | CPU `terminal-v1`; memoization often costs 3-4%, bypass sometimes ~1%. No unconditional combined enablement. |
| 117-125 | Prepare directed side operands and reuse last certified shared endpoint: `4b782a55b`. | Local counters found 320/16,668 reuse in projected union and 1,994/9,472 shared; neither needed adaptive refinement. CPU `alg-v1` mixed. |
| 126-129 | **Open:** actual reuse of the prepared native determinant heads for the crossing coordinate construction; prepared look-ahead/angular rays not found. | Existing experiment prepared side/rejection, then initialized the point separately from its original producers. A new fused path must match determinant operand order/rounding and leave the support-canonical fallbacks intact. Prefer evidence of shared cost before adding state; do not redo the rejected last-endpoint cache. |
| 131-135 | Partial: `ee4a211d2` rectangles, `451182686` / `97a96448c` bounded atom XOR with ownership and independent oracle. Landed `6cc4b014f` bypasses final merge if original boundaries cannot share an edge, still two differences. | CPU `xor-ab` independent append; bounded atom validated 1,296 sheared and 128 triangle cases but slower on concave (1.40x time) and slanted-coincident (2.03x). **General** single noding pass with arbitrary coincidence/ownership is not implemented/proved. |
| 137-141 | Implemented `27ce26e94` and distinct hook `e65316026`: x intervals, exact overlapping groups. | CPU `cluster-ab` five pairs, ~.996 GCC/.991 Clang; separate borrowed-key confirmation ~1.008x speed. Keep experimental. |
| 143-146 | Terminal support vector/output projection in `ee4a211d2`; `3687a7bd7` landed disjoint-owning-root streaming. | CPU `terminal-v1`; no blanket omission of supports for escaping/intermediate output. Untouched-source handling exists within the terminal experiment, not a global ring guarantee. |
| 148-154 | Implemented in `3566bd692` and documented in `3aac84aa0f`: retained endpoint descriptors, filters and checked capacities, no auxiliary entries. | Local dedicated 4,096 query entries to zero; CPU `alg-v1` mixed (decimal Clang .980/GCC 1.023). No need to reimplement midpoint/splice. |
| 156-160 | Landed direct order of native/exactly-convertible bounds including infinity (`7721b58d3`, `9ed1a3f04`, `16ebaf04f`); zero-error constructed coordinates in `6987622b4`/filter experiments. | CPU `section-v2` projected 4.7-5.5% vs predecessor; `section-v1` other initialization neutral/negative. Bad scalar-probe results do not disprove direct admitted native order. |
| 162-172 | `a04487baf` / `3b4225480` classify zero-error side/dot independently of topology; scalar probe, exact-class and unconditional filter comparisons exist. | Local zero-error III side 9.6 to 5.4 ns, dot 10.3 to 4.8 ns; mixed slower. **Open:** general compact shared X/Y decoding with no owning snapshots everywhere, not justified as a blanket replacement. |
| 174-181 | Landed narrow certified stored originals `ad501de5d`: role, monotone chains, chord sides and boundary proof. | Local tangent-hole full union 6.77x/14.83x (64/128 edges); no similar result for generated walks. |
| 182-195 | Ring-reduction scratch reused in midpoint experiment; no proof of generated simple cycles/global edge ownership. 128,084 observed walks don't prove it; visit IDs/biconnectivity alone admit counterexamples. | **Open with proof blocker:** remove generated cancellation or replace contact midpoints with local incidence. No change to ownership, inter-ring coincidence or contact testing recommended from the current evidence. |
| 197-199 | First-miss reserve 4/16 `0a74f7528`; precision-copy avoidance in `1ac84aaed`. Operation ownership and no initial-hit allocation preserved. | CPU `cache-ab`: 4 slots helps eager general 2.4-3.0%, approximate output mixed by compiler; no default tuning. |
| 201-203 | Landed monotone scalar, native exponent/strict-side cleanup `09060f3be` and disjoint result `6ea63648b`; structural zeros/construction and predicate analyses in the filter experiments. | CPU strict/next/profile-final. Not whole-op machinery fusion; see 126-129. |
| 205-208 | Active limbs `1ac84aaed` and type-width-proved unchecked product `b07317a40` are separate experiments. | Local arm64 unchecked: wide orientation ~5% faster with Clang; stretched public ops 1.6-3.6% faster, GCC compact orientation +13.6%, one-ULP symdiff +2.6% time. **No x86 unchecked A/B identified.** Combining with active limbs or restricting by width is a variation, not an untouched TODO. |
| 210-218 | Landed distinct bounded/nearest scheduling and more selective N<=4 retries (`9dc72f815`, `0b30b3206`, `96965a0fc`). | CPU `measured3-*`, final `measured4-*`, `one-ulp-ab`, `bounded-ab`: separate exact topology and emission policies. No blanket better correct-rounding path; large-N tradeoff distinct. |
| 220-244 | Division-free bound, L admission, axis-selective passes, emission/refinement behavior are implemented. | Bounded construction tests exercise per-axis limits against nearest rational reference. The code's explicit exponent contract still matters. |
| 246-261 | Quotient residual implemented; same-head EFT budgets with effect-based selective retry in `0b30b3206`/ `96965a0fc`; larger-N sweep `533c09081`. | CPU N=1 and 4 projected-local ~.814/.801 and .773/.766 GCC/Clang, N=100 .963/.947, other N=100 cases near parity. No unconditional tails. |
| 262-278 | Compensated/Kahan replacement, fresh quotient/candidates, rounded factor bound / zero-tail gate were tried in `3615f8d11`. | CPU `relative-v1` / `tight-v1`; some projected improvements, mixed decimal. **Open** as written: evaluate/bound the nonzero TwoDiff cross terms to certify original determinant beyond those guards. Economy not demonstrated. |
| 279-282 | Current retry picks determinant by weighted error, skips the other if optimistic improvement cannot rescue an axis. | **Open:** salvaging a failed initial positive L by separately certifying D and rebuilding L and thresholds, rather than proceeding to construction fallback. Neither the retained same-head retry nor a magnitude-only trial covers all such cases. |
| 284-296 | `4075deefa` full-range three intersections, signed radix-2^21 exponent blocks with capacities, `357434d0a` real on-demand dispatch, no record-cache growth. | Local rational/edge and time/footprint evidence; near-parallel overlay 19.42 to 13.72 us, 16/12,348 ordinary orientation calls used it. No new x86 dispatch confirmation found. Broad EFT residual/nonoverlap backend for other expressions is not established; fixed ints remain. |

## Enrichment and measurement backlog

| TODO lines | Audit |
| --- | --- |
| 298-301 | Scope is valid Cartesian polygons, not generic graph/BG buffer claims. |
| 302-313 | `a5f1cd5d4` complete experiment for cluster position ordinals (including blockers), successor copy state/full round-trip behavior, reference helpers. CPU `alg-v1` bundled enrichment; don't attribute bundle timing solely to ordinals. |
| 315-322 | `c3daf553e` antisymmetry and constrained exact spike elision. Kept legacy spike/provenance handling; measured in cluster bundle, not permission to remove all repairs. |
| 324-329 | `0ca7b2760` cyclic reverse scan, difference events, half-open intervals; implemented and exercised in CPU rank-study/algorithm cluster matrix. This is independent of numerical precision. |
| 331-339 | `36878aecb` guards actual arrival direction and only selects a unique eligible ray; retains equal-ray/ownership/reversed/round-trip logic. CPU rank evidence: basin difference 2,047/2,047 selections, union 0/2,559 (ties); paired benefit ~.999x. Broader elision lacks proof. |
| 341-348 | Stale wording. CPU refreshed sampled profiles (`profile-final`, `profile-xor`) and paired GCC 13/Clang 22.1.8 x86-64-v2, also public B2 qualification. These shares aren't per-stage wall time. This audit runs focused local tests only. |
| 349-352 | Already removed volatile, reverted blanket noinline `f66ea07fa`; not open work. |

The superseding x86 `measured4` CSVs (policy4 header snapshot) show 45 matched
independent cells per policy. The saved summary script gives GCC/Clang ratios
to same-compiler develop of .871/1.152 (stored), 1.084/1.449 (N=1),
.991/1.320 (N=4), .909/1.192 (N=100), 1.248/1.667 (eager nearest), and
1.328/1.789 (deferred nearest). Exact fingerprints were unchanged across
policies. Simple disagreed on three symmetric-difference groups; develop
disagreed on projected-local symmetric difference. These ratios are not the
effect of any new filter in this audit nor equivalent-work results for simple.

## New experiment: eager side and direction box

**Why this gap:** the primitive results already supplied a useful proof
and a possible cheap stage, but did not answer whether producer admission,
operation setup and fallback pay for it in public calls. This trial answers
that before considering invasive prepared-turn fusion or general XOR.

Experimental selector `BOOST_GEOMETRY_EXACT_EAGER_PREDICATE_BOX`: 0 (default)
existing stored-error stage; 1 local maxima then existing stage; 2
operation-wide threshold, then the same local maxima and existing stage.
The variants preserve the eager N-step representation, construction,
refinement, support routing and fixed fallback. Compile all translation
units with the same mode. Only binary64 areal native source pairs and
eager bounded-ULP policies with N<=100 are admitted. A side/direction
query must contain at least one certified proper-crossing record, and
all indirect operands must carry the certificate from the same context.
Native operands must compare within the original maxima. Native-only
orientation, clipping, queries and other producer stages stay unchanged.

The operation scan uses the original immutable polygons/rings and matches
source identities (including boundary/reverse views). Only the proper-crossing
constructor stamps a record: all four side signs are nonzero and each pair
straddles, so the true intersection lies within both closed segments and
within the original source box. The stamp survives retention/finalization;
neither range IDs after finalization nor a common Root type prove ownership.
Source rejection, not a more liberal inference, handles mixed contexts,
unavailable producers or a native operand outside the box. The extra bit
uses existing record padding on the tested arm64 ABI (72 bytes).

### Error proof

Strict binary64 RN/even, gradual underflow and no contraction. Write
`u=2^-53` and `M>=1` for one axis maximum. It encloses either
the true coordinates (operation mode) or the current binary64 coordinates
(local mode). Source maxima are at most 2^64; local maxima below 2^65.
If the true point is bounded, rounding it to nearest remains within the
representable endpoints +/-M. From either the nearest value or the current
view, each step while within +/-2M has size at most 4u*M, including in
the subnormal band and across a binade. N<=100 steps cannot leave +/-2M.
The true point is within 2u*M of its nearest value on this band. Thus
both operands are bounded by 2M, and their error is at most
`(4N+4)u*M`. Native operands have zero error. This absolute envelope
does not assume a constructed approximate zero is exact.

For a rounded difference, the two coordinate errors plus its rounding
are bounded by `D*u*M`, with `D=8N+16`: the unrounded approximate
difference has magnitude <=4M, and RN error is <=4u*M + 2^-1075 <=5u*M.
True differences are at most 4M and computed differences less than 5M.
For a product with axis scales A,B, the two difference perturbations
cost at most `(4+5)D*u*A*B` (expand by replacing one operand at a time);
rounding the product costs at most `26u*A*B`. The final sum/difference
rounding costs at most `26u*S`: each rounded product is <25*A*B because
its factors are bounded by `(4+D*u)*M`, and the subnormal absolute
allowance is <u*S. Here S is `2*M_x*M_y` for side or
`M_x^2+M_y^2` for direction. Total error is at most
`(9D+52)u*S = (72N+196)u*S`.

The nonnegative threshold evaluation has at most three rounded operations
per term, with no overflow or underflow (inputs to the threshold >=1
and <2^65). `(80N+240)(1-u)^3 > 72N+196` for 1<=N<=100, so the
implemented rounded threshold encloses the total. A strict comparison
certifies a nonzero sign; equality takes the existing path. The proof
does not require actual preliminary floating differences or products to
be normal.

### Reproduce the local A/B

```sh
boost=${BOOST_ROOT:-/Users/tinko/dev/boost}
out=$(mktemp -d)
for mode in 0 1 2; do
    clang++ -std=c++17 -O2 -DNDEBUG -ffp-contract=off -pthread \
        -Wno-constant-conversion -Iinclude -isystem "$boost" \
        -DBG_BENCH_EXACT -DBG_BENCH_ULPS=100 \
        -DBOOST_GEOMETRY_EXACT_EAGER_PREDICATE_BOX=$mode \
        test/strategies/simple_cartesian_benchmark.cpp -o "$out/bench-$mode"
done
"$out/bench-2" --workload projected --operation sym_difference --seconds .12 --repeat 1
```

Use the same build with `-DBOOST_GEOMETRY_EXACT_FILTER_STATISTICS` and
`--check 1` to print the side/dot queries, eligible calls, box passes and
local passes to stderr. Do not time that build. Run the oracle
`exact_eager_predicate_box` (mode 2 by default; use the selector for 1).

### Local evidence and decision

Apple Clang 21, macOS arm64, C++17, flags above. Three cases per group;
seven rotated rounds (120 ms per timed cell), comparing per-call CPU
time. Modes 1 and 2 both perform the operation scan; 2 vs 1 isolates
using the global threshold after shared setup. Against 0 the scan and
admission costs are included. No separate timing for setup here (the
earlier primitive box scan was ~5.4 us for 3,072 vertices).

Changes are median per-round percentage changes, negative is faster.
Each triple is intersection / symmetric difference / areal relation.

| Group | Local vs stored | Box+local vs stored | Box+local vs local |
| --- | ---: | ---: | ---: |
| Degenerate | -0.3% / +0.6% / -0.6% | +0.6% / -0.3% / -0.4% | +0.7% / -0.6% / +1.0% |
| General | -1.0% / -1.2% / -0.6% | -2.1% / -1.6% / +1.4% | -1.2% / -2.2% / +1.6% |
| Decimal | -1.8% / -1.1% / +0.3% | -1.6% / -0.6% / +2.1% | +0.3% / +0.5% / +2.5% |
| Projected | -0.8% / -0.2% / +0.9% | +0.3% / -1.1% / +0.8% | +1.4% / -0.6% / +0.1% |
| Projected-local | -0.4% / -1.2% / -0.1% | -0.9% / -1.0% / +0.2% | -0.1% / -0.1% / -0.2% |

For scale, projected symmetric difference medians were 66.65 / 66.50 /
65.86 us (stored / local / box+local). Box+local vs local ranged from
-3.0% to +3.6% per pair. Separate counters for the three-case corpora:

| Operation | Side eligible | Box hits | Subsequent local hits |
| --- | ---: | ---: | ---: |
| General and decimal intersection (each) | 19 | 19 | 0 |
| General and decimal symmetric difference (each) | 42 | 42 | 0 |
| Projected intersection | 111 | 110 | 0 |
| Projected symmetric difference | 231 | 229 | 0 |
| Projected-local symmetric difference | 231 | 231 | 0 |
| General and decimal areal relate (each) | 2 | 0 | 0 |
| Projected and projected-local areal relate (each) | 1 | 0 | 0 |

No eligible public dot calls in that corpus. All 45 check-only operation
cells (counts and quantized geometry/DE-9IM digests) and all timed digests
matched across the three modes. This is not an exact-coordinate oracle.
The separate rational test verifies threshold enclosure and every
accepted sign over 512 triples per budget at 1, 4 and 100 steps, including
ties, subnormals, outward drift at 2^64, varying scales, and cancellation.
It tests original-source retention and finalization, same-context
admission, rejection of unsupported or out-of-box operands, and analytic
rectangle overlay/relate. A disjoint distant component defeats the
global threshold while the local stage still certifies signs.
Both modes passed on Apple Clang 21, and mode 2 passed on local
GCC 16.2, Clang ASan/UBSan, and the registered CMake/CTest target.
The existing borrowed predicate, borrowed turn and construction
rounding tests also passed with mode 2 enabled.

**Recommendation:** no default change. Closing the public box gap
found plenty of side passes but no consistent operation-wide benefit.
The results here are local arm64, not evidence of an x86 win. Do not
repeat root-wide box variants absent a workload with more certified
predicate traffic. The next credible native stored-output candidate is
sharing determinant heads between rejection and crossing construction
where the support order is known; measure that overlap first. General
XOR and generated-ring simplification still need ownership proofs, not
just new timings.

### Server follow-up

On x86-64-v2, GCC 13 and Clang 22, five rotated repeats over three public
operations and five workloads gave the following geometric means of
median time ratios (smaller is faster). All output fingerprints matched.

| Compiler | Local / stored | Box + local / stored | Box + local / local |
|---|---:|---:|---:|
| GCC | .9977 | 1.0007 | 1.0030 |
| Clang | 1.0101 | 1.0043 | .9943 |

This confirms the decision not to enable the extra box stage by default.
The measured5 matrix after the borrowed-output copy and orientation changes
had 45 operation/workload pairs and three repetitions per policy.
Stored-output exact took .8610 / 1.1333 times same-compiler develop
(GCC / Clang). Absolute Clang/GCC exact time was .9887; develop was .7511.
