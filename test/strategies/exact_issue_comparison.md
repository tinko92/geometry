# Exact Issue Symptom Comparison

**Attribution: frozen exact `bc65c6d60` already includes PR1491 algorithm fixes. The 20 selected-symptom improvements versus official develop belong to the combined PR1491-plus-exact branch, NOT 20 numerical fixes attributable to the exact strategy. 15 are incremental versus PR1491. The other five (#1490, #1421, #1420, #1238, #1178) already pass with PR1491's algorithm fixes.** Even the incremental comparison is a branch comparison, not a numerical-only causal isolation of the strategy.

## Optimized Regression

Optimized live exact headers were compared with frozen-before on **all 82 in-scope groups, native and retained: 164 records across 47 unique issues**. **Zero rational-WKT, predicate, or other semantic-field differences**; 0 records with separately listed public-validator diagnostic differences. No baseline/default rebuilds; the same frozen evaluator was reused. The full live header manifest was identical before and after seven serial optimized-only builds and all runs. Header manifest SHA256: `420afc12bcd54803badc5bfcb7beb55c1f8ddb884ba6d41f0a974f1e12767eab`. Detailed per-record comparisons, compiler commands, header hashes, and logs: `/tmp/bg-exact-issue-comparison-optimized`. The frozen symptom classifications/counts therefore remain unchanged under this optimized header snapshot.

Frozen comparison of all **47 executable unique issues**, plus unavailable #653. **82 in-scope case groups**, plus the #494 representation control, produce **332 new build/mode records** (83 groups times four modes: develop default-native, PR1491 default-native, frozen exact-native, frozen exact-retained). Each group may execute multiple direction/order calls or pipeline steps. These are not old unit-test pass totals. No whole issue is declared fixed.

## Counts

| Classification | Develop | PR1491 |
| --- | ---: | ---: |
| B: Reported Public Symptom Removed | 20 | 15 |
| A: Already Passing Reported Symptom | 13 | 19 |
| C: Retained Only; Not Publicly Solved | 5 | 5 |
| D: Binary Contract Conflicts With Intended Answer | 3 | 3 |
| R: Additional Reference Improvement; Not Counted | 4 | 3 |
| F: Legacy Dimensional Failure | 1 | 1 |
| P: Unsupported Predicate Premise | 1 | 1 |
| E: Unavailable | 1 | 1 |

**Strict selected reported-public-symptom improvements: 20 for combined PR1491-plus-exact versus develop; 15 incremental versus PR1491.** Denominator 47 executed issue selections, not 46 passing exact contracts. The table totals 48 only because unavailable #653 is included separately.

Criterion B: a quoted reported public cardinality/validity/predicate or gross outside-output symptom is observed in the selected native baseline call(s), and removed for all selected variants with exact strategy AND native final output. Not retained-only, intended-answer conflicts, additional reference-only improvements, or ULP-only differences.

Develop B: #1506, #1490, #1421, #1420, #1419, #1414, #1253, #1241, #1238, #1232, #1211, #1179, #1178, #1174, #1044, #1035, #858, #831, #528, #511.

PR1491 B: #1506, #1419, #1414, #1253, #1241, #1232, #1211, #1179, #1174, #1044, #1035, #858, #831, #528, #511.

R is intentionally excluded: a newly measured reference defect can improve without reproducing the report's quoted symptom. C is also excluded even if every retained-exact contract passes. D does not turn an intended decimal answer into a binary-input fix claim. All classifications concern these selections and this toolchain, not every source model/compiler/report variant.

## Provenance

- Develop: `0b4851a77d91633cf86c3a5e679e67d8c60200e9`.
- PR1491: `72b2d832fb9596a2a23c0c164a5426ae720d0cd6`.
- Exact: frozen before `bc65c6d60fe2d891ed978d0c8cdf91bb678712b9`, never live optimization headers.
- Every exported header was checked against its commit's Git blob; complete path/blob/SHA256 manifests are in `/tmp/bg-exact-issue-comparison-final/{develop,pr1491,before}-headers.tsv`.
- Flags: `-std=c++14 -O0 -g0 -ffp-contract=off -pthread -DBOOST_GEOMETRY_DISABLE_DEPRECATED_14_WARNING`; default public API with no strategy for either baseline; `COMPARISON_EXACT` uses `cartesian_exact<NativeCoordinate>` on frozen before. No fast-math, rescaling, test-failure, or no-robustness macro.
- Compiler: Apple clang version 21.0.0 (clang-2100.4.1.1). Compilation serial, one compiler; Perl alarm bounds compile/runtime/evaluation. All 332 records executed and evaluated without exception, timeout, or process failure.
- Same native parse type, coordinates, orientation/closure correction on all builds. Initial rational-WKT input digests match across all four build/mode combinations; all original selected inputs are rational-valid. No snapping/simplification.

Data hashes and every per-issue issue-body/comment hash, baseline commit, exact commit, input digest, flags/limits, and raw record paths are in [exact_issue_comparison_results.json](exact_issue_comparison_results.json). Machine logs and full rational-coordinate output are in `/tmp/bg-exact-issue-comparison-final`; the manifest records full compiler identity and source SHA256 values.

## Reading The Evidence

Public outputs use the native coordinate type, including integral outputs. Retained outputs preserve exact constructed coordinates. Both run separately; pipelines genuinely feed rounded native output back in the native branch. Rational validation always uses frozen-before `cartesian<cpp_rational>`, with source values lifted from native bits (0.1 lifting is tested). Each build's default public validator is recorded separately. A default-validator rejection and a rationally invalid geometry are not conflated.

Pair tests preserve polygon versus multipolygon operand shapes, with native float/int parsing first; areal output uses multipolygon containers. Source-specific API groups cover line/areal tuple, segment/point, ring/line, box/line, native float crosses, bit-pattern points, downloaded multi/multi arrays, and CCW full/cropped difference. Source-model normalizations and remaining pipeline singleton wrapping are explicit below.

Independent rational segment solves/clipping support #528, #1253, #510, #511, #656 and #858; the downloaded #1044 upper-chain integral is independent of overlay. Other area/set/reference checks share BG machinery and are not independent complete-set goldens. A native rounded line normally differs from its rational clipping set: #510 is the control demonstrating why this alone is not a counted fix.

## Per Issue

A/B/C/D/E/F/P/R refer to the counts above. Every row has a separate outcome for each baseline, even when equal.

| Issue | Develop | PR1491 | Measured Outcome | Scope / Flags |
| --- | --- | --- | --- | --- |
| #1506 | B | B | Reverse union empty -> one component; both orders retain input-area sum. | full pair; both orders |
| #1490 | B | A | Develop: both invalid difference variants -> one valid native component / one hole each. PR1491 already passes both. | two local pairs; integer output; old grid area expectation 3900 is not a valid oracle |
| #1487 | A | A | Both baselines already return one component, area about 18983.885; no missing-intersection fix counted. | selected strict-FP pair; exact native default-validator regression in forward order |
| #1448 | F | F | One invalid singleton linestring in baseline, native exact and retained exact. | legacy dimensional/cardinality contract unchanged; not rounding |
| #1421 | B | A | Develop: difference 4 -> 6 components, area 16 -> 24; symmetric difference 2 -> 1, area 28 -> 36. PR1491 already passes. | both reported operations; small integer graph case |
| #1420 | B | A | Develop: forward union 3 components / area 36 -> 4 / 32; reverse also gives 4 / 32. PR1491 already passes. | small integer traversal; both orders |
| #1419 | B | B | Second native-feedback difference empty -> one valid component. | all three inputs; native feedback AND retained exact run; polygon operands singleton-wrapped in pipeline |
| #1414 | B | B | Forward union empty -> one valid component; both orders agree. | full large pair; no independent complete set golden |
| #1406 | R | A | Develop loses the contact point; exact returns (1,1.5), no line. PR1491 already does so. Neither baseline has the quoted spurious line after normalization. | two normalized variants; polygon/line two-slot tuple; original wrong orientation not endorsed |
| #1397 | C | C | Retained intersection valid; final integer output still self-intersects. | normalized CCW/open source model; native integer rounding still unsolved |
| #1395 | A | A | Baseline and exact native intersections valid under these strict-FP normalized calls. | original MSVC /fp:fast and custom CCW/open model not reproduced |
| #1394 | R | R | Reported invalidity absent, but baseline union drops a large component: area 463619.5 -> 1842454. | gross new reference improvement, NOT quoted invalidity fix; integer native area differs from retained |
| #1363 | R | R | Baseline returns bogus positive-area intersection in one order; retained/native exact empty. Quoted invalidity is not reproduced. | local regression pair; external Godbolt program not independently replayed |
| #1354 | A | A | Baseline already retains one component / one hole, same rounded result area. | full local pair; no inherited validity suppression |
| #1290 | P | P | overlaps == !touches is not a general identity; no fix inferred from that premise. | native feedback and exact retained branches run; corrected polygons replace source rings |
| #1253 | B | B | Two large outside lines -> empty in both directions, matching independent rational clipping. | downloaded binary64 inputs; intended corner segment is not the binary oracle |
| #1241 | B | B | Empty difference -> one valid component. | full pair; source polygon API shapes preserved |
| #1238 | B | A | Develop within true -> false; PR1491 already false in the report call order. | native multi-linestring predicate; exact false |
| #1236 | C | C | Native-feedback within remains false; retained exact within true. Rounded inputs really fail rational within. | all intermediates native-valid, so this is predicate persistence, not invalidity |
| #1235 | A | A | Baseline already has zero isolated points / lines and two valid polygons. | polygon/multipolygon tuple API; only rounded-coordinate differences |
| #1232 | B | B | Missing (7,4) tail and reversal mismatch removed; (5,-8) head also retained. | native line/multi-line difference; partial oracle, not a complete set golden |
| #1212 | A | A | Both reduced and full baseline results already retain two holes. | int64 adaptation; original mpq_rational API not reproduced |
| #1211 | B | B | Oversized intersection area 0.05117 -> about 6.30e-11, below the reported 0.01 bound. | native float parsed first; source CCW/open model normalized |
| #1201 | D | D | Gross million-area error disappears, but binary input has nonempty area about 1.25e-28, not the requested empty answer. | excluded from strict fix count; exact native default validator rejects tiny output |
| #1179 | B | B | All four pairs: intersecting-interior native union outputs -> valid native unions, in both orders. | all four supplied pairs; no independent complete output WKT |
| #1178 | B | A | Develop: both disconnected-interior variants -> two native-valid components. PR1491 already passes both. | all supplied variants; small integer topology |
| #1174 | B | B | Two empty quadrants -> four nonempty rational-valid rounded pieces. | single polygon/box calls; native default validator rejects the tiny-edge pieces; synthetic concatenation is NOT a valid multipolygon contract |
| #1053 | A | A | Both supplied self-coverage/within selections already true. | full multipolygon plus original polygon selection; 0.075 comment variant excluded |
| #1044 | B | B | Bug variant area 0.24 -> rational upper-chain area about 0.09; shifted/control variants pass. | multi/multi source API plus polygon pair; all four original/shifted bug/control variants; both orders |
| #1035 | B | B | Invalid integer difference -> native-valid difference. | normalized CCW/open model; actual integer output changed, not just valid retained output |
| #1034 | C | C | Retained difference has four valid components; rounded integer output still invalid. | full ticket_12752 pair; no spike cleanup |
| #900 | C | C | Native feedback changes wrong-area result but remains invalid (10 components); retained exact has four valid components and passes identity. | all four stars; native feedback != retained exact |
| #888 | A | A | All 38 baseline native intermediates already valid; final 80 components. | all 39 inputs / balanced reduction plus three reduced pairs; no independent final-set golden |
| #858 | B | B | Shared integer endpoints disagree in baseline sweep; all eight pairs agree with exact-native nearest-even output. | 16 source ring/line calls plus line reversal; native integer output policy is explicitly part of this fix |
| #852 | A | A | Baseline and exact native intersection already valid; ULP-level change is not a fix. | full pair; no blanket rounded-output guarantee |
| #831 | B | B | crosses true/false under line reversal -> true/true, both matrices 101F00212. | native float line/polygon API; rational matrix reference from audit |
| #704 | D | D | Binary overlap is real (matrix 212101212); empty/shared-decimal-edge expectation not imposed. | corrected comment pair; all native outputs already nonempty; no ULP-only fix count |
| #669 | A | A | No reported gross difference failure on baseline full/crop/rotation selections; tiny binary overlap is legitimate. | full downloaded CCW source API, rotations, native crop feedback, retained crop and default-crop feedback; not an unchanged-area oracle |
| #663 | D | D | Baseline empty; retained exact tiny nonempty overlap; native exact collapses to invalid zero-area polygon. | binary intended-empty conflict AND native representation regression; no fix |
| #656 | A | A | Baseline already equals the exact boundary chain, with no extra outside line. | native multi-line/multi-polygon API; independent rational clipping |
| #653 | E | E | No numeric coordinates available; no synthetic substitute. | not one of the 47 selected executions |
| #631 | A | A | Baseline already returns the supplied rectangle, area 214748364800. | INT_MIN case; no claim for all overflow configurations |
| #602 | C | C | Retained four-component difference valid; rounded integer difference still invalid. | full ticket_11674 pair; no validity suppression |
| #588 | A | A | Baseline native difference within red already true; frozen body-WKT within also true. | native feedback plus retained exact; later high-precision printed triangle not replayed |
| #528 | B | B | One spurious boundary intersection -> zero, matching independent rational segment solve. | original uint64 bit patterns; native polygon/line point-output API |
| #511 | B | B | Long left outside prefix removed: minimum result x about -21.856 -> -15.2803. | normalized source polygon-first API; retained output matches independent clipping; final endpoint rounding is not exact set equality |
| #510 | A | A | Baseline already returns the same native tiny single segment; exact-retained endpoints improve internally only. | normalized source polygon-first API; native rational set inequality due to rounding is not counted |
| #440 | R | R | Both baselines already direction-consistent, but wrong point y=-2690131; exact-native y=-2665585 is the correct rounded reference. | eight native integer segment directions/orders; new gross reference improvement, NOT the reported direction-asymmetry fix |

## Important Qualifications

- **#1487 is not counted against either baseline.** Both native original polygon calls already return a nonempty intersection of the expected scale in both orders here. Develop includes PR1508; no inference about its necessity follows from selected exact tests. Exact native forward output is rational-valid but frozen-default `is_valid` reports a spike. This new public-validator failure is recorded, not hidden as an issue fix.
- **#1034/#602/#1397 remain native-invalid** despite valid exact retained output. **#1236 remains false** with native feedback, even under exact predicates; rational validation confirms the rounded-set within failure. **#900 remains native-invalid**, despite an exact-retained identity pass. These five are C, not fixes.
- **#1174** corrects the quoted empty/wrong quadrants with four nonempty rational-valid native pieces. Frozen-default `is_valid` rejects the tiny-edge pieces; no default-validation fix is claimed. Concatenating adjacent quadrant polygons is only a synthetic area aggregate, not a valid output multipolygon contract, for either exact or native coordinates.
- **#858 integer rounding is a public behavior change**, not merely better internal endpoints: all eight adjacent clips now share an identical native integer endpoint. **#440** also changes native integer output to correctly rounded y=-2665585, but baseline directions already agree, so it is R, not a reported direction-symptom fix. **#1035** becomes native-valid; #1034/#602/#1397 do not. Integer rounding cannot be advertised as generally topology-preserving.
- **#1490 grid reference correction:** the frozen fixture's 3900 reference is wrong for its listed coordinates. Outer area=5550, hole area=150; subtractor area=1500 contains that hole, so the difference area is 5550-1500=**4050**. The old 3900 result subtracts the 150-area hole twice. This report leaves fixtures/audit untouched and counts only the observed invalid-result removal, supported by the corrected area calculation.
- **#1201** removes the huge wrong result but does not satisfy the report's intended empty answer; binary overlap is real and default native validation rejects the tiny result. **#704** already has tiny nonempty native overlap. **#663** changes baseline empty output into an invalid rounded tiny polygon. These three are D, not B.
- **#1394/#1363 and develop #1406** improve gross public reference behavior but not the quoted invalidity/spurious-line symptom in these normalized baseline runs. R keeps those improvements separate from the strict count.
- **#1448** stays an invalid singleton even before rounding; a dimensional contract fix is separate. **#653** has no reproducible numeric input. No fabricated data or whole-issue upper-bound count is used.

## Reproduce

From this worktree, with frozen headers at the stated paths:

```sh
node test/strategies/exact_issue_comparison_run.js /tmp/bg-exact-issue-comparison-final
node test/strategies/exact_issue_comparison_report.js /tmp/bg-exact-issue-comparison-final
```

For optimized-only verification after the frozen run exists, without rebuilding defaults or the evaluator:

```sh
node test/strategies/exact_issue_comparison_optimized.js /tmp/bg-exact-issue-comparison-optimized
node test/strategies/exact_issue_comparison_report.js /tmp/bg-exact-issue-comparison-final /tmp/bg-exact-issue-comparison-optimized
```

`COMPARISON_INCLUDE` can override the optimized header path; `COMPARISON_REFERENCE` can override the frozen reference directory. No production headers, existing cases, build registration, remote branches, or tracker state are modified by these tools.
