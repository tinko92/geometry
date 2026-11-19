# Exact Cartesian Issue Audit

Snapshot: 2026-11-11, boostorg/geometry, base commit
`72b2d832fb9596a2a23c0c164a5426ae720d0cd6` plus concurrent worktree changes.
The inventory contains all 144 open issues at capture time. This document
accounts for the 48 Cartesian relate/set-operation candidates, with nine
numeric-foundation items listed separately. Open status is a snapshot, not a
claim about future tracker state. **No issue is declared fixed.**

## Current Replay Status

- 47/48 candidates have selected-case execution evidence: 46 pass the stated
  exact contracts; #1448 is a guarded dimensional advisory under the preserved
  legacy box-intersection contract, not a rounding artifact.
- The current additional runner executes **52 typed pairs + 7 pipelines**:
  59 records: 58 exact-contract passes plus the #1448 advisory in the final
  compatible-clipping integration run.
- Fourteen API groups cover eleven geometry issues and three numeric-foundation
  examples. The unchanged original 16-case unit also passes, including
  representation-only control #494 outside the 48.
- Missing data: only #653 numeric coordinates. The #1044 ZIP and #1253/#669
  JSON sources were obtained, inspected without executing their programs,
  embedded as fixtures, and replayed.
- No stock/default baseline comparison has been established. Passing selected
  contracts is not proof that every original report variant is resolved.

Latest detailed results are `/tmp/bg-exact-additional-all-results.jsonl`
(59 records), not the older 45-record
`/tmp/bg-exact-additional-results.jsonl`. Compact machine-readable accounting:
`/tmp/bg-exact-additional-replay.json`,
`/tmp/bg-exact-issue-catalog.json`, and
`/tmp/bg-exact-issue-coverage.json`. Full retrieved issue bodies/comments:
`/tmp/bg-exact-in-scope-issue-data.json`.

## Contract And Coverage

Native float/double/integer values are parsed **before** exact lifting. Their
binary/integer values, not intended decimal collinearity, define the geometry.
Only orientation/closure is normalized; no snapping or simplification is used.
Full pair below means the available pair's complete input coordinates, **not**
the entire original program, all model variants, or a full independent output
oracle. The typed runner singleton-wraps areal/linear operands into multi
geometries where appropriate; separate API groups preserve their source API
shapes. Selected/reduced/adapted variants are labelled explicitly.

Pipeline rows replay all listed inputs and operations with exact retained
intermediates, but remain **partial report replays** where the original program
feeds native rounded output into later operations. Most checks are exact
validity, inclusion, reversal/order invariance, and area/algebraic identities;
they share the overlay engine and are not independent complete-set goldens.
The #528 segment oracle and #440/#858 rational formulas are separate arithmetic
references. The DE-9IM/linear/tuple expectations are stated per row.

Exact output validity and native output representation are separate contracts.
In particular #1034, #602, #1397 and #663 have valid exact output but invalid
nearest-even native-coordinate output in these replays. Conversely #1448's
singleton line is **already invalid with exact coordinates**, not rounding.
Existing box-intersection regressions explicitly require this legacy contact
cardinality, so it is preserved in this exact-arithmetic refactor. A dimensional
API change needs a separate targeted fix. Its empty-linear-output assertion
remains guarded by `BOOST_GEOMETRY_TEST_FAILURES`, not an active new contract.
Two exact-distinct clipped endpoints that round to the same integer or double
output point must still retain both output slots; the API unit checks both.

Suite keys: **I** = [exact_issues.cpp](exact_issues.cpp);
**A** = [exact_additional_issues.cpp](exact_additional_issues.cpp);
**P** = its full-input exact pipeline branch;
**API** = [exact_additional_issue_api.cpp](exact_additional_issue_api.cpp).
Input fixtures: [original](exact_issue_cases.hpp) and
[additional](exact_additional_issue_cases.hpp).

## All 48 Candidates

| Issue | Replay Scope | Exact Observation | Limits / Suite |
| --- | --- | --- | --- |
| [#1506](https://github.com/boostorg/geometry/issues/1506) | Full pair | Union: one valid component in both orders; exact input-area sum. | I; no stock baseline. |
| [#1490](https://github.com/boostorg/geometry/issues/1490) | Both local pairs | Original and grid differences valid, one component each. | I/A; native integer output is a separate representation. |
| [#1487](https://github.com/boostorg/geometry/issues/1487) | Linked/local pair | Intersection: one valid component; reported exact area about 18983.884971862466. | I; PR #1508 is not assumed present or unnecessary. |
| [#1448](https://github.com/boostorg/geometry/issues/1448) | Guarded dimensional advisory | Exact line/box contact retains the legacy one-point linestring. | A/API; exact singleton is invalid as a line, NOT rounding. Existing API cardinality is preserved; dimensional change is separate. |
| [#1421](https://github.com/boostorg/geometry/issues/1421) | Full pair, both operations | Difference: six components; symmetric difference: one; valid. | I/A; exact algebraic/area checks, not a stock comparison. |
| [#1420](https://github.com/boostorg/geometry/issues/1420) | Full pair | Union: four valid components. | I; small-integer traversal control, not an FP-only diagnosis. |
| [#1419](https://github.com/boostorg/geometry/issues/1419) | Full-input exact pipeline | Three operands, two differences; final one valid component. | P; exact intermediates. Original repeated native-output feedback not replayed. |
| [#1414](https://github.com/boostorg/geometry/issues/1414) | Full large pair | Union: one valid component; both orders and area partition pass. | A; no independently supplied full result WKT. |
| [#1406](https://github.com/boostorg/geometry/issues/1406) | Selected normalized variants | Original and common-vertex variants return one contact point, no line. | A; source polygon orientation corrected; not an unchanged source-model replay. |
| [#1397](https://github.com/boostorg/geometry/issues/1397) | Full integer pair | Intersection: three exact-valid components. | A; native integer rounding self-intersects. |
| [#1395](https://github.com/boostorg/geometry/issues/1395) | Selected strict-FP variant | Normalized pair produces valid exact intersection. | I; source custom CCW/open model and /fp:fast configuration not reproduced; fast-math unsupported. |
| [#1394](https://github.com/boostorg/geometry/issues/1394) | Full integer pair | Union: three valid components, operand-order/area checks pass. | A; no pure-FP cause claimed. |
| [#1363](https://github.com/boostorg/geometry/issues/1363) | Selected local regression pair | Both intersections empty; input matrix FF2FF1212. | A; local issue_1363 data used; external Godbolt program not independently retrieved. |
| [#1354](https://github.com/boostorg/geometry/issues/1354) | Full large pair | Difference: one valid component, one hole. | A; no legacy validity suppression. |
| [#1290](https://github.com/boostorg/geometry/issues/1290) | Full-input exact pipeline | Difference valid, subset of minuend, no subtractor interior overlap. | P; overlaps == !touches is not a general identity; native-output feedback differs. |
| [#1253](https://github.com/boostorg/geometry/issues/1253) | Full downloaded inputs/API | Independent rational boundary oracle finds zero intersections; exact linear output empty. | A/API; both line endpoints outside triangle x extent. Do not impose the source's intended one-segment answer; this is not #1448's point contact. |
| [#1241](https://github.com/boostorg/geometry/issues/1241) | Full pair | Difference: one valid component. | I; near-collinear coordinates retained. |
| [#1238](https://github.com/boostorg/geometry/issues/1238) | Full pair/predicate | within is false in the report's call order. | I; source operand naming/order preserved. |
| [#1236](https://github.com/boostorg/geometry/issues/1236) | Full-input exact pipeline | Three intermediates valid; final within is true, two components. | P; retained exact branch, not repeated native-output feedback. |
| [#1235](https://github.com/boostorg/geometry/issues/1235) | Full inputs/API | Tuple: zero isolated points, zero lines, two valid areal components. | A/API; areal set matches areal-only result. |
| [#1232](https://github.com/boostorg/geometry/issues/1232) | Full inputs/API, partial oracle | Difference retains (5,-8) and (7,4); reversal/validity pass. | A/API; no independent complete line-difference golden; failure-gated legacy golden not imposed. |
| [#1212](https://github.com/boostorg/geometry/issues/1212) | Adapted full pair + reduced pair | Full integer data: one valid component, two holes. | I/A; int64 input adaptation; original mpq_rational API not tested. |
| [#1211](https://github.com/boostorg/geometry/issues/1211) | Full inputs, normalized model | Native float intersection: one valid component. | I; CCW/open source model normalized; exact float bits retained. |
| [#1201](https://github.com/boostorg/geometry/issues/1201) | Selected corrected comment pair | Tiny nonempty intersection, one valid component. | I; intended empty answer conflicts with supplied binary input. |
| [#1179](https://github.com/boostorg/geometry/issues/1179) | All four supplied pairs | Each union valid and nonempty, order/area checks pass. | A; no independent full-output WKT oracle. |
| [#1178](https://github.com/boostorg/geometry/issues/1178) | Both supplied variants | Each difference: two valid components. | I/A; island-enclosing graph controls. |
| [#1174](https://github.com/boostorg/geometry/issues/1174) | Full-input quadrant pipeline | Four nonempty valid pieces; disjoint interiors and exact area sum. | P; full close-vertex oval retained; native rounded piece output not separately asserted. |
| [#1053](https://github.com/boostorg/geometry/issues/1053) | Selected self-relation variants | Full multipolygon within itself; matrix 2FFF1FFF2. Polygon covered_by self passes. | I/A; potentially invalid 0.075 comment variant not replayed. |
| [#1044](https://github.com/boostorg/geometry/issues/1044) | Full ZIP inputs/API | bug/control arrays, each original and native-y-shifted, pass both orders with one valid component. | A/API; subset checks for both operands plus independent rational upper-chain area and bounding box. Decimal 0.09 is not imposed exactly. |
| [#1035](https://github.com/boostorg/geometry/issues/1035) | Full pair | Difference: three valid exact components. | A; normalized integer model; no full native-output pipeline claim. |
| [#1034](https://github.com/boostorg/geometry/issues/1034) | Full existing pair | Difference: four valid exact components. | A; ticket_12752 reused without spike cleanup/validity suppression; rounded integer output invalid. |
| [#900](https://github.com/boostorg/geometry/issues/900) | Full-input exact pipeline | Symmetric-difference identity passes; final four valid components. | P; original four stars, not local small-grid issue_900; native feedback branch unrun. |
| [#888](https://github.com/boostorg/geometry/issues/888) | Full-input exact pipeline + three reduced pairs | All 39 inputs and 38 intermediates valid; final 80 components. | A/P; reduced results 22/15/117 components. No independent complete final-set golden or native-feedback replay. |
| [#858](https://github.com/boostorg/geometry/issues/858) | Full integer sweep/API | All 16 clips valid; shared endpoint equals (9-2*y/9,y) exactly. | A/API; endpoint checked before integral output rounding. |
| [#852](https://github.com/boostorg/geometry/issues/852) | Full pair | Intersection: one valid exact component. | I; exact working result is not a guarantee that every rounded serialization is valid. |
| [#831](https://github.com/boostorg/geometry/issues/831) | Full native-float API | crosses true both directions; matrix 101F00212 matches rational reference. | A/API; float parsed before lifting. |
| [#704](https://github.com/boostorg/geometry/issues/704) | Selected corrected comment pair | Tiny valid overlap; input matrix 212101212. | A; intended decimal shared-edge/empty result is not the binary64 contract. |
| [#669](https://github.com/boostorg/geometry/issues/669) | Full downloaded pair + crop/rotation variants | Full result: one component/one hole; three crop branches: one component/no holes. Eight start-point variants pass exact area partition and inclusions. | A/API; genuine binary-input overlap about 3.85e-12, so unchanged area is not an exact oracle. Exact, rounded-exact and source-default crop feedback all replayed. |
| [#663](https://github.com/boostorg/geometry/issues/663) | Full pair | Tiny nonzero exact-valid intersection. | A; native double rounding collapses dimension; do not impose intended empty answer. |
| [#656](https://github.com/boostorg/geometry/issues/656) | Full inputs/API | Exact linear set is MULTILINESTRING((0 0,1 0,1 1)). | A/API; reversed input agrees; isolated contact does not retain outside line. |
| [#653](https://github.com/boostorg/geometry/issues/653) | Unavailable coordinates | No replay. | Reviewed screenshots/body/comments lack copyable numeric data; no synthetic substitute. |
| [#631](https://github.com/boostorg/geometry/issues/631) | Full INT_MIN pair | Symmetric difference: one valid component. | I; rectangle shape is the supplied reference; subtraction/promotion overflow control. |
| [#602](https://github.com/boostorg/geometry/issues/602) | Full existing pair | Difference: four valid exact components. | A; ticket_11674 reused without validity suppression; rounded integer output invalid. |
| [#588](https://github.com/boostorg/geometry/issues/588) | Exact pipeline + selected frozen output | Retained exact difference within red; body-WKT blue also within under exact strategy. | P; later higher-precision printed triangle variant not replayed; do not generalize to all rounded outputs. |
| [#528](https://github.com/boostorg/geometry/issues/528) | Full bit-pattern API | Independent rational 2x2 oracle gives zero boundary intersections; output agrees. | API; original uint64 representations decoded by memcpy, not rounded decimal comments. |
| [#511](https://github.com/boostorg/geometry/issues/511) | Full pair | One valid clipped line; reversal and subset checks pass. | A; source polygon orientation normalized; no independent full-output golden. |
| [#510](https://github.com/boostorg/geometry/issues/510) | Full pair | One valid tiny clipped line; reversal and subset checks pass. | A; source polygon orientation normalized. |
| [#440](https://github.com/boostorg/geometry/issues/440) | Full integer-segment API | x=4700000; y=-2665034-1351*10043/24634 exactly. | A/API; all endpoint-direction and operand-order variants pass. |

## Numeric Foundation: Nine Related Items

These are separate from the 48 geometry candidates. Three literal example
groups were replayed under contraction fast and off; the other six are
design/policy items, not supplied standalone literal fixtures. This does not
claim to reproduce every original compiler or the linked CI pipelines.

| Issue | Relationship | Audit Status |
| --- | --- | --- |
| [#1390](https://github.com/boostorg/geometry/issues/1390) | Investigate the need for a topological predicate | Predicate-design task; no standalone geometry fixture in the retrieved body. |
| [#1389](https://github.com/boostorg/geometry/issues/1389) | Investigate the reason for wrong ix / xi turns | Wrong ix/xi-turn investigation; linked design discussion, no standalone literal regression replay. |
| [#1388](https://github.com/boostorg/geometry/issues/1388) | Add predicate involving get_distance_measure | get_distance_measure predicate design; no issue-specific geometry replay. |
| [#1387](https://github.com/boostorg/geometry/issues/1387) | Make predicate to walk through consecutive nearly-equal points | Nearly-equal consecutive-point policy; no issue-specific geometry replay. |
| [#1380](https://github.com/boostorg/geometry/issues/1380) | Change behaviour of the define `BOOST_GEOMETRY_NO_MULTIPRECISION_INTEGER` | Integer-promotion/multiprecision macro policy; not a copied set-operation repro. |
| [#1148](https://github.com/boostorg/geometry/issues/1148) | Incorrect results of boost::geometry::math::smaller | Literal 1, 1-eps, 1-eps/2 for float/double: exact ordering, strict-order laws and all six minmax permutations pass. Legacy math::smaller and linked CircleCI failures are not changed/replayed. |
| [#1094](https://github.com/boostorg/geometry/issues/1094) | Try to make threshold in get_cluster adaptive | Adaptive cluster-threshold policy task; exact no-epsilon behavior is not proof of resolving the legacy policy. |
| [#918](https://github.com/boostorg/geometry/issues/918) | side_robust causes regressions when rescaling is turned off | All seven printed triples match the rational oracle and antisymmetry/cyclic checks. Signs: -1,-1,+1,+1,-1,-1,-1. Row 2's historical old sign +1 is not an exact binary oracle. Named difference-suite configuration not independently replayed here. |
| [#702](https://github.com/boostorg/geometry/issues/702) | Minor consistency issue with modern CPUs in side_by_triangle | Original six binary64 inputs give +1 and -1 under endpoint swap, matching rational arithmetic; cyclic invariance passes. No claim about changing legacy side_by_triangle on every GCC/FMA configuration. |

## Verification And Remaining Work

C++14 final integration build and CTest passed **9/9** (10.76 seconds), after
refreshing both issue units for the restored #1448 compatibility advisory.
Source manifests and dependency inputs remained stable during the final run.
This normal CTest run does not enable BOOST_GEOMETRY_TEST_FAILURES and does
not establish that #1448 is fixed. Full integration evidence:
`/tmp/bg-exact-final-post-box/README.md` and its `LastTest.log`.
Jam registration is updated but Jam execution was not performed here.

PropertyTree is a private dependency of exact_additional_issues in CMake and
Jam, and is included in the root's test-only Boost dependency list. It was not
added to the production geometry target.

All nine exact units are registered. The eight that use implicit point or ratio
caches have Jam `<threading>multi` and private CMake `Threads::Threads`;
arithmetic alone does not need it. Older Unix toolchains may require native
thread link flags even for single-threaded callers because caches use
`std::call_once`. No thread dependency was added to the production geometry
target. Main reports `exact_segment` (13,622 pairs) and `exact_gc` passing;
their results are separate from this issue audit's execution evidence.

Next: obtain #653's numeric coordinates; replay remaining source-model and
native-feedback variants; compare identical inputs with stock/default behavior.
Keep any #1448 dimensional/cardinality change separate from this refactor.
[PR #1508](https://github.com/boostorg/geometry/pull/1508) is an open PR, not an
issue in this list; its path-out-of-cluster change for #1487 is not assumed in
the PR1491 base. No tracker comments, closures, or fix claims were posted.
