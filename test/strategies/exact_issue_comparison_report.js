// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.
// Generate an assertion-checked symptom report from NEW comparison records.
const fs=require('fs'), path=require('path'), crypto=require('crypto'), assert=require('assert');
const dir=path.resolve(process.argv[2]||'/tmp/bg-exact-issue-comparison-final');
const records=fs.readFileSync(path.join(dir,'results.jsonl'),'utf8').trim().split('\n').map(JSON.parse);
const manifest=JSON.parse(fs.readFileSync(path.join(dir,'manifest.json'),'utf8'));
const issueSource=JSON.parse(fs.readFileSync('/tmp/bg-exact-in-scope-issue-data.json','utf8'));
const hash=x=>crypto.createHash('sha256').update(x).digest('hex');
const yes=x=>x===true||x==='true';
const rows=(n,b,m='native')=>records.filter(r=>+r.issue===n&&r.build===b&&r.mode===m);
const g=(r,k='result')=>r.geometries[k];
const valid=(r,k='result')=>yes(g(r,k).rational_valid);
const publicValid=(r,k='result')=>valid(r,k)&&yes(g(r,k).public_valid);
const count=(r,k='result')=>+g(r,k).components;
const parts=r=>['result','reverse'].filter(k=>r.geometries?.[k]);
const allValid=rs=>rs.every(r=>parts(r).every(k=>publicValid(r,k)));
const allCount=(rs,c)=>rs.every(r=>parts(r).every(k=>count(r,k)===c&&valid(r,k)));
const area=(r,k='result')=>+g(r,k).area_double;
function numbers(wkt) {
  return (wkt.match(/-?\d+(?:\/\d+)?/g)||[]).map(x=>{const[n,d='1']=x.split('/');return Number(n)/Number(d);});
}
function xmin(geometry) {return Math.min(...numbers(geometry.wkt).filter((_,i)=>i%2===0));}
const checks=(rs,key)=>rs.every(r=>yes(r.checks?.[key]));
function measured(r) {
  return {status:r.status,predicate:r.predicate,reverse_predicate:r.reverse_predicate,matrix:r.matrix,
    reverse_matrix:r.reverse_matrix,checks:r.checks,reference:r.reference,result_shared_x:r.result_shared_x,
    outputs:Object.fromEntries(Object.entries(r.geometries||{}).filter(([k])=>
      k!=='a'&&k!=='b'&&!/^input\d+$/.test(k)).map(([k,v])=>
        [k,Object.fromEntries(Object.entries(v).filter(([f])=>f!=='wkt'))]))};
}
function sweep(rs) {
  let equal=true;
  for(let y=0;y<8;++y) {
    const a=rs.find(r=>r.case===`issue_858_sweep_${y}_0`),b=rs.find(r=>r.case===`issue_858_sweep_${y}_1`);
    assert(a&&b); equal=equal&&a.result_shared_x===b.result_shared_x;
  }
  return equal;
}
// These predicates encode the reported symptom, not a universal correctness
// test. In particular exact set inequality alone is never a fix criterion.
const symptomPass={
  1506:rs=>allCount(rs,1)&&rs.every(r=>parts(r).every(k=>Math.abs(area(r,k)-area(r,'a')-area(r,'b'))<1e-18)),
  1490:rs=>allCount(rs,1)&&allValid(rs),
  1420:rs=>allCount(rs,4)&&rs.every(r=>parts(r).every(k=>area(r,k)===32)),
  1421:rs=>rs.every(r=>allCount([r],r.operation==='difference'?6:1)&&area(r)===(r.operation==='difference'?24:36)),
  1419:rs=>allCount(rs,1)&&allValid(rs),
  1414:rs=>allCount(rs,1)&&allValid(rs),
  1241:rs=>allCount(rs,1)&&allValid(rs),
  1253:rs=>allCount(rs,0)&&checks(rs,'linear_oracle'),
  1238:rs=>rs.every(r=>!yes(r.predicate)),
  1232:rs=>checks(rs,'head')&&checks(rs,'tail')&&checks(rs,'reverse_set')&&allValid(rs),
  1211:rs=>allCount(rs,1)&&rs.every(r=>area(r)>0&&area(r)<0.01)&&allValid(rs),
  1179:rs=>allValid(rs)&&rs.every(r=>parts(r).every(k=>count(r,k)>0)),
  1178:rs=>allCount(rs,2)&&allValid(rs),
  1174:rs=>rs.every(r=>[1,2,3,4].every(i=>count(r,`step${i}`)>0&&valid(r,`step${i}`))),
  1044:rs=>rs.every(r=>parts(r).every(k=>valid(r,k)&&Math.abs(area(r,k)-(+r.reference.top_chain_area_double))<1e-12)),
  1035:rs=>allValid(rs),
  858:rs=>sweep(rs)&&allValid(rs),
  831:rs=>rs.every(r=>yes(r.predicate)&&yes(r.reverse_predicate)&&r.matrix==='101F00212'&&r.reverse_matrix===r.matrix),
  528:rs=>checks(rs,'boundary_oracle')&&rs.every(r=>count(r)===0),
  511:rs=>allValid(rs)&&rs.every(r=>xmin(g(r))>=xmin(g(r,'b'))-1e-12),
};
const category={A:'Already Passing Reported Symptom',B:'Reported Public Symptom Removed',
 C:'Retained Only; Not Publicly Solved',D:'Binary Contract Conflicts With Intended Answer',
 E:'Unavailable',F:'Legacy Dimensional Failure',P:'Unsupported Predicate Premise',
 R:'Additional Reference Improvement; Not Counted'};
// Categories are conservative: R does not silently substitute a different
// symptom for the report's quoted validity/direction/spurious-line failure.
const specs=[
 [1506,'B','B','Reverse union empty -> one component; both orders retain input-area sum.','full pair; both orders'],
 [1490,'B','A','Develop: both invalid difference variants -> one valid native component / one hole each. PR1491 already passes both.','two local pairs; integer output; old grid area expectation 3900 is not a valid oracle'],
 [1487,'A','A','Both baselines already return one component, area about 18983.885; no missing-intersection fix counted.','selected strict-FP pair; exact native default-validator regression in forward order'],
 [1448,'F','F','One invalid singleton linestring in baseline, native exact and retained exact.','legacy dimensional/cardinality contract unchanged; not rounding'],
 [1421,'B','A','Develop: difference 4 -> 6 components, area 16 -> 24; symmetric difference 2 -> 1, area 28 -> 36. PR1491 already passes.','both reported operations; small integer graph case'],
 [1420,'B','A','Develop: forward union 3 components / area 36 -> 4 / 32; reverse also gives 4 / 32. PR1491 already passes.','small integer traversal; both orders'],
 [1419,'B','B','Second native-feedback difference empty -> one valid component.','all three inputs; native feedback AND retained exact run; polygon operands singleton-wrapped in pipeline'],
 [1414,'B','B','Forward union empty -> one valid component; both orders agree.','full large pair; no independent complete set golden'],
 [1406,'R','A','Develop loses the contact point; exact returns (1,1.5), no line. PR1491 already does so. Neither baseline has the quoted spurious line after normalization.','two normalized variants; polygon/line two-slot tuple; original wrong orientation not endorsed'],
 [1397,'C','C','Retained intersection valid; final integer output still self-intersects.','normalized CCW/open source model; native integer rounding still unsolved'],
 [1395,'A','A','Baseline and exact native intersections valid under these strict-FP normalized calls.','original MSVC /fp:fast and custom CCW/open model not reproduced'],
 [1394,'R','R','Reported invalidity absent, but baseline union drops a large component: area 463619.5 -> 1842454.','gross new reference improvement, NOT quoted invalidity fix; integer native area differs from retained'],
 [1363,'R','R','Baseline returns bogus positive-area intersection in one order; retained/native exact empty. Quoted invalidity is not reproduced.','local regression pair; external Godbolt program not independently replayed'],
 [1354,'A','A','Baseline already retains one component / one hole, same rounded result area.','full local pair; no inherited validity suppression'],
 [1290,'P','P','overlaps == !touches is not a general identity; no fix inferred from that premise.','native feedback and exact retained branches run; corrected polygons replace source rings'],
 [1253,'B','B','Two large outside lines -> empty in both directions, matching independent rational clipping.','downloaded binary64 inputs; intended corner segment is not the binary oracle'],
 [1241,'B','B','Empty difference -> one valid component.','full pair; source polygon API shapes preserved'],
 [1238,'B','A','Develop within true -> false; PR1491 already false in the report call order.','native multi-linestring predicate; exact false'],
 [1236,'C','C','Native-feedback within remains false; retained exact within true. Rounded inputs really fail rational within.','all intermediates native-valid, so this is predicate persistence, not invalidity'],
 [1235,'A','A','Baseline already has zero isolated points / lines and two valid polygons.','polygon/multipolygon tuple API; only rounded-coordinate differences'],
 [1232,'B','B','Missing (7,4) tail and reversal mismatch removed; (5,-8) head also retained.','native line/multi-line difference; partial oracle, not a complete set golden'],
 [1212,'A','A','Both reduced and full baseline results already retain two holes.','int64 adaptation; original mpq_rational API not reproduced'],
 [1211,'B','B','Oversized intersection area 0.05117 -> about 6.30e-11, below the reported 0.01 bound.','native float parsed first; source CCW/open model normalized'],
 [1201,'D','D','Gross million-area error disappears, but binary input has nonempty area about 1.25e-28, not the requested empty answer.','excluded from strict fix count; exact native default validator rejects tiny output'],
 [1179,'B','B','All four pairs: intersecting-interior native union outputs -> valid native unions, in both orders.','all four supplied pairs; no independent complete output WKT'],
 [1178,'B','A','Develop: both disconnected-interior variants -> two native-valid components. PR1491 already passes both.','all supplied variants; small integer topology'],
 [1174,'B','B','Two empty quadrants -> four nonempty rational-valid rounded pieces.','single polygon/box calls; native default validator rejects the tiny-edge pieces; synthetic concatenation is NOT a valid multipolygon contract'],
 [1053,'A','A','Both supplied self-coverage/within selections already true.','full multipolygon plus original polygon selection; 0.075 comment variant excluded'],
 [1044,'B','B','Bug variant area 0.24 -> rational upper-chain area about 0.09; shifted/control variants pass.','multi/multi source API plus polygon pair; all four original/shifted bug/control variants; both orders'],
 [1035,'B','B','Invalid integer difference -> native-valid difference.','normalized CCW/open model; actual integer output changed, not just valid retained output'],
 [1034,'C','C','Retained difference has four valid components; rounded integer output still invalid.','full ticket_12752 pair; no spike cleanup'],
 [900,'C','C','Native feedback changes wrong-area result but remains invalid (10 components); retained exact has four valid components and passes identity.','all four stars; native feedback != retained exact'],
 [888,'A','A','All 38 baseline native intermediates already valid; final 80 components.','all 39 inputs / balanced reduction plus three reduced pairs; no independent final-set golden'],
 [858,'B','B','Shared integer endpoints disagree in baseline sweep; all eight pairs agree with exact-native nearest-even output.','16 source ring/line calls plus line reversal; native integer output policy is explicitly part of this fix'],
 [852,'A','A','Baseline and exact native intersection already valid; ULP-level change is not a fix.','full pair; no blanket rounded-output guarantee'],
 [831,'B','B','crosses true/false under line reversal -> true/true, both matrices 101F00212.','native float line/polygon API; rational matrix reference from audit'],
 [704,'D','D','Binary overlap is real (matrix 212101212); empty/shared-decimal-edge expectation not imposed.','corrected comment pair; all native outputs already nonempty; no ULP-only fix count'],
 [669,'A','A','No reported gross difference failure on baseline full/crop/rotation selections; tiny binary overlap is legitimate.','full downloaded CCW source API, rotations, native crop feedback, retained crop and default-crop feedback; not an unchanged-area oracle'],
 [663,'D','D','Baseline empty; retained exact tiny nonempty overlap; native exact collapses to invalid zero-area polygon.','binary intended-empty conflict AND native representation regression; no fix'],
 [656,'A','A','Baseline already equals the exact boundary chain, with no extra outside line.','native multi-line/multi-polygon API; independent rational clipping'],
 [653,'E','E','No numeric coordinates available; no synthetic substitute.','not one of the 47 selected executions'],
 [631,'A','A','Baseline already returns the supplied rectangle, area 214748364800.','INT_MIN case; no claim for all overflow configurations'],
 [602,'C','C','Retained four-component difference valid; rounded integer difference still invalid.','full ticket_11674 pair; no validity suppression'],
 [588,'A','A','Baseline native difference within red already true; frozen body-WKT within also true.','native feedback plus retained exact; later high-precision printed triangle not replayed'],
 [528,'B','B','One spurious boundary intersection -> zero, matching independent rational segment solve.','original uint64 bit patterns; native polygon/line point-output API'],
 [511,'B','B','Long left outside prefix removed: minimum result x about -21.856 -> -15.2803.','normalized source polygon-first API; retained output matches independent clipping; final endpoint rounding is not exact set equality'],
 [510,'A','A','Baseline already returns the same native tiny single segment; exact-retained endpoints improve internally only.','normalized source polygon-first API; native rational set inequality due to rounding is not counted'],
 [440,'R','R','Both baselines already direction-consistent, but wrong point y=-2690131; exact-native y=-2665585 is the correct rounded reference.','eight native integer segment directions/orders; new gross reference improvement, NOT the reported direction-asymmetry fix'],
];
assert.equal(specs.length,48); assert.equal(new Set(specs.map(s=>s[0])).size,48);
for(const r of records) {
  assert.equal(r.status,'ok',`${r.build}/${r.case}/${r.mode}`);
  assert.equal(r.evaluation,'ok',`${r.build}/${r.case}/${r.mode}`);
  assert.equal(r.process.status,'0');
}
// Compare initial parsed inputs, not intermediates that intentionally differ.
const inputDigest=r=>hash(JSON.stringify(Object.fromEntries(Object.entries(r.geometries||{})
  .filter(([k])=>k==='a'||k==='b'||/^input\d+$/.test(k)).map(([k,v])=>[k,v.wkt]))));
for(const r of records.filter(r=>r.build==='develop')) {
  for(const b of ['pr1491','before']) for(const m of b==='before'?['native','retained']:['native']) {
    const other=records.find(s=>s.case===r.case&&s.build===b&&s.mode===m);
    assert(other,`missing ${b}/${r.case}/${m}`); assert.equal(inputDigest(r),inputDigest(other),`input mismatch ${r.case}`);
  }
}
for(const r of records) for(const [k,v] of Object.entries(r.geometries||{}))
  if(k==='a'||k==='b'||/^input\d+$/.test(k)) assert(yes(v.rational_valid),`invalid original input: ${r.case}/${k}`);
const totals={develop:{},pr1491:{}};
const issues=specs.map(([issue,develop,pr1491,observation,scope])=>{
  const source=issueSource.issues.find(s=>s.number===issue); assert(source);
  const result={issue,url:source.url,source_sha256:hash(JSON.stringify(source)),observation,scope,
    full_issue_fixed:false,full_issue_status:'not established; selected-case comparison only',baselines:{}};
  for(const [b,c] of [['develop',develop],['pr1491',pr1491]]) {
    const base=rows(issue,b),exact=rows(issue,'before'),retained=rows(issue,'before','retained');
    if(issue!==653) {assert(base.length);assert.equal(base.length,exact.length);assert.equal(base.length,retained.length);}
    const outcome={category:c,classification:category[c],reported_public_fix:c==='B',
      baseline_source:manifest.revisions[b],exact_source:manifest.revisions.before,
      attribution:b==='develop'?'combined PR1491 algorithm fixes plus exact branch; not exact-numerics attribution'
        :'incremental exact-branch improvement over PR1491; not a numerical-only causal isolation',
      flags:{selected_case_only:true,scope_limits:scope,baseline_strategy:'public default; no strategy argument',
        exact_strategy:'cartesian_exact<NativeCoordinate>; frozen before',
        coordinate_types:[...new Set(base.map(r=>r.coordinate))],
        native_feedback_pipeline:base.some(r=>r.operation==='pipeline'||r.case.includes('crop')),
        quoted_symptom_not_reproduced:c==='A'||c==='R',
        retained_only:c==='C',binary_contract_conflict:c==='D',
        integer_final_output:base.some(r=>r.coordinate.startsWith('int'))},
      cases:base.map(r=>({case:r.case,initial_input_sha256:inputDigest(r),
        baseline_record:path.join(dir,'runs',`${b}-${r.case}-native.json`),
        exact_native_record:path.join(dir,'runs',`before-${r.case}-native.json`),
        exact_retained_record:path.join(dir,'runs',`before-${r.case}-retained.json`),
        baseline:measured(r),exact_native:measured(exact.find(e=>e.case===r.case)),
        exact_retained:measured(retained.find(e=>e.case===r.case))}))};
    if(symptomPass[issue]&&base.length) {
      outcome.baseline_reported_symptom_pass=symptomPass[issue](base);
      outcome.exact_native_reported_symptom_pass=symptomPass[issue](exact);
      if(c==='B') {assert(!outcome.baseline_reported_symptom_pass,`${issue}/${b}: no baseline symptom`);
        assert(outcome.exact_native_reported_symptom_pass,`${issue}/${b}: native symptom persists`);}
    }
    if(c==='C') {
      if(issue===1236) {assert(exact.every(r=>!yes(r.predicate)));assert(retained.every(r=>yes(r.predicate)));}
      else {assert(exact.some(r=>!valid(r)));assert(retained.every(r=>valid(r)));}
    }
    totals[b][c]=(totals[b][c]||0)+1; result.baselines[b]=outcome;
  }
  return result;
});
const provenanceFiles=['exact_issue_comparison.cpp','exact_issue_comparison_evaluate.cpp','exact_issue_comparison_run.js',
  'exact_issue_comparison_report.js','exact_issue_comparison_optimized.js','exact_issue_cases.hpp','exact_additional_issue_cases.hpp','exact_cartesian_issues.md',
  '../algorithms/overlay/overlay_cases.hpp','../algorithms/overlay/multi_overlay_cases.hpp'];
const provenance=Object.fromEntries(provenanceFiles.map(f=>[f,hash(fs.readFileSync(path.resolve(__dirname,f)))]));
for(const [f,h] of Object.entries(manifest.source_sha256)) assert.equal(provenance[f],h,`source changed: ${f}`);
const summary={criterion:'B: a quoted reported public cardinality/validity/predicate or gross outside-output symptom is observed in the selected native baseline call(s), and removed for all selected variants with exact strategy AND native final output. Not retained-only, intended-answer conflicts, additional reference-only improvements, or ULP-only differences.',
  selected_issues:47,unavailable:1,records:records.length,
  case_groups:new Set(records.map(r=>r.case)).size,
  in_scope_case_groups:new Set(records.filter(r=>+r.issue!==494).map(r=>r.case)).size,
  representation_control_issue:494,
  counts:totals,category,manifest:path.join(dir,'manifest.json'),
  source_sha256:provenance,issue_source_snapshot_sha256:hash(fs.readFileSync('/tmp/bg-exact-in-scope-issue-data.json')),
  compiler:manifest.compiler,flags:manifest.flags,revisions:manifest.revisions,headers:manifest.headers,issues};
summary.attribution={
  develop:{count:totals.develop.B,meaning:'combined PR1491 algorithm fixes plus exact branch versus official develop',
    numerical_fixes_attributable_to_exact_strategy:'not established'},
  pr1491:{count:totals.pr1491.B,meaning:'incremental selected-symptom improvements of the exact branch versus PR1491',
    numerical_only_causal_isolation:false},
  already_pr1491_fixes:issues.filter(i=>i.baselines.develop.category==='B'&&i.baselines.pr1491.category==='A').map(i=>i.issue)};
assert.deepEqual([...summary.attribution.already_pr1491_fixes].sort((a,b)=>a-b),[1178,1238,1420,1421,1490]);
const optimizedDir=process.argv[3]||process.env.COMPARISON_OPTIMIZED_RESULTS;
if(optimizedDir) {
  const location=path.resolve(optimizedDir);
  const optimized=JSON.parse(fs.readFileSync(path.join(location,'comparison.json'),'utf8'));
  const run=JSON.parse(fs.readFileSync(path.join(location,'manifest.json'),'utf8'));
  assert.equal(optimized.reference,dir,'Optimized run used a different frozen reference');
  assert.equal(run.reference_manifest_sha256,hash(fs.readFileSync(path.join(dir,'manifest.json'))));
  assert.equal(optimized.manifest_sha256,hash(fs.readFileSync(path.join(location,'manifest.json'))));
  assert.equal(optimized.records,164);assert.equal(optimized.case_groups,82);assert.equal(optimized.selected_issues,47);
  assert(optimized.headers_unchanged);assert.equal(optimized.failures.length,0);
  assert.equal(optimized.differences.length,0,'Investigate optimized semantic differences before publishing unchanged counts');
  summary.optimized_regression={status:'verified',location,comparison_sha256:hash(fs.readFileSync(path.join(location,'comparison.json'))),
    header_manifest:optimized.optimized_header_manifest,headers_unchanged:true,records:164,case_groups:82,selected_issues:47,
    semantic_differences:0,public_validator_differences:optimized.public_validator_differences,
    compiler:run.compiler,flags:run.flags,compile_jobs:run.compile_jobs,started:run.started,finished:run.finished,
    source_sha256:run.source_sha256,runner_sha256:run.runner_sha256};
} else summary.optimized_regression={status:'pending stable-header compile window; frozen counts only'};
const json=path.join(__dirname,'exact_issue_comparison_results.json');
fs.writeFileSync(json,JSON.stringify(summary,null,2)+'\n');
const list=(b,c)=>issues.filter(i=>i.baselines[b].category===c).map(i=>'#'+i.issue).join(', ');
let md=`# Exact Issue Symptom Comparison\n\n`;
md+=`**Attribution: frozen exact \`bc65c6d60\` already includes PR1491 algorithm fixes. The ${totals.develop.B} selected-symptom improvements versus official develop belong to the combined PR1491-plus-exact branch, NOT ${totals.develop.B} numerical fixes attributable to the exact strategy. ${totals.pr1491.B} are incremental versus PR1491. The other five (#1490, #1421, #1420, #1238, #1178) already pass with PR1491's algorithm fixes.** Even the incremental comparison is a branch comparison, not a numerical-only causal isolation of the strategy.\n\n`;
if(summary.optimized_regression.status==='verified') {
  const opt=summary.optimized_regression;
  md+=`## Optimized Regression\n\nOptimized live exact headers were compared with frozen-before on **all 82 in-scope groups, native and retained: 164 records across 47 unique issues**. **Zero rational-WKT, predicate, or other semantic-field differences**; ${opt.public_validator_differences.length} records with separately listed public-validator diagnostic differences. No baseline/default rebuilds; the same frozen evaluator was reused. The full live header manifest was identical before and after seven serial optimized-only builds and all runs. Header manifest SHA256: \`${opt.header_manifest.sha256}\`. Detailed per-record comparisons, compiler commands, header hashes, and logs: \`${opt.location}\`. The frozen symptom classifications/counts therefore remain unchanged under this optimized header snapshot.\n\n`;
} else md+=`Optimized live-header regression verification is pending a stable-header compile window; the counts below are frozen-before comparisons.\n\n`;
md+=`Frozen comparison of all **47 executable unique issues**, plus unavailable #653. **${summary.in_scope_case_groups} in-scope case groups**, plus the #494 representation control, produce **${records.length} new build/mode records** (${summary.case_groups} groups times four modes: develop default-native, PR1491 default-native, frozen exact-native, frozen exact-retained). Each group may execute multiple direction/order calls or pipeline steps. These are not old unit-test pass totals. No whole issue is declared fixed.\n\n`;
md+=`## Counts\n\n| Classification | Develop | PR1491 |\n| --- | ---: | ---: |\n`;
for(const c of ['B','A','C','D','R','F','P','E']) md+=`| ${c}: ${category[c]} | ${totals.develop[c]||0} | ${totals.pr1491[c]||0} |\n`;
md+=`\n**Strict selected reported-public-symptom improvements: ${totals.develop.B} for combined PR1491-plus-exact versus develop; ${totals.pr1491.B} incremental versus PR1491.** Denominator 47 executed issue selections, not 46 passing exact contracts. The table totals 48 only because unavailable #653 is included separately.\n\n`;
md+=`Criterion B: ${summary.criterion.slice(3)}\n\nDevelop B: ${list('develop','B')}.\n\nPR1491 B: ${list('pr1491','B')}.\n\n`;
md+=`R is intentionally excluded: a newly measured reference defect can improve without reproducing the report's quoted symptom. C is also excluded even if every retained-exact contract passes. D does not turn an intended decimal answer into a binary-input fix claim. All classifications concern these selections and this toolchain, not every source model/compiler/report variant.\n\n`;
md+=`## Provenance\n\n- Develop: \`${manifest.revisions.develop}\`.\n- PR1491: \`${manifest.revisions.pr1491}\`.\n- Exact: frozen before \`${manifest.revisions.before}\`, never live optimization headers.\n- Every exported header was checked against its commit's Git blob; complete path/blob/SHA256 manifests are in \`${dir}/{develop,pr1491,before}-headers.tsv\`.\n- Flags: \`${manifest.flags.join(' ')}\`; default public API with no strategy for either baseline; \`COMPARISON_EXACT\` uses \`cartesian_exact<NativeCoordinate>\` on frozen before. No fast-math, rescaling, test-failure, or no-robustness macro.\n- Compiler: ${manifest.compiler.split('\n')[0]}. Compilation serial, one compiler; Perl alarm bounds compile/runtime/evaluation. All ${records.length} records executed and evaluated without exception, timeout, or process failure.\n- Same native parse type, coordinates, orientation/closure correction on all builds. Initial rational-WKT input digests match across all four build/mode combinations; all original selected inputs are rational-valid. No snapping/simplification.\n\n`;
md+=`Data hashes and every per-issue issue-body/comment hash, baseline commit, exact commit, input digest, flags/limits, and raw record paths are in [exact_issue_comparison_results.json](exact_issue_comparison_results.json). Machine logs and full rational-coordinate output are in \`${dir}\`; the manifest records full compiler identity and source SHA256 values.\n\n`;
md+=`## Reading The Evidence\n\nPublic outputs use the native coordinate type, including integral outputs. Retained outputs preserve exact constructed coordinates. Both run separately; pipelines genuinely feed rounded native output back in the native branch. Rational validation always uses frozen-before \`cartesian<cpp_rational>\`, with source values lifted from native bits (0.1 lifting is tested). Each build's default public validator is recorded separately. A default-validator rejection and a rationally invalid geometry are not conflated.\n\n`;
md+=`Pair tests preserve polygon versus multipolygon operand shapes, with native float/int parsing first; areal output uses multipolygon containers. Source-specific API groups cover line/areal tuple, segment/point, ring/line, box/line, native float crosses, bit-pattern points, downloaded multi/multi arrays, and CCW full/cropped difference. Source-model normalizations and remaining pipeline singleton wrapping are explicit below.\n\n`;
md+=`Independent rational segment solves/clipping support #528, #1253, #510, #511, #656 and #858; the downloaded #1044 upper-chain integral is independent of overlay. Other area/set/reference checks share BG machinery and are not independent complete-set goldens. A native rounded line normally differs from its rational clipping set: #510 is the control demonstrating why this alone is not a counted fix.\n\n`;
md+=`## Per Issue\n\nA/B/C/D/E/F/P/R refer to the counts above. Every row has a separate outcome for each baseline, even when equal.\n\n| Issue | Develop | PR1491 | Measured Outcome | Scope / Flags |\n| --- | --- | --- | --- | --- |\n`;
for(const i of issues) md+=`| #${i.issue} | ${i.baselines.develop.category} | ${i.baselines.pr1491.category} | ${i.observation} | ${i.scope} |\n`;
md+=`\n## Important Qualifications\n\n- **#1487 is not counted against either baseline.** Both native original polygon calls already return a nonempty intersection of the expected scale in both orders here. Develop includes PR1508; no inference about its necessity follows from selected exact tests. Exact native forward output is rational-valid but frozen-default \`is_valid\` reports a spike. This new public-validator failure is recorded, not hidden as an issue fix.\n- **#1034/#602/#1397 remain native-invalid** despite valid exact retained output. **#1236 remains false** with native feedback, even under exact predicates; rational validation confirms the rounded-set within failure. **#900 remains native-invalid**, despite an exact-retained identity pass. These five are C, not fixes.\n- **#1174** corrects the quoted empty/wrong quadrants with four nonempty rational-valid native pieces. Frozen-default \`is_valid\` rejects the tiny-edge pieces; no default-validation fix is claimed. Concatenating adjacent quadrant polygons is only a synthetic area aggregate, not a valid output multipolygon contract, for either exact or native coordinates.\n- **#858 integer rounding is a public behavior change**, not merely better internal endpoints: all eight adjacent clips now share an identical native integer endpoint. **#440** also changes native integer output to correctly rounded y=-2665585, but baseline directions already agree, so it is R, not a reported direction-symptom fix. **#1035** becomes native-valid; #1034/#602/#1397 do not. Integer rounding cannot be advertised as generally topology-preserving.\n- **#1490 grid reference correction:** the frozen fixture's 3900 reference is wrong for its listed coordinates. Outer area=5550, hole area=150; subtractor area=1500 contains that hole, so the difference area is 5550-1500=**4050**. The old 3900 result subtracts the 150-area hole twice. This report leaves fixtures/audit untouched and counts only the observed invalid-result removal, supported by the corrected area calculation.\n- **#1201** removes the huge wrong result but does not satisfy the report's intended empty answer; binary overlap is real and default native validation rejects the tiny result. **#704** already has tiny nonempty native overlap. **#663** changes baseline empty output into an invalid rounded tiny polygon. These three are D, not B.\n- **#1394/#1363 and develop #1406** improve gross public reference behavior but not the quoted invalidity/spurious-line symptom in these normalized baseline runs. R keeps those improvements separate from the strict count.\n- **#1448** stays an invalid singleton even before rounding; a dimensional contract fix is separate. **#653** has no reproducible numeric input. No fabricated data or whole-issue upper-bound count is used.\n\n`;
md+=`## Reproduce\n\nFrom this worktree, with frozen headers at the stated paths:\n\n\`\`\`sh\nnode test/strategies/exact_issue_comparison_run.js /tmp/bg-exact-issue-comparison-final\nnode test/strategies/exact_issue_comparison_report.js /tmp/bg-exact-issue-comparison-final\n\`\`\`\n\nFor optimized-only verification after the frozen run exists, without rebuilding defaults or the evaluator:\n\n\`\`\`sh\nnode test/strategies/exact_issue_comparison_optimized.js /tmp/bg-exact-issue-comparison-optimized\nnode test/strategies/exact_issue_comparison_report.js /tmp/bg-exact-issue-comparison-final /tmp/bg-exact-issue-comparison-optimized\n\`\`\`\n\n\`COMPARISON_INCLUDE\` can override the optimized header path; \`COMPARISON_REFERENCE\` can override the frozen reference directory. No production headers, existing cases, build registration, remote branches, or tracker state are modified by these tools.\n`;
fs.writeFileSync(path.join(__dirname,'exact_issue_comparison.md'),md);
console.log(JSON.stringify({counts:totals,records:records.length,strict:{develop:totals.develop.B,pr1491:totals.pr1491.B}},null,2));
