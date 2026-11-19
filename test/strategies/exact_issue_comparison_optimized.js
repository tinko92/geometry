// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.
// Build ONLY optimized exact variants, then compare with frozen-before records.
// COMPARISON_INCLUDE overrides live include; COMPARISON_REFERENCE selects the
// completed frozen run. No baseline or evaluator compilation is performed.
const fs=require('fs'),path=require('path'),crypto=require('crypto'),cp=require('child_process'),assert=require('assert');
const repo=path.resolve(__dirname,'../..');
const out=path.resolve(process.argv[2]||'/tmp/bg-exact-issue-comparison-optimized');
const reference=path.resolve(process.env.COMPARISON_REFERENCE||'/tmp/bg-exact-issue-comparison-final');
const include=path.resolve(process.env.COMPARISON_INCLUDE||path.join(repo,'include'));
const boost='/Users/tinko/dev/boost';
const sha=x=>crypto.createHash('sha256').update(x).digest('hex');
const referenceManifest=JSON.parse(fs.readFileSync(path.join(reference,'manifest.json'),'utf8'));
const referenceRecords=fs.readFileSync(path.join(reference,'results.jsonl'),'utf8').trim().split('\n').map(JSON.parse)
  .filter(r=>r.build==='before'&&+r.issue!==494);
const flags=referenceManifest.flags;
const evaluator=path.join(reference,'evaluate');
assert(fs.existsSync(evaluator),'Frozen evaluator missing; do not rebuild defaults');
fs.mkdirSync(path.join(out,'runs'),{recursive:true});
function bounded(command,args,seconds,input) {
  return cp.spawnSync('perl',['-e','alarm shift; exec @ARGV or die $!',String(seconds),command,...args],
    {cwd:repo,encoding:'utf8',input,timeout:(seconds+3)*1000,maxBuffer:128*1024*1024});
}
function headers() {
  const files=cp.execFileSync('rg',['--files',include],{encoding:'utf8'}).trim().split('\n').sort();
  const entries=files.map(f=>`${path.relative(include,f)}\t${sha(fs.readFileSync(f))}`);
  const text=entries.join('\n')+'\n';
  return {files:files.length,sha256:sha(text),text};
}
const start=headers();
fs.writeFileSync(path.join(out,'optimized-headers.tsv'),start.text);
const sourceNames=['exact_issue_comparison.cpp','exact_issue_comparison_evaluate.cpp','exact_issue_cases.hpp',
  'exact_additional_issue_cases.hpp','exact_cartesian_issues.md'];
const sources=Object.fromEntries(sourceNames.map(s=>[s,sha(fs.readFileSync(path.join(__dirname,s)))]));
// The legacy public validator creates auxiliary exact working points. Only its
// lifetime guard differs from the frozen runner; fixtures and evaluator do not.
const scopeGuard = '#ifdef COMPARISON_EXACT\n    ec::point_cache_scope operation;\n#endif\n';
for(const [s,h] of Object.entries(referenceManifest.source_sha256)) {
  const source = fs.readFileSync(path.join(__dirname,s),'utf8');
  const original = s === 'exact_issue_comparison.cpp' ? source.replace(scopeGuard,'') : source;
  assert.equal(sha(original),h,`Frozen runner input changed beyond operation scope: ${s}`);
}
const manifest={label:'optimized',include,reference,reference_manifest_sha256:sha(fs.readFileSync(path.join(reference,'manifest.json'))),
  frozen_exact_source:referenceManifest.revisions.before,source_sha256:sources,
  runner_adaptation:'An enclosing exact-operation scope only; all other runner bytes, fixtures and evaluator are hash-verified unchanged.',
  runner_sha256:sha(fs.readFileSync(__filename)),compiler:cp.execFileSync('clang++',['--version'],{encoding:'utf8'}),
  flags,header_manifest:{files:start.files,sha256:start.sha256},evaluator_sha256:sha(fs.readFileSync(evaluator)),
  compile_jobs:1,started:new Date().toISOString(),commands:[],runs:[]};
const units=[[1,'double'],[1,'float'],[1,'int32_t'],[1,'int64_t'],[2,''],[3,''],[4,'']];
let records=[];
for(const [group,coordinate] of units) {
  const bin=path.join(out,`optimized-${group}${coordinate?'-'+coordinate:''}`);
  const args=[...flags,`-DCOMPARISON_GROUP=${group}`,'-DCOMPARISON_EXACT',
    ...(coordinate?[`-DCOMPARISON_COORDINATE=${coordinate.startsWith('int')?'std::':''}${coordinate}`,
      `-DCOMPARISON_COORDINATE_NAME="${coordinate}"`]:[]),
    `-I${include}`,`-I${boost}`,path.join(__dirname,'exact_issue_comparison.cpp'),'-o',bin];
  manifest.commands.push({command:'clang++',args});
  const signature=sha(JSON.stringify({args,sources,headers:manifest.header_manifest}));
  if(!fs.existsSync(bin+'.signature')||fs.readFileSync(bin+'.signature','utf8')!==signature) {
    console.log(`Compile ${path.basename(bin)}`);
    const p=bounded('clang++',args,240);
    fs.writeFileSync(bin+'.compile.log',(p.stdout||'')+(p.stderr||''));
    assert.equal(p.status,0,`Compile failed: ${bin} (${p.signal})`);
    fs.writeFileSync(bin+'.signature',signature);
  }
  assert.equal(bounded(bin,['--binary-check'],5).status,0);
  const list=bounded(bin,['--list'],5); assert.equal(list.status,0);
  for(const row of list.stdout.trim().split('\n')) {
    const [fixture,issue,coordinate]=row.split('\t');
    if(+issue===494) continue;
    for(const mode of ['native','retained']) {
      console.log(`Run optimized ${fixture} ${mode}`);
      const stem=path.join(out,'runs',`optimized-${fixture}-${mode}`);
      const p=bounded(bin,[fixture,mode],90);
      fs.writeFileSync(stem+'.stdout',p.stdout||''); fs.writeFileSync(stem+'.stderr',p.stderr||'');
      let r;
      try {r=JSON.parse(p.stdout);}
      catch {r={case:fixture,issue,coordinate,status:p.signal==='SIGALRM'?'timeout':'process_failure'};}
      r.build='optimized';r.mode=mode;r.process={status:p.status,signal:p.signal,error:p.error?.message};
      const ev=bounded(evaluator,[],30,JSON.stringify(r));
      fs.writeFileSync(stem+'.evaluation.stderr',ev.stderr||'');
      try {r=JSON.parse(ev.stdout);} catch {r.evaluation=ev.signal==='SIGALRM'?'timeout':'process_failure';}
      fs.writeFileSync(stem+'.json',JSON.stringify(r,null,2)+'\n');
      records.push(r);manifest.runs.push({case:fixture,mode,process:r.process,evaluation:r.evaluation});
      fs.writeFileSync(path.join(out,'results.jsonl'),records.map(r=>JSON.stringify(r)).join('\n')+'\n');
      fs.writeFileSync(path.join(out,'manifest.json'),JSON.stringify(manifest,null,2)+'\n');
    }
  }
}
const end=headers();
manifest.finished=new Date().toISOString();manifest.headers_unchanged=start.sha256===end.sha256;
manifest.records=records.length;manifest.case_groups=new Set(records.map(r=>r.case)).size;
manifest.selected_issues=new Set(records.map(r=>+r.issue)).size;
fs.writeFileSync(path.join(out,'manifest.json'),JSON.stringify(manifest,null,2)+'\n');
assert(manifest.headers_unchanged,'Live headers changed during compilation/run: results are not a stable optimized snapshot');
for(const [s,h] of Object.entries(sources)) assert.equal(sha(fs.readFileSync(path.join(__dirname,s))),h,`Runner source changed: ${s}`);
assert.equal(records.length,164);assert.equal(manifest.case_groups,82);assert.equal(manifest.selected_issues,47);
const differences=[],validatorDifferences=[],failures=[];
function diff(a,b,prefix='') {
  if(a===b) return [];
  if(a&&b&&typeof a==='object'&&typeof b==='object')
    return [...new Set([...Object.keys(a),...Object.keys(b)])].sort().flatMap(k=>diff(a[k],b[k],prefix?prefix+'.'+k:k));
  return [{path:prefix,before:a,optimized:b}];
}
const semantic=r=>Object.fromEntries(Object.entries(r).filter(([k])=>!['build','process','mode'].includes(k)));
const byKey=new Map(referenceRecords.map(r=>[r.case+'/'+r.mode,r]));
const comparisons=[];
for(const r of records) {
  const key=r.case+'/'+r.mode,old=byKey.get(key); assert(old,`Missing frozen reference: ${key}`);
  if(r.status!=='ok'||r.evaluation!=='ok'||String(r.process.status)!=='0') failures.push({case:r.case,mode:r.mode,status:r.status,evaluation:r.evaluation,process:r.process});
  const all=diff(semantic(old),semantic(r));
  const validator=all.filter(d=>/^geometries\.[^.]+\.public_(valid|reason)$/.test(d.path));
  const other=all.filter(d=>!validator.includes(d));
  if(other.length) differences.push({case:r.case,issue:r.issue,mode:r.mode,differences:other});
  if(validator.length) validatorDifferences.push({case:r.case,issue:r.issue,mode:r.mode,differences:validator});
  comparisons.push({case:r.case,issue:r.issue,mode:r.mode,semantic_identical:other.length===0,
    public_validator_identical:validator.length===0,
    frozen_record_sha256:sha(fs.readFileSync(path.join(reference,'runs',`before-${r.case}-${r.mode}.json`))),
    optimized_record_sha256:sha(fs.readFileSync(path.join(out,'runs',`optimized-${r.case}-${r.mode}.json`)))});
}
const comparison={reference,optimized:out,manifest_sha256:sha(fs.readFileSync(path.join(out,'manifest.json'))),
  frozen_exact_source:manifest.frozen_exact_source,optimized_header_manifest:manifest.header_manifest,
  headers_unchanged:manifest.headers_unchanged,records:records.length,case_groups:manifest.case_groups,
  selected_issues:manifest.selected_issues,failures,differences,public_validator_differences:validatorDifferences,
  policy:'Compare every rational WKT and semantic/predicate/evaluation field. Public validator diagnostics are listed separately, never silently treated as geometry equality.',
  comparisons};
fs.writeFileSync(path.join(out,'comparison.json'),JSON.stringify(comparison,null,2)+'\n');
console.log(JSON.stringify({records:records.length,groups:manifest.case_groups,issues:manifest.selected_issues,
  failures:failures.length,semantic_differences:differences.length,public_validator_differences:validatorDifferences.length,
  headers_unchanged:manifest.headers_unchanged},null,2));
if(failures.length||differences.length||validatorDifferences.length) process.exitCode=1;
