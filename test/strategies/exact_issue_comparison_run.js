// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.
// node test/strategies/exact_issue_comparison_run.js [output-directory]
// Deliberately serial compilation: never more than one compiler at a time.
const fs = require('fs');
const path = require('path');
const crypto = require('crypto');
const cp = require('child_process');
const repo = path.resolve(__dirname, '../..');
const out = path.resolve(process.argv[2] || '/tmp/bg-exact-issue-comparison');
const boost = '/Users/tinko/dev/boost';
const frozen = '/tmp/bg-exact-performance';
const sources = ['exact_issue_comparison.cpp', 'exact_issue_comparison_evaluate.cpp',
  'exact_issue_cases.hpp', 'exact_additional_issue_cases.hpp', 'exact_cartesian_issues.md'];
const sha = data => crypto.createHash('sha256').update(data).digest('hex');
fs.mkdirSync(out, {recursive:true});
fs.mkdirSync(path.join(out,'runs'), {recursive:true});
function bounded(command, args, seconds, input) {
  return cp.spawnSync('perl', ['-e','alarm shift; exec @ARGV or die $!',String(seconds),command,...args],
    {cwd:repo, encoding:'utf8', input, timeout:(seconds+3)*1000, maxBuffer:128*1024*1024});
}
const flags = ['-std=c++14','-O0','-g0','-ffp-contract=off','-pthread',
  '-DBOOST_GEOMETRY_DISABLE_DEPRECATED_14_WARNING'];
const hashes = Object.fromEntries(sources.map(s=>[s,sha(fs.readFileSync(path.join(__dirname,s)))]));
const revisions = {develop:'0b4851a77d91633cf86c3a5e679e67d8c60200e9',
  pr1491:'72b2d832fb9596a2a23c0c164a5426ae720d0cd6',before:'bc65c6d60fe2d891ed978d0c8cdf91bb678712b9'};
const manifest = {repo,boost,flags,compiler:cp.execFileSync('clang++',['--version'],{encoding:'utf8'}),
  revisions,source_sha256:hashes,started:new Date().toISOString(),compile_jobs:1,
  timeouts:{compile:240,native:20,exact:90,evaluate:30},headers:{},runs:[]};
// Verify the exported trees against Git blob identities, not directory names.
for(const [name,revision] of Object.entries(revisions)) {
  const tree=cp.execFileSync('git',['ls-tree','-r',revision,'include'],{cwd:repo,encoding:'utf8'});
  const entries=[];
  for(const row of tree.trim().split('\n')) {
    const [meta,file]=row.split('\t'), expected=meta.split(' ')[2];
    const data=fs.readFileSync(path.join(frozen,name,file));
    const actual=crypto.createHash('sha1').update(`blob ${data.length}\0`).update(data).digest('hex');
    if(actual!==expected) throw Error(`Frozen tree mismatch: ${name}/${file}`);
    entries.push(`${file}\t${expected}\t${sha(data)}`);
  }
  fs.writeFileSync(path.join(out,`${name}-headers.tsv`),entries.join('\n')+'\n');
  manifest.headers[name]={files:entries.length,sha256:sha(entries.join('\n')+'\n')};
}
function compile(name,group,evaluator=false,coordinate='') {
  const bin=path.join(out,evaluator?'evaluate':`${name}-${group}${coordinate?'-'+coordinate:''}`);
  const args=[...flags, ...(evaluator?[]:[`-DCOMPARISON_GROUP=${group}`]),
    ...(name==='before'&&!evaluator?['-DCOMPARISON_EXACT']:[]),
    ...(coordinate?[`-DCOMPARISON_COORDINATE=${coordinate.startsWith('int')?'std::':''}${coordinate}`,
      `-DCOMPARISON_COORDINATE_NAME="${coordinate}"`]:[]),
    `-I${frozen}/${name}/include`,`-I${boost}`,
    path.join(__dirname,evaluator?sources[1]:sources[0]),'-o',bin];
  const signature=sha(JSON.stringify({args,hashes,headers:manifest.headers[name]}));
  const stamp=bin+'.signature';
  if(!fs.existsSync(stamp)||fs.readFileSync(stamp,'utf8')!==signature) {
    console.log(`Compile ${path.basename(bin)}`);
    const p=bounded('clang++',args,240);
    fs.writeFileSync(bin+'.compile.log',(p.stdout||'')+(p.stderr||''));
    if(p.status!==0) throw Error(`Compile failed: ${bin} (${p.status}/${p.signal})`);
    fs.writeFileSync(stamp,signature);
  }
  return bin;
}
const evaluator=compile('before',0,true);
let records=[];
const units=[[1,'double'],[1,'float'],[1,'int32_t'],[1,'int64_t'],[2,''],[3,''],[4,'']];
for(const name of Object.keys(revisions)) for(const [group,coordinate] of units) {
  const bin=compile(name,group,false,coordinate);
  const check=bounded(bin,['--binary-check'],5);
  if(check.status!==0) throw Error('Binary lifting invariant failed');
  const list=bounded(bin,['--list'],5);
  if(list.status!==0) throw Error('Listing failed');
  for(const row of list.stdout.trim().split('\n')) {
    const [fixture,issue,coordinate]=row.split('\t');
    for(const mode of name==='before'?['native','retained']:['native']) {
      console.log(`Run ${name} ${fixture} ${mode}`);
      const stem=path.join(out,'runs',`${name}-${fixture}-${mode}`);
      const p=bounded(bin,[fixture,mode],name==='before'?90:20);
      fs.writeFileSync(stem+'.stdout',p.stdout||''); fs.writeFileSync(stem+'.stderr',p.stderr||'');
      let r;
      try { r=JSON.parse(p.stdout); }
      catch { r={case:fixture,issue,coordinate,status:p.signal==='SIGALRM'?'timeout':'process_failure'}; }
      r.build=name; r.mode=mode; r.process={status:p.status,signal:p.signal,error:p.error?.message};
      const ev=bounded(evaluator,[],30,JSON.stringify(r));
      fs.writeFileSync(stem+'.evaluation.stderr',ev.stderr||'');
      try { r=JSON.parse(ev.stdout); }
      catch { r.evaluation=ev.signal==='SIGALRM'?'timeout':'process_failure'; }
      fs.writeFileSync(stem+'.json',JSON.stringify(r,null,2)+'\n');
      records.push(r);
      manifest.runs.push({build:name,case:fixture,mode,process:r.process,evaluation:r.evaluation});
      fs.writeFileSync(path.join(out,'manifest.json'),JSON.stringify(manifest,null,2)+'\n');
      fs.writeFileSync(path.join(out,'results.jsonl'),records.map(r=>JSON.stringify(r)).join('\n')+'\n');
    }
  }
}
for(const s of sources) if(sha(fs.readFileSync(path.join(__dirname,s)))!==hashes[s]) throw Error(`Source changed: ${s}`);
manifest.finished=new Date().toISOString();
manifest.records=records.length;
manifest.selected_issues=[...new Set(records.map(r=>Number(r.issue)).filter(n=>n!==494))].sort((a,b)=>b-a);
if(manifest.selected_issues.length!==47) throw Error(`Expected 47 selected issues, got ${manifest.selected_issues.length}`);
fs.writeFileSync(path.join(out,'manifest.json'),JSON.stringify(manifest,null,2)+'\n');
console.log(`Completed ${records.length} records / ${manifest.selected_issues.length} selected issues`);
