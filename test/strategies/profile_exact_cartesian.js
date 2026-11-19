// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.
// node profile_exact_cartesian.js collect|summary OUTPUT [BENCHMARK_DIRECTORY]
const fs = require('fs'), path = require('path'), cp = require('child_process');
const crypto = require('crypto'), assert = require('assert');
const out = path.resolve(process.argv[3] || '/tmp/bg-exact-profile-20261113');
const bench = path.resolve(process.argv[4] || '/tmp/bg-operation-lazy-performance');
const defaultCases = [['convex','union'], ['disjoint','union'], ['convex','relate_aa'],
  ['convex','relate_ll'], ['convex','covered_by'], ['convex','clip'],
  ['star','sym_difference'], ['shared','relate_aa'], ['star','relate_aa'], ['convex','relate_la']];
const allowedCases = defaultCases.concat([['shared','union'], ['shared','intersection'],
  ['shared','difference'], ['shared','sym_difference']]);
const cases = process.env.PROFILE_CASES
  ? process.env.PROFILE_CASES.split(',').map(item=>item.split(':')) : defaultCases;
assert(cases.every(c=>c.length===2 && allowedCases.some(d=>c[0]===d[0] && c[1]===d[1])),
  'Unknown profiling workload/operation');
const hash = file => crypto.createHash('sha256').update(fs.readFileSync(file)).digest('hex');
fs.mkdirSync(out, {recursive:true});
function parse(file) {
  const text = fs.readFileSync(file, 'utf8');
  const body = text.split('Call graph:\n')[1]?.split('Total number in stack')[0];
  assert(body, 'Missing call graph: '+file);
  const roots = [], stack = [];
  for (const line of body.split('\n')) {
    const m = line.match(/^([ +!:|]*)(\d+) (.*)$/);
    if (!m) continue;
    const node = {level:m[1].length, samples:+m[2], frame:m[3], children:[]};
    while (stack.length && stack[stack.length-1].level >= node.level) stack.pop();
    (stack.length ? stack[stack.length-1].children : roots).push(node);
    stack.push(node);
  }
  const main = roots.find(n=>n.frame.includes('com.apple.main-thread'));
  assert(main, 'Missing main thread: '+file);
  const rows = [];
  function walk(n, parents) {
    const self = n.samples - n.children.reduce((a,c)=>a+c.samples,0);
    assert(self>=0, 'Overlapping child samples: '+n.frame);
    const frames = parents.concat(n.frame);
    if (self) rows.push({samples:self, frames});
    n.children.forEach(c=>walk(c,frames));
  }
  walk(main,[]);
  assert.equal(rows.reduce((a,r)=>a+r.samples,0),main.samples);
  const calls = rows.filter(r=>r.frames.some(f=>f.startsWith('execute(')));
  const samples = calls.reduce((a,r)=>a+r.samples,0);
  assert(samples>100, 'Insufficient call samples: '+file);
  const groups = {}, stages = {}, inclusive = {}, flat = new Map();
  const add=(map,key,n)=>map[key]=(map[key]||0)+n;
  for (const r of calls) {
    const leaf=r.frames[r.frames.length-1].split('  (in ')[0];
    const all=r.frames.join('\n');
    flat.set(leaf,(flat.get(leaf)||0)+r.samples);
    // Primitive groups are disjoint; stage/inclusive views are separate axes.
    let group = 'other';
    if (/nextafter/.test(leaf)) group='nextafter';
    else if (/fegetround/.test(leaf)) group='rounding_mode';
    else if (/memmove|memcpy|memset/.test(leaf)) group='memory_copy_clear';
    else if (/malloc|free|operator new|operator delete/.test(leaf)) group='allocation';
    else if (/point_detail::compare_approximations|exact_cartesian::compare\(/.test(leaf)) group='coordinate_filter';
    else if (/interval|enclose_corners|arithmetic_detail::(upward|downward)/.test(leaf)) group='interval_arithmetic';
    else if (/native_point_view.*validate/.test(leaf)) group='input_validation';
    else if (/multiprecision|fixed_(multiply|gcd|orientation|integer|coordinate|round)|divide_qr/.test(leaf)) group='fixed_exact';
    else if (/construction_detail::|orientation_detail::orientation/.test(leaf)) group='expansions_rounding';
    else if (/forward_error::|filtered_orientation/.test(leaf)) group='forward_filter';
    else if (/point_cache|operation_detail::|point_precision_detail::/.test(leaf)) group='cache';
    else if (/::exact_point<[^>]+>::exact_point\(|\bmake_point</.test(leaf)) group='working_point_setup';
    else if (/point_detail::compare|compare_points|static_compare_coordinate|compare_along_segment/.test(leaf)) group='point_comparison';
    else if (/native_predicates_detail::orientation/.test(leaf)) group='native_orientation';
    add(groups,group,r.samples);
    let stage='other';
    for(const f of r.frames) {
      if(/detail::get_turns::|get_turn_info|turn_info::/.test(f)) stage='turn_generation';
      if(/detail::sectionalize::/.test(f)) stage='sectionalize';
      if(/enrich_intersection_points|enrich_sort|enrich_assign|gather_cluster/.test(f)) stage='enrichment';
      if(/traversal::|traverse::|traverse\(/.test(f)) stage='traversal';
      if(/ring_topology_detail::|exact_ring_properties|assign_parents/.test(f)) stage='ring_topology';
      if(/append_output|assign_output|assign_round_point/.test(f)) stage='output_conversion';
    }
    add(stages,stage,r.samples);
    for (const [key,re] of Object.entries({
      bgl:/biconnected_components|depth_first_visit|depth_first_search/,
      validation:/native_point_view.*validate/,
      fixed:/multiprecision|fixed_multiply|fixed_orientation/,
      expansions:/construction_detail::|orientation_detail::orientation/,
      output:/append_output|assign_output|assign_round_point/,
      ring_topology:/ring_topology_detail::|assign_parents/,
      cache:/point_cache_detail::|point_precision_detail::|operation_detail::/}))
      if(re.test(all)) add(inclusive,key,r.samples);
  }
  return {samples,total_samples:main.samples,groups,stages,inclusive,
    top:[...flat.entries()].sort((a,b)=>b[1]-a[1]).slice(0,18).map(([frame,samples])=>({frame,samples})),
    profile_sha256:hash(file)};
}
async function collect() {
  const original=JSON.parse(fs.readFileSync(path.join(bench,'benchmark-manifest.json')));
  const manifest={benchmark:bench,benchmark_manifest_sha256:hash(path.join(bench,'benchmark-manifest.json')),
    source_revision:original.optimized_revision,source_headers:original.headers,
    source_header_status:original.optimized_header_status,develop_revision:original.revisions.develop,
    profiler_sha256:hash(__filename),
    duration:3,interval_ms:1,benchmark_seconds:8,started:new Date().toISOString(),runs:[]};
  for(const coordinate of (process.env.PROFILE_COORDINATES || 'double,int64').split(',')) for(const [workload,operation] of cases)
    for(const variant of (process.env.PROFILE_VARIANTS || 'optimized,develop').split(',')) {
      const name=`${variant}-${coordinate}-${workload}-${operation}`;
      const binary=path.join(bench,`${variant}-${coordinate}`), stem=path.join(out,name);
      const args=['--workload',workload,'--operation',operation,'--seconds','8','--repeat','1'];
      const stdout=fs.openSync(stem+'.csv','w'), stderr=fs.openSync(stem+'.stderr','w');
      const child=cp.spawn(binary,args,{stdio:['ignore',stdout,stderr]});
      const ended=new Promise((resolve,reject)=>{child.once('exit',(code,signal)=>resolve({code,signal}));child.once('error',reject);});
      let sampling,exit;
      try {
        await new Promise(resolve=>setTimeout(resolve,250));
        sampling=cp.spawnSync('sample',[String(child.pid),'3','1','-file',stem+'.txt'],
          {encoding:'utf8',timeout:30000});
        fs.writeFileSync(stem+'.sample.log',(sampling.stdout||'')+(sampling.stderr||''));
        assert.equal(sampling.status,0, 'sample failed: '+name);
        exit=await ended;
        assert.equal(exit.code,0,'Benchmark failed: '+name);
      } finally {
        if(child.exitCode===null && child.signalCode===null) { child.kill('SIGTERM'); await ended; }
        fs.closeSync(stdout);fs.closeSync(stderr);
      }
      const analysis=parse(stem+'.txt');
      manifest.runs.push({name,variant,coordinate,workload,operation,pid:child.pid,
        binary_sha256:hash(binary),args,sample_exit:sampling.status,exit,...analysis});
      fs.writeFileSync(path.join(out,'profiles.json'),JSON.stringify(manifest,null,2)+'\n');
      console.log(name+' samples='+analysis.samples);
    }
  manifest.finished=new Date().toISOString();
  fs.writeFileSync(path.join(out,'profiles.json'),JSON.stringify(manifest,null,2)+'\n');
}
function summary() {
  const m=JSON.parse(fs.readFileSync(path.join(out,'profiles.json')));
  for(const r of m.runs) {
    Object.assign(r,parse(path.join(out,r.name+'.txt')));
    const pct=o=>Object.fromEntries(Object.entries(o).map(([k,v])=>[k,+(100*v/r.samples).toFixed(1)]));
    console.log(JSON.stringify({name:r.name,samples:r.samples,groups:pct(r.groups),stages:pct(r.stages),inclusive:pct(r.inclusive)}));
  }
  m.analysis_script_sha256=hash(__filename);
  fs.writeFileSync(path.join(out,'profiles.json'),JSON.stringify(m,null,2)+'\n');
}
function experiments() {
  const variants=(process.env.PROFILE_VARIANTS || 'original,develop,compare,interval,combined,validation').split(',');
  const binary=(v,c)=>['original','develop'].includes(v)
    ? path.join(bench,`${v==='original'?'optimized':v}-${c}`) : path.join(out,`${v}-${c}`);
  const manifest={started:new Date().toISOString(),seconds:0.2,repetitions:3,
    source_revision:'ad5c0f0cb',runner_sha256:hash(__filename),binaries:{}};
  for(const c of ['double','int64']) for(const v of variants) manifest.binaries[v+'-'+c]=hash(binary(v,c));
  const records=[];
  for(let repetition=0;repetition<3;++repetition) for(const c of ['double','int64'])
    for(const [workload,operation] of cases) for(const v of variants.slice(repetition).concat(variants.slice(0,repetition))) {
      const p=cp.spawnSync(binary(v,c),['--label',v,'--workload',workload,'--operation',operation,
        '--seconds','0.2','--repeat','1'],{encoding:'utf8',timeout:30000});
      assert.equal(p.status,0,`Experiment failed: ${v}/${c}/${workload}/${operation}`);
      const [header,line]=p.stdout.trim().split('\n');
      const r=Object.fromEntries(line.split(',').map((x,i)=>[header.split(',')[i],x]));
      records.push({...r,repetition});
    }
  const groups=[];
  for(const c of ['double','int64']) for(const [workload,operation] of cases) {
    let digest;
    for(const v of variants) {
      const a=records.filter(r=>r.label===v&&r.coordinate===c&&r.workload===workload&&r.operation===operation);
      assert.equal(a.length,3);
      assert.equal(new Set(a.map(r=>r.corpus_digest)).size,1);
      if(v==='original') digest=a[0].corpus_digest;
      else if(v!=='develop') assert.equal(a[0].corpus_digest,digest,`Digest mismatch: ${v}/${c}/${workload}/${operation}`);
      const timings=a.map(r=>+r.ns_per_call).sort((a,b)=>a-b);
      groups.push({variant:v,coordinate:c,workload,operation,ns:timings[1],min:timings[0],max:timings[2],digest:a[0].corpus_digest});
    }
  }
  manifest.finished=new Date().toISOString();
  fs.writeFileSync(path.join(out,'experiments.json'),JSON.stringify({manifest,groups,records},null,2)+'\n');
  console.log('Completed '+records.length+' serial experimental timings; exact digests unchanged.');
}
function patches() {
  const repo=path.resolve(__dirname,'../..'), destination=path.join(out,'patches');
  fs.mkdirSync(destination,{recursive:true});
  for(const v of ['instrumented','compare','interval','combined','validation','provenance','endpoints','rounding','all']) {
    const root=v==='instrumented'?path.join(out,v):path.join(out,'experiments',v);
    const files=cp.execFileSync('rg',['--files',root],{encoding:'utf8'}).trim().split('\n').sort();
    let patch='';
    for(const f of files) {
      const relative=path.relative(root,f), original=path.join(repo,'include',relative);
      const frozen=cp.execFileSync('git',['show','ad5c0f0cb:include/'+relative],{cwd:repo});
      assert(frozen.equals(fs.readFileSync(original)),'Production header changed: '+relative);
      const d=cp.spawnSync('diff',['-u','-L','a/'+relative,'-L','b/'+relative,original,f],{encoding:'utf8'});
      assert(d.status===0||d.status===1); patch+=d.stdout;
    }
    fs.writeFileSync(path.join(destination,v+'.patch'),patch);
  }
}
if(process.argv[2]==='collect') collect().catch(e=>{console.error(e);process.exitCode=1;});
else if(process.argv[2]==='summary') summary();
else if(process.argv[2]==='experiments') experiments();
else if(process.argv[2]==='patches') patches();
else throw Error('Expected collect, summary, experiments or patches');
