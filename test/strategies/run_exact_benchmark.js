// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.
// node test/strategies/run_exact_benchmark.js prepare|build|run|paired|allocations|summary [directory]
// Set BOOST_ROOT for the umbrella headers; CXX optionally selects the compiler.
// BENCH_DEVELOP_REVISION pins a newer develop baseline without changing exact snapshots.
// BENCH_BORROWED_ARCHIVE selects the borrowed implementation before performance fixes.
// Whole-operation policy comparison: BENCH_VARIANTS=develop,previous,stage1,optimized,nearest
// BENCH_STAGE1_ARCHIVE selects the pre-borrowing tar; BENCH_STAGE1_INCLUDE can
// instead select an already extracted immutable include directory.
// BENCH_REQUIRE_GEOMETRY_MATCH=1 makes cross-variant spatial digest mismatches
// fail summary, after writing diagnostics (quantization boundaries can differ).
// BENCH_PROJECTED=1 BENCH_COORDINATES=double selects the projected/local pair.
const fs = require('fs');
const path = require('path');
const cp = require('child_process');
const crypto = require('crypto');
const os = require('os');
const repo = path.resolve(__dirname, '../..');
const out = path.resolve(process.argv[3] || '/tmp/bg-exact-performance');
const boost = process.env.BOOST_ROOT || path.resolve(repo, '../..');
const cxx = process.env.CXX || 'clang++';
const fractional = process.env.BENCH_FRACTIONAL === '1';
const projected = process.env.BENCH_PROJECTED === '1';
const workloads = projected ? ['projected', 'projected-local'] : ['disjoint','convex','star','holes','shared'];
const timing = {seconds:Number(process.env.BENCH_SECONDS || 0.2),
  repetitions:Number(process.env.BENCH_REPETITIONS || 3)};
if (!(timing.seconds > 0) || !Number.isInteger(timing.repetitions) || timing.repetitions < 1)
  throw Error('Invalid benchmark duration or repetition count');
const source = path.join(__dirname, 'exact_cartesian_benchmark.cpp');
const flags = ['-std=c++17', '-O3', '-DNDEBUG', '-gline-tables-only',
  '-fno-omit-frame-pointer', '-ffp-contract=off', '-pthread',
  '-DBOOST_GEOMETRY_DISABLE_DEPRECATED_14_WARNING'];
const revisions = {
  before: 'bc65c6d60fe2d891ed978d0c8cdf91bb678712b9',
  previous: '2886f2acd6228f6d739ac331b4fac1492809c625',
  pr1491: '72b2d832fb9596a2a23c0c164a5426ae720d0cd6',
  develop: process.env.BENCH_DEVELOP_REVISION || 'aaf0560e520041a833c11ad5d61017b684b46f9a'
};
fs.mkdirSync(out, {recursive: true});
function exec(command, args, options = {}) {
  return cp.execFileSync(command, args, {cwd: repo, stdio: 'inherit', ...options});
}
function hashTree(root) {
  const hash = crypto.createHash('sha256');
  function visit(directory) {
    for (const entry of fs.readdirSync(directory, {withFileTypes: true}).sort((a,b) => a.name.localeCompare(b.name))) {
      const file = path.join(directory, entry.name);
      if (entry.isDirectory()) { visit(file); }
      else { hash.update(path.relative(root, file)); hash.update(fs.readFileSync(file)); }
    }
  }
  visit(root);
  return hash.digest('hex');
}
const variants = [];
const selectedNames = (process.env.BENCH_VARIANTS || 'develop,previous,optimized').split(',')
  .map(s=>s.trim()).map(n=>n === 'eager' ? 'nearest' : n);
const variantNames = ['develop', 'pr1491', 'before', 'previous', 'stage1', 'borrowed', 'optimized', 'nearest', 'no-algopt'];
if (selectedNames.some(n => ![...variantNames, 'develop-robust'].includes(n))) throw Error('Unsupported benchmark variant');
const stage1Archive = path.resolve(process.env.BENCH_STAGE1_ARCHIVE || '/tmp/bg-exact-before-borrowing-20261116.tar');
const stage1Include = path.resolve(process.env.BENCH_STAGE1_INCLUDE || path.join(out, 'stage1/include'));
const borrowedArchive = path.resolve(process.env.BENCH_BORROWED_ARCHIVE || '/tmp/bg-exact-borrowed-before-performance-20261116.tar.gz');
const borrowedInclude = path.join(out, 'borrowed/include');
const currentNames = ['optimized', 'nearest', 'no-algopt'];
const exactNames = ['before', 'previous', 'stage1', 'borrowed', ...currentNames];
const allocationVariants = new Set(['develop', 'previous', 'stage1', 'borrowed', 'optimized', 'nearest']
  .map(n=>n+'-double'));
const coordinates = (process.env.BENCH_COORDINATES || 'double,int64').split(',');
if (coordinates.some(c => !['double', 'int64'].includes(c))) throw Error('Unsupported benchmark coordinate');
if (projected && (fractional || coordinates.some(c => c !== 'double')))
  throw Error('BENCH_PROJECTED requires BENCH_COORDINATES=double without BENCH_FRACTIONAL');
for (const coordinate of coordinates) {
  const type = coordinate === 'int64' ? ['-DBG_BENCH_INT64']
    : fractional ? ['-DBG_BENCH_FRACTIONAL_INPUTS'] : [];
  for (const name of variantNames) {
    if (!selectedNames.includes(name)) continue;
    variants.push({name: `${name}-${coordinate}`, headers: currentNames.includes(name)
      ? path.resolve(process.env.BENCH_INCLUDE || path.join(repo, 'include'))
      : name === 'stage1' ? stage1Include : name === 'borrowed' ? borrowedInclude : path.join(out, name, 'include'),
      defines: [...type, ...(exactNames.includes(name) ? ['-DBG_BENCH_EXACT'] : []),
        ...(name === 'nearest' ? ['-DBG_BENCH_NEAREST'] : []),
        ...(name === 'no-algopt' ? ['-DBOOST_GEOMETRY_EXACT_DISABLE_ALGORITHM_OPTIMIZATIONS'] : [])]});
  }
}
if (selectedNames.includes('develop-robust')) variants.push({name:'develop-robust-double', headers:path.join(out,'develop/include'),
  defines:['-DBOOST_GEOMETRY_DEFAULT_STRATEGY_SIDE_USE_SIDE_ROBUST']});
function prepare() {
  let borrowed;
  if (selectedNames.includes('borrowed')) {
    const sha256 = crypto.createHash('sha256').update(fs.readFileSync(borrowedArchive)).digest('hex');
    const directory = path.dirname(borrowedInclude), stamp = path.join(directory, 'snapshot.json');
    fs.mkdirSync(directory, {recursive:true});
    if (fs.existsSync(stamp)) {
      if (JSON.parse(fs.readFileSync(stamp,'utf8')).sha256 !== sha256)
        throw Error('Borrowed archive changed; use a fresh output directory');
    } else {
      exec('tar', ['-xf', borrowedArchive, '-C', directory, 'include']);
      fs.writeFileSync(stamp, JSON.stringify({archive:borrowedArchive,sha256},null,2)+'\n');
    }
    borrowed = {archive:borrowedArchive,sha256,headers:borrowedInclude};
  }
  for (const [name, revision] of Object.entries(revisions)) {
    if (!selectedNames.includes(name) && !(name === 'develop' && selectedNames.includes('develop-robust'))) continue;
    const directory = path.join(out, name);
    fs.mkdirSync(directory, {recursive:true});
    const archive = path.join(out, name + '.tar');
    exec('git', ['archive', '--format=tar', `--output=${archive}`, revision, 'include']);
    exec('tar', ['-xf', archive, '-C', directory]);
    fs.unlinkSync(archive);
  }
  let stage1;
  if (selectedNames.includes('stage1')) {
    stage1 = {headers:stage1Include};
    if (!process.env.BENCH_STAGE1_INCLUDE) {
      const sha256 = crypto.createHash('sha256').update(fs.readFileSync(stage1Archive)).digest('hex');
      const directory = path.dirname(stage1Include), stamp = path.join(directory, 'snapshot.json');
      if (fs.existsSync(stamp) && JSON.parse(fs.readFileSync(stamp,'utf8')).sha256 !== sha256)
        throw Error('Stage-1 archive changed; use a fresh output directory');
      fs.mkdirSync(directory,{recursive:true});
      if (!fs.existsSync(stamp)) {
        exec('tar', ['-xf', stage1Archive, '-C', directory, 'include', 'test']);
        fs.writeFileSync(stamp, JSON.stringify({archive:stage1Archive,sha256},null,2)+'\n');
      }
      if (!fs.existsSync(path.join(directory,'test')))
        exec('tar', ['-xf', stage1Archive, '-C', directory, 'test']);
      stage1.archive = stage1Archive; stage1.sha256 = sha256;
      stage1.tests = path.join(directory,'test');
    }
    if (!fs.statSync(stage1Include).isDirectory()) throw Error('Missing Stage-1 include directory');
  }
  fs.writeFileSync(path.join(out, 'benchmark-preparation.json'), JSON.stringify({revisions, flags, boost,
    timing, fractional, projected, workloads, stage1, borrowed, variants:variants.map(v=>v.name), variant_config:variants,
    source_sha256:crypto.createHash('sha256').update(fs.readFileSync(source)).digest('hex')}, null, 2)+'\n');
}
function build() {
  prepare();
  const manifest = {revisions, flags, boost, timing, fractional, projected, workloads,
    digest_contract:'corpus_digest retains component/hole/point counts and relation matrices; geometry_digest adds native exported coordinates quantized to 2^-16 units, sorted undirected edges and hole/component hierarchy. Relation/covered_by use their observed results. Spatial regression signal, not an exact equivalence oracle.',
    output_contract:'Only native output containers survive public calls; no benchmark-owned implicit Root records or precision scopes.',
    variants, stage1:JSON.parse(fs.readFileSync(path.join(out,'benchmark-preparation.json'),'utf8')).stage1,
    borrowed:JSON.parse(fs.readFileSync(path.join(out,'benchmark-preparation.json'),'utf8')).borrowed,
    compiler:exec(cxx,['--version'],{encoding:'utf8',stdio:'pipe'}),
    optimized_revision:exec('git',['rev-parse','HEAD'],{encoding:'utf8',stdio:'pipe'}).trim(),
    optimized_header_status:exec('git',['status','--porcelain','--','include'],{encoding:'utf8',stdio:'pipe'}),
    host:{platform:os.platform(), release:os.release(), cpu:os.cpus()[0]?.model, cores:os.cpus().length},
    source_sha256:crypto.createHash('sha256').update(fs.readFileSync(source)).digest('hex'), headers:{}};
  manifest.runner_sha256 = crypto.createHash('sha256').update(fs.readFileSync(__filename)).digest('hex');
  for (const variant of variants) {
    const digest = hashTree(variant.headers);
    manifest.headers[variant.name] = digest;
    for (const allocations of [false,true]) {
      if (allocations && !allocationVariants.has(variant.name)) { continue; }
      const binary = path.join(out, variant.name + (allocations ? '-alloc' : ''));
      const args = [...flags,...variant.defines,...(allocations?['-DBG_BENCH_TRACK_ALLOCATIONS']:[]),
        '-I'+variant.headers, '-I'+boost, source, '-o', binary];
      const signature = crypto.createHash('sha256').update(JSON.stringify({digest,args,source:manifest.source_sha256,
        compiler:manifest.compiler})).digest('hex');
      if (fs.existsSync(binary+'.signature') && fs.readFileSync(binary+'.signature','utf8') === signature) { continue; }
      console.log('Building '+path.basename(binary));
      exec(cxx,args);
      fs.writeFileSync(binary+'.signature', signature);
    }
  }
  fs.writeFileSync(path.join(out,'benchmark-manifest.json'), JSON.stringify(manifest,null,2)+'\n');
}
function run(allocations, paired = false) {
  const manifest=JSON.parse(fs.readFileSync(path.join(out,'benchmark-manifest.json'),'utf8'));
  if (JSON.stringify(timing)!==JSON.stringify(manifest.timing)) throw Error('Timing settings changed; rebuild manifest first');
  if (fractional !== manifest.fractional) throw Error('Input control changed; rebuild manifest first');
  if (projected !== manifest.projected) throw Error('Workload suite changed; rebuild in a separate output directory');
  const sourceHash=crypto.createHash('sha256').update(fs.readFileSync(source)).digest('hex');
  if(sourceHash!==manifest.source_sha256) { throw Error('Benchmark source changed; rebuild first'); }
  for(const variant of variants) {
    if(hashTree(variant.headers)!==manifest.headers[variant.name]) { throw Error('Headers changed; rebuild '+variant.name); }
  }
  const selected = allocations ? variants.filter(v => allocationVariants.has(v.name)) : variants;
  if (selected.length === 0) throw Error('No selected benchmark binaries');
  const cases = paired ? workloads.flatMap(workload =>
    ['union','intersection','difference','sym_difference','relate_aa','relate_la','relate_ll','covered_by','clip']
      .map(operation => ({workload,operation}))) : [{}];
  const runInfo = {mode:allocations?'allocations':paired?'case-interleaved':'batch',timing,
    workloads, started:new Date().toISOString(),source_sha256:sourceHash};
  // Never time concurrent workloads. Paired runs keep variants close in time
  // for each case, avoiding whole-binary timing drift on a busy desktop.
  for (let repetition=0; repetition<(allocations?1:timing.repetitions); ++repetition) {
    for (const {workload,operation} of cases) {
    for (const variant of selected.slice(repetition%selected.length).concat(selected.slice(0,repetition%selected.length))) {
      const name=variant.name+(allocations?'-alloc':'');
      const file=path.join(out,`${name}${paired?'-'+workload+'-'+operation:''}-${repetition}.csv`);
      const descriptor=fs.openSync(file,'w');
      try {
        exec(path.join(out,name),['--label',name,'--seconds',allocations?'0.06':String(timing.seconds),'--repeat','1',
          ...(paired?['--workload',workload,'--operation',operation]:projected?['--workload','projected-pair']:[])],
          {stdio:['ignore',descriptor,'inherit']});
      } finally { fs.closeSync(descriptor); }
    }
    }
  }
  for(const variant of selected) {
    if(hashTree(variant.headers)!==manifest.headers[variant.name]) { throw Error('Headers changed during timing'); }
  }
  runInfo.finished=new Date().toISOString();
  fs.writeFileSync(path.join(out,`benchmark-${runInfo.mode}-run.json`),JSON.stringify(runInfo,null,2)+'\n');
}
function summary() {
  const groups = new Map();
  for (const name of fs.readdirSync(out).filter(n=>/^(develop|pr1491|before|previous|stage1|borrowed|optimized|nearest|eager|deferred|no-algopt)-.*-\d+\.csv$/.test(n))) {
    const lines = fs.readFileSync(path.join(out,name),'utf8').trim().split('\n');
    const keys=lines.shift().split(',');
    for(const line of lines) {
      const record=Object.fromEntries(line.split(',').map((v,i)=>[keys[i],v]));
      const key=[record.label,record.coordinate,record.workload,record.operation].join(',');
      if(!groups.has(key)) { groups.set(key,[]); }
      groups.get(key).push(record);
    }
  }
  const rows=[];
  for(const [key,records] of groups) {
    records.sort((a,b)=>Number(a.ns_per_call)-Number(b.ns_per_call));
    const r=records[Math.floor(records.length/2)];
    if(new Set(records.map(r=>r.corpus_digest)).size!==1) { throw Error('Unstable digest: '+key); }
    if (records.some(r=>!/^\d+$/.test(r.geometry_digest || ''))) throw Error('Missing geometric digest; rebuild and rerun '+key);
    if (new Set(records.map(r=>r.geometry_digest)).size!==1) throw Error('Unstable geometric digest: '+key);
    rows.push({label:r.label,coordinate:r.coordinate,workload:r.workload,operation:r.operation,
      samples:records.length,ns_per_call:Number(r.ns_per_call),
      cpu_ns_per_call:Number(r.cpu_seconds)*1e9/Number(r.calls),
      min_ns_per_call:Number(records[0].ns_per_call),max_ns_per_call:Number(records[records.length-1].ns_per_call),
      allocations_per_call:Number(r.allocations)/Number(r.calls),bytes_per_call:Number(r.allocated_bytes)/Number(r.calls),
      corpus_digest:r.corpus_digest,geometry_digest:r.geometry_digest});
  }
  for(const r of rows.filter(r=>r.label==='optimized-'+r.coordinate)) {
    for(const reference of ['before','previous','stage1','borrowed','nearest','eager','deferred','no-algopt']) {
      if (!selectedNames.includes(reference)) continue;
      const other=rows.find(s=>s.label===reference+'-'+r.coordinate && s.workload===r.workload && s.operation===r.operation);
      if(!other || other.corpus_digest!==r.corpus_digest) { throw Error('Exact digest changed: '+JSON.stringify(r)); }
    }
  }
  if (rows.length === 0) throw Error('No benchmark CSV rows');
  for (const r of rows) {
    const allocation = r.label.endsWith('-alloc') ? '-alloc' : '';
    const sameCase = s=>s.coordinate===r.coordinate && s.workload===r.workload && s.operation===r.operation
      && s.label.endsWith('-alloc')===r.label.endsWith('-alloc');
    const reference = ['optimized','stage1','previous'].map(n=>rows.find(s=>sameCase(s)
      && s.label===n+'-'+r.coordinate+allocation)).find(Boolean) || rows.find(sameCase);
    r.geometry_reference = reference.label;
    r.geometry_matches_reference = r.geometry_digest === reference.geometry_digest;
  }
  fs.writeFileSync(path.join(out,'benchmark-summary.json'),JSON.stringify(rows,null,2)+'\n');
  const columns=Object.keys(rows[0]);
  fs.writeFileSync(path.join(out,'benchmark-summary.csv'),[columns.join(','),
    ...rows.map(r=>columns.map(k=>r[k]).join(','))].join('\n')+'\n');
  console.log('Wrote '+rows.length+' case summaries to '+path.join(out,'benchmark-summary.json'));
  const mismatches = rows.filter(r=>!r.geometry_matches_reference);
  if (mismatches.length) {
    console.warn(mismatches.length+' geometric digest mismatches; inspect quantization-boundary/output-policy differences');
    if (process.env.BENCH_REQUIRE_GEOMETRY_MATCH === '1') throw Error('Geometric benchmark digest changed');
  }
}
const action=process.argv[2];
if(action==='prepare') { prepare(); }
else if(action==='build') { build(); }
else if(action==='run') { run(false); }
else if(action==='paired') { run(false,true); }
else if(action==='allocations') { run(true); }
else if(action==='summary') { summary(); }
else { throw Error('Expected prepare, build, run, paired, allocations, or summary'); }
