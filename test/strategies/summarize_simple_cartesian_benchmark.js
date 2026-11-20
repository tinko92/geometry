// Local read-only analysis; never builds or runs benchmarks.
// node summarize_simple_cartesian_benchmark.js RESULTS_DIRECTORY [--json]
// Raw filename suffixes carry repetition IDs; CSV repeat is always zero.
const fs = require('fs');
const path = require('path');
const assert = require('assert');
const categories = ['degenerate','general','decimal'];
const operations = ['union','intersection','difference','sym_difference','relate_aa','relate_la','relate_ll','covered_by','clip'];
const controlFields = ['corpus_digest','geometry_digest','output_components','output_holes','output_points'];
// Retain eager/deferred only for reading archived results.
const variantNames = ['develop','simple','exact','nearest','before','before-eager','before-deferred','eager','deferred'];
const median = values => {
  const sorted = [...values].sort((a,b)=>a-b), middle = Math.floor(sorted.length/2);
  return sorted.length%2 ? sorted[middle] : (sorted[middle-1]+sorted[middle])/2;
};
const gm = values => Math.exp(values.reduce((sum,x)=>sum+Math.log(x),0)/values.length);

function parse(contents, filename) {
  assert(!contents.includes('"'),'Unexpected quoted CSV');
  const lines = contents.trim().split(/\r?\n/), keys = lines.shift().split(',');
  assert(lines.length && keys.includes('corpus_size'),'Incomplete or incompatible CSV: '+filename);
  const repetition = /-(\d+)\.csv$/.exec(path.basename(filename));
  return lines.map(line=>{
    const values = line.split(',');
    assert.equal(values.length,keys.length,'Malformed CSV: '+filename);
    return {...Object.fromEntries(keys.map((key,i)=>[key,values[i]])),
      repetition:repetition?Number(repetition[1]):null};
  });
}

function analyze(rows, repetitions = 3) {
  assert(Number.isInteger(repetitions) && repetitions>0,'Invalid repetition count');
  const variants = [...new Set(rows.map(r=>r.label))];
  assert(variants.includes('develop') && variants.length>1,'Need develop and at least one candidate');
  assert(variants.every(v=>variantNames.includes(v)),'Unexpected variant');
  const groups = new Map();
  let totalCalls = 0n, filenameRepetitions = 0;
  for (const r of rows) {
    assert.equal(r.coordinate,'double','Only double benchmark inputs are supported');
    assert(categories.includes(r.workload) && operations.includes(r.operation),'Unknown group');
    assert.equal(r.corpus_size,'3','Corpus size changed');
    assert(/^\d+$/.test(r.calls),'Invalid call count');
    const calls = BigInt(r.calls);
    assert(calls>0n && calls%3n===0n,'Calls not a positive complete-corpus multiple');
    totalCalls += calls;
    assert(Number.isFinite(Number(r.ns_per_call)) && Number(r.ns_per_call)>0,'Invalid timing');
    if (r.repetition!==null) {
      assert(r.repetition>=0 && r.repetition<repetitions,'Unexpected filename repetition');
      ++filenameRepetitions;
    }
    for (const field of controlFields) assert(/^\d+$/.test(r[field]),'Missing control '+field);
    const key = JSON.stringify([r.label,r.workload,r.operation]);
    if (!groups.has(key)) groups.set(key,[]);
    groups.get(key).push(r);
  }
  for (const variant of variants) for (const category of categories) for (const operation of operations) {
    const key = JSON.stringify([variant,category,operation]), samples = groups.get(key);
    assert(samples && samples.length===repetitions,'Unexpected sample count: '+key);
    if (samples.every(r=>r.repetition!==null))
      assert.equal(new Set(samples.map(r=>r.repetition)).size,repetitions,'Duplicate filename repetition: '+key);
    for (const field of controlFields)
      assert.equal(new Set(samples.map(r=>r[field])).size,1,'Unstable control '+field+': '+key);
  }
  const pairs = variantNames.filter(v=>v!=='develop' && variants.includes(v)).map(v=>[v,'develop']);
  const cases = [], aggregates = [], discrepancies = [];
  for (const [candidate,reference] of pairs) {
    const matched = [];
    for (const category of categories) for (const operation of operations) {
      const a = groups.get(JSON.stringify([candidate,category,operation]));
      const b = groups.get(JSON.stringify([reference,category,operation]));
      const candidateNs = median(a.map(r=>Number(r.ns_per_call)));
      const referenceNs = median(b.map(r=>Number(r.ns_per_call)));
      const issues = controlFields.filter(field=>a[0][field]!==b[0][field]);
      const result = {candidate,reference,category,operation,candidate_ns:candidateNs,
        reference_ns:referenceNs,ratio:candidateNs/referenceNs,issues};
      cases.push(result); matched.push(result);
      if (issues.length) discrepancies.push({candidate,reference,category,operation,
        fields:Object.fromEntries(issues.map(field=>[field,{candidate:a[0][field],reference:b[0][field]}]))});
    }
    const selections = [['overall','all',matched],
      ...categories.map(c=>['category',c,matched.filter(r=>r.category===c)]),
      ...operations.map(o=>['operation',o,matched.filter(r=>r.operation===o)])];
    for (const [group,name,selected] of selections) aggregates.push({candidate,reference,group,name,
      cases:selected.length,ratio:gm(selected.map(r=>r.ratio)),
      geometry_discrepancies:selected.filter(r=>r.issues.includes('geometry_digest')).length,
      structure_discrepancies:selected.filter(r=>r.issues.includes('corpus_digest')).length});
  }
  return {meaning:'Every factor is candidate/develop median wall time; >1 means slower. Equal-weight geometric mean of matched cases. Digest/stat controls are diagnostics, not correctness oracles; differing cases are not excluded.',
    assurances:{variants,groups_per_variant:27,repetitions,rows:rows.length,
      all_calls_positive_and_divisible_by_3:true,total_calls:String(totalCalls),
      filename_repetitions_checked:filenameRepetitions===rows.length,
      csv_repeat_field_ignored:true,within_variant_controls_stable:true},
    aggregates,cases,discrepancies};
}

try {
  const args = process.argv.slice(2);
  if (args.includes('--self-test')) {
    const rows = [];
    for (const label of ['exact','develop']) for (const workload of categories)
    for (const operation of operations) for (let repetition=0; repetition<3; ++repetition)
      rows.push({label,coordinate:'double',workload,operation,repetition,repeat:'0',corpus_size:'3',calls:'9007199254740993',
        ns_per_call:label==='exact'?'20':'10',corpus_digest:'1',geometry_digest:'18446744073709551615',
        output_components:'2',output_holes:'0',output_points:'10'});
    const result = analyze(rows);
    assert(Math.abs(result.aggregates[0].ratio-2)<1e-12);
    assert.equal(result.cases.length,27); assert.equal(result.discrepancies.length,0);
    assert(result.assurances.filename_repetitions_checked);
    assert.throws(()=>analyze(rows.slice(1)));
    assert.throws(()=>analyze(rows.map((r,i)=>i===0?{...r,calls:'4'}:r)));
    assert.equal(analyze(rows.filter(r=>r.repetition===0),1).assurances.rows,54);
    assert.equal(analyze(rows.filter(r=>r.repetition<2),2).assurances.rows,108);
    assert.equal(median([1,3]),2);
    const values = {develop:10,simple:5,exact:20,nearest:24,eager:30,deferred:25,before:22,'before-eager':32,'before-deferred':27};
    const matrix = variantNames.flatMap(label=>rows.filter(r=>r.label==='develop')
      .map(r=>({...r,label,ns_per_call:String(values[label])})));
    const full = analyze(matrix);
    assert.equal(full.assurances.rows,variantNames.length*27*3);
    assert.equal(full.cases.length,(variantNames.length-1)*27);
    assert(full.cases.every(r=>r.reference==='develop'));
    for (const label of variantNames.filter(v=>v!=='develop')) {
      const overall = full.aggregates.find(r=>r.candidate===label && r.group==='overall');
      assert(Math.abs(overall.ratio-values[label]/10)<1e-12);
    }
    assert.throws(()=>analyze(matrix.map((r,i)=>i===0?{...r,coordinate:'int64'}:r)));
    assert.throws(()=>analyze(matrix.map((r,i)=>i===0?{...r,label:'unknown'}:r)));
    const parsed = parse('label,coordinate,workload,operation,repeat,calls,ns_per_call,corpus_size,corpus_digest,geometry_digest,output_components,output_holes,output_points\n'+
      'before-eager,double,general,union,0,3,20,3,1,18446744073709551615,2,0,10\n','before-eager-general-union-2.csv');
    assert.equal(parsed[0].repetition,2);
    assert.equal(parsed[0].geometry_digest,'18446744073709551615');
    console.log('Compact analyzer self-test passed.');
  } else {
    const input = args.find(a=>a!=='--json');
    assert(input,'Usage: node summarize_simple_cartesian_benchmark.js RESULTS_DIRECTORY [--json]');
    const directory = path.resolve(input), raw = path.join(directory,'raw');
    const metadata = path.join(directory,'run-metadata.txt');
    const run = fs.existsSync(metadata) ? fs.readFileSync(metadata,'utf8') : '';
    if (run) assert(/^finished=/m.test(run),'Run not complete');
    const repetitionSetting = /^repetitions=(\d+)$/m.exec(run);
    const files = fs.existsSync(raw)?fs.readdirSync(raw).filter(n=>n.endsWith('.csv')).map(n=>path.join(raw,n))
      : [path.join(directory,'results.csv')];
    const result = analyze(files.flatMap(f=>parse(fs.readFileSync(f,'utf8'),f)),
      repetitionSetting ? Number(repetitionSetting[1]) : 3);
    const declaredVariants = /^variants=(.+)$/m.exec(run);
    if (declaredVariants) assert.deepEqual([...result.assurances.variants].sort(),
      declaredVariants[1].trim().split(/\s+/).sort(),'Run metadata/CSV variant mismatch');
    if (args.includes('--json')) console.log(JSON.stringify(result,null,2));
    else {
      console.log(result.meaning);
      for (const r of result.aggregates) console.log(`${r.candidate}/${r.reference}\t${r.group}\t${r.name}\t${r.cases}\t${r.ratio.toFixed(3)}`);
      console.log('ASSURANCES '+JSON.stringify(result.assurances));
      console.log('DISCREPANCIES '+JSON.stringify(result.discrepancies));
    }
  }
} catch (error) { console.error('Analysis failed: '+error.message); process.exitCode=1; }
