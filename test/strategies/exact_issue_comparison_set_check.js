// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.
// node test/strategies/exact_issue_comparison_set_check.js RUN_DIR [FROZEN_DIR]
// Read-only output comparison. Reuse the frozen rational evaluator; never
// rebuild, canonicalize, or modify raw run outputs or repository reports.
const fs=require('fs'),path=require('path'),cp=require('child_process'),crypto=require('crypto'),assert=require('assert');
assert(process.argv[2],'RUN_DIR is required');
const out=path.resolve(process.argv[2]);
const reference=path.resolve(process.argv[3]||'/tmp/bg-exact-issue-comparison-final');
const sha=x=>crypto.createHash('sha256').update(x).digest('hex');
const referenceText=fs.readFileSync(path.join(reference,'results.jsonl'),'utf8');
const currentText=fs.readFileSync(path.join(out,'results.jsonl'),'utf8');
const old=new Map(referenceText.trim().split('\n').map(JSON.parse)
  .filter(r=>r.build==='before').map(r=>[r.case+'/'+r.mode,r]));
const current=currentText.trim().split('\n').map(JSON.parse);
const manifest=JSON.parse(fs.readFileSync(path.join(out,'manifest.json'),'utf8'));
function differences(a,b,p='') {
  if(a===b)return [];
  if(a&&b&&typeof a==='object'&&typeof b==='object')return [...new Set([...Object.keys(a),...Object.keys(b)])]
    .flatMap(k=>differences(a[k],b[k],p?p+'.'+k:k));
  return [{path:p,before:a,adaptive:b}];
}
// Diagnose rotations/permutations only. Never rewrite, reverse, simplify,
// snap, or otherwise canonicalize the recorded output coordinates.
function polygons(wkt) {
  const tokens=wkt.match(/MULTIPOLYGON|[(),]|[+-]?\d+(?:\/\d+)?/g);
  assert(tokens&&tokens.join('')===wkt.replace(/\s/g,''));
  assert.equal(tokens[0],'MULTIPOLYGON');let i=1;
  function list(){assert.equal(tokens[i++],'(');const a=[];if(tokens[i]===')'){++i;return a;}
    for(;;){if(tokens[i]==='(')a.push(list());else a.push(tokens[i++]+' '+tokens[i++]);
      if(tokens[i]===')'){++i;return a;}assert.equal(tokens[i++],',');}}
  const result=list();assert.equal(i,tokens.length);return result;
}
function rotation(a,b){
  if(a.length!==b.length)return -1;if(!a.length)return 0;
  if(a[0]!==a[a.length-1]||b[0]!==b[b.length-1])return a.every((p,i)=>p===b[i])?0:-1;
  const n=a.length-1;if(!n)return 0;
  for(let offset=0;offset<n;++offset)if(a[0]===b[offset]&&a.slice(0,n).every((p,i)=>p===b[(i+offset)%n]))return offset;
  return -1;
}
function polygonMatch(a,b){
  if(a.length!==b.length)return null;if(!a.length)return {outer_rotation:0,hole_matches:[]};
  const outer=rotation(a[0],b[0]);if(outer<0)return null;
  const used=new Set(),holes=[];
  for(let i=1;i<a.length;++i){let match=-1,offset=-1;
    for(let j=1;j<b.length;++j)if(!used.has(j)&&(offset=rotation(a[i],b[j]))>=0){match=j;break;}
    if(match<0)return null;used.add(match);holes.push({before:i,adaptive:match,rotation:offset});}
  return {outer_rotation:outer,hole_matches:holes};
}
function permutation(a,b){
  if(!a.startsWith('MULTIPOLYGON')||!b.startsWith('MULTIPOLYGON'))return null;
  const x=polygons(a),y=polygons(b);if(x.length!==y.length)return null;
  const used=new Set(),matches=[];
  for(let i=0;i<x.length;++i){let match=-1,detail=null;
    for(let j=0;j<y.length;++j)if(!used.has(j)&&(detail=polygonMatch(x[i],y[j]))){match=j;break;}
    if(match<0)return null;used.add(match);matches.push({before:i,adaptive:match,...detail});}
  return {exact_permutations_only:true,component_order_changed:matches.some(m=>m.before!==m.adaptive),
    rotated_outer_rings:matches.filter(m=>m.outer_rotation!==0).length,
    reordered_holes:matches.reduce((n,m)=>n+m.hole_matches.filter(h=>h.before!==h.adaptive).length,0),
    rotated_holes:matches.reduce((n,m)=>n+m.hole_matches.filter(h=>h.rotation!==0).length,0)};
}
const raw=[],other=[],failures=[],cache=new Map();
fs.mkdirSync(path.join(out,'set-checks'),{recursive:true});
for(const r of current){
  const before=old.get(r.case+'/'+r.mode);assert(before);
  if(r.status!=='ok'||r.evaluation!=='ok'||String(r.process.status)!=='0')failures.push({case:r.case,mode:r.mode,status:r.status,evaluation:r.evaluation});
  const all=differences(before,r).filter(d=>!['build','mode','process'].some(k=>d.path===k||d.path.startsWith(k+'.')));
  for(const d of all){
    const m=/^geometries\.([^.]+)\.wkt$/.exec(d.path);
    if(!m){other.push({case:r.case,issue:r.issue,mode:r.mode,...d});continue;}
    const key=sha(d.before+'\n'+d.adaptive);let check=cache.get(key);
    if(!check){
      const input={issue:0,case:'frozen-versus-adaptive-exact-set',status:'ok',geometries:{result:{wkt:d.before},reverse:{wkt:d.adaptive}}};
      const p=cp.spawnSync('perl',['-e','alarm shift; exec @ARGV or die $!','30',path.join(reference,'evaluate')],
        {input:JSON.stringify(input),encoding:'utf8',timeout:33000,maxBuffer:64*1024*1024});
      const file=path.join(out,'set-checks',key+'.json');
      fs.writeFileSync(file+'.stderr',p.stderr||'');
      try {const e=JSON.parse(p.stdout);fs.writeFileSync(file,JSON.stringify(e,null,2)+'\n');
        check={process_status:p.status,signal:p.signal,evaluation:e.evaluation,
          before_valid:e.geometries?.result?.rational_valid,adaptive_valid:e.geometries?.reverse?.rational_valid,
          exact_rational_set_equal:e.checks?.reverse_set,evidence:file,permutation:permutation(d.before,d.adaptive)};
      }catch(error){check={process_status:p.status,signal:p.signal,error:String(error),evidence:file};}
      cache.set(key,check);
    }
    raw.push({case:r.case,issue:r.issue,mode:r.mode,geometry:m[1],before_wkt_sha256:sha(d.before),adaptive_wkt_sha256:sha(d.adaptive),...check});
  }
}
const unresolved=raw.filter(r=>r.process_status!==0||r.evaluation!=='ok'
  ||r.exact_rational_set_equal!=='true'||r.before_valid!=='true'||r.adaptive_valid!=='true');
const summary={source_stable:manifest.headers_unchanged,source_header_manifest:manifest.header_manifest,
  checker_sha256:sha(fs.readFileSync(__filename)),reference,reference_manifest_sha256:sha(fs.readFileSync(path.join(reference,'manifest.json'))),
  run_manifest_sha256:sha(fs.readFileSync(path.join(out,'manifest.json'))),
  raw_results_sha256:sha(currentText),reference_results_sha256:sha(referenceText),
  records:current.length,groups:new Set(current.map(r=>r.case)).size,issues:new Set(current.map(r=>+r.issue)).size,
  failures,raw_wkt_changes:raw.length,raw_wkt_changed_records:new Set(raw.map(r=>r.case+'/'+r.mode)).size,
  non_wkt_differences:other,unresolved_set_changes:unresolved,
  semantic_diff_fields:other.length+unresolved.length,
  exact_permutation_only_changes:raw.filter(r=>r.permutation?.exact_permutations_only).length,
  raw_wkt_differences:raw,policy:'Raw WKT preserved. Exact rational DE-9IM set equality on valid outputs; independent exact point-sequence matching only diagnoses rotations/component permutations. No forced canonicalization or area-only equivalence.'};
fs.writeFileSync(path.join(out,'semantic-summary.json'),JSON.stringify(summary,null,2)+'\n');
console.log(JSON.stringify({...summary,raw_wkt_differences:raw.map(r=>({case:r.case,mode:r.mode,geometry:r.geometry,equal:r.exact_rational_set_equal,permutation:r.permutation}))},null,2));
if(!manifest.headers_unchanged||current.length!==164||summary.groups!==82||summary.issues!==47
  ||new Set(current.map(r=>r.case+'/'+r.mode)).size!==164
  ||failures.length||other.length||unresolved.length)process.exitCode=1;
