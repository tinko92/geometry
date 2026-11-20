"""Extract literal features from the local Geometry Git objects."""
import argparse
import hashlib
import json
from decimal import Decimal as D
from pathlib import Path
import re
import subprocess as sp

ROOT=Path(__file__).resolve().parent
DEV="082b991818062b403a2b2079f21cf69456c15e3b"
CAND="021b24818680ff3af91c0bcd92f5609831ddcab3"
PAIR=re.compile(r"([-+]?\d+(?:\.\d+)?)\s+([-+]?\d+(?:\.\d+)?)")
def sha(data):return hashlib.sha256(data).hexdigest()
def bounds(s):
    points=[(D(x),D(y)) for x,y in PAIR.findall(s)]
    x,y=zip(*points);return (min(x),min(y),max(x),max(y)),len(points)
def window(box):
    l,b,r,t=box
    return f"POLYGON(({l} {b},{l} {t},{r} {t},{r} {b},{l} {b}))"
def component(s):
    # Copy the first polygon's complete parenthesized text, including any holes.
    assert s.startswith("MULTIPOLYGON(")
    s=s[len("MULTIPOLYGON("):];level=0
    for i,c in enumerate(s):
        level+=(c=="(")-(c==")")
        if not level:return "POLYGON"+s[:i+1]
    raise ValueError("unterminated polygon")
def main(repo, out, revision):
    global ROOT
    ROOT=Path(out); ROOT.mkdir(exist_ok=False)
    inputs=ROOT/"inputs"; inputs.mkdir()
    hashes={}; records={}
    def load(path):
        if path not in records:
            data=sp.check_output(["git","-C",repo,"show",revision+":"+path])
            hashes[path]=sha(data); records[path]=data.decode().splitlines()
        return records[path]
    cases={}
    def add(name,a,b,kind,ops,source,decimal=False):
        d=inputs/name;d.mkdir()
        for k,s in [("a",a),("b",b)]: (d/(k+".wkt")).write_text(s+"\n")
        cases[name]=dict(kind=kind,ops=ops,provenance=source,decimal=decimal,operands={
            k:dict(sha256=sha((s+"\n").encode()),bbox=[str(n) for n in bounds(s)[0]],vertices=bounds(s)[1])
            for k,s in [("a",a),("b",b)]})
    for country,lineno in [("nl",2),("it",21),("uk",11),("no",36)]:
        path=f"test/algorithms/buffer/data/{country}.wkt";a=load(path)[lineno-1]
        (x,y,X,Y),_=bounds(a)
        box=[round(x+D(".23")*(X-x)),round(y+D(".18")*(Y-y)),round(x+D(".73")*(X-x)),round(y+D(".82")*(Y-y))]
        add(country+"-window",a,window(box),"polygon",["intersection","difference","relate_aa"],
            dict(a_file=path,one_based_line=lineno,b="rectangle: fractional bbox (.23,.18) to (.73,.82), ties-even integer bounds"),
            country in ("nl","it"))
    world="example/data/world.wkt"
    usa=component(load(world)[189]);mex=component(load(world)[133])
    add("border-dissolve",usa,mex,"polygon",["union","sym_difference","relate_aa"],
        dict(a_file=world,a_line=190,a_polygon=1,b_file=world,b_line=134,b_polygon=1,
             note="One-based; first complete mainland polygons copied from separate multipolygon records."),True)
    roads="example/data/roads.wkt"; lines=load(roads)
    for name,coords in [("columbia-tracks","-122.60 45.42 -121.00 45.85"),("socal-tracks","-118.70 33.75 -118.00 34.35"),
                        ("florida-coast","-82.90 24.40 -79.90 27.20")]:
        box=list(map(D,coords.split()));l,b,r,t=box;ids=[]
        for i,s in enumerate(lines,1):
            (x,y,X,Y),_=bounds(s)
            if X>=l and x<=r and Y>=b and y<=t:ids.append(i)
        assert all(lines[i-1].startswith("LINESTRING(") for i in ids)
        tracks="MULTILINESTRING("+",".join(lines[i-1][len("LINESTRING"):] for i in ids)+")"
        mask=usa if name=="florida-coast" else window(box)
        add(name,tracks,mask,"track",["clip","relate_la"],dict(a_file=roads,one_based_lines=ids,
            selection_bbox=coords,selection="all records whose envelope intersects the selection bbox; no pre-clipping",
            b="world.wkt line 190 polygon 1" if name=="florida-coast" else "selection bbox as polygon"))
    (ROOT/"corpus.json").write_text(json.dumps(dict(cases=cases,source_file_sha256=hashes),indent=2)+"\n")
if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("repo");p.add_argument("out");p.add_argument("--revision",default=CAND)
    a=p.parse_args();main(a.repo,a.out,a.revision)
