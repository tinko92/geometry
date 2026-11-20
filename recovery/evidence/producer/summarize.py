import csv
import math
import statistics as s
import sys
from pathlib import Path

fields = "input_fingerprint corpus_digest geometry_digest output_components output_holes output_points".split()
ws = "degenerate general decimal projected projected-local".split()
ops = "union intersection difference sym_difference relate_aa relate_la relate_ll covered_by clip".split()
cells = [(w,o) for w in ws for o in ops]
def read(p):
    with p.open() as f:
        return list(csv.DictReader(f))

def summarize(name):
    root=Path("timing")/name
    samples=read(root/"samples.csv")
    assert len(samples)==2250
    times={}
    refs={}
    print(name, "cpu/wall",sum(float(r["cpu_seconds"]) for r in samples)/sum(float(r["seconds"]) for r in samples))
    for n in (1,2,4,8,100):
        for v in ("before","after"):
            label=f"{v}-{n}"
            rr=read(root/(label+".check"))
            assert len(rr)==45
            for r in rr:
                cell=r["workload"],r["operation"]
                assert r["label"]==label and cell in cells
                t=tuple(r[x] for x in fields)
                assert refs.setdefault((n,cell),t)==t, (name,label,cell)
    for r in samples:
        v,n=r["label"].split("-"); n=int(n)
        i=int(r["repeat"]); cell=r["workload"],r["operation"]
        assert tuple(r[f] for f in fields)==refs[n,cell]
        assert int(r["calls"])>0 and int(r["calls"])%3==0 and r["corpus_size"]=="3" and r["coordinate"]=="double"
        k=(n,cell,i,v); val=float(r["ns_per_call"])
        assert k not in times and math.isfinite(val) and val>0 and i in range(5);times[k]=val
    with (root/"paired.csv").open("w") as out:
        writer=csv.writer(out)
        writer.writerow(["N","workload","operation","before_ns","after_ns","after_over_before"])
        for n in (1,2,4,8,100):
            ratios=[]
            for cell in cells:
                a,b = (s.median(times[n,cell,i,v] for i in range(5)) for v in ("before","after"))
                ratios.append(b/a);writer.writerow([n,*cell,a,b,b/a])
            reps = [s.geometric_mean(times[n,c,i,"after"]/times[n,c,i,"before"] for c in cells) for i in range(5)]
            print(n,round(s.geometric_mean(ratios),5), "reps",*(round(r,4) for r in reps))

for name in sys.argv[1:] or ["gcc","clang"]:
    summarize(name)
