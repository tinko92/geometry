"""Independent regularized polygon-set oracle. Exact rational open slab sections."""
import re
from fractions import Fraction as Q
from functools import lru_cache
from collections import defaultdict
@lru_cache(None)
def segments(w):
    if not w.startswith(('POLYGON','MULTIPOLYGON')):raise ValueError(w[:18])
    seg=[]
    for ring in re.findall(r'\(([^()]*)\)',w):
        if not ring.strip():continue
        pts=[tuple(map(Q,s.split())) for s in ring.split(',')]
        if any(len(p)!=2 for p in pts) or pts[0]!=pts[-1]:raise ValueError('nonclosed ring')
        for (x,y),(u,v) in zip(pts,pts[1:]):
            if u==x:continue
            m=(v-y)/(u-x);seg.append((min(x,u),max(x,u),m,y-m*x))
    return seg
def select(op,b):
    a,c=bool(b&1),bool(b&2)
    return {'intersection':a and c,'union':a or c,'difference':a and not c,'sym_difference':a!=c}[op]
@lru_cache(None)
def compare(op,a,b,c):
    E=[(*e,1<<k) for k,w in enumerate((a,b,c)) for e in segments(w)]
    # Bound work for the report. Large unchecked cases remain explicitly unknown.
    if len(E)>450:return 'budget','more than 450 nonvertical edge instances'
    cuts={x for lo,hi,*_ in E for x in (lo,hi)}
    for i,(l,h,m,k,_) in enumerate(E):
        for L,H,M,K,_ in E[:i]:
            if m==M or min(h,H)<=max(l,L):continue
            x=(K-k)/(m-M)
            if max(l,L)<x<min(h,H): cuts.add(x)
    keys=sorted(cuts)
    for lo,hi in zip(keys,keys[1:]):
        x=(lo+hi)/2; events=defaultdict(int)
        for l,h,m,k,bit in E:
            if l<x<h:events[m*x+k]^=bit
        mask=0;last=None
        for y,event in sorted(events.items()):
            if last is not None and select(op,mask)!=bool(mask&4):
                return 'unequal','witness ('+str(x)+', '+str((y+last)/2)+'); input A/B/output memberships = '+str(tuple(bool(mask&v) for v in (1,2,4)))
            mask^=event;last=y
        assert mask==0
    return 'equal',str(len(keys)-1 if keys else 0)+' slabs; all open vertical intervals agree'
if __name__=='__main__':
    square='POLYGON((0 0,0 2,2 2,2 0,0 0))'; empty='MULTIPOLYGON()'; small='POLYGON((0 0,0 1,1 1,1 0,0 0))'
    assert compare('intersection',small,square,small)[0]=='equal'
    assert compare('difference',square,small,square)[0]=='unequal'
    assert compare('sym_difference',square,square,empty)[0]=='equal'
    print("PASS: exact open-slab controls")
