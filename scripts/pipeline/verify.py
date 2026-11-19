import json,sys,statistics,bisect
from collections import defaultdict
c=json.load(open('/mnt/user-data/uploads/data/WheatNHR_liftover/output/_ref_cache.json'))
ref=c['ref'];v21=c['v21']
by=defaultdict(list)
for g,(ch,p) in ref.items():
    t=v21.get(g)
    if t: by[ch].append((p,t[0],t[1]))
for ch in by: by[ch].sort()
keys={ch:[x[0] for x in by[ch]] for ch in by}
def local(ch,x,k=6):
    """v2.1 pos of x from its k nearest same-chr anchors (local slope, robust)"""
    L=by[ch];i=bisect.bisect_left(keys[ch],x)
    cand=L[max(0,i-k):i+k]
    cand=[c for c in cand if c[1]==ch]
    if len(cand)<3: return None,None
    cand.sort(key=lambda c:abs(c[0]-x)); near=cand[:k]
    med=statistics.median(t for _,_,t in near)
    near=[c for c in near if abs(c[2]-med)<2]
    if len(near)<3: return None,None
    near.sort()
    inc=sum(1 for a,b in zip(near,near[1:]) if b[2]>=a[2])
    sl=1 if inc>=(len(near)-1)/2 else -1
    r,_,t=min(near,key=lambda c:abs(c[0]-x))
    return t+sl*(x-r),abs(r-x)
def span(ch,lo,hi):
    """all anchors inside interval -> v2.1 span (median-clustered)"""
    L=by[ch];a=bisect.bisect_left(keys[ch],lo);b=bisect.bisect_right(keys[ch],hi)
    ins=[t for p,c2,t in L[a:b] if c2==ch]
    return ins
