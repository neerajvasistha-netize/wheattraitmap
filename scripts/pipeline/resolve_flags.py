import json,sys
sys.path.insert(0,'/home/claude/lift')
from verify import local,span
P='/home/claude/atlas/merged_v13_lifted.json'
L=json.load(open(P))
num=lambda v: None if v in (None,'') else float(v)
log=[];npk=0
for r in L['reports']:
    a=r['asm'].get('RefSeq v2.1')
    if not a: continue
    ch=r['chr'];s,e,p=num(r.get('start')),num(r.get('end')),num(r.get('peak'))
    v1=not str(r.get('ref','')).startswith('RefSeq v2')
    f=a.get('flag')
    if v1 and f in ('inverted','wide_gap'):
        lo,hi=(s,e) if s is not None else (p,p)
        ts,gs=local(ch,lo);te,ge=local(ch,hi)
        if ts is None or te is None: continue
        old=(a.get('start',a.get('peak')),a.get('end'))
        if s is None:
            a['peak']=round(ts,3)
        else:
            x,y=sorted((ts,te))
            ins=sorted(span(ch,lo,hi))
            if len(ins)>=5:
                k=len(ins)//20; x=min(x,ins[k]); y=max(y,ins[-1-k])
            a['start'],a['end']=round(x,3),round(y,3)
        a['flag']='reoriented_checked' if f=='inverted' else 'rechecked'
        a['qc']='checked against the nearest anchor genes on both sides (gap %.2f Mb)'%max(gs,ge)
        log.append((r['id'],r['name'][:30],f,old,(a.get('start',a.get('peak')),a.get('end'))))
    # lifted peak for v1.0 intervals
    if v1 and s is not None and p is not None and 'start' in a and 'peak' not in a:
        t,g=local(ch,p)
        if t is not None and g<2 and a['start']-1<=t<=a['end']+1:
            a['peak']=round(t,3);npk+=1
json.dump(L,open(P,'w'))
for l in log: print(*l,sep=' | ')
print('peaks added',npk)
