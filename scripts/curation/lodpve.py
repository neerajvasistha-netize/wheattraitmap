import json,re,collections
M=json.load(open('atlas/merged_v13f.json'));d2t=json.load(open('d2t.json'))
srcdoi={k:re.sub(r'^https?://(dx\.)?doi\.org/','',v.get('url','')).lower() for k,v in M['sources'].items()}
NUM=re.compile(r'(?<![\w.])[−-]?\d+(?:\.\d+)?(?![\w.])')
HL=re.compile(r'\bLOD\b')
HP=re.compile(r'\bPVE\b|\bR2\b|\bR²\b|\bR\^2\b|Var(?:iance)?\s*\(%\)|\bPV\b|%\s*Var|Expl\.|explained|\bR 2\b')
cache={}
def lines(f):
    if f not in cache: cache[f]=open(f,errors='ignore').read().split('\n')
    return cache[f]
def tabley(l):
    seg=[x for x in re.split(r'\s{2,}',l.strip()) if x]
    words=re.findall(r'\b[a-z]{4,}\b',l)
    return len(seg)>=4 and len(words)<=6
def colval(row,hdr,pat,lo,hi):
    m=pat.search(hdr)
    if not m: return None
    cx=(m.start()+m.end())/2
    best=None
    for t in NUM.finditer(row):
        v=float(t.group(0).replace('−','-'))
        if not (lo<=v<=hi): continue
        d=abs((t.start()+t.end())/2-cx)
        if d<=6 and (best is None or d<best[0]): best=(d,t.group(0))
    return best[1] if best else None
def filt(lod,pve,h):
    if re.search(r'thresh',h,re.I) and not re.search(r'max',h,re.I): lod=None
    if re.search(r'thresh',h,re.I) and re.search(r'max',h,re.I): lod=None
    lod=lod if lod and '.' in lod and 2<=float(lod)<=60 else None
    pve=pve if pve and 1<=float(pve)<=95 and ('.' in pve or float(pve)>=2) else None
    hasL=re.search(r'\bLOD\b(?!\s*(support|thr|threshold|interval))',h)
    if hasL and HP.search(h) and not (lod and pve and lod!=pve): return None,None
    return lod,pve
out={};stat=collections.Counter()
for r in M['reports']:
    if r['cls'] not in ('QTL','MQTL','DG'): continue
    if str(r.get('pve') or '').strip() and 'LOD' in str(r.get('stat') or ''): continue
    fs=d2t.get(srcdoi.get(r['source'],''),[])
    nm=r['name'].split(' ')[0]
    if len(nm)<5 or not fs: stat['no text/name']+=1;continue
    found=None
    for f in fs:
        L=lines(f)
        for i,l in enumerate(L):
            if nm not in l or not tabley(l) or len(NUM.findall(l))<3: continue
            # header search upwards
            for j in range(i-1,max(-1,i-45),-1):
                h=L[j]
                if (HL.search(h) or HP.search(h)) and tabley(h):
                    lod=colval(l,h,HL,1.5,150);pve=colval(l,h,HP,0.1,100)
                    lod,pve=filt(lod,pve,h)
                    if lod or pve: found=(lod,pve,f,i);break
                    # header may span two lines
                    if j>0:
                        h2=L[j-1];lod=colval(l,h2,HL,1.5,150);pve=colval(l,h2,HP,0.1,100);lod,pve=filt(lod,pve,h2+' '+h)
                        if lod or pve: found=(lod,pve,f,i);break
                    break
            if found: break
        if found: break
    if found:
        out[r['id']]=found[:2];stat['found']+=1
    else: stat['not found']+=1
print(stat)
json.dump(out,open('lodpve_found.json','w'))
import itertools
d={r['id']:r for r in M['reports']}
for k,v in itertools.islice(out.items(),15): print(k,d[k]['name'][:30],v)
