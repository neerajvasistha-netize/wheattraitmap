import json,re,os
from collections import defaultdict
M=json.load(open('/home/claude/atlas/merged_v13f.json'))
d2t=json.load(open('/home/claude/d2t.json'))
def doi(s):
    m=re.search(r'10\.\d{4,9}/[^\s"\',;]+',s or ''); return m.group(0).lower().rstrip('.') if m else ''
TXT={}
for k,v in M['sources'].items():
    d=doi(v.get('url') or ''); f=d2t.get(d) or d2t.get(d.upper())
    if isinstance(f,list): f=f[0] if f else None
    if f and os.path.exists(f): TXT[k]=f
QR=re.compile(r'\b(Q[A-Za-z][a-zA-Z]{0,7}[-.][A-Za-z0-9]{2,12}[-.]?[1-7][ABD][SL]?(?:\.\d{1,2})?[a-z]?)\b')
GN=re.compile(r'\b((?:Sr|Lr|Yr|Pm|Sb|Snn|Tsn|Bt|Cre|Wsm|Bdv|Ym|Fhb|Qfhb|Tsc|Stb|Rlnn)\d{1,3}(?:-[1-7][ABD])?)\b')
have=defaultdict(set)
for r in M['reports']: have[r['source']].add(str(r.get('name') or '').strip())
rows=[]
for s,f in TXT.items():
    try: t=re.sub(r'[ \t]+',' ',open(f,errors='ignore').read())
    except Exception: continue
    found=set(QR.findall(t))|set(GN.findall(t))
    found={x for x in found if len(x)>4}
    missed=sorted(found-have[s])
    rows.append((len(missed),len(have[s]),s,M['sources'][s]['short'],missed))
rows.sort(reverse=True)
tot=sum(r[0] for r in rows)
print('studies with local text:',len(rows))
print('locus names present in the text but absent from the atlas:',tot)
print('studies with >=5 such names:',sum(1 for r in rows if r[0]>=5))
print()
print('%-5s %-5s %-28s %s'%('miss','have','study','examples'))
for m,h,s,sh,ex in rows[:25]:
    print('%-5d %-5d %-28s %s'%(m,h,sh[:27],', '.join(ex[:5])))
json.dump([{'missed':m,'have':h,'key':s,'short':sh,'names':ex} for m,h,s,sh,ex in rows],
          open('/home/claude/audit_missed.json','w'),indent=1)
