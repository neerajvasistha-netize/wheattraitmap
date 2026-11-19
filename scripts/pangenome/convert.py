import json,glob,os,re,csv,sys,statistics
from collections import defaultdict,Counter
U='/mnt/user-data/uploads/WheatNHR_liftover'
sys.path.insert(0,'/tmp/claude-0/-home-claude/00583e43-a309-532b-9f1b-7957bcf4ea95/scratchpad/pan/lib')
import liftover_atlas as L
c=json.load(open(U+'/output/_ref_cache.json'))
ref={g:tuple(v) for g,v in c['ref'].items()}; v21={g:tuple(v) for g,v in c['v21'].items()}
NAMES={"alchemy":"Alchemy","arinalrfor":"ArinaLrFor","jagger":"Jagger","julius":"Julius","kariega":"Kariega","lancer":"Lancer","landmark":"Landmark","mace":"Mace","mattis":"SY Mattis","norin61":"Norin 61","paragon":"Paragon","renan":"Renan","stanley":"Stanley"}
T={}
for f in sorted(glob.glob('anchors_full/anchors_*.tsv')):
    k=re.sub(r'anchors_Triticum_aestivum_','',os.path.basename(f)).split('.')[0]
    if k in NAMES: T[NAMES[k]]={r['v11_gene']:(r['chr'],float(r['mid_Mb'])) for r in csv.DictReader(open(f),delimiter='\t')}
print('cultivars',len(T),{k:len(v) for k,v in T.items()})
anc=defaultdict(list);[anc[ch].append((g,(ch,p))) for g,(ch,p) in ref.items()]
rev=defaultdict(list);[rev[ch].append((g,(ch,p))) for g,(ch,p) in v21.items()]
num=lambda v: None if v in (None,'') else float(v)
CH=set('1A 1B 1D 2A 2B 2D 3A 3B 3D 4A 4B 4D 5A 5B 5D 6A 6B 6D 7A 7B 7D'.split())
need=json.load(open('need.json'));out=[];stats=defaultdict(Counter);skip=Counter()
for x in need:
    r=x['rec'];a=dict(r.get('asm') or {});ch=r.get('chr')
    if ch not in CH: skip['chromosome']+=1;continue
    isv2=str(r.get('ref','')).startswith('RefSeq v2')
    s,e,p=num(r.get('start')),num(r.get('end')),num(r.get('peak'))
    if isv2:
        b=a.get('RefSeq v1.0')
        if b and (num(b.get('start')) is not None or num(b.get('peak')) is not None):
            s1,e1,p1=num(b.get('start')),num(b.get('end')),num(b.get('peak'))
        else:
            lo,hi=(s,e) if s is not None and e is not None else (p,p)
            if lo is None: skip['no position']+=1;continue
            back=L.lift(ch,lo,hi,rev,ref)
            if 'start' not in back: skip['v2.1 to v1.0 failed']+=1;continue
            if s is not None and e is not None: s1,e1,p1=back['start'],back['end'],None
            else: s1,e1,p1=None,None,back['start']
            a['RefSeq v1.0']={k:v for k,v in {'chr':back.get('chr',ch),'start':s1,'end':e1,'peak':p1,'flag':'converted from '+r['ref']}.items() if v is not None}
    else:
        s1,e1,p1=s,e,p
        if 'RefSeq v2.1' not in a or (num(a['RefSeq v2.1'].get('start')) is None and num(a['RefSeq v2.1'].get('peak')) is None):
            lo,hi=(s1,e1) if s1 is not None and e1 is not None else (p1,p1)
            if lo is not None:
                res=L.lift(ch,lo,hi,anc,v21);o={'chr':res.get('chr'),'flag':res['flag']}
                if 'start' in res:
                    if s1 is None or e1 is None: o['peak']=res['start']
                    else: o['start'],o['end']=res['start'],res['end']
                a['RefSeq v2.1']=o
    lo,hi=(s1,e1) if s1 is not None and e1 is not None else (p1,p1)
    if lo is None: skip['no v1.0 position']+=1;continue
    for name,tg in T.items():
        if a.get(name) and (num(a[name].get('start')) is not None or num(a[name].get('peak')) is not None): continue
        res=L.lift(ch,lo,hi,anc,tg);o={'chr':res.get('chr'),'flag':res['flag']}
        if 'start' in res:
            if s1 is None or e1 is None:
                o['peak']=res['start']
            else:
                o['start'],o['end']=res['start'],res['end']
                if p1 is not None:
                    pk=L.lift(ch,p1,p1,anc,tg)
                    if 'start' in pk and pk.get('chr')==o['chr']: o['peak']=pk['start']
        a[name]=o;stats[name][res['flag']]+=1
    a.pop('_note',None)
    nr=dict(r);nr['asm']=a;out.append(nr)
json.dump(out,open('converted.json','w'))
print('converted',len(out),'skipped',dict(skip))
tot=Counter();[tot.update(v) for v in stats.values()];print('flags over all cultivars',dict(tot))
