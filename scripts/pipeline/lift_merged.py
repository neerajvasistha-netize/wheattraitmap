import json,glob,os,re,csv,sys
from collections import defaultdict
sys.path.insert(0,'/home/claude/lift')
import liftover_atlas as L
U='/mnt/user-data/uploads/data/WheatNHR_liftover'
c=json.load(open(U+'/output/_ref_cache.json'))
ref={g:tuple(v) for g,v in c['ref'].items()}
v21={g:tuple(v) for g,v in c['v21'].items()}
NAMES={"alchemy":"Alchemy","arinalrfor":"ArinaLrFor","jagger":"Jagger","julius":"Julius","kariega":"Kariega","lancer":"Lancer","landmark":"Landmark","mace":"Mace","mattis":"SY Mattis","norin61":"Norin 61","paragon":"Paragon","renan":"Renan","stanley":"Stanley"}
targets={'RefSeq v2.1':v21}
for f in sorted(glob.glob('/home/claude/pan/anchors_*.tsv')):
    k=re.sub(r'anchors_Triticum_aestivum_','',os.path.basename(f)).split('.')[0]
    d={r['v11_gene']:(r['chr'],float(r['mid_Mb'])) for r in csv.DictReader(open(f),delimiter='\t')}
    if k in NAMES and len(d)>1000: targets[NAMES[k]]=d
print('targets',list(targets))
anc=defaultdict(list)
for g,(ch,p) in ref.items(): anc[ch].append((g,(ch,p)))
rev=defaultdict(list)
for g,(ch,p) in v21.items(): rev[ch].append((g,(ch,p)))
M=json.load(open('/home/claude/atlas/merged_v13f.json'))
num=lambda v: None if v in (None,'') else float(v)
stats=defaultdict(lambda: defaultdict(int))
for r in M['reports']:
    ch=r['chr']; s,e,p=num(r.get('start')),num(r.get('end')),num(r.get('peak'))
    asm={}
    if ch not in L.CHR.__class__.__mro__ and False: pass
    if s is None and p is None or ch not in ('1A 1B 1D 2A 2B 2D 3A 3B 3D 4A 4B 4D 5A 5B 5D 6A 6B 6D 7A 7B 7D'.split()):
        r['asm']={'_note':'no physical position in source'}; continue
    lo,hi=(s,e) if s is not None else (p,p)
    v2=str(r.get('ref','')).startswith('RefSeq v2')
    if v2:
        back=L.lift(ch,lo,hi,rev,ref)
        if 'start' not in back:
            r['asm']={'RefSeq v2.1':{'chr':ch,'start':s,'end':e,'peak':p,'flag':'as_reported' if r['ref']=='RefSeq v2.1' else 'v2.0_as_reported'},'_note':'could not convert to v1.0'}; continue
        lo1,hi1=back['start'],back['end']
        asm['RefSeq v1.0']={'chr':ch,'flag':'converted from '+r['ref']}
        if s is None: asm['RefSeq v1.0']['peak']=lo1
        else: asm['RefSeq v1.0'].update(start=lo1,end=hi1)
        asm['RefSeq v2.1']={'chr':ch,'start':s,'end':e,'peak':p,'flag':'as_reported' if r['ref']=='RefSeq v2.1' else 'v2.0_as_reported'}
        asm['RefSeq v2.1']={k:v for k,v in asm['RefSeq v2.1'].items() if v is not None}
    else:
        lo1,hi1=lo,hi
    for name,tg in targets.items():
        if name=='RefSeq v2.1' and v2: continue
        res=L.lift(ch,lo1,hi1,anc,tg)
        o={'chr':res.get('chr'),'flag':res['flag']}
        if 'start' in res:
            if s is None: o['peak']=res['start']
            else: o['start'],o['end']=res['start'],res['end']
        asm[name]=o; stats[name][res['flag']]+=1
    r['asm']=asm
json.dump(M,open('/home/claude/atlas/merged_v13_lifted.json','w'))
for k,v in stats.items(): print(k,dict(v))
