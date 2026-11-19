import json,glob
from collections import defaultdict
meta=json.load(open('meta.json'))
NAMES={'cs':'RefSeq v1.0','refseqv2':'RefSeq v2.1','alchemy':'Alchemy','arinalrfor':'ArinaLrFor','jagger':'Jagger','julius':'Julius','kariega':'Kariega','lancer':'Lancer','landmark':'Landmark','mace':'Mace','mattis':'SY Mattis','norin61':'Norin 61','paragon':'Paragon','renan':'Renan','stanley':'Stanley'}
res={}
for f in sorted(glob.glob('hits_*.json')):
    g=f[5:-5];H=json.load(open(f));d=defaultdict(list)
    for x in H:d[x['id']].append(x)
    for rid,m in meta.items():
        L=len(m['trimmed']);sp=int(m['snp'][0])
        hs=sorted(d[rid],key=lambda x:x['mm'])
        if not hs: continue
        best=hs[0]['mm'];top=[x for x in hs if x['mm']==best]
        out=[]
        for x in top:
            chrom=x['chr'].replace('chr','')
            chrom=chrom[:2] if chrom[:1].isdigit() else chrom
            snp=x['pos']+sp if x['strand']=='+' else x['pos']+(L-1-sp)
            out.append({'chr':chrom,'raw':x['chr'],'start':x['pos'],'end':x['pos']+L-1,'snp':snp,'strand':x['strand'],'mm':best})
        res.setdefault(rid,{})[NAMES[g]]=out
json.dump(res,open('placements.json','w'),indent=1)
for rid in meta:
    print(rid, meta[rid]['marker'])
    for g,o in res.get(rid,{}).items():
        print('   %-12s'%g, '; '.join(f"{x['chr']}:{x['snp']:,} {x['strand']} mm{x['mm']}" for x in o[:4]), '(+%d more)'%(len(o)-4) if len(o)>4 else '')
