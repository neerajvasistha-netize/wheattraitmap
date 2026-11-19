"""Yang et al. 2022 (STB): positions from Supplementary Table 3 (DataSheet_4.xlsx, IWGSC physical regions and marker sequences)
plus in-silico mapping of the 14 representative DArT-seq tags on RefSeq v1.0 and v2.1 (end-to-end match, up to 3 mismatches)."""
import json,sys,copy
sys.path.insert(0,'.')
from load import load
R,S,base=load()
W='/home/claude/wheatnhr-atlas/web/seed_ad66.json'
D=json.load(open(W));cur={r['id']:r for r in D['reports']}
P=json.load(open('../yang/placements.json'));M=json.load(open('../yang/meta.json'))
TAG=' In-silico mapping of the marker tag and Supplementary Table 3 regions, Nov 2026 (ad66.py).'
U={ # id: (marker, start, end, use_tag_peak, extra)
 'V0670':('pav_1234699',568,593,True,''),
 'V0676':('pav_4991454',0,20,False,' The 32 bp tag itself matches only 1B 94.47 Mb (exact, v1.0), outside the region the paper gives, so no peak is set.'),
 'V0672':('snp_1112131_F_0_20',40,100,True,''),
 'V0673':('pav_1209089',736,800,False,' The presence/absence (silicoDArT) tag has no match in Chinese Spring (RefSeq v1.0 or v2.1), so only the region is given.'),
 'V0674':('pav_4990595',0,100,False,' The tag has no match on 3A; its best match is 7A 214.54 Mb with 2 mismatches, so only the region is given.'),
 'V0679':('snp_5325269_F_0_37',179,478,False,' The tag has no unique 3A match (best 1A 310.87 Mb, 1 mismatch; 3A copies with 2 mismatches), so only the region is given.'),
 'V0675':('snp_4910674_F_0_26',0,50,True,''),
 'V0677':('pav_3022794',716,730,False,' The tag has no match on 4A (best 1B 392.88 Mb, 2 mismatches), so only the region is given.'),
 'V0678':('snp_2262549_F_0_28',570,590,True,''),
 'V0671':('snp_1233403_F_0_47',0,70,False,' The tag has no match on 6A (best 6D 16.84 Mb, 3 mismatches), so only the region is given.'),
 'V0669':('snp_3026774_F_0_34',440,615,True,''),
 'V0680':('snp_2253221_F_0_65',126,131,True,''),
 'R176':('pav_9364734',705,736,True,''),
}
out=[]
for i,(mk,a,b,tp,extra) in U.items():
    r=copy.deepcopy(cur.get(i) or R[i]);r.pop('asm',None)
    m=M[mk];h=P[mk]['v1.0'];h2=P[mk]['v2.1']
    pk=None;tagtxt=''
    if tp and len(h)==1 and a<=h[0]['snp']/1e6<=b:
        x=h[0];pk=round(x['snp']/1e6,6)
        what='SNP base' if m['snp'][1] else 'tag start'
        tagtxt=f" The marker tag aligns uniquely to {x['chr']} {x['start']:,}-{x['end']:,} bp on RefSeq v1.0 ({x['mm']} mismatch{'es' if x['mm']!=1 else ''}; {what} at {x['snp']:,} bp)"+(f" and to {h2[0]['chr']} {h2[0]['snp']:,} bp on RefSeq v2.1" if len(h2)==1 else '')+'.'
    r.update(start=float(a),end=float(b),peak=pk,ref='RefSeq v1.0',grade='B',marker=mk,curation='hand')
    if 'snpseq' not in r: r['snpseq']=[{'m':mk,'s':m['tag']}]
    r['evidence']=f"Supplementary Table 3: {m['qtl'].lstrip('p')} representative marker {mk}, {m['chr_paper']} {m['mb_paper']} Mb (IWGSC)"
    r['note']=f"GWAS of 273 Australian bread wheat accessions for Zymoseptoria tritici resistance (DArTseq SNP and silicoDArT). Supplementary Table 3 gives the physical region {m['chr_paper']} {a}-{b} Mb and the marker sequence.{tagtxt}{extra}{TAG}"
    out.append(r)
json.dump(out,open('yang_upd.json','w'),indent=0)
for r in out: print(r['id'],r['name'],r['start'],r['end'],r['peak'])
