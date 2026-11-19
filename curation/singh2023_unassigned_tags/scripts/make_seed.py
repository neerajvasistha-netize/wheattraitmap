"""Build web/seed_ad65.json: chromosome and physical position for the 11 unassigned Singh et al. 2023 MTAs,
from whole-genome alignment of their trimmed DArT tags (placements.json from analyse.py)."""
import json,re,collections,csv
W='/home/claude/wheatnhr-atlas/web/';C='/home/claude/wheatnhr-atlas/curation/singh2023_unassigned_tags/'
meta=json.load(open('meta.json'));P=json.load(open('placements.json'))
src=json.load(open(W+'seed_ad59.json'));recs={r['id']:r for r in src['reports']}
ORDER=['RefSeq v1.0','RefSeq v2.1','Alchemy','ArinaLrFor','Jagger','Julius','Kariega','Lancer','Landmark','Mace','Norin 61','Paragon','Renan','SY Mattis','Stanley']
out=[];rows=[]
for rid,m in meta.items():
    r=json.loads(json.dumps(recs[rid]));pl=P.get(rid,{})
    chrs=collections.Counter(h[0]['chr'] for g,h in pl.items() if len(h)==1 and h[0]['mm']<=2 and h[0]['chr']!='Un')
    cons=chrs.most_common(1)[0][0]
    L=len(m['trimmed']);asm={}
    for g in ORDER:
        h=pl.get(g)
        if not h: e={'chr':cons,'flag':'unmapped','qc':'tag not found (more than 3 mismatches)'}
        else:
            on=[x for x in h if x['chr']==cons]; extra=''
            if len(h)>1:
                others=[x for x in h if x not in on[:1]]
                if len(on)==1 or (len(on)>1 and max(x['snp'] for x in on)-min(x['snp'] for x in on)<1_000_000):
                    extra='; identical-quality copy also at '+', '.join(f"{x['chr']} {x['snp']:,} bp" for x in others)
                    h=on[:1]
                else:
                    h=None
            if h is None:
                e={'chr':cons,'flag':'unmapped','qc':'tag matches several sites equally well ('+', '.join(f"{x['chr']} {x['snp']:,} bp" for x in pl[g])+f', {pl[g][0]["mm"]} mismatches); not placed'}
            else:
                x=h[0];desc=f'DArT tag {x["start"]:,}-{x["end"]:,} bp ({x["strand"]} strand), {x["mm"]} mismatch{"es" if x["mm"]!=1 else ""} over {L} bp; SNP at {x["snp"]:,} bp'
                if x['mm']>2: e={'chr':cons,'flag':'unmapped','qc':f'tag not found on {cons}; nearest match {x["chr"]} {x["snp"]:,} bp with {x["mm"]} mismatches (homoeologous copy), not used'}
                elif x['chr']==cons or x['chr']=='Un': e={'chr':x['chr'],'peak':round(x['snp']/1e6,6),'flag':'tag_aligned','qc':desc+extra+('; unanchored chromosome (chrUn) in this assembly' if x['chr']=='Un' else '')}
                elif x['mm']==0: e={'chr':x['chr'],'peak':round(x['snp']/1e6,6),'flag':'other_chromosome','qc':desc+f'; exact match on {x["chr"]} instead of {cons} in this assembly'+extra}
                else: e={'chr':cons,'flag':'unmapped','qc':f'tag not found on {cons}; nearest match {x["chr"]} {x["snp"]:,} bp with {x["mm"]} mismatch{"es" if x["mm"]!=1 else ""} (probably the homoeologous copy), not used'}
        asm[g]=e
        hh=pl.get(g) or []
        rows.append([rid,m['marker'],m['clone'],g,'; '.join(x['chr'] for x in hh),'; '.join(str(x['start']) for x in hh),'; '.join(str(x['end']) for x in hh),'; '.join(x['strand'] for x in hh),hh[0]['mm'] if hh else '','; '.join(str(x['snp']) for x in hh),e['flag'],e.get('qc','')])
    v1,v2=asm['RefSeq v1.0'],asm['RefSeq v2.1']
    if v1['flag']=='tag_aligned' and v1['chr']!='Un': base,ref='RefSeq v1.0','RefSeq v1.0'
    else: base,ref='RefSeq v2.1','RefSeq v2.1'
    b=asm[base];r['chr']=b['chr'];r['peak']=b['peak'];r['start']=None;r['end']=None;r['ref']=ref;r['grade']='B'
    r['name']=re.sub(r', unassigned SNP\)',')',r['name'])
    if base=='RefSeq v1.0': del asm['RefSeq v1.0']
    else: pass
    r['asm']=asm
    n=r['note']
    n=re.sub(r' Chromosome not assigned by the vendor \(NA in Table 2\)\.','',n)
    n=n.replace(' Not placed: DArT tag with no match in RefSeq v1.0.','');n=re.sub(r' Not placed: DArT tag seed hits? [^;]*?mismatch(es)?\.','',n)
    v1t=(f"RefSeq v1.0 {v1['chr']} {int(round(v1['peak']*1e6)):,} bp" if v1['chr']!='Un' else f"RefSeq v1.0 places it on the unanchored chrUn at {int(round(v1['peak']*1e6)):,} bp, so the RefSeq v2.1 position {cons} {int(round(v2['peak']*1e6)):,} bp is the record's position") if 'peak' in v1 else 'not found in RefSeq v1.0'
    placed=sum(1 for g in ORDER if (asm.get(g,{}).get('flag')=='tag_aligned' and asm[g]['chr']==cons) or (g=='RefSeq v1.0' and base=='RefSeq v1.0'))
    n+=(f" Chromosome NA in Table 2 (not assigned by the vendor). Placed by aligning the DArT tag to the genome assemblies: the 69 bp tag was trimmed at the sequencing-adapter read-through"
        f" (AGATCGGAAGAGC) to its {L} bp genomic part (PstI to MspI site), and matched end to end (up to 3 mismatches, no gaps) against IWGSC RefSeq v1.0, v2.1 and 13 pangenome cultivars"
        f" (Ensembl Plants). Placed on {cons} in {placed} of 15 assemblies; position given is the SNP base ({v1t}). Each assembly's position comes from its own alignment, not from liftover."
        f" Grade B: position derived from the marker sequence, not printed in the paper. Tag alignment Nov 2026 (ad65.py).")
    r['note']=n
    r['evidence']=r['evidence']+f"; DArT tag aligned to {ref} {r['chr']} {int(round(r['peak']*1e6)):,} bp (SNP base)"
    out.append(r)
    print(rid,m['marker'],ref,r['chr'],r['peak'],placed,'/15')
json.dump({'reports':out,'sources':{}},open(W+'seed_ad65.json','w'),indent=0,ensure_ascii=False)
import os;os.makedirs(C,exist_ok=True)
with open(C+'tag_placements.tsv','w',newline='') as fh:
    w=csv.writer(fh,delimiter='\t');w.writerow(['record','marker','dart_clone','assembly','chr','tag_start','tag_end','strand','mismatches','snp_position','flag','note']);w.writerows(rows)
with open(C+'tags_trimmed.tsv','w',newline='') as fh:
    w=csv.writer(fh,delimiter='\t');w.writerow(['record','marker','dart_clone','snp_offset','ref_allele','alt_allele','tag_69bp','trimmed_tag'])
    for rid,m in meta.items(): w.writerow([rid,m['marker'],m['clone']]+m['snp']+[m['tag'],m['trimmed']])
