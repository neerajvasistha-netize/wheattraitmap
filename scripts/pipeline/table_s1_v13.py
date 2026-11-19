import json
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment
from openpyxl.utils import get_column_letter as L
import re
from openpyxl.cell.cell import ILLEGAL_CHARACTERS_RE as ILL
clean=lambda row:[ILL.sub('',v) if isinstance(v,str) else v for v in row]
M=json.load(open('/home/claude/atlas/merged_v13_lifted.json'))
ASM=['RefSeq v2.1','RefSeq v1.0','Alchemy','ArinaLrFor','Jagger','Julius','Kariega','Lancer','Landmark','Mace','Norin 61','Paragon','Renan','Stanley','SY Mattis']
F='Arial'; H=Font(name=F,bold=True,color='FFFFFF',size=10); HF=PatternFill('solid',fgColor='1F4E78'); base=Font(name=F,size=10)
wb=Workbook(); ws=wb.active; ws.title='README'
txt=['WheatNHR-Atlas Table S1, release 1.3 (draft, not released)',
 'One row per report (gene, QTL, MTA, QTN). Published position and assembly are kept as given in the source; positions are also given on IWGSC RefSeq v2.1 (atlas frame), RefSeq v1.0 and 13 chromosome-level pangenome assemblies (Ensembl Plants release 63).',
 'Conversion: gene-anchor liftover (IWGSC v1.1 to v2.1 gene mapping; pangenome by CDS matching, minimap2, identity >= 95%, coverage >= 60%). Quality flags: ok; as_reported = position given on v2.1 in the paper; v2.0_as_reported = v2.0 position used directly (v2.1 keeps the v2.0 layout); reoriented_checked = segment reversed between v1.0 and v2.1, both ends confirmed from the nearest anchor genes; rechecked = recalculated from the nearest anchor genes on both sides; inverted = gene order reversed relative to Chinese Spring; other_chromosome = on another chromosome in that cultivar; wide_gap/few_anchors = approximate; unmapped.',
 'Curation: "atlas v1.2" = previous release; "hand" = curated from the paper; "auto" = automatically extracted from paper text (grade B), must be checked against the paper (column "text in paper").',
 'Grades: A sequence-anchored position from the source; B rounded, flanking-marker, second-hand or auto-extracted; C pre-RefSeq; D arm or bin only.']
for i,t in enumerate(txt,1):
    c=ws.cell(i,1,t); c.font=Font(name=F,size=11 if i>1 else 14,bold=(i==1)); c.alignment=Alignment(wrap_text=True)
ws.column_dimensions['A'].width=140
r2=wb.create_sheet('Reports')
cols=['id','name','cls','disease','chr','start','end','peak','ref','marker','pve','stat','effect','grade','curation','source','note','evidence']
hdr=cols[:]
for a in ASM: hdr+=[a+' chr',a+' start',a+' end',a+' peak',a+' quality']
hdr.insert(hdr.index('source')+1,'citation')
r2.append(hdr)
for x in M['reports']:
    row=[x.get(c) for c in cols]
    row.insert(cols.index('source')+1,M['sources'].get(x['source'],{}).get('citation','') if isinstance(M['sources'].get(x['source']),dict) else str(M['sources'].get(x['source'],'')))
    A=x.get('asm',{})
    for a in ASM:
        if a=='RefSeq v1.0' and a not in A and not str(x.get('ref','')).startswith('RefSeq v2') and x.get('chr'):
            o={'chr':x['chr'],'start':x.get('start'),'end':x.get('end'),'peak':x.get('peak'),'flag':'as reported' if (x.get('start') is not None or x.get('peak') is not None) else ''}
        else: o=A.get(a,{})
        row+=[o.get('chr'),o.get('start'),o.get('end'),o.get('peak'),o.get('flag')]
    r2.append(clean(row))
pm={}
for x in M['reports']:
    for p in x.get('primers') or []:
        e=pm.setdefault(p['m'],[p.get('type',''),p['F'],p['R'],p.get('motif',''),p.get('tm',''),p['src'],[]]);e[6].append(x['id'])
s4=wb.create_sheet('Marker primers');s4.append(['marker','type','forward (5-3)','reverse (5-3)','motif','annealing temp (C)','primer source','atlas records'])
for k in sorted(pm): v=pm[k];s4.append(clean([k]+v[:6]+[', '.join(v[6])]))
for c in s4[1]: c.font=H; c.fill=HF
for i,wd in enumerate([14,10,30,30,12,10,40,40],1): s4.column_dimensions[L(i)].width=wd
s3=wb.create_sheet('Sources'); s3.append(['key','citation','url','file'])
for k,v in M['sources'].items():
    s3.append(clean([k,v.get('citation',''),v.get('url',''),v.get('file','')] if isinstance(v,dict) else [k,str(v),'','']))
for w in (r2,s3):
    for c in w[1]: c.font=H; c.fill=HF; c.alignment=Alignment(wrap_text=True,vertical='top')
    for row in w.iter_rows(min_row=2):
        for c in row: c.font=base
    w.freeze_panes='B2'
    for i in range(1,w.max_column+1): w.column_dimensions[L(i)].width=14
r2.column_dimensions['B'].width=28; r2.column_dimensions[L(hdr.index('citation')+1)].width=50
s3.column_dimensions['B'].width=90
r2.auto_filter.ref=r2.dimensions
wb.save('/mnt/user-data/outputs/Table_S1_WheatNHR_Atlas_v1.3.xlsx'); print('saved',r2.max_row-1,'rows')
