import json,re,openpyxl
from openpyxl.styles import Font
S=json.load(open('/home/claude/sr/sr_v21.json'));su=S['summary'];m=round(0.05/S['bonf'])
M=json.load(open('merged_v13_lifted.json'));nrep=len(M['reports'])
GC=[{'id':s['region'],'chr':s['chr'],'lo':s['start_Mb_v2.1'],'hi':s['end_Mb_v2.1'],'g':s['HC_genes'],'n':s['NLR'],'nf':s['NLR_fold_zone'],'np':s['NLR_p_zone'],'ns':bool(s['NLR_sig_zone_bonf']),'w':s['WAK/WAKL'],'wf':s['WAK/WAKL_fold_zone'],'wp':s['WAK/WAKL_p_zone'],'ws':bool(s['WAK/WAKL_sig_zone_bonf']),'pr':s['mean_%present_cultivars']} for s in su]
a=open('web/app.js').read()
a=re.sub(r'const GC=\[.*?\];','const GC='+json.dumps(GC,separators=(',',':'))+';',a,count=1,flags=re.S)
a=a.replace('correction for 200 tests','correction for %d tests'%m)
a=re.sub(r'correction for \d+ tests','correction for %d tests'%m,a)
open('web/app.js','w').write(a)
old=openpyxl.load_workbook('/mnt/user-data/outputs/Shared_regions_gene_content_v1.3_RefSeqv2.1.xlsx')
readme=[r[0] for r in old['README'].iter_rows(values_only=True)]
wb=openpyxl.Workbook();ws=wb.active;ws.title='README'
for t in readme:
    t=re.sub(r'\(\d+ regions from \d+ merged loci, \d+ reports\)','(%d regions, %d reports)'%(len(su),nrep),t)
    t=re.sub(r'Bonferroni threshold for \d+ tests: p < [0-9.]+','Bonferroni threshold for %d tests: p < %.5f'%(m,S['bonf']),t)
    ws.append([t])
ws['A1'].font=Font(bold=True,name='Arial')
cols=[k for k in su[0] if k not in ('NLR/WAK absent in >=1 cultivar',)]
for title,rows in (('Regions',su),('Enriched regions',[s for s in su if s['NLR_sig_zone_bonf'] or s['WAK/WAKL_sig_zone_bonf']])):
    w=wb.create_sheet(title);w.append(cols)
    for s in rows: w.append([s.get(c) if not isinstance(s.get(c),bool) else ('yes' if s[c] else '') for c in cols])
    for c in w[1]: c.font=Font(bold=True,name='Arial')
w=wb.create_sheet('Genes');g=S['genes'];gc=list(dict.fromkeys(k for x in g for k in x));w.append(gc)
for x in g: w.append([x.get(c,"") for c in gc])
for c in w[1]: c.font=Font(bold=True,name='Arial')
wb.save('/mnt/user-data/outputs/Shared_regions_gene_content_v1.3_RefSeqv2.1.xlsx')
print(len(GC),m,sum(1 for s in GC if s['ns'] or s['ws']),'enriched')
