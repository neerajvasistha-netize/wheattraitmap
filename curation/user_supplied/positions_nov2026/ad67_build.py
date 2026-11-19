"""Second batch of uploaded papers (19 Nov 2026): Dhariwal et al. 2020 IJMS (FHB, AAC Tenacious) checked and completed;
three stripe/leaf/stem rust QTL papers added. Writes upd67.json (records) and src67.json (new sources)."""
import sys,json,copy
sys.path.insert(0,'.')
from load import load
from mp import look
R,S,base=load()
out={};TAG=' Physical position from the uploaded paper, Nov 2026 (ad67.py).'
def mk_peak(marker,ch,ref,lo,hi):
    L=[x for x in look(marker,ch) if x[2]==ref]
    if L and lo-0.5<=L[0][1]<=hi+0.5: return L[0][1]
    return None
# --- Dhariwal et al. 2020 IJMS 21:4497, Table 3 (closest marker positions, RefSeq v2.0)
T3=[('QFhi.lrdc-2B','2B','DI','BS00012078_51',265366661,265366634,68.8,'64.2-71.8',3.2,4.4,'T'),
('QFhi.lrdc-2D','2D','DI','Ppd-D1',36205433,36204772,40.9,'40.2-63.7',15.5,25.0,'T'),
('QFhi.lrdc-4A','4A','DI','BobWhite_c10610_1096',688409042,688408942,118.0,'114.8-121.6',8.0,11.9,'I'),
('QFhi.lrdc-4D','4D','DI','BobWhite_c28101_376',515801764,515801671,146.0,'145.1-146.0',2.8,4.0,'I'),
('QFhi.lrdc-5D','5D','DI','RAC875_c92929_92',279491136,279491036,0.1,'0.0-0.6',2.6,3.9,'T'),
('QFhi.lrdc-7A','7A','DI','IACX2471',673715328,673715447,206.4,'205.9-206.8',2.6,2.1,'T'),
('QFhs.lrdc-2D.1','2D','DS','Ppd-D1',36205433,36204772,40.9,'40.2-63.4',12.7,19.0,'T'),
('QFhs.lrdc-2D.2','2D','DS','BobWhite_c17782_194',555098707,555098805,115.6,'105.2-126.4',5.4,15.2,'T'),
('QFhs.lrdc-4A','4A','DS','BobWhite_c10610_1096',688409042,688408942,118.0,'114.8-121.3',5.4,10.1,'I'),
('QFhs.lrdc-7B','7B','DS','wsnp_CAP7_c44_26549',17091249,17091049,32.3,'24.3-35.5',3.6,6.3,'I'),
('QFhb.lrdc-2D.1','2D','VRI','Ppd-D1',36205433,36204772,40.9,'40.2-63.4',12.4,18.0,'T'),
('QFhb.lrdc-2D.2','2D','VRI','BobWhite_c17782_194',555098707,555098805,114.6,'106.4-120.8',3.6,12.2,'T'),
('QFhb.lrdc-4A','4A','VRI','BobWhite_c10610_1096',688409042,688408942,118.0,'114.7-121.3',7.1,11.5,'I'),
('QFhb.lrdc-5D','5D','VRI','RAC875_c92929_92',279491136,279491036,0.2,'0.0-0.4',2.5,3.0,'T'),
('QFhb.lrdc-7B','7B','VRI','wsnp_CAP7_c44_26549',17091249,17091049,31.6,'24.1-35.2',3.9,7.3,'I'),
('QDon.lrdc-2D.1','2D','DON','Ppd-D1',36205433,36204772,40.6,'37.9-41.9',3.9,12.3,'T'),
('QDon.lrdc-2D.2','2D','DON','BobWhite_c17782_194',555098707,555098805,112.6,'105.1-126.6',14.4,34.5,'T'),
('QDon.lrdc-3B.1','3B','DON','Kukri_c13345_481',753508932,753508832,132.8,'132.2-133.2',2.7,9.3,'T'),
('QDon.lrdc-3B.2','3B','DON','Excalibur_rep_c97324_623',771943609,771943709,147.3,'144.5-150.3',2.7,4.0,'T'),
('QDon.lrdc-4B','4B','DON','Ex_c101685_705',31725418,31725518,49.5,'46.4-52.4',6.6,8.0,'T'),
('QDon.lrdc-5A.1','5A','DON','IACX2540',621553818,621553699,164.8,'156.4-169.0',3.1,3.0,'T'),
('QDon.lrdc-7B','7B','DON','Kukri_c19823_491',60918162,60918262,52.9,'52.2-57.7',3.4,5.0,'I')]
TR={'DI':'disease incidence','DS':'disease severity','VRI':'visual rating index','DON':'deoxynivalenol content'}
srcI='X5_1919';by={R[i]['name']:i for i in R if R[i].get('source')==srcI}
k=[int(i[2:]) for i in R if i.startswith('NK')];k=max(k+[51])
import re
seq=[int(r['id'][2:]) for r in json.load(open('/home/claude/wheatnhr-atlas/web/seed_ad66.json'))['reports'] if r['id'].startswith('NK')];k=max([k]+seq)
def nid():
    global k;k+=1;return f'NK{k:04d}'
for nm,ch,tr,m,a,b,cm,ci,lod,r2,don in T3:
    lo,hi=sorted([a,b])
    rec=dict(name=nm,cls='QTL',disease='Fusarium head blight',chr=ch,start=None,end=None,peak=round(lo/1e6,6),ref='RefSeq v2.0',marker=m,stat=f'LOD {lod}',pve=str(r2),source=srcI,grade='A',curation='hand',
      evidence=f'Table 3: {nm} ({TR[tr]}) {ch}, peak {cm} cM (interval {ci} cM), LOD {lod}, R2 {r2}%, closest marker {m} at {lo:,}-{hi:,} bp (IWGSC RefSeq v2.0)',
      note=f'AAC Innova x AAC Tenacious doubled haploid population, pooled data; QTL for FHB {TR[tr]}, favourable allele from {"AAC Tenacious" if don=="T" else "AAC Innova"}. Position is that of the closest marker (Table 3, RefSeq v2.0).'+TAG)
    if nm in by:
        r=copy.deepcopy(R[by[nm]]);r.pop('asm',None);r.update(rec);out[r['id']]=r
    else: i=nid();out[i]=dict(id=i,**rec)
# QTL named only in the epistasis text (individual environments, Supplementary Table S6): flanking markers of the epistasis interval
EP={'QFhs.lrdc-2B':('2B',98.6478,107.7114,'RefSeq v2.1','Tdurum_contig53156_111; Tdurum_contig1653_190','72.6-82.5 cM on 2B'),
    'QFhb.lrdc-2B':('2B',98.6478,107.7114,'RefSeq v2.1','Tdurum_contig53156_111; Tdurum_contig1653_190','72.6-82.5 cM on 2B'),
    'QDon.lrdc-5A.2':('5A',645.4687,664.2744,'RefSeq v1.0','Kukri_c108256_381; BS00090847_51','184.7-196.9 cM on 5A')}
for nm,(ch,a,b,ref,mm,txt) in EP.items():
    r=copy.deepcopy(R[by[nm]]);r.pop('asm',None)
    r.update(chr=ch,start=a,end=b,peak=None,ref=ref,marker=mm,grade='B',curation='hand',evidence=f'Results: epistasis region {txt} overlaps {nm} (Table 4, Supplementary Table S6)',
      note=f'QTL detected in individual environments (Supplementary Table S6, not in the paper). Placed at the epistasis interval the paper says overlaps it ({txt}), using the 90K positions of its flanking markers {mm.replace("; "," and ")}.'+TAG)
    out[r['id']]=r
# the old X5 records that held cM values read as Mb and duplicate names: QFhs.lrdc-2B handled above; others matched by name
# --- new sources
SRC={'NK_Ciechanowska2022':{'citation':'Ciechanowska I, Semagn K, McCallum B, Randhawa H, Strenzke K, Dhariwal R, Virginillo M, Iqbal M, Spaner D (2022) Quantitative trait locus mapping of rust resistance and agronomic traits in spring wheat. Canadian Journal of Plant Science 102:1139-1150','url':'https://doi.org/10.1139/CJPS-2022-0023','short':'Ciechanowska et al. (2022)'},
     'NK_Farzand2021':{'citation':"Farzand M, Dhariwal R, Hiebert CW, Spaner D, Randhawa HS (2021) Mapping quantitative trait loci associated with stripe rust resistance from the Canadian wheat cultivar 'AAC Innova'. Canadian Journal of Plant Pathology 43:S227-S241",'url':'https://doi.org/10.1080/07060661.2021.1982011','short':'Farzand et al. (2021)'},
     'NK_Farzand2022':{'citation':'Farzand M, Dhariwal R, Hiebert CW, Spaner D, Randhawa HS (2022) QTL mapping for adult plant field resistance to stripe rust in the AAC Cameron/P2711 spring wheat population. Crop Science 62:1088-1106','url':'https://doi.org/10.1002/csc2.20741','short':'Farzand et al. (2022)'}}
# Ciechanowska 2022 Table 3 (RefSeq v2.0)
for nm,dis,ch,cm,a,b,l,rr,lod,pve in [('QLr.dms-4A','Leaf rust','4A',3127.0,646.4,648.4,'Ra_c63534_581','RAC875_c6939_1042',12.1,9.0),('QSr.dms-1A','Stem rust','1A',1305.0,536.8,543.6,'Excalibur_rep_c103592_955','RAC875_rep_c69334_132',5.6,6.0),('QSr.dms-2B','Stem rust','2B',3143.0,694.9,695.2,'BobWhite_c3871_428','BS00065914_51',9.5,22.3)]:
    i=nid();out[i]=dict(id=i,name=nm,cls='QTL',disease=dis,chr=ch,start=a,end=b,peak=None,ref='RefSeq v2.0',marker=f'{l}; {rr}',stat=f'LOD {lod}',pve=str(pve),source='NK_Ciechanowska2022',grade='A',curation='hand',
      evidence=f'Table 3: {nm} {ch} {cm} cM, flanking {l} ({a} Mb) and {rr} ({b} Mb), LOD {lod}, PVE {pve}% (IWGSC RefSeq v2.0)',
      note='HYAYT12-10 x GP146 doubled haploid population (167 lines, 90K SNP array), least-squares means across environments; interval between the flanking markers on RefSeq v2.0.'+TAG)
# Farzand 2021 (AAC Innova / AAC Proclaim) Table 3, RefSeq v1.0; combined-environment row where present
for nm,ch,env,cm,ci,lod,add,r2,m,a,b,don in [('QYr.lrdc-2A','2A','Combined',238.4,'218.1-239.9',21.25,14.90,33,'Excalibur_c10461_520',717,765,'I'),('QYr.lrdc-2B.1','2B','LET20',108.8,'103.9-112.8',3.41,7.78,5,'CAP12_rep_c3980_87',115,174,'I'),('QYr.lrdc-2B.2','2B','LET17',145.4,'139.3-149.3',5.73,9.48,7,'Tdurum_contig10219_706',694,744,'I'),('QYr.lrdc-2B.3','2B','CRE20',168.9,'163.3-172.5',3.90,7.44,6,'RAC875_c65882_668',750,763,'I'),('QYr.lrdc-2B.4','2B','Combined',184.2,'182.6-189.3',5.92,6.74,7,'BS00074091_51',766,779,'I'),('QYr.lrdc-2D','2D','Combined',6.3,'4.3-8.4',3.70,-5.43,4,'wsnp_JD_rep_c63957_40798083',21,22,'P'),('QYr.lrdc-3B','3B','LET20',26.3,'25.3-27.4',3.27,7.21,5,'RAC875_c36432_197',16,18,'I'),('QYr.lrdc-5A','5A','Combined',64.7,'63.6-65.7',5.99,-6.88,7,'Tdurum_contig17216_310',700,702,'P'),('QYr.lrdc-5B','5B','LET17',88.4,'87.2-88.7',2.61,14.22,3,'RAC875_c19099_308',519,534,'I'),('QYr.lrdc-5D','5D','LET17',5.4,'5.4-5.4',2.51,6.01,3,'Excalibur_c50584_358',549,550,'I'),('QYr.lrdc-7D','7D','CRE20',0.0,'0.0-1.2',4.41,-7.55,6,'wsnp_Ra_c6894_11980338',55,59,'P')]:
    pk=mk_peak(m,ch,'RefSeq v1.0',a,b);i=nid()
    out[i]=dict(id=i,name=nm,cls='QTL',disease='Stripe rust',chr=ch,start=float(a),end=float(b),peak=pk,ref='RefSeq v1.0',marker=m,stat=f'LOD {lod}',pve=str(r2),source='NK_Farzand2021',grade='A',curation='hand',
      evidence=f'Table 3: {nm} {ch} ({env}) {cm} cM (CI {ci} cM), LOD {lod}, R2 {r2}%, closest marker {m}, physical position {a}-{b} Mb (IWGSC RefSeq v1.0)',
      note=f"AAC Innova x AAC Proclaim doubled haploid population, field stripe rust severity at Lethbridge and Creston (2016, 2017, 2020); resistance allele from {'AAC Innova' if don=='I' else 'AAC Proclaim'}."+(f' Peak at the 90K position of the closest marker {m}.' if pk else '')+TAG)
# Farzand 2022 (AAC Cameron / P2711) Table 3, RefSeq v1.0
for nm,ch,env,cm,ci,lod,r2,m,a,b,don,extra in [('QYr.lrdc-1A.1','1A','CRE20',63.6,'62.9-64.3',4.5,7.2,'Tdurum_contig47183_205',41.9,47.2,'Cam',''),('QYr.lrdc-1A.2','1A','Combined',153.8,'153.4-153.7',3.7,6.0,'Kukri_c52420_112',560.5,572.0,'Cam',''),('QYr.lrdc-2A.1','2A','Combined',2.0,'0.0-15.0',9.5,20.6,'VPM',4.0,36.9,'RP',' Corresponds to Yr17 on the 2NS/2AS (VPM) translocation.'),('QYr.lrdc-2A.2','2A','CRE20',79.6,'78.4-79.8',3.6,5.5,'BS00022301_51',597.5,675.0,'RP',''),('QYr.lrdc-2B','2B','Combined',94.8,'91.2-98.0',8.0,14.0,'Tdurum_contig54704_176',47.4,68.2,'Cam',''),('QYr.lrdc-3B','3B','CRE19',20.1,'19.2-20.3',4.0,6.0,'RFL_Contig4531_1195',6.3,13.8,'RP',''),('QYr.lrdc-5A','5A','Combined',169.6,'166.9-176.1',5.7,8.4,'RFL_Contig316_572',593.3,611.6,'RP','')]:
    pk=mk_peak(m,ch,'RefSeq v1.0',a,b);i=nid()
    out[i]=dict(id=i,name=nm,cls='QTL',disease='Stripe rust',chr=ch,start=a,end=b,peak=pk,ref='RefSeq v1.0',marker=m,stat=f'LOD {lod}',pve=str(r2),source='NK_Farzand2022',grade='A',curation='hand',
      evidence=f'Table 3: {nm} {ch} ({env}) {cm} cM (CI {ci} cM), LOD {lod}, R2 {r2}%, closest marker {m}, physical position {a}-{b} Mb (IWGSC RefSeq v1.0)',
      note=f"AAC Cameron x P2711 recombinant inbred lines, adult-plant field stripe rust severity at Creston and Lethbridge (2018-2020); resistance allele from {'P2711' if don=='RP' else 'AAC Cameron'}."+extra+(f' Peak at the 90K position of the closest marker {m}.' if pk else '')+TAG)
json.dump(list(out.values()),open('upd67.json','w'),indent=0,ensure_ascii=False);json.dump(SRC,open('src67.json','w'),indent=0)
from collections import Counter
print(len(out),Counter(r['source'] for r in out.values()),Counter('new' if i not in R else 'upd' for i in out))
for r in out.values(): print(r['id'],r['name'],r['chr'],r['start'],r['end'],r['peak'],r['ref'])
