"""Physical positions for ten uploaded papers (Nov 2026). Writes upd.json: list of full records (changed or new)."""
import sys,json,copy,re
sys.path.insert(0,'.')
from load import load
from mp import look
R,S,base=load()
out={};TAG=' Physical position added from the uploaded paper, Nov 2026 (ad66.py).'
def upd(i,**kw):
    r=copy.deepcopy(out.get(i) or R[i]);r.pop('asm',None)
    note=kw.pop('note',None)
    for k,v in kw.items(): r[k]=v
    if note is not None:
        old=(r.get('note') or '').strip()
        if R.get(i,{}).get('curation')=='auto' or 'auto-extracted' in old: old=''
        placed=r.get('peak') is not None or r.get('start') is not None
        r['note']=(old+(' ' if old else '')+note+(TAG if placed else ' Checked against the uploaded paper, Nov 2026 (ad66.py).')).strip()
    if r.get('curation')=='auto': r['curation']='hand'
    out[i]=r
def new(i,**kw):
    r=dict(id=i,cls='QTL',curation='hand',marker='',pve='',stat='',start=None,end=None,peak=None);r.update(kw);r['note']=r.get('note','')+TAG;out[i]=r
MB=lambda x:round(x/1e6,6)
# 1. Xu et al. 2024 TAG, Fusarium crown rot (Table 2 flanking SNPs named chr_position, RefSeq v1.0)
for i,(a,b) in {'V1363':(27199242,29978751),'V1364':(33084970,41351613),'V1365':(48783034,51143792),'V1366':(759449813,772606231),'V1367':(43599467,48783034),'V1368':(704324449,719292783)}.items():
    upd(i,start=MB(a),end=MB(b),peak=None,ref='RefSeq v1.0',grade='A',marker=f"{R[i]['chr']}_{a}; {R[i]['chr']}_{b}",note=f"Table 2: flanking SNPs {R[i]['chr']}_{a} and {R[i]['chr']}_{b} (SNP names give the RefSeq v1.0 position).")
fix6=" Table 2 names the 6B and 7B flanking SNPs with one digit fewer (e.g. 6B_28163872), but the Discussion gives Qfcr.cau-6B as 281.6-327.0 Mb (containing barc198) and Qfcr.cau-7B.1 as 315.6-432.5 Mb (overlapping the Collard et al. 2005 7B QTL at 342.9-359.0 Mb), so the SNP names were read as 281,638,72x etc."
upd('V1369',start=281.6,end=327.0,peak=None,ref='RefSeq v1.0',grade='B',marker='6B_28163872; 6B_32696670; barc198',note='Discussion: Qfcr.cau-6B at 281.6-327.0 Mb.'+fix6)
upd('V1370',start=315.6,end=432.5,peak=None,ref='RefSeq v1.0',grade='B',marker='7B_31556403; 7B_43246561',note='Discussion: Qfcr.cau-7B.1 at 315.6-432.5 Mb.'+fix6)
upd('V1371',start=431.27282,end=458.44119,peak=None,ref='RefSeq v1.0',grade='B',marker='7B_43127282; 7B_45844119',note='Flanking SNPs 7B_43127282 and 7B_45844119 read as 431.27 and 458.44 Mb, as for Qfcr.cau-7B.1.'+fix6)
for i,p in {'N0015':281.63872,'N0016':326.9667,'N0017':315.56403,'N0018':432.46561,'N0019':431.27282,'N0020':458.44119}.items():
    upd(i,peak=p,grade='B',note='Position read from the SNP name with the missing digit restored, as given in the Discussion (see Qfcr.cau-6B and Qfcr.cau-7B.1).')
# 2. ElDoliefy et al. 2020 Euphytica (Glenn x MN00261-4 DArT map, Fig. 2), positions of map markers from DArT tag BLAST (RefSeq v1.0)
upd('V0242',start=2.2372,end=9.3781,peak=None,grade='B',marker='wPt-664778; wPt-729972; wPt-731843; wPt-671698',note='QTL interval 23-36.6 cM (Table 2). Fig. 2 markers: wPt-664778 (21.4 cM) at 2.24 Mb and wPt-729972, wPt-731843, wPt-671698 (27.2-30.1 cM, inside the interval) at 9.34-9.38 Mb on RefSeq v1.0.')
upd('V0243',start=8.2215,end=10.1165,peak=None,grade='B',marker='wPt-744675; wPt-744265; wPt-744388; wPt-664286; wPt-664309; wPt-743380',note='QTL interval 19.9-21.4 cM (Table 2). Fig. 2 markers at 19.6-21.1 cM map to 8.22-10.12 Mb on RefSeq v1.0.')
upd('V1274',start=9.1957,end=23.1551,peak=None,grade='B',marker='wPt-740561; wPt-0002',note='QTL interval 53.1-53.7 cM (Table 2), between wPt-740561 (42.3 cM, 23.16 Mb) and wPt-0002 (57.8 cM, 9.20 Mb) on RefSeq v1.0 (map orientation reversed).')
upd('V1275',start=151.8128,end=466.7733,peak=None,grade='B',marker='wPt-741530; wPt-2564',note='QTL peak 52.7 cM (Table 2), between wPt-741530 (49.8 cM, 151.81 Mb) and wPt-2564 (55.7 cM, 466.77 Mb) on RefSeq v1.0; the nearest markers (wPt-6522, wPt-7935, rPt-1040) have no physical position, so the interval is wide and pericentromeric although the paper names it 6BS.')
# 3. Farhad et al. 2026 PLOS One, wheat blast (Table 3, RefSeq v2.1)
F=json.load(open('farhad.json'))
F.append(dict(n=40,name='Qwb.cim-7A.03',chr='7A',lead=80160295,beta='8.67',start=80160295,end=80160465,length=170,ov='2',r2='1'))
F.sort(key=lambda x:x['n'])
have={R[i]['name']:i for i in ['V1352','V1353','V1354','V1355']}
k=3
for f in F:
    pk=MB(f['lead']) if f['start']<=f['lead']<=f['end'] else None
    ev=f"Table 3 row {f['n']}: {f['name']} {f['chr']}, leading SNP {f['lead']:,} bp, region {f['start']:,}-{f['end']:,} bp ({f['ov']} {'QTL' if f['ov']=='Single' else 'SNPs, overlapping QTLs'}), beta {f['beta']}, R2 (LD) {f['r2']}"
    nt='GWAS of a Bangladeshi wheat panel for wheat blast (WB index) in Jashore (Bangladesh), Okinawa (Japan) and Quirusillas (Bolivia); positions on IWGSC RefSeq v2.1. Negative beta lowers WB index (resistance).'+(' The leading SNP printed in Table 3 lies outside the printed region, so no peak is given.' if pk is None else '')+(' Corresponds to the 2NS/2AS translocation.' if f['name']=='Qwb.cim-2A.01' else '')
    kw=dict(name=f['name'],chr=f['chr'],start=MB(f['start']),end=MB(f['end']),peak=pk,ref='RefSeq v2.1',grade='A',marker=f"SNP {f['chr']}_{f['lead']}",stat=f"beta {f['beta']}",evidence=ev)
    if f['name'] in have: upd(have[f['name']],cls='MTA',**kw,note=nt)
    else:
        k+=1;new(f'NK{k:04d}',cls='MTA',disease='Wheat blast',source='V_Farhad2026',**kw,note=nt)
# 4. Zhao et al. 2023 Mol Breeding: the two spike-length QTL with effects on FHB (Table 4, RefSeq v2.1)
for nm,ch,a,b,mk,pve,eff in [('QSl.yaas-5A','5A',448.90,479.78,'AX-111153208; AX-111102760','9.11-15.85','PSS -9.78%'),('QSl.yaas-5B','5B',642.52,658.83,'AX-110581004; AX-110942912','9.01-12.85','PSS -10.40%')]:
    k+=1;new(f'NK{k:04d}',name=nm,disease='Fusarium head blight',chr=ch,start=a,end=b,ref='RefSeq v2.1',grade='A',marker=mk,pve=pve,stat=eff,source='V_Zhao2023',
        evidence=f'Table 4: {nm} {a:.2f}-{b:.2f} Mb (RefSeq v2.1), flanking {mk.replace("; ","-")}; Table 5 / Fig. 4C: positive allele lowers FHB severity (PSS) in 180 cultivars',
        note=f'Spike length QTL (Yangmai 12 x Yanzhan 1 RILs) with a pleiotropic effect on FHB resistance: in the validation panel the positive allele reduced the percentage of symptomatic spikelets after spray inoculation ({eff}); PVE for spike length {pve}%.')
upd('V0360',note='Cited from Hu et al. (2023) in the Discussion; this paper gives no position for it.')
# 5. Semagn et al. 2023 Front Plant Sci (Table 2, RefSeq v2.0)
SEM={'Z0076':(414.2,522.7,None),'Z0081':(12.4,48.8,None),'Z0083':(None,None,640.0),'V0322':(10.2,23.5,None),'V0323':(683.3,710.1,None),'V0324':(597.6,609.4,None),'V0325':(None,None,170.8),'V0326':(None,None,159.1),'V0327':(None,None,123.9),'V0328':(3.1,11.7,None),'V0329':(None,None,477.8),'V0330':(43.0,89.2,None)}
for i,(a,b,p) in SEM.items():
    pos=f'{a}-{b}' if a is not None else f'{p}'
    upd(i,start=a,end=b,peak=p,ref='RefSeq v2.0',grade='A',evidence=f"Table 2: {R[i]['name']} {R[i]['chr']} {pos} Mb (IWGSC RefSeq v2.0)",note=f"Checked against Table 2 of the paper: {R[i]['chr']} {pos} Mb on RefSeq v2.0.")
# 6. Yang et al. 2022 Front Plant Sci, STB (Chinese Spring RefSeq v1.0)
upd('V0669',start=440.0,end=615.0,peak=454.0,ref='RefSeq v1.0',grade='B',marker='snp_3026774_F_0_34',note='Results/Discussion: tightly linked snp_3026774_F_0_34 peaks at 454 Mb; probable region 440-615 Mb on the Chinese Spring reference (Supplementary Table 3).')
upd('V0676',start=0.0,end=20.0,peak=None,ref='RefSeq v1.0',grade='B',marker='pav_4991454',note='Results: genomic region 0-20 Mb (Supplementary Table 3).')
upd('V0672',start=40.0,end=100.0,peak=None,ref='RefSeq v1.0',grade='B',marker='snp_1112131_F_0_20',note='Results: genomic region 40-100 Mb (Supplementary Table 3).')
upd('V0678',start=570.0,end=590.0,peak=586.6694,ref='RefSeq v1.0',grade='B',marker='snp_2262549_F_0_28',note='Discussion: region 570-590 Mb (Supplementary Table 3); the tagged DArT-seq SNP 2262549 aligns at 586.67 Mb on RefSeq v1.0.')
upd('V0675',peak=18.5681,start=None,end=None,ref='RefSeq v1.0',grade='B',marker='snp_4910674_F_0_26',note='Representative DArT-seq SNP 4910674 (Table 2, 17.4 cM) aligns at 18.57 Mb on RefSeq v1.0.')
upd('V0680',peak=126.5932,start=None,end=None,ref='RefSeq v1.0',grade='B',marker='snp_2253221_F_0_65',note='Representative DArT-seq SNP 2253221 (Table 2, 71.9 cM) aligns at 126.59 Mb on RefSeq v1.0.')
for i in ['V0670','V0671','V0673','V0674','V0677','V0679']:
    upd(i,note='Not placed: the paper gives only a cM position and a DArT-seq marker that could not be located on the reference.')
# 7. Genievskaya et al. 2020 PeerJ (Table 3; 90K SNP positions from Sun et al. 2020 / Hu et al. 2022)
G=json.load(open('gen_rows.json'))
byname={R[i]['name']:i for i in R if R[i].get('source')=='OA_Genievskaya2020genome'}
for n,tr,mk,q,ch,L in G:
    if not q or q not in byname: continue
    i=byname[q];v10=[x for x in L if x[2]=='RefSeq v1.0'];v21=[x for x in L if x[2]=='RefSeq v2.1']
    if v10: upd(i,peak=v10[0][1],start=None,end=None,ref='RefSeq v1.0',grade='B',marker=mk,note=f'GWAS of spring wheat accessions under field conditions in South Kazakhstan (2018, 2019; 90K array). Table 3: peak SNP {mk} ({ch}, consensus map). Physical position of the SNP from {v10[0][3]}.')
    elif v21: upd(i,peak=v21[0][1],start=None,end=None,ref='RefSeq v2.1',grade='B',marker=mk,note=f'GWAS of spring wheat accessions under field conditions in South Kazakhstan (2018, 2019; 90K array). Table 3: peak SNP {mk} ({ch}, consensus map). Physical position of the SNP from {v21[0][3]}.')
    else:
        o=look(mk)
        why=f'its physical position is on the homoeologous chromosome {o[0][0]} ({o[0][1]} Mb)' if o else 'it has no known physical position'
        upd(i,marker=mk,note=f'GWAS of spring wheat accessions in South Kazakhstan (90K array). Table 3: peak SNP {mk}; not placed because {why}, not on {ch}.')
# 8. Running et al. 2022: QTs.zhl-1A corresponds to Tsc1
t=R['R103'];upd('V0480',start=t['start'],end=t['end'],peak=t['peak'],ref=t['ref'],grade='B',note='The paper states that QTs.zhl-1A corresponds to the Tsc1 region (QTsc.ndsu-1A, 1AS); placed at the Tsc1 candidate interval delineated in the paper (record R103).')
# 9. Ren et al. 2019: QFhb.cau-1AS linked to IWB63682, which has no known physical position
upd('V0320',marker='IWB63682',note='Table 2: closest marker IWB63682 (1AS); not placed because IWB63682 has no physical position in the 90K maps used by the atlas.')
# 10. Shi et al. 2023 Agronomy (Table 1, RefSeq v1.1)
upd('V1321',start=3.1,end=13.3,peak=None,ref='RefSeq v1.1',grade='B',note='Discussion: Qfhb.nmbu.1A.1 (Nannuru et al.) at 3.1-13.3 Mb on 1A, cited as overlapping gFHB-1A.1.')
SH=[('gFHB-2B','2B',694.1,704.0,'Tdurum_contig96648_102~BS00065136_51','7.3-18.2','27','II'),('gFHB-3A.1','3A',1.3,14.9,'Tdurum_contig76679_307~IAAV1155','8.1-19.7','57','II'),('gFHB-3A.2','3A',539.7,570.8,'AX-110433978~AX-109645516','5.0-9.7','100','I'),('gFHB-3D','3D',0.8,9.7,'Excalibur_c9485_707~Excalibur_c60530_113','7.9-19.7','21','II'),('gFHB-6B','6B',132.1,151.1,'BS00003897_51~Kukri_c38732_546','8.8-16.9','13','II'),('gFHB-7A','7A',8.3,68.0,'BobWhite_c17904_259~Tdurum_contig30621_328','7.1-16.6','29','II'),('gFHB-7B','7B',128.9,199.8,'AX-109345488~Tdurum_contig10980_1229','5.8-10.4','114','I & II')]
for nm,ch,a,b,mk,pve,ns,pan in SH:
    k+=1;new(f'NK{k:04d}',cls='MTA',name=nm,disease='Fusarium head blight',chr=ch,start=a,end=b,ref='RefSeq v1.1',grade='A',marker=mk.replace('~','; '),pve=pve,source='Agron2023FHB',
        evidence=f'Table 1: {nm} {mk} {a}~{b} Mb, {ns} significant SNPs, PVE {pve}%, panel {pan}',note=f'GWAS of Chinese wheat panels (55K and 90K SNP arrays); LD block of {ns} significant SNPs, panel {pan}; positions on RefSeq v1.1.')
json.dump(list(out.values()),open('upd.json','w'),indent=0,ensure_ascii=False)
from collections import Counter
print(len(out),'records;',Counter('new' if i not in R else 'changed' for i in out),Counter(r['source'] for r in out.values()))
print('placed now:',sum(1 for r in out.values() if r.get('peak') is not None or r.get('start') is not None))
