#!/usr/bin/env python3
"""Rebuild Lhamo et al. 2024 (TAG 137:193) from Table 6."""
import json, re, shutil, sys

P = '/home/claude/atlas/merged_v13f.json'
M = json.load(open(P))
SRC = 'Z_1228'
TS, SNB = 'Tan spot', 'Septoria nodorum blotch'
ISO = {'Pti2': (TS, 'Pyrenophora tritici-repentis isolate Pti2, race 1'),
       '86-124': (TS, 'Pyrenophora tritici-repentis isolate 86-124, race 2'),
       '331-9': (TS, 'Pyrenophora tritici-repentis isolate 331-9, race 3'),
       'DW5': (TS, 'Pyrenophora tritici-repentis isolate DW5, race 5'),
       'Sn4': (SNB, 'Parastagonospora nodorum isolate Sn4'),
       'Sn2000': (SNB, 'Parastagonospora nodorum isolate Sn2000')}

# isolate, SNP, chr, position Mbp, 9K marker code or '', cM or None, -log10(p), PVE %
T6 = [
('Pti2','wsnp_Ex_c8429_14198883','1B',332.83,'IWA4806',37.19,13.1,7.0),
('Pti2','S1B_343953587','1B',343.95,'',None,14.9,6.5),
('Pti2','S1B_581234893','1B',581.23,'',None,7.8,2.1),
('Pti2','wsnp_Ra_c13959_21943006','1B',640.55,'IWA7619',101.09,9.2,11.9),
('Pti2','wsnp_BG274687B_Ta_2_1','1B',686.84,'IWA545',136.89,5.7,2.8),
('Pti2','S3A_667161197','3A',667.16,'',None,8.5,3.0),
('Pti2','S3B_733335913','3B',733.33,'',None,12.3,6.9),
('Pti2','S4A_611044981','4A',611.04,'',None,5.0,4.2),
('Pti2','S4B_658553895','4B',658.55,'',None,4.4,4.8),
('Pti2','S5B_38353962','5B',38.35,'',None,4.3,1.3),
('Pti2','S6A_610529874','6A',610.53,'',None,10.4,5.5),
('Pti2','wsnp_Ex_c8614_14453388','7A',514.89,'IWA4817',100.51,7.1,46.3),
('Pti2','S7B_749390045','7B',749.39,'',None,5.4,2.3),
('86-124','S1A_4206845','1A',4.20,'',None,5.0,2.1),
('86-124','S1A_101243220','1A',101.24,'',None,4.0,4.9),
('86-124','S1B_426502100','1B',426.50,'',None,8.9,1.8),
('86-124','S1B_675363619','1B',675.36,'',None,4.2,4.1),
('86-124','S3A_689112769','3A',689.11,'',None,9.6,8.4),
('86-124','S3B_44395858','3B',44.39,'',None,5.6,2.1),
('86-124','S4A_610623925','4A',610.62,'',None,7.0,4.3),
('86-124','S4A_628969385','4A',628.96,'',None,4.6,0.3),
('86-124','wsnp_Ex_c30581_39482788','4B',563.80,'IWA3287',71.64,6.4,5.8),
('86-124','wsnp_Ex_c6100_10676217','5B',415.64,'IWA4422',159.37,9.8,4.4),
('86-124','wsnp_Ex_c6548_11355524','5B',439.72,'IWA4533',138.45,13.7,22.3),
('86-124','S5B_709562067','5B',709.56,'',None,7.7,6.9),
('86-124','wsnp_Ex_c24301_33540742','6B',123.83,'IWA2888',38.51,5.5,6.3),
('86-124','S7B_662906335','7B',662.90,'',None,4.4,2.6),
('86-124','S7B_663741929','7B',663.74,'',None,6.4,1.3),
('331-9','S1A_23408829','1A',23.40,'',None,4.9,3.8),
('331-9','S1A_542693048','1A',542.69,'',None,4.8,2.7),
('331-9','wsnp_Ku_c11987_19473555','1B',56.31,'IWA6449',40.42,12.6,8.8),
('331-9','S3B_798111300','3B',798.11,'',None,4.2,2.4),
('331-9','S4A_689082898','4A',689.08,'',None,6.9,4.9),
('331-9','wsnp_RFL_Contig2187_1525458','5A',472.34,'IWA8308',123.52,6.2,6.9),
('331-9','S6B_702389600','6B',702.39,'',None,5.7,5.5),
('331-9','wsnp_Ex_c8614_14453388','7A',514.89,'IWA4817',100.51,5.5,44.3),
('DW5','wsnp_Ex_rep_c88749_79992654','1B',676.26,'IWA5758',None,6.1,20.1),
('DW5','wsnp_Ku_c27759_37702654','2B',184.49,'IWA6830',113.52,6.9,3.6),
('DW5','wsnp_Ex_rep_c68173_66946365','2B',632.17,'IWA5506',119.90,5.4,2.8),
('DW5','S3A_727403936','3A',727.40,'',None,5.2,7.3),
('DW5','wsnp_RFL_Contig3599_3786631','3B',725.70,'IWA8490',127.10,12.0,6.3),
('DW5','wsnp_Ex_c4286_7734046','4A',114.74,'IWA3845',151.32,6.3,13.7),
('DW5','S4A_689130908','4A',689.13,'',None,5.0,5.7),
('DW5','S4B_30865852','4B',30.86,'',None,4.2,5.4),
('DW5','wsnp_Ex_c31955_40681185','7A',714.59,'IWA3371',175.63,8.4,8.4),
('DW5','S7B_105890653','7B',105.89,'',None,5.7,3.6),
('Sn4','S1A_272986948','1A',272.98,'',None,7.4,0.5),
('Sn4','S2B_576564612','2B',576.56,'',None,15.4,7.7),
('Sn4','S2B_672777969','2B',672.77,'',None,5.8,3.5),
('Sn4','S2B_703712740','2B',703.71,'',None,4.8,1.3),
('Sn4','wsnp_Ra_c18164_27178535','3B',726.48,'IWA7680',127.10,15.1,13.2),
('Sn4','S4A_602735878','4A',602.73,'',None,7.1,4.7),
('Sn4','S4A_668404073','4A',668.40,'',None,5.1,4.0),
('Sn4','wsnp_Ex_c662_1301994','4A',723.51,'IWA4558',207.07,5.0,20.5),
('Sn4','wsnp_Ex_c18107_26909127','5A',644.13,'IWA2282',44.95,6.3,12.4),
('Sn4','S6B_117509726','6B',117.51,'',None,4.2,3.1),
('Sn4','wsnp_Ex_c3248_5987129','7B',201.26,'IWA3393',54.08,7.5,6.1),
('Sn2000','wsnp_Ex_rep_c102909_87976721','1B',532.17,'IWA5076',64.07,5.5,6.1),
('Sn2000','S1B_585328785','1B',585.33,'',None,6.8,20.8),
('Sn2000','wsnp_Ku_c19631_29148397','3B',469.79,'IWA6677',80.93,8.3,11.1),
('Sn2000','wsnp_Ku_c12503_20174234','4B',637.39,'IWA6461',85.23,8.6,6.4),
('Sn2000','wsnp_Ex_c7713_13153321','6B',684.27,'IWA4717',127.53,5.2,2.9),
('Sn2000','wsnp_Ku_rep_c70721_70356213','6B',692.78,'IWA7479',138.48,4.7,3.6),
('Sn2000','S7A_30905269','7A',30.91,'',None,6.3,14.2),
('Sn2000','wsnp_Ex_c2277_4267788','7A',567.58,'IWA2752',98.42,4.7,3.3),
]

old = [r for r in M['reports'] if r['source'] == SRC]
# keep the previously reported loci the paper discusses (meta-QTL and QSnb.nmbu-2B.1)
KEEP = {'MQTL-3B.1', 'MQTL-5B.1', 'MQTL-5A.1', 'QSnb.nmbu-2B.1'}
kept = [r for r in old if r['name'] in KEEP]
print('existing records:', len(old), '| keeping previously reported loci:', [r['name'] for r in kept])

# compare before dropping
MPOS = {}
for r in old:
    if r['name'] not in KEEP:
        MPOS[r['name']] = r.get('peak')
wrong = []
for iso, snp, ch, mb, code, cm, lp, pve in T6:
    if snp in MPOS and MPOS[snp] is not None and abs(MPOS[snp] - mb) > 0.5:
        wrong.append((snp, ch, MPOS[snp], mb))
    elif snp in MPOS and MPOS[snp] is None:
        wrong.append((snp, ch, None, mb))
print('\nstored position wrong or missing vs the paper (%d):' % len(wrong))
for snp, ch, was, now in sorted(wrong, key=lambda x: -(abs((x[2] or 0) - x[3]))):
    print('   %-30s %-3s stored %-10s paper %-9s %s' %
          (snp, ch, was if was is not None else 'none', now,
           'off by %.1f Mb' % abs(was - now) if was is not None else 'was unplaced'))

M['reports'] = [r for r in M['reports'] if r['source'] != SRC] + kept
have = {r['id'] for r in M['reports']}
nxt = max(int(r['id'][1:]) for r in M['reports'] if re.fullmatch(r'V\d+', r['id'])) + 1
NEW = []
for iso, snp, ch, mb, code, cm, lp, pve in T6:
    i = 'V%04d' % nxt
    while i in have:
        nxt += 1; i = 'V%04d' % nxt
    have.add(i); nxt += 1
    dis, isod = ISO[iso]
    plat = '9K SNP array' if snp.startswith('wsnp') else 'GBS'
    NEW.append({'id': i, 'name': snp, 'cls': 'MTA', 'disease': dis, 'chr': ch,
      'start': None, 'end': None, 'peak': mb, 'ref': 'RefSeq v1.0',
      'marker': snp + ((' (%s)' % code) if code else ''),
      'pve': '%.1f' % pve, 'stat': '-log10(p) %.1f' % lp, 'source': SRC, 'grade': 'A',
      'curation': 'auto',
      'note': ('Marker-trait association for %s in an association-mapping panel of 180 cultivated '
               'emmer wheat (Triticum turgidum ssp. dicoccum) accessions, %s (Lhamo et al. 2024, '
               'Table 6). Marker platform: %s. Position in Mb on IWGSC RefSeq v1.0 as reported.%s'
               % (dis.lower(), isod, plat,
                  ' Genetic position %.2f cM.' % cm if cm is not None else '')),
      'evidence': 'Table 6: %s %s %.2f Mb; -log10(p) %.1f; PVE %.1f%%%s'
                  % (snp, ch, mb, lp, pve, '; 9K code %s' % code if code else '')})
M['reports'].extend(NEW)
for r in kept:
    r['note'] = (r.get('note', '').rstrip('. ') + '. ' if r.get('note') else '') + \
      'This is a previously reported locus discussed by Lhamo et al. (2024) for comparison, not a ' \
      'result of that study.'
print('\nrebuilt:', len(NEW), 'associations (%d tan spot, %d SNB), all placed'
      % (sum(1 for x in NEW if x['disease'] == TS), sum(1 for x in NEW if x['disease'] == SNB)))
if '--write' in sys.argv:
    shutil.copy(P, '/home/claude/atlas/merged_v13f_preLhamo.json')
    json.dump(M, open(P, 'w'))
    print('written. total reports', len(M['reports']))
