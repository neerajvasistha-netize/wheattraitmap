#!/usr/bin/env python3
"""Rebuild Martin et al. 2025 (Crop Science 65:e70023) from Table 4 of the paper."""
import json, shutil, sys

P = '/home/claude/atlas/merged_v13f.json'
M = json.load(open(P))
SRC = 'A_crown_rot_015'
CIT = 'Martin et al. (2025) Table 4'

# QTL name, trial year, LOD, PVE, favourable-allele parent, peak marker,
# flanking marker 1 + its Chinese Spring position, flanking marker 2 + its CS position
T4 = [
 ('QFcr.cg-3D',  2020, 3.3,  5.38,  'PAN 3471', '5411532', '5411532', 495342360, 'cfd152',  538177121),
 ('QFcr.cg-7A',  2020, 5.0,  8.23,  'SST 087',  '1128642', '987757',  575259106, '3532952', 690860873),
 ('QFcr.cg-2A',  2022, 3.6,  6.76,  'PAN 3471', '3954358', '1105077',  42194565, '3034278', 698281086),
 ('QFcr.cg-5B',  2022, 4.3,  7.67,  'SST 087',  '3952912', '3028104', 492093143, '1233307', 545244564),
 ('QFcr.cg-7A',  2022, 7.8, 14.58,  'SST 087',  '2261308', '1083319',  57742411, '1112611', 670246544),
 ('QFcr.usq-1D', 2020, 3.7,  7.19,  'SST 087',  '5358480', 'gdm126b', 455803257, '1360603', 486211404),
 ('QFcr.usq-7A', 2020, 3.3,  6.07,  'SST 087',  '1068343', '1083319',  57742411, '1112611', 670246545),
]
# named in the text as carrying FCR QTL but reported only in Table S3 (not held)
TEXT_ONLY = [('QFcr.cg-1A', '1A'), ('QFcr.cg-6D', '6D')]
CHRS_TEXT = '1A, 1D, 2A, 2B, 3D, 5B, 6D, 7A and 7D'

old = [r for r in M['reports'] if r['source'] == SRC]
keep = [r for r in M['reports'] if r['source'] != SRC]
print('existing records for this source:', len(old))
for r in old:
    print('  dropped', r['id'], r['name'], r['chr'], r['start'], r['end'], r['peak'])

ids = sorted(r['id'] for r in old)
import re
nxt = max(int(r['id'][1:]) for r in M['reports'] if re.fullmatch(r'V\d+', r['id'])) + 1
have = {r['id'] for r in M['reports']}
NEW = []

def nid():
    global nxt
    i = 'V%04d' % nxt
    while i in have:
        nxt += 1; i = 'V%04d' % nxt
    have.add(i); nxt += 1
    return i

for nm, yr, lod, pve, parent, peakm, f1, p1, f2, p2 in T4:
    ch = nm.split('-')[1]
    a, b = sorted((p1 / 1e6, p2 / 1e6))
    NEW.append({
      'id': nid(), 'name': nm, 'cls': 'QTL', 'disease': 'Fusarium crown rot', 'chr': ch,
      'start': round(a, 6), 'end': round(b, 6), 'peak': round((a + b) / 2, 6),
      'ref': 'RefSeq v1.0', 'marker': '%s (peak); %s-%s (flanking)' % (peakm, f1, f2),
      'pve': '%.2f' % pve, 'stat': 'LOD %.1f' % lod, 'source': SRC, 'grade': 'A',
      'curation': 'auto',
      'note': ('QTL for Fusarium crown rot field resistance in the PAN 3471 / SST 087 doubled-haploid '
               'population, trial year %d (%s). Favourable allele from %s. Physical interval is the '
               'Chinese Spring (IWGSC RefSeq v1.0) span between the flanking DArT markers %s '
               '(%s:%s) and %s (%s:%s), taken from %s. "cg" denotes a QTL from the South African '
               'trials, "usq" from the Australian trial.'
               % (yr, 'South Africa' if '.cg-' in nm else 'Australia', parent,
                  f1, ch, format(p1, ','), f2, ch, format(p2, ','), CIT)),
      'evidence': 'Table 4: peak marker %s; flanking %s %s:%s and %s %s:%s; LOD %.1f; PVE %.2f%%'
                  % (peakm, f1, ch, format(p1, ','), f2, ch, format(p2, ','), lod, pve)})

for nm, ch in TEXT_ONLY:
    NEW.append({
      'id': nid(), 'name': nm, 'cls': 'QTL', 'disease': 'Fusarium crown rot', 'chr': ch,
      'start': None, 'end': None, 'peak': None, 'ref': 'RefSeq v1.0', 'marker': '',
      'pve': '', 'stat': '', 'source': SRC, 'grade': 'D', 'curation': 'auto',
      'note': ('QTL for Fusarium crown rot named by Martin et al. (2025) in the PAN 3471 / SST 087 '
               'doubled-haploid population. The paper states that FCR QTL were found on chromosomes '
               '%s; only those detected across multiple environments are given in Table 4, and this '
               'one appears solely in Supplementary Table S3, which is not held, so no position, LOD '
               'or PVE can be recorded.' % CHRS_TEXT)})

M['reports'] = keep + NEW
print('rebuilt:', len(NEW), 'records (%d placed)' % sum(1 for r in NEW if r['start'] is not None))
for r in NEW:
    print('  %-12s %-3s %-12s %-12s LOD %-6s PVE %-6s' %
          (r['name'], r['chr'], r['start'], r['end'], r['stat'].replace('LOD ', '') or '-', r['pve'] or '-'))
if '--write' in sys.argv:
    shutil.copy(P, '/home/claude/atlas/merged_v13f_preMartin.json')
    json.dump(M, open(P, 'w'))
    print('written. total reports', len(M['reports']))
