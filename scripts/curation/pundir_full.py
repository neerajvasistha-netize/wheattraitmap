#!/usr/bin/env python3
"""Full extraction of Pundir et al. 2023 (BMC Plant Biol 23:529) supplementary:
Table S1 digenic epistatic QTL (all named loci) and Table S2 previously reported QTL."""
import openpyxl, json, re, sys, shutil
from collections import defaultdict

XL = '/mnt/user-data/uploads/9744cba3-12870_2023_4526_MOESM1_ESM.xlsx'
P = '/home/claude/atlas/merged_v13f.json'
M = json.load(open(P))
wb = openpyxl.load_workbook(XL, data_only=True)

have = {r['id'] for r in M['reports']}
nxt = max(int(r['id'][1:]) for r in M['reports'] if re.fullmatch(r'V\d+', r['id'])) + 1
exist = {(r['name'], r['source']) for r in M['reports']}
MP = json.load(open('/home/claude/markerpos_all.json'))
def mk(s): return re.sub(r'[^a-z0-9]', '', str(s).lower())

NEW = []
def add(**kw):
    global nxt
    rid = 'V%04d' % nxt
    while rid in have:
        nxt += 1; rid = 'V%04d' % nxt
    have.add(rid); nxt += 1
    r = {'id': rid, 'start': None, 'end': None, 'peak': None, 'ref': 'RefSeq v1.0',
         'marker': '', 'pve': '', 'stat': '', 'grade': 'D', 'curation': 'auto'}
    r.update(kw)
    NEW.append(r); return rid

# ---------------- Table S1: every named QTL in the digenic epistatic analysis
ws = wb['Table S1']
QN = re.compile(r'^(QCcnr\.ccsu-([1-7][ABD])\.\d+)\s*\((\d+)\)$')
loci = defaultdict(lambda: {'cm': set(), 'lod': [], 'pve': [], 'env': set(), 'part': []})
for row in ws.iter_rows(min_row=3, values_only=True):
    if not row or not row[0]:
        continue
    m1, m2, lp, ab, env = row[0], row[1], row[2], row[3], row[4]
    g = re.match(r'LOD', str(lp or ''))
    mlp = re.match(r'^([\d.]+)\s*\(([\d.]+)\)$', str(lp or '').strip())
    if not mlp:
        continue
    lod, pve = float(mlp.group(1)), float(mlp.group(2))
    pair = []
    for cell in (m1, m2):
        q = QN.match(str(cell or '').strip())
        if q:
            pair.append((q.group(1), q.group(2), int(q.group(3))))
    for nm, ch, cm in pair:
        d = loci[(nm, ch)]
        d['cm'].add(cm); d['lod'].append(lod); d['pve'].append(pve)
        d['env'].add(str(env or '').strip())
        other = [p[0] for p in pair if p[0] != nm]
        if other:
            d['part'].append('%s (%s, LOD %.2f, PVE %.2f%%)' % (other[0], env, lod, pve))

ENVN = {'R1': 'repeat 1', 'R2': 'repeat 2', 'CD': 'combined data'}
added = skipped = 0
for (nm, ch), d in sorted(loci.items()):
    if (nm, 'Pundir2023BMC') in exist:
        skipped += 1
        continue
    cm = sorted(d['cm'])
    add(name=nm, cls='QTL', disease='Cereal cyst nematode', chr=ch,
        pve='%.2f-%.2f' % (min(d['pve']), max(d['pve'])) if len(set(d['pve'])) > 1
            else '%.2f' % d['pve'][0],
        stat='LOD %.2f-%.2f' % (min(d['lod']), max(d['lod'])) if len(set(d['lod'])) > 1
            else 'LOD %.2f' % d['lod'][0],
        source='Pundir2023BMC', grade='D',
        note=('QTL for Heterodera avenae resistance detected in the digenic epistatic analysis of '
              'the HUW 468 x C 306 RIL population (Pundir et al. 2023, Table S1), in %s. '
              'Genetic position %s cM on the study linkage map; the study reports no physical '
              'position for this locus and its flanking markers are study-internal WSNP codes, '
              'so no coordinate is assigned. Interacting partner(s): %s.'
              % (', '.join(ENVN.get(e, e) for e in sorted(d['env']) if e),
                 '/'.join(str(c) for c in cm), '; '.join(d['part'][:4]))),
        evidence='epistatic QTL, %s cM' % '/'.join(str(c) for c in cm))
    added += 1

# ---------------- Table S2: previously reported QTL not already in the atlas
REF = {
 1: 'Williams KJ, Willsmore KL, Olson S, Matic M, Kuchel H (2006) Theor Appl Genet 112:1480-1486',
 2: 'Singh K, Chhuneja P, Singh I, Sharma SK, Garg T, Garg M, Keller B, Dhaliwal HS (2010) Euphytica 176:213-222',
 3: 'Jayatilake DV, Tucker EJ, Brueggemann J, Lewis J, Garcia M, Dreisigacker S, Hayden MJ, Chalmers K, Mather DE (2015) Mol Breed 35:1-12',
 4: 'Cui L, Qiu D, Sun L, Sun Y, Ren Y, Zhang H, Li J, et al. (2020) Phytopathology 110:472-482',
 5: 'Al-Ateeq TK, Al-Doss AA, Al-Hazmi AS, Ghazy AI, Dawabah AM, Motawei MI (2021) Cereal Research Communications 49:1-7',
}
ws = wb['Table S2']
s2 = 0
for row in ws.iter_rows(min_row=3, values_only=True):
    if not row or not row[2] or str(row[0] or '').startswith('Reference'):
        continue
    pop, meth, nm, flank, closest, pve, ref = [str(x or '').strip() for x in row[:7]]
    nm = nm.replace(' (2A)', '').strip()
    if not re.match(r'^(Q|Cre)', nm):
        continue
    if (nm, 'Pundir2023BMC') in exist or any(x['name'] == nm for x in NEW):
        continue
    ch = (re.search(r'-([1-7][ABD])', nm) or re.search(r'([1-7][ABD])', row[2] or '') or [None, None])[1]
    if not ch:
        continue
    pos = None; src = ''
    for cand in re.split(r'[,–\-\s]+', (closest or flank)):
        cand = cand.strip().lstrip('X')
        e = MP.get(mk(cand))
        if e:
            hit = [x for x in e if x[0] == ch]
            if hit:
                pos = round(hit[0][1], 4); src = cand; break
    rnum = int(re.search(r'\d', ref).group(0)) if re.search(r'\d', ref) else None
    add(name=nm, cls='QTL', disease='Cereal cyst nematode', chr=ch, peak=pos,
        marker=closest or flank, pve=pve,
        stat='%s; %s' % (meth, pop[:60]), source='Pundir2023BMC',
        grade='C' if pos else 'D',
        note=('Previously reported QTL for Heterodera avenae resistance, compiled by Pundir et al. '
              '(2023) Table S2. Population: %s. Mapping method: %s. Flanking/associated markers: %s. '
              'Original study: %s.%s'
              % (pop, meth, flank, REF.get(rnum, 'see source'),
                 '' if pos else ' No physical position could be assigned from the reported markers.')),
        evidence=('physical position of %s from the atlas marker dictionary' % src) if pos else '')
    s2 += 1

print('Table S1 named QTL added:', added, '| already present, skipped:', skipped)
print('Table S2 previously reported QTL added:', s2)
print('total new records:', len(NEW))
ch = defaultdict(int)
for r in NEW:
    ch[r['chr']] += 1
print('by chromosome:', dict(sorted(ch.items())))
print('placed:', sum(1 for r in NEW if r['peak'] is not None))
for r in NEW[:3] + NEW[-2:]:
    print('  ', r['id'], r['name'], r['chr'], r['peak'], r['pve'], r['stat'][:30])

if '--write' in sys.argv:
    shutil.copy(P, '/home/claude/atlas/merged_v13f_prePundirFull.json')
    M['reports'].extend(NEW)
    json.dump(M, open(P, 'w'))
    print('written. total reports', len(M['reports']))
