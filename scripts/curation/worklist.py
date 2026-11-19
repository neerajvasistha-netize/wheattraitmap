#!/usr/bin/env python3
"""Build the download/curation worklist for WheatNHR-Atlas."""
import json, re, csv, os
import openpyxl
from openpyxl.styles import Font, PatternFill, Alignment
from collections import Counter, defaultdict

M = json.load(open('/home/claude/atlas/merged_v13f.json'))
U = json.load(open('/home/claude/atlas/papers_uploaded.json'))
AUD = json.load(open('/home/claude/audit_disease_missing.json'))

def doi(s):
    m = re.search(r'10\.\d{4,9}/[^\s"\',;]+', s or '')
    return m.group(0).lower().rstrip('.') if m else ''

# RGU library links
RGU = {}
with open('/mnt/user-data/uploads/a29d7657-WheatNHR_RGU_Library_Links.csv') as f:
    for row in csv.DictReader(f):
        d = doi(row.get('DOI'))
        if d:
            RGU[d] = (row.get('Direct DOI Link', ''), row.get('RGU Library Proxy Link', ''),
                      row.get('Disease Category', ''), row.get('Short Reference', ''))

def links(d):
    if d in RGU:
        return RGU[d][0], RGU[d][1]
    if d:
        return 'https://doi.org/' + d, 'https://ezproxy.rgu.ac.in/login?url=https://doi.org/' + d
    return '', ''

# ---------------- rows
rows = {}   # doi or key -> dict
def put(k, **kw):
    r = rows.setdefault(k, {'reasons': [], 'expected': 0, 'need': set(), 'disease': set(),
                            'short': '', 'doi': '', 'in_atlas': 0})
    for a, b in kw.items():
        if a == 'reason': r['reasons'].append(b)
        elif a == 'need': r['need'].add(b)
        elif a == 'disease':
            for x in (b if isinstance(b, (set, list)) else [b]):
                if x: r['disease'].add(x)
        elif a == 'expected': r['expected'] = max(r['expected'], b)
        else: r[a] = b or r.get(a) or ''

# A. under-extracted studies already in the atlas
nrec = Counter(r['source'] for r in M['reports'])
for k, v in AUD.items():
    key, short = k.split('|', 1)
    d = doi(M['sources'][key].get('url') or '')
    put(d or key, short=short, doi=d, in_atlas=nrec.get(key, 0),
        reason='under-extracted: %d disease loci named in the paper are not in the atlas' % len(v),
        need='supplementary tables', expected=len(v),
        disease={x.split(' (')[0] for _, x in v})

# B. papers reviewed whose loci sit in figures/supplementary, and unreadable PDFs
for p in U:
    if p.get('in_atlas'):
        continue
    nt = (p.get('notes') or '')
    d = doi(p.get('link'))
    if 'only in figures or supplementary' in nt:
        put(d or p['id'], short=(p.get('authors') or p.get('title') or '')[:80], doi=d,
            reason='loci present but reported only in figures or supplementary files',
            need='supplementary tables', disease=p.get('disease'), expected=0)
    elif 'not readable' in nt:
        put(d or p['id'], short=(p.get('authors') or p.get('title') or '')[:80], doi=d,
            reason='PDF text could not be read; a clean copy is needed',
            need='full text (clean PDF)', disease=p.get('disease'), expected=0)

# ---------------- priority
PRI = lambda r: (0 if r['expected'] >= 10 else 1 if r['expected'] > 0 else 2, -r['expected'])
out = sorted(rows.values(), key=PRI)

wb = openpyxl.Workbook()
HDR = Font(name='Arial', bold=True, color='FFFFFF'); FILL = PatternFill('solid', fgColor='2A5D8F')
F = Font(name='Arial'); WRAP = Alignment(wrap_text=True, vertical='top')

ws = wb.active; ws.title = 'Download worklist'
LINK = Font(name='Arial', color='0563C1', underline='single')
ws.append(['Priority', 'Open', 'Study', 'Disease(s)',
           'Loci already in atlas', 'Expected additional loci', 'Why it is needed',
           'What to download', 'DOI'])
for i, r in enumerate(out, 1):
    dl, pl = links(r['doi'])
    ws.append([i, 'Open' if dl else '', r['short'],
               '; '.join(sorted(r['disease']))[:80], r['in_atlas'], r['expected'] or '',
               ' | '.join(r['reasons']), '; '.join(sorted(r['need'])), r['doi']])
    row = ws.max_row
    if dl:
        for col in (2, 3):
            c = ws.cell(row=row, column=col); c.hyperlink = dl; c.font = LINK

# ---------------- critical analysis by disease
miss = Counter()
for k, v in AUD.items():
    for n, d in v:
        miss[d.split(' (')[0]] += 1
st = defaultdict(set); rec = Counter(); pl = Counter(); mkc = Counter(); pve = Counter()
for r in M['reports']:
    d = r['disease']; rec[d] += 1; st[d].add(r['source'])
    if r.get('peak') is not None or r.get('start') is not None: pl[d] += 1
    if (r.get('marker') or '').strip(): mkc[d] += 1
    if str(r.get('pve') or '').strip(): pve[d] += 1

def verdict(d):
    n, s = rec[d], len(st[d])
    v = []
    if s == 1 and n > 100:
        v.append('all records come from a single compilation; no primary-source attribution')
    if s >= 20 and n / s < 7:
        v.append('loci per study (%.1f) is low for this many studies - likely under-extracted' % (n / s))
    if n < 20:
        v.append('thin coverage; more mapping studies needed')
    if 100 * pl[d] / n < 80:
        v.append('only %.0f%% of records carry a physical position' % (100 * pl[d] / n))
    if 100 * mkc[d] / n < 75:
        v.append('only %.0f%% of records carry a marker name' % (100 * mkc[d] / n))
    if miss.get(d, 0) >= 25:
        v.append('%d named loci in held papers are still missing' % miss[d])
    return '; '.join(v) if v else 'coverage and completeness look sound'

ws2 = wb.create_sheet('Critical analysis by disease')
ws2.append(['Disease', 'Records', 'Source studies', 'Loci per study', 'With position %',
            'With marker %', 'With PVE %', 'Named loci still missing', 'Assessment'])
for d in sorted(rec, key=lambda x: -rec[x]):
    n = rec[d]
    ws2.append([d, n, len(st[d]), round(n / len(st[d]), 1), round(100 * pl[d] / n, 1),
                round(100 * mkc[d] / n, 1), round(100 * pve[d] / n, 1), miss.get(d, 0), verdict(d)])

ws3 = wb.create_sheet('Diseases not yet in atlas')
ws3.append(['Disease', 'Evidence that mappable loci exist', 'What to download'])
for d, note, what in [
 ('Root-lesion nematode (Pratylenchus thornei / P. neglectus)',
  '7 QRlnn loci named in papers already held (Zwart et al. 2010 and others); Rlnn1 is a designated gene',
  'Zwart et al. 2010 and the primary Rlnn1 mapping papers, with supplementary tables'),
 ('Rhizoctonia root rot / bare patch',
  '2 QRc loci named in papers already held', 'the primary mapping paper and its supplementary'),
 ('Flag smut (Urocystis agropyri)', 'no mapping study with physical positions found',
  'search for any QTL/GWAS study; none located so far'),
 ('Cephalosporium stripe', 'no mapping study with physical positions found',
  'search; none located so far'),
 ('Take-all (Gaeumannomyces tritici)',
  'two QTL reported in Moughan 2017 PhD thesis but chromosome and position withheld as Syngenta confidential',
  'nothing available; the positions are not published'),
 ('Snow mould, black point, ergot', 'no mapping study with physical positions found',
  'search; none located so far'),
 ('Soil-borne wheat mosaic virus (SBWMV)',
  'Sbm1 is a designated gene; a BYDV/SBWMV GWAS exists (The Plant Genome 2025) but is not open access',
  'the BYDV/SBWMV GWAS paper and the primary Sbm1 mapping paper'),
]:
    ws3.append([d, note, what])

ws4 = wb.create_sheet('Missing loci detail')
ws4.append(['Source key', 'Study', 'Locus name', 'Disease (from locus name)', 'Open', 'DOI'])
for k, v in sorted(AUD.items(), key=lambda x: -len(x[1])):
    key, short = k.split('|', 1)
    d = doi(M['sources'][key].get('url') or '')
    dl, pl = links(d)
    for n, dis in sorted(v):
        ws4.append([key, short, n, dis, 'Open' if dl else '', d])
        if dl:
            c = ws4.cell(row=ws4.max_row, column=5)
            c.hyperlink = dl; c.font = Font(name='Arial', color='0563C1', underline='single')

ws5 = wb.create_sheet('Method and limits')
for line in [
 'WheatNHR-Atlas - download and curation worklist',
 'Atlas state at generation: %d records, %d source studies with loci, %d diseases, %d records without a physical position.'
   % (len(M['reports']), len({r['source'] for r in M['reports']}), len(rec),
      sum(1 for r in M['reports'] if r.get('peak') is None and r.get('start') is None)),
 '',
 'Sheet 1 - Download worklist. One row per paper to obtain, highest expected yield first.',
 '  Priority 1 onwards: studies where 10 or more disease loci named in the paper are absent from the atlas.',
 '  Then: studies with 1-9 such loci. Then: papers already reviewed whose loci appear only in figures or',
 '  supplementary files, and the few whose PDF text could not be read.',
 '  "Expected additional loci" is blank where the number cannot be estimated from the text held.',
 '',
 'Sheet 2 - Critical analysis by disease. The assessment column applies these tests:',
 '  - a single source study behind more than 100 records means the records come from a compilation,',
 '    not from the primary literature, so no original study is credited;',
 '  - fewer than 7 loci per study across 20 or more studies suggests systematic under-extraction',
 '    (for comparison, tan spot averages 13.5 and spot blotch 11.3);',
 '  - fewer than 20 records means coverage is thin and more mapping studies are needed;',
 '  - thresholds of 80% for physical position, 75% for marker name.',
 '',
 'Sheet 3 - Diseases with no records yet, and whether any mappable loci exist to obtain.',
 '',
 'Sheet 4 - The individual locus names behind the counts in sheet 1.',
 '',
 'How the missing-loci estimate was produced, and its limits:',
 '  Each held full text was scanned for locus designations matching standard wheat QTL/gene naming,',
 '  names already in the atlas for that study were subtracted, and the remainder classified by trait',
 '  code; agronomic codes (Ht, Hd, Ph, Da, Flt and similar) were excluded as outside scope.',
 '  Meta-analysis papers are excluded: the QTL they name are their input loci from other studies.',
 '  This is an estimate, not a verified list. A name may appear because the authors cite it from',
 '  another study. The disease is inferred from the name prefix, not read from the paper. Some names',
 '  will be parsing artefacts. Every row needs the source table opened and read before anything is',
 '  added to the atlas.',
 '',
 'Generated by Claude for N.K. Vasistha, Rajiv Gandhi University, 6 November 2026.']:
    ws5.append([line])

for s in (ws, ws2, ws3, ws4):
    for c in s[1]:
        c.font = HDR; c.fill = FILL; c.alignment = WRAP
    for row in s.iter_rows(min_row=2):
        for c in row:
            if c.hyperlink is None: c.font = F
    s.freeze_panes = 'A2'
for s, w in [(ws, [8, 10, 40, 28, 12, 13, 56, 24, 30]),
             (ws2, [30, 9, 13, 13, 14, 13, 11, 14, 70]),
             (ws3, [44, 72, 56]), (ws4, [20, 34, 24, 28, 19, 30])]:
    for i, x in enumerate(w, 1):
        s.column_dimensions[openpyxl.utils.get_column_letter(i)].width = x
for s in (ws, ws2, ws3):
    for row in s.iter_rows(min_row=2):
        row[-1].alignment = WRAP
ws5.column_dimensions['A'].width = 108
for r_ in ws5.iter_rows():
    for c in r_: c.font = F

os.makedirs('/mnt/user-data/outputs', exist_ok=True)
wb.save('/mnt/user-data/outputs/WheatNHR_Atlas_papers_to_download.xlsx')
print('download worklist rows:', len(out))
print('  with an estimated yield:', sum(1 for r in out if r['expected']))
print('  total estimated additional loci:', sum(r['expected'] for r in out))
print('diseases analysed:', len(rec))
