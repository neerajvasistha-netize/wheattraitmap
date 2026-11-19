#!/usr/bin/env python3
"""For every atlas record whose paper text we hold, look for a physical position
stated in the paper next to that locus/marker name and compare it with what the
atlas stores. Emits candidates for correction - it does not change the atlas."""
import json, re, os, sys
from collections import defaultdict

M = json.load(open('/home/claude/atlas/merged_v13f.json'))
d2t = json.load(open('/home/claude/d2t.json'))

def doi(s):
    m = re.search(r'10\.\d{4,9}/[^\s"\',;]+', s or '')
    return m.group(0).lower().rstrip('.') if m else ''

TXT = {}
for k, v in M['sources'].items():
    d = doi(v.get('url') or '')
    f = d2t.get(d) or d2t.get(d.upper())
    if isinstance(f, list):
        f = f[0] if f else None
    if f and os.path.exists(f):
        TXT[k] = f

CACHE = {}
def text(src):
    if src not in CACHE:
        try:
            CACHE[src] = re.sub(r'[ \t]+', ' ', open(TXT[src], errors='ignore').read())
        except Exception:
            CACHE[src] = ''
    return CACHE[src]

CH = r'(?:[1-7][ABD])'
# "7A:717,089,588-717,089,788"  |  "717.09 Mb" | "717,089,588 bp"
P_RANGE = re.compile(r'\b(%s)\s*[:\s]\s*([\d,]{6,15})\s*[-–to]{1,2}\s*([\d,]{6,15})' % CH)
P_MB    = re.compile(r'\b(\d{1,3}(?:\.\d{1,3})?)\s*(Mb|Mbp)\b', re.I)
P_BP    = re.compile(r'\b([\d,]{7,15})\s*bp\b', re.I)

def norm(s): return re.sub(r'[^a-z0-9]', '', str(s).lower())

def candidates(win):
    out = []
    for m in P_RANGE.finditer(win):
        a = float(m.group(2).replace(',', '')) / 1e6
        b = float(m.group(3).replace(',', '')) / 1e6
        if 0 < a < 1000 and 0 < b < 1000:
            out.append((m.group(1), a, b, 'chr:bp-bp'))
    for m in P_MB.finditer(win):
        v = float(m.group(1))
        if 0 < v < 1000:
            out.append((None, v, v, 'Mb'))
    for m in P_BP.finditer(win):
        v = float(m.group(1).replace(',', '')) / 1e6
        if 0 < v < 1000:
            out.append((None, v, v, 'bp'))
    return out

TOL = 0.5          # Mb; smaller differences are rounding, not an error
hits = []
scanned = 0
for r in M['reports']:
    src = r['source']
    if src not in TXT:
        continue
    t = text(src)
    if not t:
        continue
    keys = [r.get('name'), r.get('marker')]
    keys = [k for k in keys if k and len(str(k)) > 3]
    found = None
    for key in keys:
        for m in re.finditer(re.escape(str(key)), t):
            win = t[max(0, m.start() - 160): m.end() + 260]
            for ch, a, b, kind in candidates(win):
                if ch and r.get('chr') and ch != r['chr']:
                    continue
                if ch is None and kind != 'chr:bp-bp' and r.get('chr') is None:
                    continue
                found = (key, ch, a, b, kind)
                break
            if found:
                break
        if found:
            break
    scanned += 1
    if not found:
        continue
    key, ch, a, b, kind = found
    paper = (a + b) / 2
    cur = r.get('peak')
    if cur is None and r.get('start') is not None and r.get('end') is not None:
        cur = (r['start'] + r['end']) / 2
    if cur is None:
        hits.append(('UNPLACED', r['id'], src, r['chr'], key, kind, round(paper, 4), None))
    elif abs(cur - paper) > TOL:
        hits.append(('MISMATCH', r['id'], src, r['chr'], key, kind, round(paper, 4), round(cur, 4)))
    elif r.get('start') is None and kind == 'chr:bp-bp' and b - a > 0:
        hits.append(('NO_INTERVAL', r['id'], src, r['chr'], key, kind,
                     '%.4f-%.4f' % (a, b), round(cur, 4)))

print('records scanned against their paper text:', scanned)
by = defaultdict(int)
for h in hits:
    by[h[0]] += 1
for k, v in sorted(by.items()):
    print('  %-12s %5d' % (k, v))
bysrc = defaultdict(int)
for h in hits:
    bysrc[h[2]] += 1
print('\nacross %d source studies; worst 15:' % len(bysrc))
for k, v in sorted(bysrc.items(), key=lambda x: -x[1])[:15]:
    print('   %-30s %4d  %s' % (k, v, M['sources'][k]['short']))
json.dump(hits, open('/home/claude/reconcile_hits.json', 'w'))
print('\nwritten /home/claude/reconcile_hits.json')
