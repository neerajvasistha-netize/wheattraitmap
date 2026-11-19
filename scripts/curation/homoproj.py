#!/usr/bin/env python3
"""Project a physical position from one wheat subgenome onto its homoeologue,
using array probes that map to two subgenomes as synteny anchors."""
import json, re, bisect, statistics, sys, shutil
from collections import defaultdict

MP = json.load(open('/home/claude/markerpos_all.json'))

# ---- build anchor sets: probes with a v1.0 position on two homoeologues
raw = defaultdict(list)
for k, v in MP.items():
    byc = {}
    for ch, pos, asm, src in v:
        if asm != 'RefSeq v1.0':
            continue
        byc.setdefault(ch, []).append(pos)
    if len(byc) < 2:
        continue
    ks = sorted(byc)
    for a in ks:
        for b in ks:
            if a != b and a[0] == b[0] and a[1] != b[1]:
                raw[(a, b)].append((min(byc[a]), min(byc[b])))

def monotone(pts):
    """longest non-decreasing-in-target run, greedily; removes rearranged outliers"""
    pts = sorted(set(pts))
    out = []
    for s, t in pts:
        while out and out[-1][1] > t:
            out.pop()
        out.append((s, t))
    return out

ANCH = {k: monotone(v) for k, v in raw.items() if len(v) >= 30}

CHRLEN = {}
for k_, v_ in MP.items():
    for ch, pos, asm, src in v_:
        if asm == 'RefSeq v1.0':
            CHRLEN[ch] = max(CHRLEN.get(ch, 0), pos)

def project(src_chr, pos, tgt_chr, k=7, gate=True):
    """-> (projected_Mb, n_anchors_used, span_Mb) or None"""
    A = ANCH.get((src_chr, tgt_chr))
    if not A:
        return None
    xs = [a[0] for a in A]
    i = bisect.bisect_left(xs, pos)
    if gate and not (0 < i < len(A)):
        return None                      # extrapolation beyond the anchor range
    lo, hi = max(0, i - k // 2 - 1), min(len(A), i + k // 2 + 1)
    w = A[lo:hi]
    if len(w) < 2:
        return None
    if gate and not (w[0][0] <= pos <= w[-1][0]):
        return None
    # robust local slope = median pairwise slope, intercept via medians
    sl = []
    for a in range(len(w)):
        for b in range(a + 1, len(w)):
            dx = w[b][0] - w[a][0]
            if abs(dx) > 1e-6:
                sl.append((w[b][1] - w[a][1]) / dx)
    if not sl:
        return None
    m = statistics.median(sl)
    if m <= 0:
        return None
    mx = statistics.median([p[0] for p in w]); my = statistics.median([p[1] for p in w])
    est = my + m * (pos - mx)
    span = w[-1][0] - w[0][0]
    if gate:
        L = CHRLEN.get(tgt_chr, 1e9)
        if est < 0 or est > L:
            return None                  # off the end of the target chromosome
        prop = pos * L / max(CHRLEN.get(src_chr, L), 1e-6)
        if abs(est - prop) > 60:
            return None                  # grossly inconsistent with collinearity
        # local anchors must bracket closely, else the fit is unsupported
        if min(abs(pos - w[0][0]), abs(w[-1][0] - pos)) > 25:
            return None
    return round(est, 3), len(w), round(span, 1)

def selftest():
    errs = defaultdict(list)
    for (s, t), A in ANCH.items():
        if len(A) < 60:
            continue
        for idx in range(0, len(A), max(1, len(A) // 40)):
            pos, truth = A[idx]
            B = A[:idx] + A[idx + 1:]
            ANCH[(s, t)] = B
            r = project(s, pos, t, gate=False)
            ANCH[(s, t)] = A
            if r:
                errs[(s, t)].append(abs(r[0] - truth))
    allе = [e for v in errs.values() for e in v]
    allе.sort()
    n = len(allе)
    print('leave-one-out validation on %d anchor points across %d chromosome pairs'
          % (n, len(errs)))
    print('  median error %.2f Mb | 75th %.2f Mb | 90th %.2f Mb'
          % (allе[n // 2], allе[int(n * .75)], allе[int(n * .90)]))
    return allе[n // 2], allе[int(n * .90)]

if __name__ == '__main__':
    med, p90 = selftest()
    P = '/home/claude/atlas/merged_v13f.json'
    M = json.load(open(P))
    done = skip = 0
    rows = []
    for r in M['reports']:
        note = r.get('note') or ''
        if 'maps physically to' not in note:
            continue
        if r.get('peak') is not None or r.get('start') is not None:
            continue
        m = re.search(r'Marker (\S+) maps physically to ([^,]+), not to (\S+?) as reported', note)
        if not m:
            skip += 1; continue
        marker = m.group(1); tgt = r['chr']
        srcs = re.findall(r'(\d[ABD]) ([\d.]+) Mb', m.group(2))
        best = None
        for sc, sp in srcs:
            pr = project(sc, float(sp), tgt)
            if pr and (best is None or pr[1] > best[0][1]):
                best = (pr, sc, float(sp))
        if not best:
            skip += 1; continue
        (est, na, span), sc, sp = best
        r['peak'] = est
        r['ref'] = 'RefSeq v1.0'
        r['grade'] = 'D'
        r['note'] = note.rstrip('. ') + (
            ' Position shown is the homoeologous projection of that site onto %s '
            '(%.3f Mb on %s projected to %.3f Mb on %s), computed by local robust '
            'regression on %d array probes that map to both %s and %s. '
            'Inferred, not measured: leave-one-out error of this projection is '
            '%.1f Mb (median) and %.1f Mb (90th percentile).'
            % (tgt, sp, sc, est, tgt, na, sc, tgt, med, p90))
        r['evidence'] = (r.get('evidence', '') and r['evidence'] + '; ') + \
            'homoeologous projection %s %.3f Mb -> %s %.3f Mb' % (sc, sp, tgt, est)
        rows.append((r['id'], marker, sc, sp, tgt, est, r['disease']))
        done += 1
    print('projected:', done, '| could not project:', skip)
    for x in rows:
        print('  %-6s %-30s %s %8.2f -> %s %8.2f  %s' % x)
    if '--write' in sys.argv:
        shutil.copy(P, '/home/claude/atlas/merged_v13f_preHomoProj.json')
        json.dump(M, open(P, 'w'))
        un = [r for r in M['reports'] if r.get('peak') is None and r.get('start') is None]
        print('written. unplaced now', len(un))
