#!/usr/bin/env python3
"""
WheatNHR-Atlas: convert locus positions from IWGSC RefSeq v1.0 to RefSeq v2.1
and to every pangenome assembly, using conserved genes as anchors.

Method (gene-anchor liftover)
  1. IWGSC v1.1 HC genes give anchor positions on RefSeq v1.0.
  2. Each target assembly supplies the position of the same genes:
       - RefSeq v2.1: the IWGSC v1.1 -> v2.1 gene ID mapping + v2.1 GFF3
       - pangenome cultivars: Ensembl Plants GFF3, whose genes carry the
         projected IWGSC v1.1 gene ID (TraesCS..02G..) in their attributes.
  3. For each locus, anchors within the interval +/- FLANK Mb that land on the
     same chromosome are kept; anchors off the local collinear trend are
     dropped (median offset +/- OUTLIER Mb).
  4. Locus start/end are moved with the offset of the nearest anchor; if the
     local anchor order is reversed the locus is flagged "inverted".
  5. Each result gets a quality flag: ok, few_anchors, wide_gap, inverted,
     other_chromosome, unmapped.

Usage
  python3 liftover_atlas.py --atlas input/atlas_core.tsv --refs refs --pan pangenome --out output
Only Python 3 standard library is needed.
"""
import argparse, csv, glob, gzip, io, os, re, statistics, sys, zipfile
from collections import defaultdict

FLANK = 5.0        # Mb searched either side of the locus for anchors
OUTLIER = 10.0     # Mb allowed from the median offset
MIN_ANCHORS = 3
MAX_GAP = 2.0      # Mb: distance from a locus end to its nearest anchor before flagging

V11 = re.compile(r"TraesCS(\d[ABD]|U)02G\d{6}")
V21 = re.compile(r"TraesCS(\d[ABD]|U)03G\d{7}")
CHR = re.compile(r"(?:chr)?([1-7][ABD]|Un)$", re.I)


def open_any(path):
    """Yield text lines from .gff3, .gff3.gz, .csv, .txt, or every matching member of a .zip."""
    if path.endswith(".zip"):
        with zipfile.ZipFile(path) as z:
            for name in z.namelist():
                low = name.lower()
                if low.endswith((".gff3", ".gff", ".csv", ".txt", ".tsv", ".gff3.gz")):
                    data = z.read(name)
                    if low.endswith(".gz"):
                        data = gzip.decompress(data)
                    for line in io.TextIOWrapper(io.BytesIO(data), encoding="utf-8", errors="replace"):
                        yield name, line
        return
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            yield path, line


def norm_chr(c):
    c = c.strip()
    m = CHR.search(c.split("_")[-1]) or CHR.search(c)
    return m.group(1).upper().replace("UN", "Un") if m else None


def read_gff_genes(path, idre):
    """Gene ID -> (chr, midpoint Mb). Works for GFF3 and GTF: the span of every feature
    carrying the ID (gene, mRNA, exon, CDS) is merged, so GTF files without gene lines work."""
    span = {}
    for fname, line in open_any(path):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        ch = norm_chr(f[0])
        if not ch:
            continue
        attr = f[8]
        if "previous_id=" in attr:          # IWGSC v2.1: ignore the old ID field
            attr = attr.split("previous_id=")[0]
        m = idre.search(attr)
        if not m or "LC" in attr[m.end():m.end() + 2]:
            continue
        g = m.group(0)
        a, b = int(f[3]), int(f[4])
        if g in span:
            c, lo, hi = span[g]
            if c == ch:
                span[g] = (c, min(lo, a), max(hi, b))
        else:
            span[g] = (ch, a, b)
    return {g: (c, (lo + hi) / 2e6) for g, (c, lo, hi) in span.items()}


def read_v11_to_v21(path):
    """v1.1 gene ID -> v2.1 gene ID from any IWGSC mapping table inside path."""
    mp = {}
    for fname, line in open_any(path):
        if "\tgene\t" in line:
            m11 = re.search(r"previous_id=(TraesCS\w+?02G\d{6})(?![\dL])", line)
            m21 = V21.search(line)
            if m11 and m21 and "primconf=HC" in line:
                mp[m11.group(1)] = m21.group(0)
            continue
        if "LC" in line or line.startswith(">"):
            continue
        a, b = V11.search(line), V21.search(line)
        if a and b:
            mp.setdefault(a.group(0), b.group(0))
    return mp


def find(pattern_list, folder):
    hits = []
    for p in pattern_list:
        hits += glob.glob(os.path.join(folder, p))
    return sorted(set(hits))


def lift(ch, s, e, ref_anchor, tgt):
    """Return dict with target chr/start/end and QC for one locus."""
    near = [(pos, g) for g, (c, pos) in ref_anchor.get(ch, []) if s - FLANK <= pos <= e + FLANK]
    pairs, other = [], defaultdict(list)
    for pos, g in near:
        t = tgt.get(g)
        if not t:
            continue
        if t[0] == ch:
            pairs.append((pos, t[1]))
        else:
            other[t[0]].append((pos, t[1]))
    moved = None
    if len(pairs) < MIN_ANCHORS and other:
        best = max(other, key=lambda c: len(other[c]))
        if len(other[best]) >= MIN_ANCHORS and len(other[best]) > 3 * len(pairs):
            moved, pairs, ch = best, other[best], best   # e.g. translocated segment
    if not pairs:
        return {"flag": "unmapped", "n": 0}
    pairs.sort()
    # orientation from rank correlation sign of anchor order
    slope = 1
    if len(pairs) >= 3:
        inc = sum(1 for i in range(1, len(pairs)) if pairs[i][1] >= pairs[i - 1][1])
        slope = 1 if inc >= (len(pairs) - 1) / 2 else -1
    offs = [t - slope * r for r, t in pairs]
    med = statistics.median(offs)
    keep = [(r, t) for (r, t), o in zip(pairs, offs) if abs(o - med) <= OUTLIER]
    if not keep:
        keep = pairs

    def move(x):
        r, t = min(keep, key=lambda p: abs(p[0] - x))
        return t + slope * (x - r), abs(r - x)

    ts, gs = move(s)
    te, ge = move(e)
    lo, hi = sorted((ts, te))
    flag = "ok"
    if moved:
        return {"chr": ch, "start": round(min(ts, te), 3), "end": round(max(ts, te), 3), "flag": "other_chromosome",
                "n": len(keep), "gap": round(max(gs, ge), 3)}
    if slope < 0:
        flag = "inverted"
    if len(keep) < MIN_ANCHORS:
        flag = "few_anchors"
    elif max(gs, ge) > MAX_GAP:
        flag = "wide_gap"
    return {"chr": ch, "start": round(lo, 3), "end": round(hi, 3), "flag": flag,
            "n": len(keep), "gap": round(max(gs, ge), 3)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--atlas", required=True)
    ap.add_argument("--refs", required=True)
    ap.add_argument("--pan", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    v11_files = find(["IWGSC_v1.1_HC*.gtf", "IWGSC_v1.1_HC*.gff3"], a.refs) or find(["*refseqv1.1*.zip"], a.refs)
    v21_files = find(["*refseqv2.1*_HC.gff3", "*refseqv2.1*IDmapping.csv"], a.refs) or find(["*refseqv2.1*.zip"], a.refs)
    if not v11_files or not v21_files:
        sys.exit("Missing IWGSC files in %s (need v1.1 genes and v2.1 annotation zips)" % a.refs)

    print("reading v1.1 genes:", v11_files)
    ref = {}
    for f in v11_files:
        ref.update(read_gff_genes(f, V11))
    ref_anchor = defaultdict(list)
    for g, (c, pos) in ref.items():
        ref_anchor[c].append((g, (c, pos)))
    print("  v1.1 anchor genes:", len(ref))

    print("reading v2.1 annotation:", v21_files)
    v21pos, v11to21 = {}, {}
    for f in v21_files:
        v21pos.update(read_gff_genes(f, V21))
        v11to21.update(read_v11_to_v21(f))
    tgt_v21 = {g: v21pos[n] for g, n in v11to21.items() if n in v21pos}
    print("  v1.1 genes placed on v2.1:", len(tgt_v21))

    targets = {"RefSeq_v2.1": tgt_v21}
    for f in find(["anchors_*.tsv"], a.pan):
        name = os.path.basename(f)[len("anchors_"):-4]
        name = re.sub(r"^Triticum_aestivum_", "", name)
        d = {}
        for row in csv.DictReader(open(f), delimiter="\t"):
            d[row["v11_gene"]] = (row["chr"], float(row["mid_Mb"]))
        print("  pangenome", name, "anchors:", len(d))
        if len(d) > 500:
            targets[name] = d

    # reverse map for reports given on v2.x: v2.1 -> v1.0 frame
    rev_anchor = defaultdict(list)
    for g, (c, pos) in tgt_v21.items():
        rev_anchor[c].append((g, (c, pos)))

    rows = list(csv.DictReader(open(a.atlas), delimiter="\t"))
    long_rows = []
    wide = []
    for r in rows:
        ch = r["chromosome"].strip()
        s = r["start_Mb"] or r["peak_Mb"]
        e = r["end_Mb"] or r["peak_Mb"]
        base = {"locus_name": r["locus_name"], "disease": r["disease"], "source_key": r["source_key"],
                "chromosome": ch, "ref_version": r["ref_version"], "start_Mb": s, "end_Mb": e}
        if not s or not ch:
            wide.append(dict(base, note="no physical position in source"))
            continue
        s, e = float(s), float(e)
        rv = r["ref_version"]
        if rv.startswith("RefSeq v2"):
            back = lift(ch, s, e, rev_anchor, ref)
            if "start" not in back:
                wide.append(dict(base, note="v2.x report could not be placed on v1.0: " + back["flag"]))
                continue
            s1, e1 = back["start"], back["end"]
            base.update(v10_start=s1, v10_end=e1, v10_flag=back["flag"])
        else:
            s1, e1 = s, e
        rec = dict(base)
        for tname, tg in targets.items():
            if tname == "RefSeq_v2.1" and rv.startswith("RefSeq v2"):
                res = {"chr": ch, "start": s, "end": e, "flag": "as_reported" if rv == "RefSeq v2.1" else "v2.0_as_reported", "n": ""}
            else:
                res = lift(ch, s1, e1, ref_anchor, tg)
            long_rows.append({"locus_name": r["locus_name"], "assembly": tname, "chr": res.get("chr", ""),
                              "start_Mb": res.get("start", ""), "end_Mb": res.get("end", ""),
                              "flag": res["flag"], "anchors": res.get("n", ""), "nearest_anchor_Mb": res.get("gap", "")})
            rec[tname + "_chr"] = res.get("chr", "")
            rec[tname + "_start"] = res.get("start", "")
            rec[tname + "_end"] = res.get("end", "")
            rec[tname + "_flag"] = res["flag"]
        wide.append(rec)

    with open(os.path.join(a.out, "liftover_long.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(long_rows[0].keys()), delimiter="\t")
        w.writeheader(); w.writerows(long_rows)
    keys = []
    for r in wide:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(os.path.join(a.out, "liftover_wide.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, delimiter="\t")
        w.writeheader(); w.writerows(wide)
    flags = defaultdict(lambda: defaultdict(int))
    for r in long_rows:
        flags[r["assembly"]][r["flag"]] += 1
    with open(os.path.join(a.out, "liftover_summary.tsv"), "w") as fh:
        allf = sorted({f for d in flags.values() for f in d})
        fh.write("assembly\t" + "\t".join(allf) + "\n")
        for t in targets:
            fh.write(t + "\t" + "\t".join(str(flags[t][f]) for f in allf) + "\n")
    print("done:", len(long_rows), "locus x assembly placements written to", a.out)


if __name__ == "__main__":
    main()
