"""Download each paper's supplementary tables from the PMC open-access dataset, find wheat gene IDs
with fold-change columns, keep the DEG rows, delete the download. Output: out/<PMCID>.json"""
import csv, io, json, os, re, sys, zipfile, gzip, tempfile, urllib.request, warnings, concurrent.futures as cf
import pandas as pd
warnings.filterwarnings("ignore")
B = "https://pmc-oa-opendata.s3.amazonaws.com/"
GID = re.compile(r"(TraesCS(?:[1-7][ABD]|U)0[123]G\d{6,7}(?:LC)?)(?![0-9])", re.I)
FC = re.compile(r"log\s*2?\s*[\(_ .-]*\s*(fold[\s_.-]*change|fc|ratio)|^log[\s_]*fc$|log2fc|logfc|log2\(fc\)|log2_?fold|fold[\s_.-]*change|^fc$|log2ratio", re.I)
PADJ = re.compile(r"padj|p\.?adj|fdr|q[\s_.-]*val|adj[\s_.]*p|corrected", re.I)
P = re.compile(r"^p[\s_.-]*val|^pvalue$|^p$", re.I)
SUP = re.compile(r"\.(xlsx|xlsm|xls|csv|tsv|txt|zip|gz|docx)$", re.I)
os.makedirs("out", exist_ok=True)

def tables_from(name, data):
    n = name.lower()
    try:
        if n.endswith((".xlsx", ".xlsm")):
            for sh, df in pd.read_excel(io.BytesIO(data), sheet_name=None, header=None, engine="openpyxl").items(): yield sh, df
        elif n.endswith(".xls"):
            for sh, df in pd.read_excel(io.BytesIO(data), sheet_name=None, header=None, engine="xlrd").items(): yield sh, df
        elif n.endswith((".csv", ".tsv", ".txt")):
            txt = data.decode("utf-8", "ignore")
            if not GID.search(txt[:5_000_000]): return
            sep = "\t" if txt.count("\t") > txt.count(",") else ","
            yield os.path.basename(name), pd.read_csv(io.StringIO(txt), sep=sep, header=None, dtype=str, on_bad_lines="skip", engine="python")
        elif n.endswith(".zip"):
            z = zipfile.ZipFile(io.BytesIO(data))
            for m in z.namelist():
                if SUP.search(m) and not m.startswith("__MACOSX") and z.getinfo(m).file_size < 300e6:
                    for sh, df in tables_from(m, z.read(m)): yield f"{os.path.basename(m)}:{sh}", df
        elif n.endswith(".gz"):
            yield from tables_from(n[:-3], gzip.decompress(data))
        elif n.endswith(".docx"):
            import docx
            d = docx.Document(io.BytesIO(data))
            for i, t in enumerate(d.tables):
                rows = [[c.text for c in r.cells] for r in t.rows]
                if rows: yield f"docx table {i+1}", pd.DataFrame(rows)
    except Exception as e:
        if n.endswith((".xls", ".xlsx", ".xlsm")) and not n.endswith(".txt"):
            try:
                yield from tables_from(name + ".txt", data); return
            except Exception: pass
        print("  read fail", name, str(e)[:80], file=sys.stderr)

def num(v):
    try:
        x = float(str(v).replace(",", "").strip())
        return x if x == x and abs(x) != float("inf") else None
    except: return None

def parse(sheet, df):
    """Find header row, gene column, fold-change columns. Returns list of contrasts."""
    if df.empty or df.shape[0] < 3: return []
    df = df.astype(str)
    # gene column = the column with most gene IDs
    best, bc = 0, None
    sample = df.head(3000)
    for c in df.columns:
        k = sample[c].str.contains(GID, na=False).sum()
        if k > best: best, bc = k, c
    if best < 5: return []
    # header row = last row above first gene row that has text
    first = sample.index[sample[bc].str.contains(GID, na=False)][0]
    hdr_i = None
    for i in range(first - 1, max(-1, first - 8), -1):
        row = df.loc[i]
        if (row.str.strip().replace("nan", "") != "").sum() >= 2: hdr_i = i; break
    hdr = [str(x).strip() for x in df.loc[hdr_i]] if hdr_i is not None else [str(c) for c in df.columns]
    # merged two-line headers: prepend the row above when it has a group label
    grp = [str(x).strip() for x in df.loc[hdr_i - 1]] if hdr_i is not None and hdr_i - 1 >= 0 else [""] * len(hdr)
    g_ff, last = [], ""
    for x in grp:
        last = x if x and x != "nan" else last; g_ff.append(last)
    body = df.loc[first:]
    out = []
    fc_cols = [j for j, h in enumerate(hdr) if FC.search(h) and not PADJ.search(h)]
    for j in fc_cols:
        col = body.iloc[:, j]
        vals = col.map(num)
        if vals.notna().sum() < 5: continue
        # nearest padj / p column to the right (within 4 columns, before the next FC column)
        pj = None
        nxt = [k for k in fc_cols if k > j]
        lim = min(nxt[0] if nxt else len(hdr), j + 5)
        for k in range(j + 1, lim):
            if PADJ.search(hdr[k]): pj = k; break
        if pj is None:
            for k in range(j + 1, lim):
                if P.search(hdr[k]): pj = k; break
        if pj is None:  # significance column just left of the fold change
            prv = [k for k in fc_cols if k < j]
            for k in range(j - 1, max(prv[-1] if prv else -1, j - 4), -1):
                if PADJ.search(hdr[k]): pj = k; break
            if pj is None:
                for k in range(j - 1, max(prv[-1] if prv else -1, j - 4), -1):
                    if P.search(hdr[k]): pj = k; break
        hdr_txt = hdr[j]
        if "fold" in hdr_txt.lower() and "log" not in hdr_txt.lower():
            # plain fold change -> log2
            v2 = vals.map(lambda x: None if x is None or x <= 0 else __import__("math").log2(x))
            if v2.notna().sum() >= 5 and (vals.dropna() > 0).all(): vals = v2; hdr_txt += " (converted to log2)"
            else: continue
        genes = body[bc].str.extract(GID, expand=False)
        recs = {}
        for g, v, pv in zip(genes, vals, body.iloc[:, pj] if pj is not None else [None] * len(genes)):
            if not isinstance(g, str) or v is None: continue
            g = g.upper().replace("TRAESCS", "TraesCS")
            if g.endswith("LC"): continue
            p = num(pv) if pv is not None else None
            if g not in recs or abs(v) > abs(recs[g][0]): recs[g] = (round(v, 3), p)
        if len(recs) >= 5:
            grp_names = []
            for h in hdr:
                m = re.match(r"^(.+?)[-_ .]?\d{0,2}[-_ .](count|fpkm|tpm|rpkm|readcount|reads)$", h.strip(), re.I) or re.match(r"^(count|fpkm|tpm|rpkm)[-_ .](.+?)[-_ .]?\d{0,2}$", h.strip(), re.I)
                if m:
                    gname = m.group(1) if m.lastindex and m.group(2).lower() in ("count","fpkm","tpm","rpkm","readcount","reads") else m.group(2)
                    gname = re.sub(r"[-_ .]?(rep)?\d+$", "", gname, flags=re.I)
                    if gname and gname not in grp_names: grp_names.append(gname)
            label = " | ".join(x for x in [sheet, g_ff[j] if g_ff[j] not in ("", "nan", hdr[j]) else "", hdr[j]] if x)
            if 2 <= len(grp_names) <= 6: label += " | samples: " + ", ".join(grp_names)
            out.append({"label": label, "fc_col": hdr_txt, "p_col": hdr[pj] if pj is not None else "", "n": len(recs), "genes": recs})
    # direction-only lists (IDs, no fold change) named up/down
    if not out:
        m = re.search(r"\b(up|down)[\s_-]*(regulated)?", sheet, re.I)
        if m:
            genes = sorted({g.upper().replace("TRAESCS", "TraesCS") for g in body[bc].str.extract(GID, expand=False).dropna() if not g.upper().endswith("LC")})
            if len(genes) >= 5:
                s = 1 if m.group(1).lower() == "up" else -1
                out.append({"label": sheet + " (direction only)", "fc_col": "", "p_col": "", "n": len(genes), "genes": {g: (s, None) for g in genes}, "direction_only": True})
    return out

def work(pmc):
    dst = f"out/{pmc}.json"
    if os.path.exists(dst): return pmc, "cached"
    keys = [l.split("\t") for l in open(f"lists/{pmc}.tsv").read().splitlines()]
    files = [(k, int(s)) for k, s in keys if SUP.search(k) and not re.search(r"/PMC\d+\.\d+\.txt$", k) and int(s) < 300e6]
    res = {"pmcid": pmc, "files": [], "contrasts": []}
    for k, s in files:
        try: data = urllib.request.urlopen(B + k, timeout=300).read()
        except Exception as e: res["files"].append({"file": k, "error": str(e)[:80]}); continue
        found = 0
        for sh, df in tables_from(k, data):
            for c in parse(sh, df):
                c["file"] = os.path.basename(k); res["contrasts"].append(c); found += 1
        res["files"].append({"file": os.path.basename(k), "size": s, "contrasts": found})
        del data
    json.dump(res, open(dst, "w"))
    return pmc, f"{len(res['contrasts'])} contrasts"

if __name__ == "__main__":
    rows = [r for r in csv.DictReader(open("screen.tsv"), delimiter="\t") if r["deg"] == "True" and r["seq"] == "True" and r["host"] == "True" and int(r["nsup"]) > 0]
    ids = [r["pmcid"] for r in rows] if len(sys.argv) < 2 else sys.argv[1:]
    with cf.ThreadPoolExecutor(10) as ex:
        for i, (p, st) in enumerate(ex.map(work, ids)): print(i, p, st, flush=True)
