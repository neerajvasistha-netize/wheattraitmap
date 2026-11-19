"""Build web/expr_degs.json from the extracted tables (out/*.json) and the paper metadata (screen.tsv)."""
import csv, json, glob, re, sys, datetime
from collections import Counter, defaultdict
ATLAS = "/home/claude/wheatnhr-atlas/web/"
KW = {
 "Fusarium head blight": r"fusarium head blight|\bFHB\b|head blight|F(usarium|\.) graminearum|scab",
 "Fusarium crown rot": r"crown rot|F(usarium|\.) pseudograminearum|F(usarium|\.) culmorum",
 "Spot blotch": r"spot blotch|Bipolaris sorokiniana|B\. sorokiniana|Cochliobolus sativus",
 "Tan spot": r"tan spot|Pyrenophora tritici-repentis|ToxA|Ptr\b",
 "Septoria nodorum blotch": r"(Para)?stagonospora nodorum|Septoria nodorum|P\. nodorum",
 "Septoria tritici blotch": r"Zymoseptoria|Septoria tritici|Mycosphaerella graminicola|Z\. tritici|\bSTB\b",
 "Wheat blast": r"wheat blast|Magnaporthe oryzae|Pyricularia oryzae|\bMoT\b",
 "Sharp eyespot": r"sharp eyespot|Rhizoctonia cerealis|R\. cerealis",
 "Eyespot": r"(?<!sharp )eyespot|Oculimacula|Pseudocercosporella",
 "Take-all": r"take-all|Gaeumannomyces",
 "Leaf rust": r"leaf rust|Puccinia triticina|P\. triticina|\bLr\d",
 "Stripe rust": r"stripe rust|yellow rust|Puccinia striiformis|P\. striiformis|\bPst\b|\bYr\d",
 "Stem rust": r"stem rust|Puccinia graminis|P\. graminis|\bPgt\b|\bSr\d",
 "Powdery mildew": r"powdery mildew|Blumeria|B\. graminis|\bBgt\b|\bPm\d",
 "Common bunt": r"common bunt|Tilletia caries|Tilletia foetida|T\. caries",
 "Dwarf bunt": r"dwarf bunt|Tilletia controversa|T\. controversa",
 "Karnal bunt": r"karnal bunt|Tilletia indica|T\. indica",
 "Wheat streak mosaic virus": r"wheat streak mosaic|\bWSMV\b",
 "Barley yellow dwarf virus": r"barley yellow dwarf|\bBYDV\b",
 "Wheat dwarf virus": r"wheat dwarf virus|\bWDV\b",
 "Wheat yellow mosaic virus": r"wheat yellow mosaic|\bWYMV\b",
 "Soil-borne cereal mosaic virus": r"soil-borne (cereal|wheat) mosaic|\bSBCMV\b|\bSBWMV\b",
 "Bacterial leaf streak": r"bacterial leaf streak|Xanthomonas translucens",
 "Cereal cyst nematode": r"cereal cyst nematode|Heterodera (avenae|filipjevi)",
 "Root-knot nematode": r"root-knot|Meloidogyne",
}
ABIOTIC = re.compile(r"\b(drought|heat|cold|freez|salt|salin|osmotic|PEG|ABA|nitrogen|phosph|cadmium|waterlog|submerg|UV)\b", re.I)
BIOTIC = re.compile(r"inoculat|infect|pathogen|fung|rust|mildew|virus|blight|blotch|spot|nematode|bacteri|mock|\bPst\b|\bBgt\b|\bFg\b", re.I)
SKIP_LAB = re.compile(r"CUT&Tag|ChIP|H3K\d|methyl|ATAC|\bpeaks?\b|lncRNA|circRNA|miRNA|\bGO\b|KEGG|enrich|pathway|primer|qRT|qPCR|SNP|InDel|variant|marker|haplotype|BSA|QTL|homolog|ortholog|miRNA target|effector|secret|fungal gene|pathogen gene", re.I)
STAT = re.compile(r"padj|p\.?adj|fdr|q[\s_.-]*val|adj|corrected", re.I)

def meta_rows():
    return {r["pmcid"]: r for r in csv.DictReader(open("screen.tsv"), delimiter="\t")}

def disease_of(r, txt):
    sc = Counter()
    for d, rx in KW.items():
        t = len(re.findall(rx, r["title"], re.I)); a = len(re.findall(rx, txt[:6000], re.I))
        if t or a: sc[d] = t * 100 + a
    return sc.most_common(1)[0][0] if sc else None

def cite(r):
    a = (r["first_author"] or "").split(" ")[0] or "Anon"
    return f"{a} et al. ({r['year']})"

def main():
    G = json.load(open(ATLAS + "genes_all.json"))
    hc = set(); v21 = {}
    for c, rows in G["g"].items():
        for g in rows:
            hc.add(f"{c}:{g[1]}")
            v = g[5] if len(g) > 5 else None
            if isinstance(v, int) and v >= 0: v21[f"TraesCS{c}03G{str(v).zfill(7)}"] = f"{c}:{g[1]}"
            elif isinstance(v, str): cc, n = v.split(":"); v21[f"TraesCS{cc}03G{n.zfill(7)}"] = f"{c}:{g[1]}"
    def key(gid):
        m = re.match(r"TraesCS([1-7][ABD]|U)0([12])G(\d{6})$", gid)
        if m:
            k = f"{m.group(1)}:{int(m.group(3))}"; return k if k in hc else None
        m = re.match(r"TraesCS([1-7][ABD]|U)03G(\d{7})$", gid)
        if m: return v21.get(gid)
        return None
    M = meta_rows(); S = []; C = []; GENES = defaultdict(dict); report = []
    for f in sorted(glob.glob("out/*.json")):
        d = json.load(open(f)); p = d["pmcid"]; r = M.get(p)
        if not r or not d["contrasts"]: continue
        if re.search(r"meta-?analys|meta-?transcriptom|integrat\w* analysis of (public|multiple)", r["title"], re.I): report.append((p, "skipped: meta-analysis")); continue
        if re.search(r"biocontrol|bio-control|Trichoderma|Bacillus|Pseudomonas|Streptomyces|endophyt|fungicide|nanoparticle|elicitor|priming|chitosan|silicon|biostimulant|antagonis|mycovirus|deoxynivalenol|\bDON\b|zearalenone|mycotoxin|exogenous|inducer|inhibitor|abscisic|gibberellic|proteom|iTRAQ|branched-spike|colou?ration|meta-QTL|consensus genomic regions", r["title"], re.I): report.append((p, "skipped: biocontrol, chemical or toxin treatment study")); continue
        txt = open(f"txt/{p}.txt", errors="ignore").read()
        dis = disease_of(r, txt)
        if not dis: report.append((p, "skipped: no disease identified")); continue
        kept = []
        for c in d["contrasts"]:
            lab = c["label"]
            if SKIP_LAB.search(lab) or re.search(r"\bmeta\b", lab, re.I): continue
            if ABIOTIC.search(lab) and not BIOTIC.search(lab): continue
            has_p = bool(c.get("p_col")); padj = bool(STAT.search(c.get("p_col", "")))
            if not has_p and not c.get("direction_only") and c["n"] > 30000: continue  # all-gene table without a significance test
            genes = {}
            for gid, (fc, pv) in c["genes"].items():
                k = key(gid)
                if not k or fc is None or fc != fc: continue
                if c.get("direction_only"): genes[k] = fc; continue
                if abs(fc) < 1: continue
                if has_p and (pv is None or pv >= 0.05): continue
                genes[k] = fc
            if len(genes) >= 10: kept.append((lab, genes, int(bool(c.get("direction_only"))), c.get("p_col", ""), c.get("file", "")))
        def _clean(lab):
            parts = [x.strip() for x in lab.split("|")]
            parts = [re.sub(r"\.(xlsx?|xlsm|csv|tsv|txt|zip|gz|docx)(\.txt)?$", "", x.split(":")[-1].strip(), flags=re.I) for x in parts]
            parts = [x for x in parts if x and not re.fullmatch(r"(log2?[\s_.()-]*(fold[\s_.-]*change|fc|ratio)[\s_.()a-z/-]*|logfc|fold[\s_.-]*change|fc|unshrunk\.logfc)", x, re.I)]
            return " · ".join(dict.fromkeys(parts))
        merged = []
        for t in kept:
            k = (_clean(t[0]), t[4]); hit = None
            for i, u in enumerate(merged):
                if u[0] != k: continue
                a_, b_ = set(t[1]), set(u[1][1]); j = len(a_ & b_) / max(1, len(a_ | b_))
                if j >= 0.6 or (bool(t[3]) != bool(u[1][3]) and len(a_ & b_) >= 0.8 * min(len(a_), len(b_))): hit = i; break
            if hit is None: merged.append((k, t)); continue
            cur = merged[hit][1]
            if (bool(t[3]) and not bool(cur[3])) or (bool(t[3]) == bool(cur[3]) and len(t[1]) > len(cur[1])): merged[hit] = (k, t)
        kept = [t for _, t in merged]
        # drop contrasts that are subsets of a larger contrast of the same paper (e.g. a top-50 table)
        kept.sort(key=lambda x: -len(x[1])); final = []
        for lab, g, dr, pc, fcc in kept:
            if any(len(set(g) & set(h[1])) >= 0.9 * len(g) and len(h[1]) > len(g) for h in final): continue
            if any(set(g) == set(h[1]) for h in final): continue
            final.append((lab, g, dr, pc, fcc))
        if not final: report.append((p, "skipped: no usable contrast after filters")); continue
        si = len(S)
        S.append({"id": f"E{si+1:03d}", "pmcid": p, "doi": r["doi"], "cite": cite(r), "title": r["title"].rstrip("."), "journal": r["journal"], "year": int(r["year"] or 0), "disease": dis})
        for lab, g, dr, pc, fcc in final:
            ci = len(C)
            parts = [x.strip() for x in lab.split("|")]
            parts = [re.sub(r"\.(xlsx?|xlsm|csv|tsv|txt|zip|gz|docx)(\.txt)?$", "", x.split(":")[-1].strip(), flags=re.I) for x in parts]
            parts = [x for x in parts if x and not re.fullmatch(r"(log2?[\s_.()-]*(fold[\s_.-]*change|fc|ratio)[\s_.()a-z/-]*|logfc|fold[\s_.-]*change|fc)", x, re.I)]
            dd = []
            for x in parts:
                if x not in dd: dd.append(x)
            if dd and re.fullmatch(r"sheet\s*\d*( \(direction only\))?", dd[0], re.I):
                fn = re.sub(r"\.(xlsx?|xlsm|csv|tsv|txt|zip|gz|docx)$", "", fcc or "", flags=re.I)
                dd[0] = (fn + (" (direction only)" if "direction" in dd[0] else "")) if fn else dd[0]
            clean = re.sub(r"\s+", " ", " · ".join(dd)).strip() or "DEGs"
            C.append([si, clean[:120], len(g), dr, pc[:40]])
            for k, fc in g.items():
                GENES[k][ci] = max(-99, min(99, round(fc * 10))) if not dr else (10 if fc > 0 else -10)
        report.append((p, f"{dis}: {len(final)} contrasts, {sum(len(x[1]) for x in final)} DEGs"))
    out = {"built": datetime.date.today().isoformat(), "S": S, "C": C,
           "G": {k: [x for ci, v in sorted(v.items()) for x in (ci, v)] for k, v in sorted(GENES.items())},
           "note": "log2 fold changes x10 (integers); direction-only lists stored as +/-10. Filter: |log2FC| >= 1 and p < 0.05 where the study reports a p or adjusted p column."}
    json.dump(out, open(ATLAS + "expr_degs.json", "w"), separators=(",", ":"))
    with open("build_report.tsv", "w") as fh:
        for p, m in report: fh.write(f"{p}\t{m}\n")
    print(len(S), "studies", len(C), "contrasts", len(GENES), "genes", Counter(s["disease"] for s in S).most_common())

if __name__ == "__main__": main()
