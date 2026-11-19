"""Write the expression-module curation files into the repository (curation/expression_degs/)."""
import csv, json, glob, os, datetime
D = "/home/claude/wheatnhr-atlas/curation/expression_degs/"
os.makedirs(D, exist_ok=True)
screen = list(csv.DictReader(open("screen.tsv"), delimiter="\t"))
rep = dict(l.rstrip("\n").split("\t", 1) for l in open("build_report.tsv"))
outs = {os.path.basename(f)[:-5]: json.load(open(f)) for f in glob.glob("out/*.json")}
E = json.load(open("/home/claude/wheatnhr-atlas/web/expr_degs.json"))
used = {s["pmcid"]: s for s in E["S"]}
with open(D + "papers_screened.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["pmcid", "year", "first_author", "title", "journal", "doi", "diseases_searched", "host_DEG_study_with_supplement", "supplementary_files", "tables_with_gene_IDs_and_fold_change", "decision", "atlas_study_id", "disease_assigned"])
    for r in screen:
        p = r["pmcid"]; o = outs.get(p)
        cand = r["deg"] == "True" and r["seq"] == "True" and r["host"] == "True" and int(r["nsup"]) > 0
        if not cand: dec = "not screened further: no host DEG study with supplementary files (pathogen-side, review, or no supplement)"
        elif o is None: dec = "not processed"
        elif not o["contrasts"]: dec = "no machine-readable DEG table (DEGs given as figures, PDF tables or without fold changes)"
        else: dec = rep.get(p, "")
        s = used.get(p)
        w.writerow([p, r["year"], r["first_author"], r["title"], r["journal"], r["doi"], r["diseases"], "yes" if cand else "no", r["sup"], len(o["contrasts"]) if o else "", dec, s["id"] if s else "", s["disease"] if s else ""])
with open(D + "contrasts_included.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["study_id", "citation", "disease", "pmcid", "doi", "contrast", "n_DEGs", "direction_only", "significance_column"])
    for c in E["C"]:
        s = E["S"][c[0]]
        w.writerow([s["id"], s["cite"], s["disease"], s["pmcid"], s["doi"], c[1], c[2], "yes" if c[3] else "no", c[4]])
open(D + "README.md", "w").write(f"""# Expression module (in-silico DEG analysis): curation record

Built {datetime.date.today().isoformat()}. Data file: `web/expr_degs.json`, used by Genomic tools > In-silico expression (DEGs).

## How the data were assembled
1. **Search.** PubMed Central, open-access papers, one query per disease: wheat AND (RNA-seq OR transcriptome OR differentially expressed OR gene expression) AND the pathogen or disease names. {len(screen)} unique papers (`papers_screened.tsv`).
2. **Screen.** Full text checked automatically for a wheat host transcriptome study that reports DEGs and has supplementary files.
3. **Extract.** Supplementary spreadsheets and text tables from the PMC open-access dataset were read. A table was used when it had IWGSC RefSeq gene IDs (v1.0/v1.1 `02G`/`01G`, or v2.1 `03G` converted to v1.1) and a fold-change column; the significance column next to it was used when present. Low-confidence (`LC`) genes were excluded.
4. **Filter.** |log2 fold change| >= 1, and p (or adjusted p, as reported) < 0.05 where the table gives one. Lists that give only up/down direction are kept as direction-only. Tables that repeat a larger table of the same paper were dropped.
5. **Excluded.** Meta-analyses; biocontrol, chemical or mycotoxin-treatment studies; GO/KEGG, epigenetic (CUT&Tag, ChIP, methylation), non-coding RNA and abiotic-stress tables.
6. **Disease** of each study assigned from its title and abstract.

## Files
- `papers_screened.tsv`: every paper found, and the decision for each.
- `contrasts_included.tsv`: every contrast in the atlas, with its study, DEG count and significance column.

Papers whose DEGs are only in figures or PDF tables are not yet included; they can be added by entering their tables by hand.
""")
print("written", len(screen), "papers;", len(E["S"]), "studies;", len(E["C"]), "contrasts")
