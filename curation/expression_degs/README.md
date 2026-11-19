# Expression module (in-silico DEG analysis): curation record

Built 2026-11-19. Data file: `web/expr_degs.json`, used by Genomic tools > In-silico expression (DEGs).

## How the data were assembled
1. **Search.** PubMed Central, open-access papers, one query per disease: wheat AND (RNA-seq OR transcriptome OR differentially expressed OR gene expression) AND the pathogen or disease names. 574 unique papers (`papers_screened.tsv`).
2. **Screen.** Full text checked automatically for a wheat host transcriptome study that reports DEGs and has supplementary files.
3. **Extract.** Supplementary spreadsheets and text tables from the PMC open-access dataset were read. A table was used when it had IWGSC RefSeq gene IDs (v1.0/v1.1 `02G`/`01G`, or v2.1 `03G` converted to v1.1) and a fold-change column; the significance column next to it was used when present. Low-confidence (`LC`) genes were excluded.
4. **Filter.** |log2 fold change| >= 1, and p (or adjusted p, as reported) < 0.05 where the table gives one. Lists that give only up/down direction are kept as direction-only. Tables that repeat a larger table of the same paper were dropped.
5. **Excluded.** Meta-analyses; biocontrol, chemical or mycotoxin-treatment studies; GO/KEGG, epigenetic (CUT&Tag, ChIP, methylation), non-coding RNA and abiotic-stress tables.
6. **Disease** of each study assigned from its title and abstract.

## Files
- `papers_screened.tsv`: every paper found, and the decision for each.
- `contrasts_included.tsv`: every contrast in the atlas, with its study, DEG count and significance column.

Papers whose DEGs are only in figures or PDF tables are not yet included; they can be added by entering their tables by hand.
