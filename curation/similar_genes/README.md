# Highly similar genes (gene viewer and gene prediction tool)

Files: `web/sim/similar_<chr>.json` (one per chromosome, 22 files). Keys are the six-digit gene number of a
TraesCS<chr>02G gene; each entry is `[gene ID without the TraesCS prefix, CDS identity %, aligned % of the query]`.

Method: all 106,795 high-confidence CDS from IWGSC RefSeq annotation were aligned with minimap2 (mappy, preset asm20)
against separate A, B, D and Un subgenome indices, so a gene's own copy does not hide its homoeologs.
Hits were kept at >= 85% identity over >= 50% of the query length, top 10 per gene, self-hits removed.
91,343 genes have at least one hit. Script: `sim2.py`.
