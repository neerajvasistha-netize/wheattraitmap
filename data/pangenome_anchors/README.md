# Genome-wide pangenome anchor tables (Nov 2026)

One table per cultivar assembly: every IWGSC RefSeq v1.1 high-confidence gene (106,795 with a CDS) matched to that
cultivar's genes with minimap2 (mappy, asm10): identity >= 95%, coverage of the Chinese Spring CDS >= 60%, best hit,
same-chromosome hit preferred. Columns: v11_gene, chr, mid_Mb (cultivar), identity, cultivar_gene.

78,000-89,000 genes anchored per cultivar (Kariega 51,632); 100% of the genome covered in 10-Mb bins.
Agreement with the earlier locus-restricted tables (WheatNHR_liftover, Oct 2026): 98.2-99.8% of shared genes.

Built by `scripts/pangenome/anchors_all.py`; records converted with `scripts/pangenome/convert.py`
(conversion rules in `liftover_atlas.py`: anchors within +/-5 Mb, orientation from anchor order, outliers > 10 Mb removed;
flags ok, inverted, wide_gap, few_anchors, other_chromosome, unmapped). Results in `web/seed_ad64.json`.
