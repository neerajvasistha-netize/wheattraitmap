# In-silico placement of loci from their molecular markers (19 Nov 2026)

Scope: every atlas record without a physical position that names markers (630 records).

1. **Markers with a known physical position** (atlas marker table `data/markerpos_all.json.gz`: 90K/35K/55K/660K SNP maps,
   DArT tag BLAST, Zhao et al. 2019 SSR positions) or with the position written in the marker name (chr_position SNPs).
   Earlier passes missed these because of name variants (Xiwa/Xiwb, numeric DArTseq clone IDs, wmc120.2-type suffixes).
2. **Marker sequences found and mapped on RefSeq v1.0** (`scanq.py`, every chromosome streamed from Ensembl Plants):
   - BARC SSR primer pairs (USDA-ARS BARC list, scabusa.org, `barc_primers.json`): in-silico PCR, product 60-800 bp,
     up to 2 mismatches per primer outside the 3' 12 bases. 7 of 18 needed pairs amplify.
   - RFLP probe sequences from GenBank (BCD/CDO/ABG/MWG/WG/PSR clones, `rflp_seqs.json`): 16-mer seed clusters in 5 kb
     windows; the copy on the record's chromosome is used (the homoeologous copies are listed in `scan_results.json`).
3. Rules: several markers within 30 Mb give an interval, otherwise the first-listed marker gives the position; grade B;
   placements that contradict the chromosome arm stated by the source were rejected (Xwmc120.2 1AS, QLr.fcu-3AL, QKb.cim-3DL);
   names that are marker IDs rather than positions (Navathe et al. chr:ID; Kazakh DArTseq IDs) and a locus earlier rejected as
   paralogous (S19-081) were not used.

Result: 75 records placed (`records.tsv`); 51 in `web/seed_ad68.json`, 24 live-database records updated in the live database
and in `web/seed_live.json`. Not placed: wPt DArT markers with no public sequence (about 100), SSRs without an accessible primer
sequence (GWM 8xx/9xx, many WMC, CFA, CFD), AFLP and KSU/FCP RFLP probes without sequence, and SNPs whose physical position
lies on a homoeologous chromosome.

## BARC re-run with relaxed primer matching (ad69, same day)
The BARC primer list uploaded by the owner is identical to the USDA-ARS list already used. Many BARC primers carry a 5' GCG tail
that is not in the genome, so the 11 pairs that failed were re-run on their target chromosomes (`relax.py`): each primer located by
its 3' 15 bases, product 50-1500 bp. Six more records placed (`web/seed_ad69.json`, S163-005 in the live database):
Yrxy2 (barc5, 2A), QLr.osu-7BL and YrC591 (barc182, 7B), QTLBvr5AS (barc10, 5A), FHB association 2B (barc101), QYr.caas-2BS (barc13).
QPm.caas-7AL was rejected (barc174 product on 7AS). barc121, barc147, barc173, barc45, barc98 and the 3S/6B targets give no product.

## GWM primers and the SSR sequence workbook (ad70, same day)
GWM primer sequences (224 pairs, WMS052-WMS588) were taken from US patent 6720137 (IPK Gatersleben; the Röder et al. 1998 set,
`gwm_primers.json`) and run as in-silico PCR on the target chromosomes of the 29 unplaced GWM markers (`relax2.py`). Products of the
expected size (e.g. GWM255 148 bp, GWM257 193 bp, GWM263 136 bp) placed QSb.bhu-7B (3 records), QStb-7B (History/Rubens), the
4BS FHB QTL of Lu et al. 2012, YrS2199 and the 5A QTL of Szabo-Hever et al. 2014 (`web/seed_ad70.json`, live records updated).
PmY39 was not placed (alien 2U segment). `SSR_marker_sequences.xlsx` lists every unplaced SSR marker in the atlas with its primers
where public (BARC, GWM) and the in-silico PCR product, plus the full GWM and BARC primer lists.

WMC primer table (Gupta et al. 2002; wmc24-wmc331, supplied by the owner as an image, transcribed to `wmc_table.tsv`): only wmc105, wmc169 and wmc215 are used by unplaced atlas loci. wmc169 amplifies 3A 721.2 Mb but its locus is on an alien 3S chromosome; wmc105 (6B) and wmc215 (3A) give no product. No records changed.

WMC genome-wide in-silico PCR (`epcr.py`, `wmc_epcr.json`): 42 of the 58 Gupta et al. 2002 primer pairs amplify in Chinese Spring RefSeq v1.0; 38 agree with the Zhao et al. 2019 positions (validation), 3 differ and 1 is new (WMC97, chrUn). New unique products were added to `data/markerpos_all.json.gz`. Sheet 'WMC physical map' in `SSR_marker_sequences.xlsx`.

## Röder et al. 1998 Table 3 and Abeysekara et al. 2009 (ad71, same day; papers supplied by the owner)
GWM primers were re-read from the Röder et al. 1998 table (`gwm_roder1998.json`, 223 markers; it corrects several pairs mis-read from
the patent text, e.g. GWM260, GWM264) and merged with the patent set (`gwm_merged.json`). Re-running the unplaced GWM markers placed
QLr.pser-1BL (GWM264, 1BL copy 328.94 Mb) and QLr.fcu-3AL (GWM666, 3AL copy 723.89 Mb). From Abeysekara et al. 2009 Table 1, the EST
marker XBG262975 amplifies 1A 2.64 Mb; with Xcfd58.1 next to Xcfa2153 (7.23 Mb), Snn4 and its QTL records were placed at 2.64-7.23 Mb
(V0095 had been placed at Xcfa2263, a 2A flanking marker). QSnb.fcu-3A placed by BARC57 (3A 10.31 Mb).
