# Singh et al. (2023) spot blotch MTAs with unassigned chromosomes: DArT tag placement

Source: Singh S, Gaurav SS, Vasistha NK, et al. (2023) Frontiers in Plant Science 13:1036064 (atlas source `Singh2023FPS`).
Eleven MTAs in Table 2 have chromosome "NA" (DArT-seq SNPs not assigned by the vendor). Their DArT clone IDs and 69 bp tag
sequences come from the authors' DArT-seq key. Atlas records S205-019 to S205-034, updated in `web/seed_ad65.json` (Nov 2026).

## Method
1. **Trim.** Each 69 bp tag starts at the PstI remnant (TGCAG) and, for these short fragments, runs into the sequencing adapter
   (AGATCGGAAGAGC). The adapter part was removed, leaving 37 to 61 bp of genomic sequence ending at the MspI site (CCG).
   The SNP offset in the clone ID (e.g. `0--15:G>C`) is 0-based from the first base of the tag; the reference allele was checked
   against the tag for all 11 markers. (`tags_trimmed.tsv`)
2. **Align.** Every chromosome FASTA of 15 assemblies (IWGSC RefSeq v1.0, RefSeq v2.1, and Alchemy, ArinaLrFor, Jagger, Julius,
   Kariega, Lancer, Landmark, Mace, Norin 61, Paragon, Renan, SY Mattis, Stanley) was streamed from Ensembl Plants
   (`genome_urls.txt`) and searched for each tag and its reverse complement, end to end, with up to 3 mismatches and no gaps
   (four-seed Aho-Corasick search, then full-length check). (`scripts/scan2.py`, `scripts/analyse.py`)
3. **Place.** The consensus chromosome is the one most assemblies agree on. In each assembly the best match (fewest mismatches)
   on that chromosome gives the position of the SNP base. Matches with more than 2 mismatches, or on a homoeologous chromosome
   with mismatches, are not used (flag `unmapped`); an exact match on another chromosome is kept with flag `other_chromosome`.
   Each assembly's position comes from its own alignment, not from liftover. (`scripts/make_seed.py`, `tag_placements.tsv`)
4. **Grade B**: position derived from the marker sequence, not printed in the paper.

## Result (SNP base)
| Record | Marker | Chr | RefSeq v1.0 (bp) | RefSeq v2.1 (bp) | Assemblies placed |
|---|---|---|---|---|---|
| S205-019 | M10592 | 1B | 32,775,893 | 34,951,779 | 14 of 15 |
| S205-021 | M11870 | 7A | 736,378,850 | 744,170,356 | 12 of 15 |
| S205-022 | M976 | 1A | 323,192,097 | 324,794,490 | 15 of 15 |
| S205-023 | M7025 | 7D | 5,376,644 | 5,734,786 | 15 of 15 |
| S205-025 | M11338 | 7A | chrUn 83,577,041 | 728,555,135 | 13 of 15 |
| S205-026 | M8531 | 6B | 229,522,401 | 235,600,384 | 15 of 15 |
| S205-027 | M4837 | 1A | 559,780,854 | 561,856,018 | 14 of 15 |
| S205-029 | M9264 | 6A | 4,773,855 | 6,076,788 | 13 of 15 |
| S205-030 | M12439 | 1D | 21,544,688 | 23,694,135 | 15 of 15 |
| S205-032 | M9994 | 2A | 769,631,212 | 777,461,293 | 14 of 15 |
| S205-034 | M9263 | 1B | 588,284,402 | 595,196,353 | 13 of 15 |

M11338 lies on the unanchored chrUn in RefSeq v1.0; its record uses the RefSeq v2.1 position (7A). Tags not found in a cultivar (more than 3 mismatches) may reflect presence/absence variation.
