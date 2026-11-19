# Physical positions from ten uploaded papers (Nov 2026)

Papers supplied by the atlas owner on 19 Nov 2026; records written to `web/seed_ad66.json` by `ad66_build.py`
(full record list with status in `records.tsv`). Positions were converted to RefSeq v1.0/v2.1 and the 13 pangenome
cultivars with the genome-wide gene anchors (`scripts/pangenome/convert.py`).

| Paper | Atlas source | What was done |
|---|---|---|
| Xu et al. 2024 TAG 137:49 (crown rot, C549/3642) | A_crown_rot_010 | 9 QTL placed from Table 2 flanking SNPs (named chr_position, v1.0). The 6B/7B SNP names in Table 2 lack one digit; the Discussion intervals (6B 281.6-327.0, 7B.1 315.6-432.5 Mb) were used and the six 6B/7B marker records corrected. |
| ElDoliefy et al. 2020 Euphytica 216:71 (FHB, Glenn) | V_ElDoliefy2020 | 4 QTL placed from the Fig. 2 DArT map markers (DArT tag positions, v1.0). 6BS interval is wide (151.8-466.8 Mb). |
| Farhad et al. 2026 PLOS One 21:e0349201 (wheat blast) | V_Farhad2026 | All 45 QTL of Table 3 with regions and leading SNPs (v2.1): 4 existing records placed, 41 added. |
| Zhao et al. 2023 Mol Breed 43:82 (spike length, FHB) | V_Zhao2023 | QSl.yaas-5A and QSl.yaas-5B added (Table 4, v2.1; positive allele lowers FHB severity). QFhb.yas-1AL is only cited; no position. |
| Semagn et al. 2023 Front Plant Sci 14:1190358 (FHB) | Z_1368 | 12 QTL checked against Table 2 (v2.0); 3 missing positions added, 3 wrong positions corrected, reference set to v2.0. |
| Yang et al. 2022 Front Plant Sci 13:990915 (STB) | Wai7A | 6 QTL placed from text regions or DArT-seq tag positions; 6 have only cM positions and stay unplaced. |
| Genievskaya et al. 2020 PeerJ 8:e9820 (rusts) | OA_Genievskaya2020genome | 29 of 44 MTAs placed from 90K SNP positions; 15 not placed (SNP maps to the homoeologous chromosome or has no position). |
| Running et al. 2022 Front Plant Sci 13:793925 (tan spot, Tsc1) | Running2022 | QTs.zhl-1A placed at the Tsc1 interval (the paper equates them). |
| Ren et al. 2019 TAG 132:241 (FHB, 7DL) | Z_1349 | QFhb.cau-1AS not placed: linked marker IWB63682 has no known physical position. |
| Shi et al. 2023 Agronomy 13:1712 (FHB GWAS) | Agron2023FHB / A_FHB_038 | 7 missing gFHB loci of Table 1 added (v1.1); cited Qfhb.nmbu.1A.1 placed at 3.1-13.3 Mb. |

## Yang et al. 2022 (STB): in-silico marker mapping (added later the same day)
Supplementary Table 3 of the paper (PMC open-access DataSheet_4.xlsx) gives the IWGSC physical region and the sequence of the
representative DArTseq marker of every QTL. The 14 marker tags were aligned end to end (up to 3 mismatches, no gaps) to every
chromosome of RefSeq v1.0 and v2.1 (`yang2022_tag_placements.json`). All 13 QTL now carry the paper's region; 7 also carry the
tag position as peak (1A.1, 1B.2, 3B.1, 5A.1, 6A.2, 7A.1, 7A.2; 7 QTL). Tags of 2B.1, 3A.1, 3A.2, 4A.1 and 6A.1 have no match on the
QTL chromosome in Chinese Spring (presence/absence markers or homoeologous best hits), and the 1B.1 tag maps outside its region.

## Second batch (19 Nov 2026), `web/seed_ad67.json` (`ad67_build.py`, `records_ad67.tsv`)
| Paper | Atlas source | What was done |
|---|---|---|
| Dhariwal et al. 2020 IJMS 21:4497 (FHB, AAC Tenacious) | X5_1919 | All 22 FHB QTL of Table 3 at their closest-marker positions (RefSeq v2.0); 15 existing records corrected (cM values had been read as Mb, wrong reference), 7 added. QFhs.lrdc-2B, QFhb.lrdc-2B and QDon.lrdc-5A.2 (individual environments only) placed at the overlapping epistasis intervals via flanking 90K markers. |
| Ciechanowska et al. 2022 Can J Plant Sci 102:1139 (rusts) | NK_Ciechanowska2022 (new) | QLr.dms-4A, QSr.dms-1A, QSr.dms-2B added (Table 3, RefSeq v2.0). |
| Farzand et al. 2021 Can J Plant Pathol 43:S227 (stripe rust, AAC Innova) | NK_Farzand2021 (new) | 11 QTL added (Table 3, RefSeq v1.0); peak at the closest marker where its 90K position falls in the interval. |
| Farzand et al. 2022 Crop Sci 62:1088 (stripe rust, AAC Cameron/P2711) | NK_Farzand2022 (new) | 7 QTL added (Table 3, RefSeq v1.0), including Yr17 (QYr.lrdc-2A.1). |
The other papers in this upload were the same files as the first batch and were already done.
