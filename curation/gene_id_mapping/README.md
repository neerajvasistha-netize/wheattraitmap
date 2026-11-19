# Gene ID conversion tables for GO/KEGG enrichment

Files: `web/idmap/*.gz.txt` (gzip JSON stored as base64 text, because the Claude artifact does not serve .gz files; `{input ID: "chr:gene number" or "a|b" when several}`), loaded by the
enrichment tool only when a list contains that ID type. Keys map to IWGSC RefSeq v1.1 high-confidence genes
(TraesCS<chr>02G<number>), which carry the atlas GO and KEGG annotation.

- `entrez`, `refseq`, `uniprot`: Ensembl Plants release 63 cross-reference dumps
  (ftp.ebi.ac.uk/ensemblgenomes/pub/plants/current/tsv/triticum_aestivum/Triticum_aestivum.IWGSC.63.{entrez,refseq,uniprot}.tsv.gz).
  RefSeq and UniProt accessions are stored without version or isoform suffix.
- `pan_Traes<code>`: pangenome cultivar genes matched to Chinese Spring v1.1 genes with the genome-wide anchors
  used for the pangenome coordinate conversion (minimap2 asm10, CDS identity >= 95%, coverage >= 60%):
  ARI ArinaLrFor, JAG Jagger, JUL Julius, KAR Kariega, LAC Lancer, LDM Landmark, MAC Mace, SYM Mattis,
  NOR Norin 61, PARA Paragon, RN Renan, STA Stanley.
- IWGSC v1.0 (01G) and v1.1 (02G) IDs are read directly; v2.1 (03G) IDs through the v1.1 to v2.1 projection in genes_all.json.
- Not supported: TGACv1 and low-confidence genes (no annotation).

Build: `build_idmap.py` (paths point to the working copies of the downloads and anchors).
