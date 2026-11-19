# Gene network and functional response map (Genomic tools)

Input: any gene list (all ID types of the enrichment tool), optional log2 fold change column.

Functional response map: genes are assigned to functional categories by keyword rules on the IWGSC gene
descriptions (the rules are the NMCAT list in web/index.html), grouped by the cell compartment where each
category mainly acts (apoplast and cell wall, plasma membrane, cytoplasm, nucleus, chloroplast). One cell per
gene, coloured by log2 fold change, drawn on a schematic plant cell (apoplast, cell wall, membrane, cytoplasm,
nucleus, chloroplast) with arrows for the main flow of the defence response. The layout is the atlas's own design.

Gene network: genes are linked by
- shared function: 2 or more specific GO terms (each annotated to at most 300 genes) or a KEGG pathway;
- co-expression: Pearson r >= 0.8 of log2 fold changes across the atlas DEG contrasts (expr_degs.json),
  at least 3 shared contrasts, each gene keeping its 5 strongest partners;
- sequence similarity: CDS identity >= 90% (web/sim/similar_<chr>.json).
Layout: Fruchterman-Reingold per connected component, components shelf-packed by size. Edge list exports as CSV
for Cytoscape.

## Network analysis (added Nov 2026)
- Topology per node and for the network, as defined in NetworkAnalyzer (Assenov et al. 2008, Bioinformatics 24:282):
  degree, betweenness (Brandes algorithm, normalised within the component), closeness, clustering coefficient,
  eccentricity, average shortest path; network density, diameter, radius, characteristic path length,
  centralization, heterogeneity, components; degree distribution with a power-law fit.
- Hub ranking as in cytoHubba (Chin et al. 2014, BMC Syst Biol 8:S11): MCC (from all maximal cliques, Bron-Kerbosch),
  degree, MNC, DMNC (epsilon 1.7), betweenness, closeness; top 10 drawn as a ranked circle figure.
- Modules: MCODE (Bader & Hogue 2003, BMC Bioinformatics 4:2) with degree cutoff 2, node score cutoff 0.2,
  haircut and k-core 2; each module labelled with its most enriched GO term or KEGG pathway (hypergeometric p < 0.01
  against all annotated genes).
- Node and edge tables export as CSV (edge table imports into Cytoscape as a network, node table as attributes).
