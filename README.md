# WheatNHR-Atlas

A physical genomic atlas and web database of wheat disease resistance loci (WheatRLA, Wheat Resistance Loci Atlas).

**Website:** https://neerajvasistha-netize.github.io/WheatNHR-Atlas/  
**Mirror (identical content, auto-updated from this repository):** https://wheattraitmap.neeraj-vasistha.workers.dev/

Maintained by **Neeraj Kumar Vasistha**, Molecular Plant Biology Lab, Department of Genetics and
Plant Breeding, Rajiv Gandhi University, Rono Hills, Itanagar, Arunachal Pradesh 791112, India.
ORCID [0000-0001-5527-340X](https://orcid.org/0000-0001-5527-340X).

> **Status: public preview, under curation.** The website mirrors the WheatRLA database on Claude
> (7,826 records as of 18 Nov 2026; live edits exported to `web/seed_live.json`). The contents of this repository have not been
> hand-verified end to end. See *Curation status* below before using any record.

## What this is

A single curated table of wheat disease resistance loci placed on the IWGSC RefSeq assemblies,
with a self-contained web interface for browsing and analysing them. Loci are drawn from the
published QTL mapping, GWAS, meta-analysis and gene-cloning literature.

As of the current commit:

| | |
|---|---|
| Records | 5,555 |
| With a physical position | 4,798 |
| Source studies with loci | 485 |
| Papers reviewed | 856 |
| Diseases | 24 |

Diseases covered: spot blotch, tan spot, Septoria nodorum blotch, Septoria tritici blotch,
Fusarium head blight, Fusarium crown rot, wheat blast, sharp eyespot, eyespot, leaf rust,
stripe rust, stem rust, powdery mildew, multiple rusts, common bunt, dwarf bunt, Karnal bunt,
loose smut, cereal cyst nematode, bacterial leaf streak, wheat streak mosaic virus,
barley yellow dwarf virus, wheat dwarf virus and wheat yellow mosaic virus.

## Repository layout

```
data/
  merged_v13f.json        the master record set: every locus and every source citation
  papers_uploaded.json    the paper review log (what was read, what was used, what was not)
  Table_S1_...xlsx        the master set as a spreadsheet, one row per locus
  markerpos_all.json.gz   830,150 marker positions used to place loci that lack coordinates
  sr_v21.json.gz          RefSeq v1.1 to v2.1 gene anchors for the liftover
web/
  app.js, body.html       the application source
  index.html              the built single-file page (generated; do not edit by hand)
  *_lite.json             marker, primer, gene and GO/KEGG tables the page loads at startup
scripts/pipeline/         build the page and the spreadsheet from the master record set
scripts/curation/         one script per curation task, kept so each change is reproducible
docs/                     completeness audit and the paper download worklist
```

## Record format

Each record in `data/merged_v13f.json` carries:

`id`, `name`, `cls`, `disease`, `chr`, `start`, `end`, `peak`, `ref`, `marker`, `pve`, `stat`,
`source`, `grade`, `curation`, `note`, `evidence`.

* `cls` is the kind of locus: `CG` cloned gene, `DG` designated gene, `SG` susceptibility or
  sensitivity gene, `QTL`, `MTA`, `QTN`, `BSA`, `MQTL`.
* `start`, `end` and `peak` are in Mb. `peak` is present **only where the source reports a peak**.
  Where a source gives an interval and no peak, `peak` is null and the interface says so rather
  than showing the midpoint as if it were a position.
* `ref` names the assembly the source used, normally `RefSeq v1.0` or `RefSeq v2.1`. The pipeline
  lifts every record to v2.1 for display.
* `grade` is the confidence in the physical position: **A** read directly from the source,
  **B** derived from a marker or interval the source gives, **C** inferred by marker matching in
  this atlas, **D** uncertain or unplaced.
* `curation` is `hand` or `auto`. **Every `auto` record still needs checking against its paper.**
* `note` states where the position came from, and any correction that has been applied.
* `evidence` quotes the source row the record was built from, so a reader can check it.

## Rebuilding

Run in this order from the repository root after any edit to `data/merged_v13f.json`:

```sh
python scripts/pipeline/lift_merged.py      # lift v1.x records to RefSeq v2.1
python scripts/pipeline/resolve_flags.py
python scripts/pipeline/sr.py               # shared regions between diseases
python scripts/pipeline/gc_update.py        # candidate gene counts
python scripts/pipeline/build_v13.py        # writes web/index.html
python scripts/pipeline/table_s1_v13.py     # writes the spreadsheet
```

The built page is a single HTML file that loads four JSON files from the same folder. It needs no
server beyond static file hosting, and performs all of its analysis in the browser: physical
meta-QTL, physical meta-GWAS, genes in a region, gene search, gene-family and GO/KEGG enrichment,
susceptibility gene panel, marker lookup, PCR and KASP primer design, CAPS enzyme finder, GWAS
and QTL annotation, coordinate conversion and a sequence toolkit.

## Data sources

Positions come from the IWGSC RefSeq v1.0 and v2.1 assemblies. Marker positions are compiled from
Sun et al. (2020, *Plant Biotechnol J*) for the Axiom 35K, 55K and 660K arrays, Hu et al. (2022)
for the 90K array on v2.1, Zhao et al. (2019) for SSRs, and in-house BLAST and in-silico PCR of
DArT tag and primer sequences. Every locus record names its source study, with a full citation in
the `sources` block of the master file.

**No publisher PDFs or supplementary files are included in this repository.** The review log
records what was read; the papers themselves are not redistributed.

## Curation status

This is an honest account of what is and is not verified.

* 4,798 of 5,555 records have a physical position. Of the 757 without one, most name no marker in
  the source paper at all, so no marker map can place them.
* Records marked `curation: auto` were built by script from a paper's tables and have not been
  checked by hand. They are the majority.
* `docs/WheatNHR_Atlas_completeness_audit.xlsx` lists about 1,080 locus names that appear in
  papers already read but are absent from the atlas, across 260 studies. It is a name-matching
  estimate, not a verified list of omissions; read its Method sheet before using it.
* `docs/WheatNHR_Atlas_papers_to_download.xlsx` is the resulting worklist of papers to obtain,
  highest expected yield first.
* Known systematic weakness: the rust and powdery mildew records trace to two published
  compilations rather than to the primary mapping studies, so the original study is not credited
  for those loci. This needs fixing before publication.
* Several sources have been rebuilt after errors were found, including Lhamo et al. (2024),
  Serajazari et al. (2023), Martin et al. (2025), Navathe et al. (2023), Chaturvedi et al. (2023)
  and Pundir et al. (2023). Each correction is recorded in the affected records' `note` field and
  in the matching script under `scripts/curation/`.

## Citing

The atlas itself is not yet published. Until it is, cite the original study for any locus; the
`sources` block of the master file gives the full citation and DOI for each.

## Licence

Not yet decided. All rights reserved until a licence is added. The underlying locus information
is drawn from the published literature and belongs to the original authors.
