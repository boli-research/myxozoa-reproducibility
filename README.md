# Myxozoa Reproducibility

Code, data and documentation for reproducing the analyses and figures of the manuscript:

> **Integrated mitochondrial and nuclear genomics reveals distinct evolutionary patterns in Myxozoa**

The study integrates complete mitochondrial genomes and draft nuclear genomes of two *Myxobolus* species (*M. episquamalis* isolates QD and XM; *M. pronini* isolates TZ and XY) to compare mitochondrial and nuclear evolutionary patterns across myxozoans.

## What is in this repository

- Analysis scripts for every computational step described in the manuscript, organised as small command-line programs (`scripts/`).
- Curated metadata: isolates, NCBI accessions, sampling sites, software versions and the manuscript-to-script mapping (`metadata/`).
- Small key datasets (`data/`): 18S rDNA sequences, assembled mitochondrial genomes, the mitochondrial and nuclear supermatrices with partition files, the 76 trimmed BUSCO locus alignments, final phylogenetic trees (including gene/site concordance factors), gene-order GenBank files, BUSCO and QUAST summaries, and the OrthoFinder species tree.
- Final manuscript figures (Figs. 1–7) and supplementary figures (Figs. S1–S3) as submitted (`figures/`).
- Supplementary Tables S1–S4 (`supplementary/tables/`).
- Representative light-micrograph images of the two species (`images/`).
- Unit/smoke tests for the core scripts (`tests/`).

Raw sequencing reads, draft nuclear genome assemblies and other large datasets are **not** stored here; they are available from the NCBI records listed below and in `metadata/accessions.tsv`.

## Data availability

- NCBI BioProjects: [PRJNA1377137](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1377137), [PRJNA1453518](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1453518)
- Complete mitochondrial genomes (GenBank): `PX674009.1`, `PX674010.1`, `PZ514093.1`, `PZ514094.1`
- SSU rDNA sequences (GenBank): `OR621296.1`, `PZ580198.1`, `PZ564097.1`, `PZ575193.1`

The sample-level mapping is in `metadata/accessions.tsv`.

## Repository structure

```text
metadata/          Samples, accessions, sampling sites, software versions, code map
scripts/           Analysis and figure-generation code (see scripts/README.md)
  01_genome_assembly_annotation/
  02_mitogenome/
  03_phylogenetics/
  04_functional_genes/
  05_figures/
  checks/
  lib/
data/              Small key datasets (see below)
  18S_rDNA/        SSU rDNA sequences and the inter-species alignment
  mitogenomes/     Four assembled mitochondrial genomes (FASTA)
  alignments/      Supermatrices, partition files and the 76 trimmed BUSCO loci
  trees/           Final ML trees, concordance factors and partition schemes
  gene_order/      Clinker-ready annotated GenBank files for the gene-order comparison
  busco/           BUSCO summary table, 17 short summaries and QUAST report
  orthofinder/     OrthoFinder rooted species tree
figures/           Final main and supplementary figures as submitted
supplementary/     Supplementary Tables S1–S4
images/            Representative light-micrograph images
environment/       Python dependency list
tests/             Deterministic unit and smoke tests
```

## Python environment

Python 3.11 or later is recommended.

```bash
python3 -m venv .venv
.venv/bin/pip install -r environment/python-requirements.txt
```

External bioinformatics programs (fastp, SPAdes, BUSCO, NOVOPlasty, MAFFT, trimAl, IQ-TREE, HMMER, clinker and others) are listed with versions and settings in `metadata/software_versions.tsv`. Scripts that need an external program stop with a clear error when it is missing.

## Reproducing the analyses

Run the pipeline stages in numeric order; every script supports `-h`:

1. `scripts/01_genome_assembly_annotation/` — read QC, host/microbial decontamination, nuclear assembly, repeat annotation and structural gene annotation.
2. `scripts/02_mitogenome/` — NOVOPlasty mitochondrial assembly, circularity validation, feature standardisation and clinker input preparation.
3. `scripts/03_phylogenetics/` — locus extraction, alignment/trimming/concatenation and the mitochondrial, nuclear and 18S phylogenies.
4. `scripts/04_functional_genes/` — mitochondrial-metabolism gene classification and invasion/adhesion domain searches.
5. `scripts/05_figures/` — regeneration of manuscript Figs. 1–7, Supplementary Fig. S1 and figure source data.

The exact correspondence between manuscript methods and scripts is recorded in `metadata/manuscript_code_map.tsv`. Inputs are described by `metadata/samples.tsv` (paths are placeholders — point them at your local copies of the NCBI data). Large or licensed inputs (reference databases, host genomes) are not redistributed.

Notes:

- `scripts/05_figures/Supplementary_Figure_S1.py` rebuilds Supplementary Fig. S1 from the original microscopy PPTX following the image-integrity rules in `scripts/05_figures/figure_contract.md`. The 65 MB source PPTX is not tracked in git; it is available from the corresponding author on request.
- The OrthoFinder run used for gene-family analysis followed the standard command `orthofinder -f <protein_fasta_dir>`; the rooted species tree it produced is in `data/orthofinder/`.
- Tree files in `data/trees/`: `mito_AA_partitioned.*` is the final mitochondrial ML analysis (83 taxa, five protein-coding genes; Figs. 4–5) and `nuclear_partitioned.*` the final nuclear ML analysis (71 taxa, 76 BUSCO loci; Fig. 6), each with its IQ-TREE report and merged partition scheme. `nuclear_concordance.cf.*` holds the gCF/sCF results computed from `gene_trees.nwk` (76 single-locus trees) and the nuclear supermatrix. `partitions_pruned_keep18_fixed.treefile` is the Fig. 7 ordering tree, pruned from `nuclear_partitioned.treefile` to the 18 surveyed taxa. `busco_core52.*` is an earlier, superseded nuclear run (63 taxa, 52 loci) retained only for provenance — it does not correspond to any published figure.
- Script-to-figure numbering: the figure scripts were named before the figures were reordered during revision. `fig1_sampling_map.py` → Fig. 1, `fig5_genome_busco.py` → Fig. 2, `fig2_mitogenome_maps.py` → Fig. 3, `fig3_gene_order_tree.py` → Fig. 4, `fig4_cnidarian_architecture.py` → Fig. 5, `fig6_nuclear_tree.py` → Fig. 6, `fig7_coulson.py` → Fig. 7.

## Running the tests

```bash
python -m unittest discover tests -v
```

The tests use only synthetic inputs and take well under a minute.

## License and citation

Code in this repository is released under the [MIT License](LICENSE). Figures and tables may be reused under the terms of the manuscript's license with attribution. Citation metadata are provided in `CITATION.cff`.
