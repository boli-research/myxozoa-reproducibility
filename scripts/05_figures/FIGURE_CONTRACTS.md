# Figure contracts

All figures use Python/Matplotlib exclusively and export editable SVG/PDF plus 600-dpi TIFF and a PNG preview.

## Figure 1

- Core conclusion: The study compares two marine and two freshwater *Myxobolus* isolates sampled from distinct hosts and tissues.
- Archetype: schematic-led geographic composite.
- Evidence: national distribution map followed by four locality panels.
- Source data: `metadata/sampling_sites.tsv`, public administrative boundaries and author-approved host/spore artwork.
- Reviewer risk: map licensing, coordinate accuracy and artwork attribution.

## Figure 2

- Core conclusion: Conspecific mitochondrial genomes are similar, whereas the two species differ in size, gene organization and GC-skew structure.
- Archetype: quantitative grid of circular genome maps.
- Evidence: four identically scaled annotated circular maps with GC content and GC skew.
- Source data: annotated GenBank records and computed sliding-window tables.
- Reviewer risk: circularization boundary, ORF filtering, gene-name standardization and sign of GC skew.

## Figure 3

- Core conclusion: Mitochondrial gene order is conserved within species but highly labile across Myxozoa.
- Archetype: asymmetric tree plus aligned gene-order tracks.
- Evidence: mitochondrial phylogeny and standardized gene-order blocks.
- Source data: Newick tree and feature TSVs created by `standardize_mito_features.py`.
- Reviewer risk: inconsistent ORF annotation and feature overlap criteria.

## Figure 4

- Core conclusion: Myxozoans occupy Endocnidozoa within Cnidaria and display diverse mitochondrial chromosome architectures.
- Archetype: phylogeny with categorical architecture annotations.
- Evidence: mitochondrial amino-acid tree, lineage labels and circular/linear/multipartite markers.
- Reviewer risk: long branches, rooting and sparse taxon representation.

## Figure 5

- Core conclusion: Myxozoan nuclear assemblies differ markedly in size, GC content, continuity and conserved-gene recovery.
- Archetype: quantitative grid.
- Evidence: four assembly metrics plus stacked BUSCO proportions.
- Source data: one audited TSV produced by `summarize_assemblies.py`.
- Reviewer risk: technology/assembly-version differences and BUSCO lineage consistency.

## Figure 6

- Core conclusion: Nuclear single-copy orthologs recover conspecific clustering and a coherent myxozoan phylogenomic signal.
- Archetype: single hero phylogeny.
- Evidence: partitioned maximum-likelihood tree with bootstrap and concordance annotations.
- Source data: Newick tree and taxon metadata.
- Reviewer risk: taxon-count mismatch, missing loci, model merging and long-branch effects.

## Figure 7

- Core conclusion: Nuclear functional-gene recovery is pathway-specific, with broad OXPHOS reduction but recurrent recovery of several host-interaction categories.
- Archetype: asymmetric tree plus Coulson-style state matrices.
- Evidence: independent 18S tree, 48-gene mitochondrial-metabolism matrix and 30-sector invasion/adhesion matrix.
- Source data: complete evidence and state matrices from the similarity/HMM scripts.
- Reviewer risk: fragmented assemblies, unequal transcript support and interpretation of “not recovered” as absence.
