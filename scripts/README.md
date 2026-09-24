# Scripts

The repository uses small command-line programs instead of a workflow manager so
that each released step can be inspected and run independently.

| Directory | Purpose |
|---|---|
| `01_genome_assembly_annotation/` | read QC, contaminant filtering, assembly, repeats and structural annotation |
| `02_mitogenome/` | NOVOPlasty setup, circularity support, annotation standardisation and clinker input preparation |
| `03_phylogenetics/` | mitochondrial and BUSCO alignments, concatenation and IQ-TREE commands |
| `04_functional_genes/` | similarity/HMM searches and three-state recovery classification |
| `05_figures/` | Python/Matplotlib source for manuscript Figs. 1–7, Supplementary Fig. S1 and source-data export |
| `lib/` | shared parsers and validation helpers |
| `checks/` | public-release audit |

Run every script with `-h` or with no arguments to see its required inputs. Shell
wrappers stop immediately when an external program is missing. Python figures
export SVG, PDF, 600-dpi TIFF and a preview PNG from the same source.

The exact correspondence between manuscript methods and scripts is recorded in
`metadata/manuscript_code_map.tsv`. Inputs containing unpublished sequence data
belong under ignored local directories (`source_data/private/`, `work/` or
`results/`) and must not be committed.
