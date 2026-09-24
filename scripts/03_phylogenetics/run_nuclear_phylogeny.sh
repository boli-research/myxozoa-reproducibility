#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 5 ]]; then
  echo "Usage: $0 SUPERMATRIX PARTITIONS_NEX TRIMMED_LOCI_DIR OUTDIR THREADS" >&2
  exit 2
fi

matrix=$1
partitions=$2
loci_dir=$3
outdir=$4
threads=$5
command -v iqtree2 >/dev/null 2>&1 || { echo "Missing required executable: iqtree2" >&2; exit 127; }
mkdir -p "$outdir/gene_trees"

iqtree2 -s "$matrix" -p "$partitions" -m MFP+MERGE -alrt 1000 -B 1000 \
  -T "$threads" --prefix "$outdir/nuclear_partitioned"

for alignment in "$loci_dir"/*.faa; do
  locus=$(basename "$alignment" .trimmed.faa)
  iqtree2 -s "$alignment" -m MFP -T "$threads" --prefix "$outdir/gene_trees/$locus"
done
cat "$outdir"/gene_trees/*.treefile > "$outdir/gene_trees.nwk"

iqtree2 -t "$outdir/nuclear_partitioned.treefile" \
  --gcf "$outdir/gene_trees.nwk" \
  -s "$matrix" --scf 100 \
  -T "$threads" --prefix "$outdir/nuclear_concordance"
