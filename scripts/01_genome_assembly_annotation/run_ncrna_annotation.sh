#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 5 ]]; then
  echo "Usage: $0 GENOME RFAM_BLAST_DB OUTDIR SAMPLE THREADS" >&2
  exit 2
fi

genome=$1
rfam_db=$2
outdir=$3
sample=$4
threads=$5
for command in tRNAscan-SE blastn; do
  command -v "$command" >/dev/null 2>&1 || { echo "Missing required executable: $command" >&2; exit 127; }
done
mkdir -p "$outdir"

tRNAscan-SE --thread "$threads" -o "$outdir/${sample}.trnascan.tsv" "$genome"
blastn -query "$genome" -db "$rfam_db" -num_threads "$threads" -evalue 1e-5 \
  -outfmt '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore' \
  -out "$outdir/${sample}.rfam_blastn.tsv"
