#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 6 ]]; then
  echo "Usage: $0 TAXON PROTEINS_FAA GENOME_FNA REFERENCE_FAA OUTDIR THREADS" >&2
  exit 2
fi

taxon=$1
proteins=$2
genome=$3
references=$4
outdir=$5
threads=$6
for command in makeblastdb blastp tblastn; do
  command -v "$command" >/dev/null 2>&1 || { echo "Missing required executable: $command" >&2; exit 127; }
done
mkdir -p "$outdir/db" "$outdir/hits"

makeblastdb -in "$proteins" -dbtype prot -out "$outdir/db/${taxon}.proteins" >/dev/null
makeblastdb -in "$genome" -dbtype nucl -out "$outdir/db/${taxon}.genome" >/dev/null

blastp -query "$references" -db "$outdir/db/${taxon}.proteins" \
  -evalue 1e-10 -num_threads "$threads" -max_target_seqs 50 \
  -outfmt '6 qseqid sseqid pident length qlen evalue bitscore qcovs' \
  -out "$outdir/hits/${taxon}.blastp.tsv"

tblastn -query "$references" -db "$outdir/db/${taxon}.genome" \
  -evalue 1e-5 -num_threads "$threads" -max_target_seqs 50 \
  -outfmt '6 qseqid sseqid pident length qlen evalue bitscore qstart qend sstart send' \
  -out "$outdir/hits/${taxon}.tblastn.tsv"
