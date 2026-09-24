#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 || $# -gt 3 ]]; then
  echo "Usage: $0 INPUT_FASTA OUT_PREFIX [THREADS]" >&2
  exit 2
fi

input=$1
prefix=$2
threads=${3:-AUTO}
for command in mafft trimal iqtree2; do
  command -v "$command" >/dev/null 2>&1 || { echo "Missing required executable: $command" >&2; exit 127; }
done

mafft --auto "$input" > "${prefix}.aligned.fasta"
trimal -in "${prefix}.aligned.fasta" -out "${prefix}.trimmed.fasta" -automated1
iqtree2 -s "${prefix}.trimmed.fasta" -m MFP -alrt 1000 -B 1000 -T "$threads" --prefix "$prefix"
