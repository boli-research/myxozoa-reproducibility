#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 5 ]]; then
  echo "Usage: $0 TAXON PROTEINS_FAA TARGET_HMM OUTDIR THREADS" >&2
  exit 2
fi

taxon=$1
proteins=$2
hmm=$3
outdir=$4
threads=$5
for command in hmmpress hmmsearch; do
  command -v "$command" >/dev/null 2>&1 || { echo "Missing required executable: $command" >&2; exit 127; }
done
mkdir -p "$outdir"

[[ -s "${hmm}.h3f" ]] || hmmpress "$hmm"
hmmsearch --cut_ga --cpu "$threads" \
  --domtblout "$outdir/${taxon}.domtblout" \
  -o "$outdir/${taxon}.hmmsearch.txt" \
  "$hmm" "$proteins"
