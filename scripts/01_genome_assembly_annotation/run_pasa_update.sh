#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 5 ]]; then
  echo "Usage: $0 PASA_CONFIG GENOME TRANSCRIPT_FASTA EVM_GFF3 OUTDIR" >&2
  exit 2
fi

config=$1
genome=$2
transcripts=$3
evm_gff=$4
outdir=$5
command -v Launch_PASA_pipeline.pl >/dev/null 2>&1 || { echo "Missing required executable: Launch_PASA_pipeline.pl" >&2; exit 127; }
mkdir -p "$outdir"

config=$(cd "$(dirname "$config")" && pwd -P)/$(basename "$config")
genome=$(cd "$(dirname "$genome")" && pwd -P)/$(basename "$genome")
transcripts=$(cd "$(dirname "$transcripts")" && pwd -P)/$(basename "$transcripts")
evm_gff=$(cd "$(dirname "$evm_gff")" && pwd -P)/$(basename "$evm_gff")
outdir=$(cd "$outdir" && pwd -P)
cd "$outdir"

Launch_PASA_pipeline.pl -c "$config" -C -R -g "$genome" -t "$transcripts" --ALIGNERS blat,gmap --CPU 1
Launch_PASA_pipeline.pl -c "$config" -A -g "$genome" -t "$transcripts" --annots "$evm_gff"
