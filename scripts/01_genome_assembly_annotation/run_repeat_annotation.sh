#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 4 || $# -gt 5 ]]; then
  echo "Usage: $0 SAMPLE ASSEMBLY OUTDIR THREADS [REPBASE_FASTA]" >&2
  exit 2
fi

sample=$1
assembly=$2
outdir=$3
threads=$4
repbase=${5:-}
for command in BuildDatabase RepeatModeler RepeatMasker; do
  command -v "$command" >/dev/null 2>&1 || { echo "Missing required executable: $command" >&2; exit 127; }
done
mkdir -p "$outdir/modeler" "$outdir/masker"

assembly=$(cd "$(dirname "$assembly")" && pwd -P)/$(basename "$assembly")
outdir=$(cd "$outdir" && pwd -P)

cd "$outdir/modeler"
BuildDatabase -name "$sample" "$assembly"
RepeatModeler -database "$sample" -pa "$threads" -LTRStruct
consensi=$(find "$outdir/modeler" -name 'consensi.fa.classified' -print -quit)
[[ -s "$consensi" ]] || { echo "RepeatModeler classified consensus library not found" >&2; exit 3; }

library="$outdir/${sample}.repeat_library.fasta"
cp "$consensi" "$library"
if [[ -n "$repbase" ]]; then
  [[ -s "$repbase" ]] || { echo "Repbase file not found: $repbase" >&2; exit 2; }
  cat "$repbase" >> "$library"
fi

RepeatMasker -pa "$threads" -lib "$library" -dir "$outdir/masker" -gff "$assembly"
