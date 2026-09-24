#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 6 || $# -gt 8 ]]; then
  echo "Usage: $0 SAMPLE MASKED_GENOME REFERENCE_PROTEINS AUGUSTUS_SPECIES OUTDIR THREADS [RNA_R1 RNA_R2]" >&2
  exit 2
fi

sample=$1
genome=$2
proteins=$3
augustus_species=$4
outdir=$5
threads=$6
rna1=${7:-}
rna2=${8:-}
for command in augustus exonerate; do
  command -v "$command" >/dev/null 2>&1 || { echo "Missing required executable: $command" >&2; exit 127; }
done
mkdir -p "$outdir/ab_initio" "$outdir/homology" "$outdir/transcript"

augustus --species="$augustus_species" --gff3=on "$genome" > "$outdir/ab_initio/${sample}.augustus.gff3"

exonerate --model protein2genome --showalignment false --showtargetgff true --showquerygff true \
  "$proteins" "$genome" > "$outdir/homology/${sample}.exonerate.gff"

if [[ -n "$rna1" && -n "$rna2" && "$rna1" != "NA" && "$rna2" != "NA" ]]; then
  for command in hisat2-build hisat2 samtools stringtie gffread TransDecoder.LongOrfs TransDecoder.Predict; do
    command -v "$command" >/dev/null 2>&1 || { echo "Missing transcript-evidence executable: $command" >&2; exit 127; }
  done
  hisat2-build "$genome" "$outdir/transcript/${sample}.genome"
  hisat2 -x "$outdir/transcript/${sample}.genome" -1 "$rna1" -2 "$rna2" -p "$threads" \
    | samtools sort -@ "$threads" -o "$outdir/transcript/${sample}.rna.bam"
  samtools index "$outdir/transcript/${sample}.rna.bam"
  stringtie "$outdir/transcript/${sample}.rna.bam" -p "$threads" -o "$outdir/transcript/${sample}.stringtie.gtf"
  gffread "$outdir/transcript/${sample}.stringtie.gtf" -g "$genome" -w "$outdir/transcript/${sample}.transcripts.fasta"
  TransDecoder.LongOrfs -t "$outdir/transcript/${sample}.transcripts.fasta"
  TransDecoder.Predict -t "$outdir/transcript/${sample}.transcripts.fasta"
fi
