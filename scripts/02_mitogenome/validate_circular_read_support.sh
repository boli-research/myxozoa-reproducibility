#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage: $0 SAMPLE MITOGENOME_FASTA READ1 READ2 OUTDIR THREADS"
}

if [[ $# -ne 6 ]]; then
  usage >&2
  exit 2
fi

sample=$1
mitogenome=$2
read1=$3
read2=$4
outdir=$5
threads=$6
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
python_bin=${PYTHON_BIN:-python3}

for command in bowtie2 bowtie2-build samtools "$python_bin"; do
  command -v "$command" >/dev/null 2>&1 || { echo "Missing required executable: $command" >&2; exit 127; }
done

mkdir -p "$outdir"
rotated="$outdir/${sample}.rotated.fasta"
"$python_bin" "$script_dir/rotate_fasta.py" "$mitogenome" "$rotated"
bowtie2-build "$rotated" "$outdir/${sample}.rotated" >/dev/null
bowtie2 --sensitive-local -x "$outdir/${sample}.rotated" -1 "$read1" -2 "$read2" \
  -p "$threads" 2> "$outdir/${sample}.mapping.log" \
  | samtools sort -@ "$threads" -o "$outdir/${sample}.rotated.bam"
samtools index "$outdir/${sample}.rotated.bam"
samtools depth -aa "$outdir/${sample}.rotated.bam" > "$outdir/${sample}.rotated.depth.tsv"
samtools coverage "$outdir/${sample}.rotated.bam" > "$outdir/${sample}.rotated.coverage.tsv"

echo "Inspect the rotated depth profile and read pairs spanning the new sequence boundary."
