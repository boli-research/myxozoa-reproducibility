#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage: $0 SAMPLE READ1 READ2 HOST_INDEX MICROBIAL_INDEX OUTDIR THREADS"
}

if [[ $# -ne 7 ]]; then
  usage >&2
  exit 2
fi

sample=$1
read1=$2
read2=$3
host_index=$4
microbial_index=$5
outdir=$6
threads=$7

for command in fastp bowtie2 spades.py; do
  command -v "$command" >/dev/null 2>&1 || { echo "Missing required executable: $command" >&2; exit 127; }
done
for file in "$read1" "$read2"; do
  [[ -s "$file" ]] || { echo "Missing or empty read file: $file" >&2; exit 2; }
done

mkdir -p "$outdir/qc" "$outdir/decontamination" "$outdir/spades"

fastp \
  -i "$read1" -I "$read2" \
  -o "$outdir/qc/${sample}.R1.clean.fastq.gz" \
  -O "$outdir/qc/${sample}.R2.clean.fastq.gz" \
  --qualified_quality_phred 15 \
  --unqualified_percent_limit 40 \
  --n_base_limit 5 \
  --length_required 15 \
  --thread "$threads" \
  --json "$outdir/qc/${sample}.fastp.json" \
  --html "$outdir/qc/${sample}.fastp.html"

bowtie2 --sensitive --end-to-end --threads "$threads" \
  -x "$host_index" \
  -1 "$outdir/qc/${sample}.R1.clean.fastq.gz" \
  -2 "$outdir/qc/${sample}.R2.clean.fastq.gz" \
  --un-conc-gz "$outdir/decontamination/${sample}.host_unmapped.R%.fastq.gz" \
  -S /dev/null \
  2> "$outdir/decontamination/${sample}.host_bowtie2.log"

bowtie2 --sensitive --end-to-end --threads "$threads" \
  -x "$microbial_index" \
  -1 "$outdir/decontamination/${sample}.host_unmapped.R1.fastq.gz" \
  -2 "$outdir/decontamination/${sample}.host_unmapped.R2.fastq.gz" \
  --un-conc-gz "$outdir/decontamination/${sample}.clean.R%.fastq.gz" \
  -S /dev/null \
  2> "$outdir/decontamination/${sample}.microbial_bowtie2.log"

spades.py \
  -1 "$outdir/decontamination/${sample}.clean.R1.fastq.gz" \
  -2 "$outdir/decontamination/${sample}.clean.R2.fastq.gz" \
  -t "$threads" \
  -o "$outdir/spades"

echo "$outdir/spades/scaffolds.fasta"
