#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage: $0 SAMPLE SCAFFOLDS READ1 READ2 OUTDIR THREADS [BUSCO_LINEAGE]"
}

if [[ $# -lt 6 || $# -gt 7 ]]; then
  usage >&2
  exit 2
fi

sample=$1
scaffolds=$2
read1=$3
read2=$4
outdir=$5
threads=$6
lineage=${7:-eukaryota_odb10}
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
python_bin=${PYTHON_BIN:-python3}

for command in GapFiller.pl cd-hit-est quast.py busco "$python_bin"; do
  command -v "$command" >/dev/null 2>&1 || { echo "Missing required executable: $command" >&2; exit 127; }
done

mkdir -p "$outdir/gapfiller" "$outdir/cdhit" "$outdir/qc"
library_file="$outdir/gapfiller/${sample}.libraries.txt"
printf '%s\tbowtie\t%s\t%s\t367\t0.25\tFR\n' "$sample" "$read1" "$read2" > "$library_file"

GapFiller.pl -s "$library_file" -b "$outdir/gapfiller/${sample}" -l 149 -t 4
gapfilled="$outdir/gapfiller/${sample}/${sample}.gapfilled.final.fa"
[[ -s "$gapfilled" ]] || { echo "GapFiller output not found: $gapfilled" >&2; exit 3; }

cd-hit-est -i "$gapfilled" -o "$outdir/cdhit/${sample}.nonredundant.fasta" \
  -c 0.95 -n 10 -d 0 -M 16000 -T "$threads"

"$python_bin" "$script_dir/filter_contigs.py" \
  "$outdir/cdhit/${sample}.nonredundant.fasta" \
  "$outdir/${sample}.assembly.min300.fasta" \
  --min-length 300

quast.py "$outdir/${sample}.assembly.min300.fasta" -o "$outdir/qc/quast" -t "$threads"
busco -i "$outdir/${sample}.assembly.min300.fasta" -o "${sample}_busco" \
  -m genome -l "$lineage" -c "$threads" --out_path "$outdir/qc"
