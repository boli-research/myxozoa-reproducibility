#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 7 ]]; then
  echo "Usage: $0 SAMPLE READ1 READ2 SEED OUTDIR PYTHON_BIN NOVOPLASTY_PL" >&2
  exit 2
fi

sample=$1
read1=$2
read2=$3
seed=$4
outdir=$5
python_bin=$6
novoplasty=$7
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
[[ -x "$novoplasty" || -s "$novoplasty" ]] || { echo "NOVOPlasty script not found: $novoplasty" >&2; exit 2; }

for replicate in 1 2 3; do
  run_dir="$outdir/replicate_${replicate}"
  mkdir -p "$run_dir"
  config="$run_dir/${sample}.config.txt"
  "$python_bin" "$script_dir/build_novoplasty_config.py" \
    --project "${sample}_replicate_${replicate}" --read1 "$read1" --read2 "$read2" --seed "$seed" --output "$config"
  (cd "$run_dir" && perl "$novoplasty" -c "$config" > novoplasty.log 2>&1)
done

echo "Compare the three assemblies, circularization boundaries and terminal overlaps before creating a consensus."
