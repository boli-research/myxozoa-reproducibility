#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 3 ]]; then
  echo "Usage: $0 SPECIES_NAME MASKED_GENOME TRAINING_GFF" >&2
  exit 2
fi

species=$1
genome=$2
training_gff=$3
command -v autoAugTrain.pl >/dev/null 2>&1 || { echo "Missing required executable: autoAugTrain.pl" >&2; exit 127; }
[[ -n "${AUGUSTUS_CONFIG_PATH:-}" ]] || { echo "AUGUSTUS_CONFIG_PATH must point to a writable configuration directory" >&2; exit 2; }

autoAugTrain.pl --species="$species" --genome="$genome" --trainingset="$training_gff"
