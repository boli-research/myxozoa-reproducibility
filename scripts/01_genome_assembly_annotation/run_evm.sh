#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 7 ]]; then
  echo "Usage: $0 GENOME GENE_PREDICTIONS_GFF3 PROTEIN_ALIGNMENTS_GFF3 TRANSCRIPT_ALIGNMENTS_GFF3 WEIGHTS OUTDIR THREADS" >&2
  exit 2
fi

genome=$1
genes=$2
proteins=$3
transcripts=$4
weights=$5
outdir=$6
threads=$7
for command in partition_EVM_inputs.pl write_EVM_commands.pl execute_EVM_commands.pl recombine_EVM_partial_outputs.pl convert_EVM_outputs_to_GFF3.pl; do
  command -v "$command" >/dev/null 2>&1 || { echo "Missing EVM executable: $command" >&2; exit 127; }
done
mkdir -p "$outdir"

genome=$(cd "$(dirname "$genome")" && pwd -P)/$(basename "$genome")
genes=$(cd "$(dirname "$genes")" && pwd -P)/$(basename "$genes")
proteins=$(cd "$(dirname "$proteins")" && pwd -P)/$(basename "$proteins")
transcripts=$(cd "$(dirname "$transcripts")" && pwd -P)/$(basename "$transcripts")
weights=$(cd "$(dirname "$weights")" && pwd -P)/$(basename "$weights")
outdir=$(cd "$outdir" && pwd -P)
cd "$outdir"

partition_EVM_inputs.pl --genome "$genome" --gene_predictions "$genes" \
  --protein_alignments "$proteins" --transcript_alignments "$transcripts" \
  --segmentSize 100000 --overlapSize 10000 --partition_listing partitions_list.out

write_EVM_commands.pl --genome "$genome" --weights "$weights" --gene_predictions "$genes" \
  --protein_alignments "$proteins" --transcript_alignments "$transcripts" \
  --output_file_name evm.out --partitions partitions_list.out > commands.list

execute_EVM_commands.pl commands.list --CPU "$threads"
recombine_EVM_partial_outputs.pl --partitions partitions_list.out --output_file_name evm.out
convert_EVM_outputs_to_GFF3.pl --partitions partitions_list.out --output_file_name evm.out --genome "$genome"
