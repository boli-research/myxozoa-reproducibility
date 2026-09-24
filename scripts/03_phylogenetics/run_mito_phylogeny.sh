#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 || $# -gt 3 ]]; then
  echo "Usage: $0 SUPERMATRIX OUT_PREFIX [THREADS]" >&2
  exit 2
fi

matrix=$1
prefix=$2
threads=${3:-AUTO}
command -v iqtree3 >/dev/null 2>&1 || { echo "Missing required executable: iqtree3" >&2; exit 127; }

iqtree3 -s "$matrix" -m MTZOA+F+R6 -alrt 1000 -B 1000 --bnni -T "$threads" --prefix "$prefix"
