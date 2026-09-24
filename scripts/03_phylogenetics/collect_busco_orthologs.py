#!/usr/bin/env python3
"""Collect single-copy BUSCO proteins present in a chosen fraction of taxa."""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path, help="TSV with taxon and single_copy_dir columns")
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--minimum-fraction", type=float, default=0.70)
    args = parser.parse_args()
    if not 0 < args.minimum_fraction <= 1:
        parser.error("--minimum-fraction must be in (0, 1]")

    with args.manifest.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        parser.error("Manifest contains no taxa")

    loci: dict[str, dict[str, Path]] = defaultdict(dict)
    for row in rows:
        directory = Path(row["single_copy_dir"])
        if not directory.is_dir():
            raise FileNotFoundError(directory)
        for path in directory.glob("*.faa"):
            loci[path.stem][row["taxon"]] = path

    required = math.ceil(len(rows) * args.minimum_fraction)
    selected = sorted(locus for locus, taxa in loci.items() if len(taxa) >= required)
    raw_dir = args.output_dir / "raw_loci"
    raw_dir.mkdir(parents=True, exist_ok=True)
    for locus in selected:
        output_records = []
        for taxon, path in sorted(loci[locus].items()):
            records = list(SeqIO.parse(path, "fasta"))
            if not records:
                continue
            record = records[0]
            output_records.append(SeqRecord(record.seq, id=taxon, description=locus))
        SeqIO.write(output_records, raw_dir / f"{locus}.faa", "fasta")

    with (args.output_dir / "selected_buscos.tsv").open("w", encoding="utf-8") as handle:
        handle.write("busco_id\ttaxa_present\ttaxa_total\tpresence_fraction\n")
        for locus in selected:
            count = len(loci[locus])
            handle.write(f"{locus}\t{count}\t{len(rows)}\t{count / len(rows):.6f}\n")


if __name__ == "__main__":
    main()
