#!/usr/bin/env python3
"""Remove short contigs and write an auditable length summary."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from Bio import SeqIO


def n50(lengths: list[int]) -> int:
    if not lengths:
        return 0
    target = sum(lengths) / 2
    cumulative = 0
    for length in sorted(lengths, reverse=True):
        cumulative += length
        if cumulative >= target:
            return length
    return 0


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_fasta", type=Path)
    parser.add_argument("output_fasta", type=Path)
    parser.add_argument("--min-length", type=int, default=300)
    parser.add_argument("--summary", type=Path)
    args = parser.parse_args()

    if args.min_length < 1:
        parser.error("--min-length must be positive")
    args.output_fasta.parent.mkdir(parents=True, exist_ok=True)

    input_lengths: list[int] = []
    kept_lengths: list[int] = []
    kept = []
    for record in SeqIO.parse(args.input_fasta, "fasta"):
        length = len(record.seq)
        input_lengths.append(length)
        if length >= args.min_length:
            kept.append(record)
            kept_lengths.append(length)

    SeqIO.write(kept, args.output_fasta, "fasta")
    summary = {
        "input_contigs": len(input_lengths),
        "input_bases": sum(input_lengths),
        "input_n50": n50(input_lengths),
        "minimum_length_bp": args.min_length,
        "kept_contigs": len(kept_lengths),
        "kept_bases": sum(kept_lengths),
        "kept_n50": n50(kept_lengths),
        "removed_contigs": len(input_lengths) - len(kept_lengths),
    }
    summary_path = args.summary or args.output_fasta.with_suffix(".summary.json")
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
