#!/usr/bin/env python3
"""Calculate circular sliding-window GC content and GC skew."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from Bio import SeqIO


def circular_window(sequence: str, start: int, size: int) -> str:
    length = len(sequence)
    return "".join(sequence[(start + index) % length] for index in range(size))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_sequence", type=Path, help="Single-record FASTA or GenBank file")
    parser.add_argument("output_tsv", type=Path)
    parser.add_argument("--window", type=int, default=1000)
    parser.add_argument("--step", type=int, default=100)
    args = parser.parse_args()

    fmt = "genbank" if args.input_sequence.suffix.lower() in {".gb", ".gbk", ".gbff"} else "fasta"
    records = list(SeqIO.parse(args.input_sequence, fmt))
    if len(records) != 1:
        parser.error("The input must contain exactly one sequence")
    sequence = str(records[0].seq).upper()
    if args.window < 1 or args.step < 1 or args.window > len(sequence):
        parser.error("Require 1 <= window <= sequence length and step >= 1")

    args.output_tsv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_tsv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["sequence", "position_bp", "window_bp", "gc_content", "gc_skew"], delimiter="\t")
        writer.writeheader()
        for start in range(0, len(sequence), args.step):
            window = circular_window(sequence, start, args.window)
            g = window.count("G")
            c = window.count("C")
            valid = sum(window.count(base) for base in "ACGT")
            writer.writerow({
                "sequence": records[0].id,
                "position_bp": start + 1,
                "window_bp": args.window,
                "gc_content": (g + c) / valid if valid else 0.0,
                "gc_skew": (g - c) / (g + c) if g + c else 0.0,
            })


if __name__ == "__main__":
    main()
