#!/usr/bin/env python3
"""Rotate a single circular FASTA record without changing sequence content."""

from __future__ import annotations

import argparse
from pathlib import Path

from Bio import SeqIO


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_fasta", type=Path)
    parser.add_argument("output_fasta", type=Path)
    parser.add_argument("--offset", type=int, help="0-based rotation offset; default is half the sequence length")
    args = parser.parse_args()

    records = list(SeqIO.parse(args.input_fasta, "fasta"))
    if len(records) != 1:
        parser.error("The input must contain exactly one FASTA record")
    record = records[0]
    offset = args.offset if args.offset is not None else len(record.seq) // 2
    if not 0 <= offset < len(record.seq):
        parser.error("--offset must be within the sequence")
    record.seq = record.seq[offset:] + record.seq[:offset]
    record.id = f"{record.id}_rot{offset}"
    record.description = f"circular rotation offset={offset} original_length={len(record.seq)}"
    args.output_fasta.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write([record], args.output_fasta, "fasta")


if __name__ == "__main__":
    main()
