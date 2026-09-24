#!/usr/bin/env python3
"""Combine standardized mitochondrial feature TSV files with schema validation."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_tsv", type=Path)
    parser.add_argument("input_tsv", nargs="+", type=Path)
    args = parser.parse_args()
    tables = [pd.read_csv(path, sep="\t") for path in args.input_tsv]
    columns = list(tables[0].columns)
    for path, table in zip(args.input_tsv, tables):
        if list(table.columns) != columns:
            parser.error(f"Column mismatch in {path}")
    args.output_tsv.parent.mkdir(parents=True, exist_ok=True)
    pd.concat(tables, ignore_index=True).to_csv(args.output_tsv, sep="\t", index=False)


if __name__ == "__main__":
    main()
