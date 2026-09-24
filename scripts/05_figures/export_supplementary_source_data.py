#!/usr/bin/env python3
"""Export the approved Table S3/S4 sheets as tab-separated Fig. 7 source data."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def export_sheet(workbook: Path, sheet: str, output: Path, rename_species: bool = True) -> None:
    data = pd.read_excel(workbook, sheet_name=sheet)
    if rename_species and "species" in data.columns:
        data = data.rename(columns={"species": "taxon"})
    output.parent.mkdir(parents=True, exist_ok=True)
    data.to_csv(output, sep="\t", index=False, na_rep="")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("table_s3", type=Path)
    parser.add_argument("table_s4", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()

    for path in (args.table_s3, args.table_s4):
        if not path.is_file():
            parser.error(f"Workbook not found: {path}")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    export_sheet(args.table_s3, "Gene metadata", args.output_dir / "Fig7_mito_targets.tsv", rename_species=False)
    export_sheet(args.table_s3, "State matrix", args.output_dir / "Fig7_mito_state_matrix.tsv")
    export_sheet(args.table_s3, "Detection evidence", args.output_dir / "Fig7_mito_detection_evidence.tsv")
    export_sheet(args.table_s4, "Domain mapping", args.output_dir / "Fig7_invasion_targets.tsv", rename_species=False)
    export_sheet(args.table_s4, "Presence matrix", args.output_dir / "Fig7_invasion_state_matrix.tsv")
    export_sheet(args.table_s4, "HMM best evidence", args.output_dir / "Fig7_invasion_detection_evidence.tsv")


if __name__ == "__main__":
    main()
