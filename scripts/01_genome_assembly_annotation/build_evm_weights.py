#!/usr/bin/env python3
"""Write the EVidenceModeler weights stated in the manuscript."""

from __future__ import annotations

import argparse
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--transcript-source", default="stringtie")
    parser.add_argument("--protein-source", default="exonerate")
    parser.add_argument("--ab-initio-source", default="augustus")
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        f"TRANSCRIPT\t{args.transcript_source}\t10\n"
        f"PROTEIN\t{args.protein_source}\t5\n"
        f"ABINITIO_PREDICTION\t{args.ab_initio_source}\t1\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
