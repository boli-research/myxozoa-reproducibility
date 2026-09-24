#!/usr/bin/env python3
"""Create a NOVOPlasty 4.x mitochondrial assembly configuration file."""

from __future__ import annotations

import argparse
from pathlib import Path


TEMPLATE = """Project:
-----------------------
Project name          = {project}
Type                  = mito
Genome Range          = {minimum}-{maximum}
K-mer                 = {kmer}
Max memory            =
Extended log          = 0
Save assembled reads  = no
Seed Input            = {seed}
Reference sequence    =
Variance detection    =
Chloroplast sequence  =

Dataset 1:
-----------------------
Read Length           = {read_length}
Insert size           = {insert_size}
Platform              = illumina
Single/Paired         = PE
Combined reads        =
Forward reads         = {read1}
Reverse reads         = {read2}
Store Hash            =

Optional:
-----------------------
Insert size auto      = yes
Use Quality Scores    = no
"""


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", required=True)
    parser.add_argument("--read1", type=Path, required=True)
    parser.add_argument("--read2", type=Path, required=True)
    parser.add_argument("--seed", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--insert-size", type=int, default=367)
    parser.add_argument("--read-length", type=int, default=150)
    parser.add_argument("--genome-min", type=int, default=12000)
    parser.add_argument("--genome-max", type=int, default=50000)
    parser.add_argument("--kmer", type=int, default=39)
    args = parser.parse_args()

    for path in (args.read1, args.read2, args.seed):
        if not path.is_file():
            parser.error(f"Input file does not exist: {path}")
    if args.genome_min >= args.genome_max:
        parser.error("--genome-min must be smaller than --genome-max")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(TEMPLATE.format(
        project=args.project,
        minimum=args.genome_min,
        maximum=args.genome_max,
        kmer=args.kmer,
        seed=args.seed.resolve(),
        read_length=args.read_length,
        insert_size=args.insert_size,
        read1=args.read1.resolve(),
        read2=args.read2.resolve(),
    ), encoding="utf-8")


if __name__ == "__main__":
    main()
