#!/usr/bin/env python3
"""Align loci, trim them with the selected manuscript method and concatenate taxa."""

from __future__ import annotations

import argparse
import math
import shlex
import shutil
import subprocess
from collections import OrderedDict
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def run_alignment(input_path: Path, aligned_path: Path) -> None:
    with aligned_path.open("w", encoding="utf-8") as handle:
        subprocess.run(["mafft", "--auto", str(input_path)], check=True, stdout=handle)


def run_trimmer(aligned: Path, trimmed: Path, method: str, trimal_args: str) -> None:
    if method == "none":
        shutil.copyfile(aligned, trimmed)
    elif method == "trimal":
        command = ["trimal", "-in", str(aligned), "-out", str(trimmed)] + shlex.split(trimal_args)
        subprocess.run(command, check=True)
    elif method == "gblocks":
        subprocess.run(["Gblocks", str(aligned), "-t=p"], check=True)
        result = Path(str(aligned) + "-gb")
        if not result.is_file():
            raise RuntimeError(f"Gblocks output was not created: {result}")
        shutil.move(result, trimmed)


def read_alignment(path: Path) -> OrderedDict[str, str]:
    records = OrderedDict((record.id, str(record.seq)) for record in SeqIO.parse(path, "fasta"))
    lengths = {len(sequence) for sequence in records.values()}
    if len(lengths) != 1:
        raise ValueError(f"Alignment contains inconsistent sequence lengths: {path}")
    return records


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_dir", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--genes", nargs="+", help="Ordered locus names; defaults to all *.faa stems")
    parser.add_argument("--trimmer", choices=["trimal", "gblocks", "none"], default="trimal")
    parser.add_argument("--trimal-args", default="-automated1")
    args = parser.parse_args()

    required = ["mafft"]
    if args.trimmer == "trimal":
        required.append("trimal")
    if args.trimmer == "gblocks":
        required.append("Gblocks")
    missing = [command for command in required if shutil.which(command) is None]
    if missing:
        raise SystemExit("Missing required executable(s): " + ", ".join(missing))

    genes = args.genes or sorted(path.stem for path in args.input_dir.glob("*.faa"))
    if not genes:
        parser.error("No loci were found")
    aligned_dir = args.output_dir / "aligned_loci"
    trimmed_dir = args.output_dir / "trimmed_loci"
    aligned_dir.mkdir(parents=True, exist_ok=True)
    trimmed_dir.mkdir(parents=True, exist_ok=True)

    loci: list[tuple[str, OrderedDict[str, str]]] = []
    taxa: set[str] = set()
    for gene in genes:
        source = args.input_dir / f"{gene}.faa"
        if not source.is_file():
            raise FileNotFoundError(source)
        aligned = aligned_dir / f"{gene}.aligned.faa"
        trimmed = trimmed_dir / f"{gene}.trimmed.faa"
        run_alignment(source, aligned)
        run_trimmer(aligned, trimmed, args.trimmer, args.trimal_args)
        alignment = read_alignment(trimmed)
        loci.append((gene, alignment))
        taxa.update(alignment)

    taxa_order = sorted(taxa)
    concatenated = {taxon: [] for taxon in taxa_order}
    partitions: list[tuple[str, int, int]] = []
    position = 1
    for gene, alignment in loci:
        length = len(next(iter(alignment.values())))
        partitions.append((gene, position, position + length - 1))
        for taxon in taxa_order:
            concatenated[taxon].append(alignment.get(taxon, "-" * length))
        position += length

    matrix = args.output_dir / "supermatrix.faa"
    SeqIO.write(
        [SeqRecord(Seq("".join(concatenated[taxon])), id=taxon, description="") for taxon in taxa_order],
        matrix,
        "fasta",
    )
    with (args.output_dir / "partitions.nex").open("w", encoding="utf-8") as handle:
        handle.write("#nexus\nbegin sets;\n")
        for gene, start, end in partitions:
            handle.write(f"  charset {gene} = {start}-{end};\n")
        handle.write("end;\n")
    with (args.output_dir / "matrix_summary.tsv").open("w", encoding="utf-8") as handle:
        handle.write("taxa\tloci\talignment_positions\ttrimmer\ttrimmer_parameters\n")
        handle.write(f"{len(taxa_order)}\t{len(loci)}\t{position - 1}\t{args.trimmer}\t{args.trimal_args if args.trimmer == 'trimal' else 'default'}\n")


if __name__ == "__main__":
    main()
