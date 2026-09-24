#!/usr/bin/env python3
"""Extract the five conserved mitochondrial proteins from annotated GenBank files."""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


GENES = ["cox1", "nad1", "cox2", "nad5", "cytb"]
ALIASES = {
    "coi": "cox1", "co1": "cox1", "cox i": "cox1",
    "coii": "cox2", "co2": "cox2", "cox ii": "cox2",
    "cob": "cytb", "cytochrome b": "cytb",
}


def gene_name(feature) -> str | None:
    values = []
    for key in ("gene", "product", "locus_tag", "label"):
        values.extend(feature.qualifiers.get(key, []))
    text = re.sub(r"[_-]+", " ", " ".join(values).lower()).strip()
    if text in ALIASES:
        return ALIASES[text]
    for gene in GENES:
        if re.search(rf"\b{gene}\b", text):
            return gene
    return None


def protein_sequence(feature, record) -> str:
    translations = feature.qualifiers.get("translation", [])
    if translations:
        return re.sub(r"\s+", "", translations[0]).rstrip("*")
    nucleotide = feature.extract(record.seq)
    codon_start = int(feature.qualifiers.get("codon_start", ["1"])[0]) - 1
    return str(nucleotide[codon_start:].translate(table=4, to_stop=True))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path, help="TSV with taxon and genbank columns")
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    sequences: dict[str, list[SeqRecord]] = defaultdict(list)
    missing: list[dict[str, str]] = []

    with args.manifest.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    for row in rows:
        taxon = row["taxon"]
        records = list(SeqIO.parse(row["genbank"], "genbank"))
        if len(records) != 1:
            raise ValueError(f"Expected one record for {taxon}: {row['genbank']}")
        record = records[0]
        found: dict[str, str] = {}
        for feature in record.features:
            if feature.type != "CDS":
                continue
            name = gene_name(feature)
            if name and name not in found:
                found[name] = protein_sequence(feature, record)
        for gene in GENES:
            if gene in found:
                sequences[gene].append(SeqRecord(Seq(found[gene]), id=taxon, description=gene))
            else:
                missing.append({"taxon": taxon, "gene": gene})

    for gene in GENES:
        SeqIO.write(sequences[gene], args.output_dir / f"{gene}.faa", "fasta")
    with (args.output_dir / "missing_genes.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["taxon", "gene"], delimiter="\t")
        writer.writeheader()
        writer.writerows(missing)


if __name__ == "__main__":
    main()
