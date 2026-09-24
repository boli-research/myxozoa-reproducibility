#!/usr/bin/env python3
"""Classify 48 mitochondrial-metabolism genes as recovered, partial or not recovered."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


BLASTP_FIELDS = ["qseqid", "sseqid", "pident", "length", "qlen", "evalue", "bitscore", "qcovs"]
TBLASTN_FIELDS = ["qseqid", "sseqid", "pident", "length", "qlen", "evalue", "bitscore", "qstart", "qend", "sstart", "send"]


def gene_from_query(query: str) -> str:
    return query.split("|", 1)[0].upper()


def read_hits(path: Path, fields: list[str]) -> list[dict]:
    hits = []
    with path.open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t", fieldnames=fields):
            try:
                for key in ("pident", "length", "qlen", "evalue", "bitscore"):
                    row[key] = float(row[key])
                if "qcovs" in row:
                    row["qcovs"] = float(row["qcovs"])
            except (TypeError, ValueError):
                continue
            row["gene"] = gene_from_query(row["qseqid"])
            hits.append(row)
    return hits


def best(rows: list[dict]) -> dict | None:
    return max(rows, key=lambda row: row["bitscore"], default=None)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path, help="TSV with taxon, blastp and tblastn columns")
    parser.add_argument("targets", type=Path, help="metadata/mitochondrial_metabolism_targets.tsv")
    parser.add_argument("output_prefix", type=Path)
    parser.add_argument("--blastp-evalue", type=float, default=1e-10)
    parser.add_argument("--minimum-query-coverage", type=float, default=0.40)
    parser.add_argument("--tblastn-evalue", type=float, default=1e-5)
    parser.add_argument("--minimum-rescue-aa", type=int, default=80)
    args = parser.parse_args()

    with args.targets.open(newline="", encoding="utf-8") as handle:
        targets = list(csv.DictReader(handle, delimiter="\t"))
    target_lookup = {row["gene"].upper(): row for row in targets}
    with args.manifest.open(newline="", encoding="utf-8") as handle:
        taxa = list(csv.DictReader(handle, delimiter="\t"))

    evidence: list[dict] = []
    wide: dict[str, dict[str, str]] = {}
    for item in taxa:
        taxon = item["taxon"]
        protein_hits = defaultdict(list)
        genome_hits = defaultdict(list)
        for row in read_hits(Path(item["blastp"]), BLASTP_FIELDS):
            if row["evalue"] <= args.blastp_evalue and row["qcovs"] / 100 >= args.minimum_query_coverage:
                protein_hits[row["gene"]].append(row)
        for row in read_hits(Path(item["tblastn"]), TBLASTN_FIELDS):
            if row["evalue"] <= args.tblastn_evalue and row["length"] >= args.minimum_rescue_aa:
                genome_hits[row["gene"]].append(row)

        wide[taxon] = {"taxon": taxon}
        for target in targets:
            gene = target["gene"].upper()
            protein = best(protein_hits[gene])
            genome = best(genome_hits[gene])
            if protein:
                state, code, basis = "Recovered", 1.0, "BLASTP"
            elif genome:
                state, code, basis = "Partially retained", 0.5, "TBLASTN rescue"
            else:
                state, code, basis = "Not recovered", 0.0, "No qualifying hit"
            wide[taxon][gene] = state
            evidence.append({
                "taxon": taxon,
                "module": target["module"],
                "gene": gene,
                "state": state,
                "state_code": code,
                "evidence_basis": basis,
                "best_protein_target": protein["sseqid"] if protein else "",
                "best_protein_evalue": protein["evalue"] if protein else "",
                "best_protein_query_coverage_percent": protein["qcovs"] if protein else "",
                "best_genome_contig": genome["sseqid"] if genome else "",
                "best_genome_evalue": genome["evalue"] if genome else "",
                "best_genome_aligned_aa": genome["length"] if genome else "",
            })

    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    evidence_path = Path(str(args.output_prefix) + ".evidence.tsv")
    with evidence_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(evidence[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(evidence)
    matrix_path = Path(str(args.output_prefix) + ".state_matrix.tsv")
    with matrix_path.open("w", newline="", encoding="utf-8") as handle:
        fields = ["taxon"] + [row["gene"].upper() for row in targets]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(wide.values())


if __name__ == "__main__":
    main()
