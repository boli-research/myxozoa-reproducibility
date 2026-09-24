#!/usr/bin/env python3
"""Convert HMMER domtblout files into sector counts and recovered/not-recovered states."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def read_domtblout(path: Path) -> list[dict]:
    hits = []
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split(maxsplit=22)
            if len(fields) < 22:
                continue
            hits.append({
                "target": fields[0],
                "hmm_name": fields[3],
                "hmm_accession": fields[4],
                "full_evalue": float(fields[6]),
                "full_score": float(fields[7]),
                "domain_i_evalue": float(fields[12]),
                "domain_score": float(fields[13]),
                "alignment_accuracy": float(fields[21]),
            })
    return hits


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path, help="TSV with taxon and domtblout columns")
    parser.add_argument("profiles_tsv", type=Path)
    parser.add_argument("output_prefix", type=Path)
    args = parser.parse_args()

    with args.profiles_tsv.open(newline="", encoding="utf-8") as handle:
        profiles = list(csv.DictReader(handle, delimiter="\t"))
    lookup: dict[str, list[dict[str, str]]] = defaultdict(list)
    sectors: list[tuple[str, str]] = []
    for row in profiles:
        lookup[row["hmm_name"]].append(row)
        pair = (row["module"], row["sector"])
        if pair not in sectors:
            sectors.append(pair)
    with args.manifest.open(newline="", encoding="utf-8") as handle:
        taxa = list(csv.DictReader(handle, delimiter="\t"))

    evidence = []
    matrix = []
    for item in taxa:
        taxon = item["taxon"]
        grouped: dict[tuple[str, str], list[dict]] = defaultdict(list)
        for hit in read_domtblout(Path(item["domtblout"])):
            for profile in lookup.get(hit["hmm_name"], []):
                grouped[(profile["module"], profile["sector"])].append({**hit, **profile})
        row = {"taxon": taxon}
        for module, sector in sectors:
            hits = grouped[(module, sector)]
            targets = {hit["target"] for hit in hits}
            state = "Recovered" if targets else "Not recovered"
            row[f"{module}:{sector}"] = state
            best = min(hits, key=lambda hit: hit["domain_i_evalue"], default=None)
            evidence.append({
                "taxon": taxon,
                "module": module,
                "sector": sector,
                "state": state,
                "unique_target_count": len(targets),
                "domtblout_rows": len(hits),
                "best_target": best["target"] if best else "",
                "best_hmm_name": best["hmm_name"] if best else "",
                "best_hmm_accession": best["hmm_accession"] if best else "",
                "best_domain_i_evalue": best["domain_i_evalue"] if best else "",
                "best_domain_score": best["domain_score"] if best else "",
            })
        matrix.append(row)

    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    evidence_path = Path(str(args.output_prefix) + ".evidence.tsv")
    with evidence_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(evidence[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(evidence)
    matrix_path = Path(str(args.output_prefix) + ".presence_matrix.tsv")
    with matrix_path.open("w", newline="", encoding="utf-8") as handle:
        fields = ["taxon"] + [f"{module}:{sector}" for module, sector in sectors]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(matrix)


if __name__ == "__main__":
    main()
