#!/usr/bin/env python3
"""Filter contigs using tabular BLASTN hits and explicit taxon rules."""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO


FIELDS = ["qseqid", "sseqid", "pident", "length", "evalue", "bitscore", "sscinames"]


def read_hits(path: Path, max_evalue: float, min_identity: float) -> dict[str, list[dict]]:
    hits: dict[str, list[dict]] = defaultdict(list)
    with path.open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t", fieldnames=FIELDS)
        for row in reader:
            try:
                row["pident"] = float(row["pident"])
                row["evalue"] = float(row["evalue"])
                row["bitscore"] = float(row["bitscore"])
            except (TypeError, ValueError):
                continue
            if row["evalue"] <= max_evalue and row["pident"] >= min_identity:
                hits[row["qseqid"]].append(row)
    return hits


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("assembly", type=Path)
    parser.add_argument("blast_hits", type=Path, help="BLAST outfmt 6 with qseqid sseqid pident length evalue bitscore sscinames")
    parser.add_argument("output_fasta", type=Path)
    parser.add_argument("--max-evalue", type=float, default=1e-75)
    parser.add_argument("--min-identity", type=float, default=85.0)
    parser.add_argument("--target-regex", default=r"Myxozoa|Myxobol|Cnidaria")
    parser.add_argument("--decision", choices=["best-hit", "any-nontarget"], default="best-hit")
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()

    target = re.compile(args.target_regex, re.IGNORECASE)
    hits = read_hits(args.blast_hits, args.max_evalue, args.min_identity)
    remove: set[str] = set()
    decisions: list[dict] = []

    for contig, significant in hits.items():
        ordered = sorted(significant, key=lambda row: row["bitscore"], reverse=True)
        considered = ordered[:1] if args.decision == "best-hit" else ordered
        nontarget = [row for row in considered if not target.search(row.get("sscinames") or "")]
        if nontarget:
            remove.add(contig)
        best = ordered[0]
        decisions.append({
            "contig": contig,
            "decision": "remove" if contig in remove else "keep",
            "best_subject": best["sseqid"],
            "best_taxon": best.get("sscinames", ""),
            "best_identity": best["pident"],
            "best_evalue": best["evalue"],
            "best_bitscore": best["bitscore"],
            "significant_hits": len(ordered),
        })

    args.output_fasta.parent.mkdir(parents=True, exist_ok=True)
    records = list(SeqIO.parse(args.assembly, "fasta"))
    kept = [record for record in records if record.id not in remove]
    SeqIO.write(kept, args.output_fasta, "fasta")

    report = args.report or args.output_fasta.with_suffix(".contaminant_report.tsv")
    report.parent.mkdir(parents=True, exist_ok=True)
    with report.open("w", newline="", encoding="utf-8") as handle:
        fields = list(decisions[0]) if decisions else ["contig", "decision"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(decisions)

    summary = {
        "input_contigs": len(records),
        "removed_contigs": len(remove),
        "kept_contigs": len(kept),
        "max_evalue": args.max_evalue,
        "min_identity_percent": args.min_identity,
        "target_regex": args.target_regex,
        "decision_rule": args.decision,
    }
    report.with_suffix(".json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
