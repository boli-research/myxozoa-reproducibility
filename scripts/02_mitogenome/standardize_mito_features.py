#!/usr/bin/env python3
"""Standardize conserved genes and apply the manuscript ORF filter to GenBank records."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

from Bio import SeqIO


CORE = {"cox1", "cox2", "nad1", "nad5", "cytb", "rrns", "rrnl"}
ALIASES = {
    "coi": "cox1", "co1": "cox1", "cox i": "cox1",
    "coii": "cox2", "co2": "cox2", "cox ii": "cox2",
    "cob": "cytb", "cytochrome b": "cytb",
    "12s": "rrnS", "12s rrna": "rrnS", "rrns": "rrnS",
    "16s": "rrnL", "16s rrna": "rrnL", "rrnl": "rrnL",
}


def label(feature) -> str:
    values = []
    for key in ("gene", "locus_tag", "product", "label", "note"):
        values.extend(feature.qualifiers.get(key, []))
    text = " ".join(values).strip()
    cleaned = re.sub(r"[_-]+", " ", text.lower()).strip()
    if cleaned in ALIASES:
        return ALIASES[cleaned]
    for name in CORE:
        if re.search(rf"\b{re.escape(name)}\b", cleaned, re.IGNORECASE):
            return "rrnS" if name == "rrns" else "rrnL" if name == "rrnl" else name
    if feature.type == "rRNA":
        if "small" in cleaned or "12s" in cleaned:
            return "rrnS"
        if "large" in cleaned or "16s" in cleaned:
            return "rrnL"
    if feature.type == "CDS" and ("orf" in cleaned or not text):
        return "ORF"
    return text or feature.type


def positions(feature) -> set[int]:
    return set(int(position) for position in feature.location)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_genbank", type=Path)
    parser.add_argument("output_tsv", type=Path)
    parser.add_argument("--taxon", help="Override the GenBank record name")
    parser.add_argument("--min-orf-bp", type=int, default=150, help="ORFs must be strictly longer than this value")
    parser.add_argument("--max-overlap-fraction", type=float, default=0.10)
    args = parser.parse_args()

    records = list(SeqIO.parse(args.input_genbank, "genbank"))
    if len(records) != 1:
        parser.error("The input must contain exactly one GenBank record")
    record = records[0]
    taxon = args.taxon or record.annotations.get("organism") or record.id

    candidates = []
    core_features = []
    for feature in record.features:
        name = label(feature)
        if name.lower() in CORE:
            core_features.append((feature, name))
        elif feature.type == "CDS" and name == "ORF":
            candidates.append((feature, name))

    kept = list(core_features)
    for feature, name in candidates:
        feature_positions = positions(feature)
        if len(feature_positions) <= args.min_orf_bp:
            continue
        reject = False
        for core_feature, _ in core_features:
            core_positions = positions(core_feature)
            overlap = len(feature_positions & core_positions)
            denominator = min(len(feature_positions), len(core_positions))
            if denominator and overlap / denominator > args.max_overlap_fraction:
                reject = True
                break
        if not reject:
            kept.append((feature, name))

    kept.sort(key=lambda item: int(item[0].location.start))
    args.output_tsv.parent.mkdir(parents=True, exist_ok=True)
    fields = ["taxon", "record", "feature", "feature_type", "start", "end", "strand", "length_bp", "order"]
    with args.output_tsv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for order, (feature, name) in enumerate(kept, start=1):
            writer.writerow({
                "taxon": taxon,
                "record": record.id,
                "feature": name,
                "feature_type": feature.type,
                "start": int(feature.location.start) + 1,
                "end": int(feature.location.end),
                "strand": feature.location.strand or 0,
                "length_bp": len(feature.location),
                "order": order,
            })


if __name__ == "__main__":
    main()
