#!/usr/bin/env python3
"""Combine QUAST and BUSCO summaries into the source table for Fig. 5."""

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path


BUSCO_RE = re.compile(
    r"C:(?P<C>[\d.]+)%\[S:(?P<S>[\d.]+)%,D:(?P<D>[\d.]+)%\],F:(?P<F>[\d.]+)%,M:(?P<M>[\d.]+)%,n:(?P<n>\d+)"
)


def parse_quast(path: Path) -> dict[str, float]:
    rows = {}
    with path.open(encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if len(row) >= 2:
                rows[row[0].strip()] = row[1].strip()

    def number(*keys: str) -> float:
        for key in keys:
            if key in rows:
                return float(str(rows[key]).replace(",", ""))
        raise KeyError(f"None of the QUAST fields were found: {keys}")

    return {
        "assembly_size_mb": number("Total length", "Total length (>= 0 bp)") / 1_000_000,
        "contigs": number("# contigs", "# contigs (>= 0 bp)"),
        "n50_bp": number("N50"),
        "gc_percent": number("GC (%)"),
    }


def parse_busco(path: Path) -> dict[str, float]:
    if path.suffix.lower() == ".json":
        data = json.loads(path.read_text(encoding="utf-8"))
        results = data.get("results", data)
        aliases = {
            "busco_complete": ["Complete percentage", "C"],
            "busco_single": ["Single copy percentage", "S"],
            "busco_duplicated": ["Multi copy percentage", "D"],
            "busco_fragmented": ["Fragmented percentage", "F"],
            "busco_missing": ["Missing percentage", "M"],
            "busco_n": ["n_markers", "n"],
        }
        parsed = {}
        for output, keys in aliases.items():
            for key in keys:
                if key in results:
                    parsed[output] = float(results[key])
                    break
        if len(parsed) == len(aliases):
            return parsed

    text = path.read_text(encoding="utf-8", errors="replace")
    match = BUSCO_RE.search(text.replace(" ", ""))
    if not match:
        raise ValueError(f"Could not find BUSCO summary line in {path}")
    values = {key: float(value) for key, value in match.groupdict().items()}
    return {
        "busco_complete": values["C"],
        "busco_single": values["S"],
        "busco_duplicated": values["D"],
        "busco_fragmented": values["F"],
        "busco_missing": values["M"],
        "busco_n": int(values["n"]),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path, help="TSV: taxon, group, quast_report, busco_summary")
    parser.add_argument("output_tsv", type=Path)
    args = parser.parse_args()

    output: list[dict] = []
    with args.manifest.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            result = {"taxon": row["taxon"], "group": row.get("group", "")}
            result.update(parse_quast(Path(row["quast_report"])))
            result.update(parse_busco(Path(row["busco_summary"])))
            output.append(result)

    fields = [
        "taxon", "group", "assembly_size_mb", "contigs", "n50_bp", "gc_percent",
        "busco_complete", "busco_single", "busco_duplicated", "busco_fragmented",
        "busco_missing", "busco_n",
    ]
    args.output_tsv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_tsv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(output)


if __name__ == "__main__":
    main()
