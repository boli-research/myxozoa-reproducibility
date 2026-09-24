#!/usr/bin/env python3
"""Extract Pfam HMM profiles matching the 30 invasion/adhesion sectors."""

from __future__ import annotations

import argparse
import csv
import re
import shutil
import subprocess
from pathlib import Path


def parse_block(block: str) -> tuple[str, str]:
    name = ""
    accession = ""
    for line in block.splitlines():
        if line.startswith("NAME"):
            name = line.split(maxsplit=1)[1].strip()
        elif line.startswith("ACC"):
            accession = line.split(maxsplit=1)[1].strip()
    return name, accession


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pfam_hmm", type=Path)
    parser.add_argument("mapping_tsv", type=Path)
    parser.add_argument("output_hmm", type=Path)
    parser.add_argument("--profiles-tsv", type=Path)
    parser.add_argument("--press", action="store_true")
    args = parser.parse_args()

    with args.mapping_tsv.open(newline="", encoding="utf-8") as handle:
        mapping = list(csv.DictReader(handle, delimiter="\t"))
    blocks = args.pfam_hmm.read_text(encoding="utf-8", errors="replace").split("//")
    selected: list[str] = []
    profiles: list[dict[str, str]] = []
    seen: set[str] = set()
    for block in blocks:
        name, accession = parse_block(block)
        if not name:
            continue
        for row in mapping:
            if re.search(row["pattern"], name, flags=re.IGNORECASE):
                if name not in seen:
                    selected.append(block.rstrip() + "\n//\n")
                    seen.add(name)
                profiles.append({
                    "hmm_name": name,
                    "accession": accession,
                    "module": row["module"],
                    "sector": row["sector"],
                    "matched_pattern": row["pattern"],
                })

    args.output_hmm.parent.mkdir(parents=True, exist_ok=True)
    args.output_hmm.write_text("".join(selected), encoding="utf-8")
    profile_path = args.profiles_tsv or args.output_hmm.with_suffix(".profiles.tsv")
    with profile_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["hmm_name", "accession", "module", "sector", "matched_pattern"], delimiter="\t")
        writer.writeheader()
        writer.writerows(profiles)
    if args.press:
        if shutil.which("hmmpress") is None:
            raise SystemExit("Missing required executable: hmmpress")
        subprocess.run(["hmmpress", "-f", str(args.output_hmm)], check=True)


if __name__ == "__main__":
    main()
