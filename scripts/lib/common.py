from __future__ import annotations

import csv
import shutil
from pathlib import Path
from typing import Iterable


def existing_file(value: str | Path) -> Path:
    path = Path(value)
    if not path.is_file():
        raise FileNotFoundError(f"Input file does not exist: {path}")
    return path


def ensure_parent(value: str | Path) -> Path:
    path = Path(value)
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def require_commands(names: Iterable[str]) -> None:
    missing = [name for name in names if shutil.which(name) is None]
    if missing:
        raise SystemExit("Missing required executable(s): " + ", ".join(missing))


def delimiter_for(path: str | Path) -> str:
    return "," if Path(path).suffix.lower() == ".csv" else "\t"


def read_records(path: str | Path) -> list[dict[str, str]]:
    with existing_file(path).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter_for(path)))


def write_records(path: str | Path, fieldnames: list[str], rows: Iterable[dict]) -> None:
    output = ensure_parent(path)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter=delimiter_for(path))
        writer.writeheader()
        writer.writerows(rows)
