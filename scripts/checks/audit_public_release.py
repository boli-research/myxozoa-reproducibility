#!/usr/bin/env python3
"""Fail closed when a proposed public release contains common private artefacts."""

from __future__ import annotations

import argparse
import re
from pathlib import Path


SKIP_DIRS = {".git", ".venv", "work", "results", "__pycache__"}
TEXT_SUFFIXES = {
    "", ".bib", ".cff", ".csv", ".json", ".md", ".py", ".r", ".sh",
    ".tsv", ".txt", ".yaml", ".yml",
}
PATTERNS = {
    "macOS personal path": re.compile("/" + r"Users/[^/\s]+/"),
    "Linux personal path": re.compile("/" + r"home/[^/\s]+/"),
    "private key": re.compile(r"-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----"),
    "GitHub token": re.compile(r"\b(?:ghp|github_pat)_[A-Za-z0-9_]{20,}\b"),
    "generic assigned secret": re.compile(
        r"(?i)\b(?:api[_-]?key|access[_-]?token|password|secret)\s*[:=]\s*['\"]?[A-Za-z0-9_./+\-=]{12,}"
    ),
}


def iter_release_files(root: Path):
    for path in root.rglob("*"):
        if not path.is_file() or any(part in SKIP_DIRS for part in path.parts):
            continue
        yield path


def audit(root: Path, max_mb: float) -> list[str]:
    findings: list[str] = []
    max_bytes = int(max_mb * 1024 * 1024)
    for path in iter_release_files(root):
        rel = path.relative_to(root)
        size = path.stat().st_size
        if size > max_bytes:
            findings.append(f"large file ({size / 1024 / 1024:.1f} MiB): {rel}")
        if path.suffix.lower() not in TEXT_SUFFIXES:
            continue
        try:
            text = path.read_text(encoding="utf-8")
        except UnicodeDecodeError:
            findings.append(f"binary content with text-like suffix: {rel}")
            continue
        for label, pattern in PATTERNS.items():
            for match in pattern.finditer(text):
                line = text.count("\n", 0, match.start()) + 1
                findings.append(f"{label}: {rel}:{line}")
    return findings


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", nargs="?", type=Path, default=Path.cwd())
    parser.add_argument("--max-mb", type=float, default=10.0)
    args = parser.parse_args()

    findings = audit(args.root.resolve(), args.max_mb)
    if findings:
        print("Public-release audit failed:")
        for item in findings:
            print(f"- {item}")
        raise SystemExit(1)
    print("Public-release audit passed: no personal paths, credential patterns, or oversized files found.")


if __name__ == "__main__":
    main()
