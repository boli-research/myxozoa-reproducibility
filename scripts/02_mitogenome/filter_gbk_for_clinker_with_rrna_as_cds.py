#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import glob
from pathlib import Path
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature

def norm(s: str) -> str:
    s = (s or "").strip().lower()
    s = s.replace(" ", "").replace("_", "").replace("-", "")
    s = s.replace("cytochromeb", "cytb")
    return s

def all_qual_text(feat) -> str:
    chunks = []
    for k, vals in feat.qualifiers.items():
        if isinstance(vals, list):
            chunks.extend([str(v) for v in vals if v is not None])
        else:
            chunks.append(str(vals))
    return " ".join(chunks).lower()

def is_known_cds(feat, keep_genes_norm: set) -> bool:
    if feat.type != "CDS":
        return False

    text = all_qual_text(feat)
    gene = feat.qualifiers.get("gene", [""])[0]
    g = norm(gene)

    if g in ("cytb", "cob") or "cytb" in text:
        g = "cob"

    return g in keep_genes_norm

def location_len_bp(loc) -> int:
    parts = getattr(loc, "parts", [loc])
    total = 0
    for p in parts:
        total += int(p.end) - int(p.start)
    return max(0, total)

def rrna_name_by_text_or_length(rrna_feats):
    """
    Input: list of rRNA features
    Output: dict {id(feat): 'rrnL' or 'rrnS'}
    Strategy:
      1) If text contains rrnS/rrnL/12S/16S -> use it
      2) Else fallback to length: longer = rrnL, shorter = rrnS
    """
    names = {}

    # First pass: text-based
    undecided = []
    for feat in rrna_feats:
        text = all_qual_text(feat)
        if any(k in text for k in ["rrnl", "16s", "large subunit"]):
            names[id(feat)] = "rrnL"
        elif any(k in text for k in ["rrns", "12s", "small subunit"]):
            names[id(feat)] = "rrnS"
        else:
            undecided.append(feat)

    # Fallback: length-based
    if undecided:
        by_len = sorted(undecided, key=lambda f: location_len_bp(f.location), reverse=True)
        if len(by_len) >= 1:
            names[id(by_len[0])] = "rrnL"
        if len(by_len) >= 2:
            names[id(by_len[1])] = "rrnS"
        # any extra (rare) still label by size
        for f in by_len[2:]:
            names[id(f)] = "rrnS"

    return names

def make_pseudo_cds_from_rrna(rrna_feat, rrna_name: str):
    L = location_len_bp(rrna_feat.location)
    aa_len = max(1, L // 3)
    fake_tr = "M" * aa_len

    quals = {
        "gene": [rrna_name],
        "product": [("16S ribosomal RNA (pseudo-CDS)" if rrna_name == "rrnL"
                     else "12S ribosomal RNA (pseudo-CDS)")],
        "translation": [fake_tr],
        "note": ["pseudo-CDS generated for clinker visualisation"],
    }
    return SeqFeature(location=rrna_feat.location, type="CDS", qualifiers=quals)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in_glob", default="*.gb*", help="Input GenBank files")
    ap.add_argument("--outdir", default="gbk_filtered", help="Output directory")
    ap.add_argument("--genes", default="cox1,cox2,nad1,nad5,cob",
                    help="Known protein CDS genes to keep")
    ap.add_argument("--keep-original-rrna", action="store_true",
                    help="Also keep original rRNA features")
    args = ap.parse_args()

    keep_norm = set(norm(x) for x in args.genes.split(","))
    if "cytb" in keep_norm:
        keep_norm.add("cob")

    outdir = Path(args.outdir)
    outdir.mkdir(exist_ok=True)

    for fp in sorted(glob.glob(args.in_glob)):
        rec = SeqIO.read(fp, "genbank")

        rrna_feats = [f for f in rec.features if f.type == "rRNA"]
        rrna_map = rrna_name_by_text_or_length(rrna_feats)

        new_feats = []
        pseudo_count = 0

        for feat in rec.features:
            if feat.type == "source":
                new_feats.append(feat)
                continue

            if is_known_cds(feat, keep_norm):
                new_feats.append(feat)
                continue

            if feat.type == "rRNA":
                name = rrna_map.get(id(feat))
                if name:
                    new_feats.append(make_pseudo_cds_from_rrna(feat, name))
                    pseudo_count += 1
                    if args.keep_original_rrna:
                        new_feats.append(feat)
                continue

        rec.features = new_feats
        out_fp = outdir / (Path(fp).stem + ".filtered.gbk")
        SeqIO.write(rec, out_fp, "genbank")

        print(f"[OK] {fp} -> {out_fp}  pseudo_rRNA_CDS={pseudo_count}")

if __name__ == "__main__":
    main()

