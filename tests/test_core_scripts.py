from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord


ROOT = Path(__file__).resolve().parents[1]
PYTHON = sys.executable


def run_script(relative: str, *arguments: object) -> None:
    subprocess.run([PYTHON, str(ROOT / relative), *map(str, arguments)], check=True, cwd=ROOT)


class CoreScriptTests(unittest.TestCase):
    def test_filter_contigs(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            temp = Path(temp)
            source = temp / "input.fasta"
            SeqIO.write([
                SeqRecord(Seq("A" * 299), id="short", description=""),
                SeqRecord(Seq("C" * 300), id="boundary", description=""),
                SeqRecord(Seq("G" * 500), id="long", description=""),
            ], source, "fasta")
            output = temp / "filtered.fasta"
            run_script("scripts/01_genome_assembly_annotation/filter_contigs.py", source, output, "--min-length", 300)
            with output.open(encoding="utf-8") as handle:
                self.assertEqual([record.id for record in SeqIO.parse(handle, "fasta")], ["boundary", "long"])

    def test_standardize_mito_orfs(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            temp = Path(temp)
            record = SeqRecord(Seq("ATG" * 400), id="mito", name="mito", description="test")
            record.annotations["molecule_type"] = "DNA"
            record.features = [
                SeqFeature(FeatureLocation(0, 300, strand=1), type="CDS", qualifiers={"gene": ["cox1"]}),
                SeqFeature(FeatureLocation(330, 510, strand=1), type="CDS", qualifiers={"gene": ["ORF1"]}),
                SeqFeature(FeatureLocation(20, 200, strand=1), type="CDS", qualifiers={"gene": ["ORF_nested"]}),
                SeqFeature(FeatureLocation(600, 750, strand=1), type="CDS", qualifiers={"gene": ["ORF_short"]}),
            ]
            genbank = temp / "mito.gb"
            SeqIO.write([record], genbank, "genbank")
            output = temp / "features.tsv"
            run_script("scripts/02_mitogenome/standardize_mito_features.py", genbank, output, "--taxon", "Test_taxon")
            data = pd.read_csv(output, sep="\t")
            self.assertEqual(data["feature"].tolist(), ["cox1", "ORF"])

    def test_classify_mito_states(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            temp = Path(temp)
            targets = temp / "targets.tsv"
            targets.write_text("gene\tmodule\tsector_order\toverall_order\nGENEA\tCI\t1\t1\nGENEB\tCI\t2\t2\nGENEC\tCI\t3\t3\n", encoding="utf-8")
            blastp = temp / "taxon.blastp.tsv"
            blastp.write_text("GENEA|ref|sp\tprotein1\t75\t100\t120\t1e-40\t200\t83.3\n", encoding="utf-8")
            tblastn = temp / "taxon.tblastn.tsv"
            tblastn.write_text("GENEB|ref|sp\tcontig1\t50\t90\t120\t1e-12\t100\t1\t90\t100\t370\n", encoding="utf-8")
            manifest = temp / "manifest.tsv"
            manifest.write_text(f"taxon\tblastp\ttblastn\nTaxon_A\t{blastp}\t{tblastn}\n", encoding="utf-8")
            prefix = temp / "states"
            run_script("scripts/04_functional_genes/classify_mito_states.py", manifest, targets, prefix)
            matrix = pd.read_csv(temp / "states.state_matrix.tsv", sep="\t")
            self.assertEqual(matrix.loc[0, "GENEA"], "Recovered")
            self.assertEqual(matrix.loc[0, "GENEB"], "Partially retained")
            self.assertEqual(matrix.loc[0, "GENEC"], "Not recovered")

    def test_fig4_through_fig7_smoke(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            temp = Path(temp)
            tree = temp / "tree.nwk"
            tree.write_text("((Taxon_A:0.1,Taxon_B:0.1)95:0.2,(Taxon_C:0.1,Taxon_D:0.1)90:0.2);\n", encoding="utf-8")

            fig4_metadata = pd.DataFrame({
                "taxon": ["Taxon_A", "Taxon_B", "Taxon_C", "Taxon_D"],
                "lineage": ["Myxozoa", "Myxozoa", "Hydrozoa", "Porifera"],
                "architecture": ["absent", "multiple_linear", "single_linear", "single_circular"],
                "new_isolate": [1, 0, 0, 0],
            })
            fig4_tsv = temp / "fig4.tsv"
            fig4_metadata.to_csv(fig4_tsv, sep="\t", index=False)
            run_script("scripts/05_figures/fig4_cnidarian_architecture.py", tree, fig4_tsv, temp / "Fig4")

            fig5_data = pd.DataFrame({
                "taxon": ["Taxon_A", "Taxon_B", "Taxon_C", "Taxon_D"],
                "group": ["new", "new", "Other Myxozoa", "Other Myxozoa"],
                "assembly_size_mb": [150, 160, 70, 30],
                "contigs": [50000, 60000, 10000, 300],
                "n50_bp": [3000, 5000, 9000, 700000],
                "gc_percent": [31, 39, 24, 30],
                "busco_single": [20, 18, 25, 40],
                "busco_duplicated": [2, 1, 3, 2],
                "busco_fragmented": [20, 21, 15, 10],
                "busco_missing": [58, 60, 57, 48],
            })
            fig5_tsv = temp / "fig5.tsv"
            fig5_data.to_csv(fig5_tsv, sep="\t", index=False)
            run_script("scripts/05_figures/fig5_genome_busco.py", fig5_tsv, temp / "Fig5")

            fig6_metadata = fig4_metadata[["taxon", "lineage", "new_isolate"]]
            fig6_tsv = temp / "fig6.tsv"
            fig6_metadata.to_csv(fig6_tsv, sep="\t", index=False)
            run_script("scripts/05_figures/fig6_nuclear_tree.py", tree, fig6_tsv, temp / "Fig6")

            mito_targets = pd.DataFrame({
                "gene": [f"G{i}" for i in range(1, 7)], "module": ["CI"] * 6,
                "sector_order": list(range(1, 7)), "overall_order": list(range(1, 7)),
            })
            mito_targets_path = temp / "mito_targets.tsv"
            mito_targets.to_csv(mito_targets_path, sep="\t", index=False)
            mito_matrix = pd.DataFrame({
                "taxon": ["Taxon_A", "Taxon_B", "Taxon_C", "Taxon_D"],
                **{f"G{i}": ["Recovered", "Partially retained", "Not recovered", "Recovered"] for i in range(1, 7)},
            })
            mito_matrix_path = temp / "mito.tsv"
            mito_matrix.to_csv(mito_matrix_path, sep="\t", index=False)
            inv_targets = pd.DataFrame({"module": ["Adhesion"] * 6, "sector": [f"D{i}" for i in range(1, 7)], "pattern": [f"D{i}" for i in range(1, 7)]})
            inv_targets_path = temp / "inv_targets.tsv"
            inv_targets.to_csv(inv_targets_path, sep="\t", index=False)
            inv_matrix = pd.DataFrame({
                "taxon": ["Taxon_A", "Taxon_B", "Taxon_C", "Taxon_D"],
                **{f"Adhesion:D{i}": ["Recovered", "Not recovered", "Recovered", "Not recovered"] for i in range(1, 7)},
            })
            inv_matrix_path = temp / "inv.tsv"
            inv_matrix.to_csv(inv_matrix_path, sep="\t", index=False)
            run_script("scripts/05_figures/fig7_coulson.py", tree, mito_matrix_path, mito_targets_path, inv_matrix_path, inv_targets_path, temp / "Fig7")

            for prefix in ("Fig4", "Fig5", "Fig6", "Fig7"):
                for suffix in ("svg", "pdf", "tiff", "png"):
                    output = temp / f"{prefix}.{suffix}"
                    self.assertTrue(output.is_file() and output.stat().st_size > 0, output)


if __name__ == "__main__":
    unittest.main()
