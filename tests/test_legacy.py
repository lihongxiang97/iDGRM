import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import pandas as pd

from idgrm.config import IDGRMConfig
from idgrm.io import read_expression, read_pairs, read_samples
from idgrm.legacy import load_differential_expression_results
from idgrm.pipeline import run_primary_classification_from_differential_expression


class LegacyTests(unittest.TestCase):
    def test_imports_deprecated_deseq2_shape_without_splitting_gene_ids(self):
        with TemporaryDirectory() as temp_name:
            temp = Path(temp_name)
            pair_path = temp / "pairs.tsv"
            pair_path.write_text("gene1\tgene2\nG-1\tG-2\n", encoding="utf-8")
            header = '"",baseMean,log2FoldChange,lfcSE,stat,pvalue,padj\n'
            (temp / "Ath.td.pairs.root.DESeq2.csv").write_text(
                header + '"G-1-G-2",10,-2,0.1,-20,1e-8,1e-7\n', encoding="utf-8"
            )
            (temp / "Ath.td.pairs.leaf.DESeq2.csv").write_text(
                header + '"G-1-G-2",10,2,0.1,20,1e-8,1e-7\n', encoding="utf-8"
            )
            expression = temp / "strict_tpm.tsv"
            expression.write_text(
                "gene_id\tr1\tr2\tl1\tl2\n"
                "G-1\t10\t10\t10\t10\n"
                "G-2\t10\t10\t10\t10\n", encoding="utf-8"
            )
            samples = temp / "samples.tsv"
            samples.write_text(
                "sample_id\ttissue\nr1\troot\nr2\troot\nl1\tleaf\nl2\tleaf\n",
                encoding="utf-8",
            )
            pairs = read_pairs(pair_path)
            matrix = read_expression(expression)
            evidence, warnings = load_differential_expression_results(
                temp, pairs, IDGRMConfig(),
                activity_expression=matrix,
                activity_samples=read_samples(samples, matrix.columns),
            )
            self.assertFalse(warnings)
            self.assertEqual(set(evidence["direction"]), {"gene1_high", "gene2_high"})
            self.assertEqual(set(evidence["pair_id"]), {"G-1__G-2"})

            result = run_primary_classification_from_differential_expression(
                differential_expression_dir=temp,
                duplicate_pairs_path=pair_path,
                output_dir=temp / "out",
                config=IDGRMConfig(),
                activity_expression_matrix_path=expression,
                activity_sample_metadata_path=samples,
            )
            self.assertEqual(
                result.classifications.iloc[0]["primary_class"],
                "SUB_OR_NEOFUNCTIONALIZED",
            )
            for output_path in result.output_paths.values():
                self.assertTrue(output_path.exists(), output_path)
            metadata = json.loads(result.output_paths["run_metadata"].read_text(encoding="utf-8"))
            self.assertEqual(len(metadata["inputs"]["differential_expression_files"]), 2)
            self.assertTrue(all(item["sha256"] for item in metadata["inputs"]["differential_expression_files"]))

    def test_measured_tpm_controls_activity_not_deseq2_base_mean(self):
        with TemporaryDirectory() as temp_name:
            temp = Path(temp_name)
            pair_path = temp / "pairs.tsv"
            pair_path.write_text("gene1\tgene2\nG-1\tG-2\n", encoding="utf-8")
            header = '"",baseMean,log2FoldChange,lfcSE,stat,pvalue,padj\n'
            for tissue, lfc in (("root", -2), ("leaf", 2)):
                (temp / f"Ath.td.pairs.{tissue}.DESeq2.csv").write_text(
                    header + f'"G-1-G-2",100,{lfc},0.1,20,1e-8,1e-7\n', encoding="utf-8"
                )
            expression = temp / "strict_tpm.tsv"
            expression.write_text(
                "gene_id\tr1\tr2\tl1\tl2\n"
                "G-1\t10\t12\t0\t0\n"
                "G-2\t0\t0\t0\t0\n", encoding="utf-8"
            )
            samples = temp / "samples.tsv"
            samples.write_text(
                "sample_id\ttissue\nr1\troot\nr2\troot\nl1\tleaf\nl2\tleaf\n",
                encoding="utf-8",
            )
            result = run_primary_classification_from_differential_expression(
                differential_expression_dir=temp,
                duplicate_pairs_path=pair_path,
                output_dir=temp / "out_tpm",
                activity_expression_matrix_path=expression,
                activity_sample_metadata_path=samples,
            )
            self.assertEqual(result.classifications.iloc[0]["primary_class"], "UNMAPPED")
            root = result.evidence.set_index("tissue").loc["root"]
            leaf = result.evidence.set_index("tissue").loc["leaf"]
            self.assertTrue(bool(root["active1"]))
            self.assertFalse(bool(root["active2"]))
            self.assertFalse(bool(leaf["active1"]))
            self.assertFalse(bool(leaf["active2"]))
            self.assertEqual(root["method"], "precomputed_differential_expression_with_tpm_activity")

    def test_de_input_rejects_missing_measured_activity(self):
        with TemporaryDirectory() as temp_name:
            temp = Path(temp_name)
            pair_path = temp / "pairs.tsv"
            pair_path.write_text("gene1\tgene2\nG-1\tG-2\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "measured activity expression"):
                run_primary_classification_from_differential_expression(
                    differential_expression_dir=temp,
                    duplicate_pairs_path=pair_path,
                    output_dir=temp / "out",
                )

    def test_missing_supplied_adjusted_p_is_not_recomputed_from_raw_p(self):
        with TemporaryDirectory() as temp_name:
            temp = Path(temp_name)
            pair_path = temp / "pairs.tsv"
            pair_path.write_text("gene1\tgene2\nG-1\tG-2\n", encoding="utf-8")
            header = '"",baseMean,log2FoldChange,pvalue,padj\n'
            (temp / "Ath.td.pairs.root.DESeq2.csv").write_text(
                header + '"G-1-G-2",100,-2,0.001,\n', encoding="utf-8"
            )
            (temp / "Ath.td.pairs.leaf.DESeq2.csv").write_text(
                header + '"G-1-G-2",100,2,0.001,0.01\n', encoding="utf-8"
            )
            expression = temp / "strict_tpm.tsv"
            expression.write_text(
                "gene_id\tr1\tr2\tl1\tl2\n"
                "G-1\t10\t10\t10\t10\nG-2\t10\t10\t10\t10\n",
                encoding="utf-8",
            )
            samples = temp / "samples.tsv"
            samples.write_text(
                "sample_id\ttissue\nr1\troot\nr2\troot\nl1\tleaf\nl2\tleaf\n",
                encoding="utf-8",
            )
            result = run_primary_classification_from_differential_expression(
                differential_expression_dir=temp,
                duplicate_pairs_path=pair_path,
                output_dir=temp / "out",
                activity_expression_matrix_path=expression,
                activity_sample_metadata_path=samples,
            )
            self.assertTrue(pd.isna(result.evidence.set_index("tissue").loc["root", "qvalue"]))
            self.assertEqual(int(result.classifications.iloc[0]["n_evaluable_tissues"]), 1)
            self.assertEqual(result.classifications.iloc[0]["primary_class"], "UNMAPPED")


if __name__ == "__main__":
    unittest.main()
