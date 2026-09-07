from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from idgrm.config import IDGRMConfig
from idgrm.pipeline import run_analysis


ROOT = Path(__file__).resolve().parent.parent


class PipelineTests(unittest.TestCase):
    def test_demo_with_ancestor(self):
        with TemporaryDirectory() as temp:
            result = run_analysis(
                expression_path=ROOT / "examples" / "expression.tsv",
                samples_path=ROOT / "examples" / "samples.tsv",
                pairs_path=ROOT / "examples" / "pairs.tsv",
                ancestor_expression_path=ROOT / "examples" / "ancestor_expression.tsv",
                ancestor_samples_path=ROOT / "examples" / "ancestor_samples.tsv",
                output_dir=Path(temp) / "out",
                config=IDGRMConfig(),
            )
            fates = result.classifications.set_index("pair_id")["extended_fate"].to_dict()
            self.assertEqual(fates["pair_sub"], "SUBFUNCTIONALIZATION")
            self.assertEqual(fates["pair_neo"], "NEOFUNCTIONALIZATION")
            self.assertEqual(fates["pair_aed"], "AED")
            self.assertEqual(fates["pair_conserved"], "NO_DIFFERENCE")
            self.assertEqual(fates["pair_loss"], "EXPRESSION_LOSS")
            self.assertEqual(fates["pair_ambiguous"], "AMBIGUOUS_SUB_NEO")
            for path in result.output_paths.values():
                self.assertTrue(path.exists(), path)

    def test_demo_without_ancestor_is_marked_proxy(self):
        with TemporaryDirectory() as temp:
            result = run_analysis(
                expression_path=ROOT / "examples" / "expression.tsv",
                samples_path=ROOT / "examples" / "samples.tsv",
                pairs_path=ROOT / "examples" / "pairs.tsv",
                output_dir=Path(temp) / "out",
            )
            indexed = result.classifications.set_index("pair_id")
            self.assertEqual(indexed.loc["pair_sub", "evidence_basis"], "expression_only_proxy")
            self.assertEqual(indexed.loc["pair_neo", "evidence_basis"], "expression_only_proxy")

    def test_insufficient_replicates_do_not_dilute_aed_denominator(self):
        with TemporaryDirectory() as temp_name:
            temp = Path(temp_name)
            expression_path = temp / "expression.tsv"
            samples_path = temp / "samples.tsv"
            pairs_path = temp / "pairs.tsv"
            expression_path.write_text(
                "gene_id\tT1a\tT1b\tT2a\tT2b\tT3a\tT4a\tT5a\tT6a\n"
                "G1\t100\t100\t10\t10\t10\t10\t10\t10\n"
                "G2\t1\t1\t10\t10\t10\t10\t10\t10\n",
                encoding="utf-8",
            )
            samples_path.write_text(
                "sample_id\ttissue\n"
                "T1a\tT1\nT1b\tT1\nT2a\tT2\nT2b\tT2\n"
                "T3a\tT3\nT4a\tT4\nT5a\tT5\nT6a\tT6\n",
                encoding="utf-8",
            )
            pairs_path.write_text("pair_id\tgene1\tgene2\np1\tG1\tG2\n", encoding="utf-8")
            result = run_analysis(
                expression_path=expression_path,
                samples_path=samples_path,
                pairs_path=pairs_path,
                output_dir=temp / "out",
            )
            row = result.classifications.iloc[0]
            self.assertEqual(row["science_class"], "AED")
            self.assertEqual(row["n_active_tissues"], 6)
            self.assertEqual(row["n_evaluable_tissues"], 2)
            self.assertEqual(row["aed_min_tissues"], 1)


if __name__ == "__main__":
    unittest.main()
