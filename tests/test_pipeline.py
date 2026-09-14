from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from idgrm.pipeline import run_primary_classification, run_subtype_refinement


ROOT = Path(__file__).resolve().parent.parent


class PipelineTests(unittest.TestCase):
    def test_hierarchical_workflow(self):
        with TemporaryDirectory() as temp_name:
            temp = Path(temp_name)
            primary = run_primary_classification(
                expression_matrix_path=ROOT / "examples" / "expression.tsv",
                sample_metadata_path=ROOT / "examples" / "samples.tsv",
                duplicate_pairs_path=ROOT / "examples" / "pairs.tsv",
                output_dir=temp / "primary",
            )
            allowed = {
                "UNMAPPED", "NO_DIFFERENCE", "ASYMMETRICALLY_EXPRESSED",
                "SUB_OR_NEOFUNCTIONALIZED",
            }
            self.assertTrue(set(primary.classifications.primary_class).issubset(allowed))
            self.assertNotIn("refined_class", primary.classifications.columns)

            refined = run_subtype_refinement(
                primary_classifications_path=primary.output_paths["primary_classifications"],
                tissue_evidence_path=primary.output_paths["tissue_evidence"],
                duplicate_pairs_path=ROOT / "examples" / "pairs.tsv",
                outgroup_expression_matrix_path=ROOT / "examples" / "ancestor_expression.tsv",
                outgroup_sample_metadata_path=ROOT / "examples" / "ancestor_samples.tsv",
                output_dir=temp / "refined",
            )
            self.assertEqual(
                set(refined.classifications.parent_primary_class),
                {"ASYMMETRICALLY_EXPRESSED", "SUB_OR_NEOFUNCTIONALIZED"},
            )
            subtypes = refined.classifications.set_index("pair_id").refined_class.to_dict()
            self.assertEqual(subtypes["pair_sub"], "SUBFUNCTIONALIZATION")
            self.assertEqual(subtypes["pair_neo"], "NEOFUNCTIONALIZATION")
            self.assertEqual(
                subtypes["pair_loss"], "EXPRESSION_SILENCING_ASSOCIATED_ASYMMETRY"
            )
            self.assertEqual(
                subtypes["pair_ambiguous"], "UNRESOLVED_SUB_OR_NEOFUNCTIONALIZATION"
            )

    def test_insufficient_replicates_do_not_dilute_asymmetry_denominator(self):
        with TemporaryDirectory() as temp_name:
            temp = Path(temp_name)
            expression_path = temp / "expression.tsv"
            samples_path = temp / "samples.tsv"
            pairs_path = temp / "pairs.tsv"
            expression_path.write_text(
                "gene_id\tT1a\tT1b\tT2a\tT2b\tT3a\tT4a\tT5a\tT6a\n"
                "G1\t100\t100\t10\t10\t10\t10\t10\t10\n"
                "G2\t1\t1\t10\t10\t10\t10\t10\t10\n", encoding="utf-8",
            )
            samples_path.write_text(
                "sample_id\ttissue\nT1a\tT1\nT1b\tT1\nT2a\tT2\nT2b\tT2\n"
                "T3a\tT3\nT4a\tT4\nT5a\tT5\nT6a\tT6\n", encoding="utf-8",
            )
            pairs_path.write_text("pair_id\tgene1\tgene2\np1\tG1\tG2\n", encoding="utf-8")
            result = run_primary_classification(
                expression_matrix_path=expression_path,
                sample_metadata_path=samples_path,
                duplicate_pairs_path=pairs_path,
                output_dir=temp / "out",
            )
            row = result.classifications.iloc[0]
            self.assertEqual(row.primary_class, "ASYMMETRICALLY_EXPRESSED")
            self.assertEqual(row.n_active_tissues, 6)
            self.assertEqual(row.n_evaluable_tissues, 2)
            self.assertEqual(row.asymmetry_min_tissues, 1)


if __name__ == "__main__":
    unittest.main()
