"""Rule boundaries, honest uncertainty and independent public-module integration."""
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import pandas as pd

from idgrm import PolarizationConfig, run_outgroup_polarization, run_primary_classification
from idgrm.polarization import _state, polarize_patterns
from idgrm.statistics import benjamini_hochberg

ROOT = Path(__file__).resolve().parent.parent


def states(rows):
    return pd.DataFrame(rows, columns=["copy1_state", "copy2_state", "outgroup_state", "direction"])


class PolarizationTests(unittest.TestCase):
    def test_complementary_partition(self):
        label, _, _ = polarize_patterns(states([
            ("ON", "OFF", "ON", "gene1_high"), ("OFF", "ON", "ON", "gene2_high")]), PolarizationConfig())
        self.assertEqual(label, "OUTGROUP_SUPPORTED_SUB_CANDIDATE")

    def test_one_copy_gain(self):
        label, _, metrics = polarize_patterns(states([
            ("ON", "OFF", "ON", "gene1_high"), ("OFF", "ON", "OFF", "gene2_high")]), PolarizationConfig())
        self.assertEqual(label, "OUTGROUP_SUPPORTED_NEO_CANDIDATE")
        self.assertEqual(metrics["candidate_gain_copy"], "copy2")

    def test_unknown_is_not_absence(self):
        label, _, _ = polarize_patterns(states([
            ("ON", "OFF", "ON", "gene1_high"), ("OFF", "ON", "UNKNOWN", "gene2_high")]), PolarizationConfig())
        self.assertEqual(label, "UNRESOLVED")

    def test_both_copy_gain_is_not_one_copy_neo(self):
        label, _, _ = polarize_patterns(states([
            ("ON", "ON", "ON", "no_difference"), ("ON", "ON", "OFF", "gene2_high")]), PolarizationConfig())
        self.assertEqual(label, "UNRESOLVED")

    def test_mixed_signal_unresolved(self):
        label, _, _ = polarize_patterns(states([
            ("ON", "OFF", "ON", "gene1_high"), ("OFF", "ON", "ON", "gene2_high"),
            ("ON", "OFF", "OFF", "gene1_high")]), PolarizationConfig())
        self.assertEqual(label, "UNRESOLVED")

    def test_inadequate_partner_retention(self):
        label, _, _ = polarize_patterns(states([
            ("ON", "OFF", "ON", "gene1_high"), ("OFF", "OFF", "ON", "both_inactive"),
            ("OFF", "ON", "OFF", "gene2_high")]), PolarizationConfig())
        self.assertEqual(label, "UNRESOLVED")

    def test_state_thresholds_and_missing_values(self):
        c = PolarizationConfig()
        for values, expected in [([0., .1], "OFF"), ([1., 2.], "ON"), ([.2, .4], "UNKNOWN"),
                                 ([0., float("nan")], "UNKNOWN"), ([0., .2], "UNKNOWN")]:
            self.assertEqual(_state(values, c), expected)

    def test_invalid_config(self):
        for c in [PolarizationConfig(off_threshold=2), PolarizationConfig(minimum_replicates=2.5),
                  PolarizationConfig(retained_coverage=float("nan"))]:
            with self.assertRaises(ValueError):
                c.validate()

    def test_bh_rejects_invalid_probabilities(self):
        with self.assertRaises(ValueError):
            benjamini_hochberg(pd.Series([.1, -1.]))

    def run_example(self, directory, **changes):
        kwargs = dict(primary_results_dir=directory / "primary",
            outgroup_expression_path=ROOT / "examples/ancestor_expression.tsv",
            outgroup_samples_path=ROOT / "examples/ancestor_samples.tsv",
            orthologs_path=ROOT / "examples/orthologs.tsv", tissue_map_path=ROOT / "examples/tissue_map.tsv",
            detection_qc_path=ROOT / "examples/detection_qc.tsv", output_dir=directory / "polarization",
            outgroup_species="Synthetic example", reference_version="fixture-v1", orthology_source="Synthetic fixture")
        kwargs.update(changes)
        return run_outgroup_polarization(**kwargs)

    def primary(self, directory):
        return run_primary_classification(ROOT / "examples/expression.tsv", ROOT / "examples/pairs.tsv",
            directory / "primary", ROOT / "examples/samples.tsv")

    def test_independent_workflow_and_no_overwrite(self):
        with TemporaryDirectory() as folder:
            directory = Path(folder)
            primary = self.primary(directory)
            before = primary.output_paths["primary_classifications"].read_bytes()
            result = self.run_example(directory)
            labels = result.classifications.set_index("pair_id").polarization_class
            self.assertEqual(labels["pair_sub"], "OUTGROUP_SUPPORTED_SUB_CANDIDATE")
            self.assertEqual(labels["pair_neo"], "OUTGROUP_SUPPORTED_NEO_CANDIDATE")
            self.assertEqual(labels["pair_aed"], "NOT_ELIGIBLE")
            self.assertEqual(before, primary.output_paths["primary_classifications"].read_bytes())
            self.assertNotIn("rule_support_score", result.classifications)
            with self.assertRaises(FileExistsError):
                self.run_example(directory)
            with self.assertRaises(FileExistsError):
                self.primary(directory)
            self.assertEqual(json.loads(result.output_paths["metadata"].read_text(encoding="utf-8"))["version"], "1.0.0")

    def test_missing_qc_cannot_generate_neo(self):
        with TemporaryDirectory() as folder:
            directory = Path(folder)
            self.primary(directory)
            qc = directory / "qc.tsv"
            qc.write_text("outgroup_gene_id\toutgroup_tissue\tassessment\nANEO\troot\tadequate\n")
            result = self.run_example(directory, detection_qc_path=qc)
            labels = result.classifications.set_index("pair_id").polarization_class
            self.assertEqual(labels["pair_neo"], "UNRESOLVED")

    def test_modified_primary_is_rejected(self):
        with TemporaryDirectory() as folder:
            directory = Path(folder)
            primary = self.primary(directory)
            with primary.output_paths["tissue_evidence"].open("a") as handle:
                handle.write("\n")
            with self.assertRaisesRegex(ValueError, "integrity"):
                self.run_example(directory)

    def test_no_effect_only_public_api(self):
        with TemporaryDirectory() as folder:
            with self.assertRaisesRegex(ValueError, "metadata"):
                run_primary_classification(ROOT / "examples/expression.tsv", ROOT / "examples/pairs.tsv", folder)

    def test_copy_exchange_and_row_order(self):
        original = states([("ON", "OFF", "ON", "gene1_high"), ("OFF", "ON", "OFF", "gene2_high")])
        exchanged = original.rename(columns={"copy1_state": "copy2_state", "copy2_state": "copy1_state"})
        exchanged["direction"] = exchanged.direction.map({"gene1_high": "gene2_high", "gene2_high": "gene1_high"})
        a = polarize_patterns(original, PolarizationConfig())
        b = polarize_patterns(exchanged.iloc[::-1], PolarizationConfig())
        self.assertEqual(a[0], b[0])
        self.assertEqual(b[2]["candidate_gain_copy"], "copy1")

    def test_sparse_coverage_not_hidden(self):
        rows = [("ON", "OFF", "ON", "gene1_high"), ("OFF", "ON", "OFF", "gene2_high")]
        rows += [("UNKNOWN", "UNKNOWN", "UNKNOWN", "no_difference")] * 5
        self.assertEqual(polarize_patterns(states(rows), PolarizationConfig())[0], "UNRESOLVED")

    def test_approximate_tissues_require_opt_in(self):
        with TemporaryDirectory() as folder:
            directory = Path(folder)
            self.primary(directory)
            mapping = directory / "map.tsv"
            mapping.write_text((ROOT / "examples/tissue_map.tsv").read_text().replace("exact", "approximate"))
            result = self.run_example(directory, tissue_map_path=mapping)
            self.assertEqual(result.classifications.set_index("pair_id").loc["pair_sub", "polarization_class"], "UNRESOLVED")
            result2 = self.run_example(directory, tissue_map_path=mapping,
                config=PolarizationConfig(allow_approximate_tissues=True), output_dir=directory / "opt_in")
            self.assertEqual(result2.classifications.set_index("pair_id").loc["pair_sub", "polarization_class"], "OUTGROUP_SUPPORTED_SUB_CANDIDATE")

    def test_ambiguous_ortholog_unresolved(self):
        with TemporaryDirectory() as folder:
            directory = Path(folder)
            self.primary(directory)
            orthologs = directory / "orthologs.tsv"
            orthologs.write_text((ROOT / "examples/orthologs.tsv").read_text().replace("single_copy", "ambiguous"))
            result = self.run_example(directory, orthologs_path=orthologs)
            self.assertEqual(result.classifications.set_index("pair_id").loc["pair_sub", "polarization_class"], "UNRESOLVED")

    def test_duplicate_tissue_mapping_rejected(self):
        with TemporaryDirectory() as folder:
            directory = Path(folder)
            self.primary(directory)
            mapping = directory / "map.tsv"
            mapping.write_text("target_tissue\toutgroup_tissue\tmatch_quality\nroot\troot\texact\nleaf\troot\texact\n")
            with self.assertRaisesRegex(ValueError, "one-to-one"):
                self.run_example(directory, tissue_map_path=mapping)

    def test_mismatched_threshold_rejected(self):
        with TemporaryDirectory() as folder:
            directory = Path(folder)
            self.primary(directory)
            with self.assertRaisesRegex(ValueError, "threshold"):
                self.run_example(directory, config=PolarizationConfig(on_threshold=2.))
