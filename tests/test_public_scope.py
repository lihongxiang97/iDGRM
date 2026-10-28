"""Public interface is restricted to four expression-pattern classes."""
import contextlib
import io
import importlib.util
from dataclasses import fields
import unittest

import idgrm
from idgrm.cli import build_parser, _configuration
from idgrm.classifier import PRIMARY_LABELS_CN


class PublicScopeTests(unittest.TestCase):
    def test_only_four_class_workflows_are_public(self):
        self.assertEqual(set(idgrm.__all__), {
            "IDGRMConfig", "AnalysisResult", "run_primary_classification",
            "run_primary_classification_from_differential_expression",
        })
        self.assertIsNone(importlib.util.find_spec("idgrm.polarization"))
        self.assertFalse(hasattr(idgrm, "PolarizationConfig"))
        self.assertFalse(hasattr(idgrm, "run_outgroup_polarization"))

    def test_removed_outgroup_command_is_rejected(self):
        with contextlib.redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit) as raised:
                build_parser().parse_args(["polarize-outgroup"])
        self.assertEqual(raised.exception.code, 2)

    def test_config_contains_only_primary_thresholds(self):
        self.assertEqual({field.name for field in fields(idgrm.IDGRMConfig)}, {
            "log2fc_threshold", "alpha", "p_adjust", "min_expression",
            "min_active_fraction", "pseudocount", "min_replicates",
            "min_evaluable_tissues", "asymmetry_tissue_fraction",
        })

    def test_refinement_not_publicly_exported(self):
        self.assertNotIn("run_subtype_refinement", idgrm.__all__)
        self.assertFalse(hasattr(idgrm, "run_subtype_refinement"))

    def test_refinement_command_rejected(self):
        with contextlib.redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit) as raised:
                build_parser().parse_args(["refine-subtypes"])
        self.assertEqual(raised.exception.code, 2)

    def test_primary_options_do_not_require_removed_parameters(self):
        args = build_parser().parse_args([
            "classify-primary", "--expression-matrix", "expression.tsv",
            "--sample-metadata", "samples.tsv", "--duplicate-pairs", "pairs.tsv",
            "--output-dir", "results",
        ])
        config = _configuration(args)
        self.assertEqual(config.alpha, .05)
        self.assertFalse(hasattr(args, "outgroup_gain_fold_change"))

    def test_exactly_four_primary_classes(self):
        self.assertEqual(set(PRIMARY_LABELS_CN), {
            "UNMAPPED", "NO_DIFFERENCE", "ASYMMETRICALLY_EXPRESSED",
            "SUB_OR_NEOFUNCTIONALIZED",
        })
