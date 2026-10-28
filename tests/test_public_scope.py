"""Public interface is restricted to four expression-pattern classes."""
import contextlib
import io
import unittest

import idgrm
from idgrm.cli import build_parser, _configuration
from idgrm.classifier import PRIMARY_LABELS_CN


class PublicScopeTests(unittest.TestCase):
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
