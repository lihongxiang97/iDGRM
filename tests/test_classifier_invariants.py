"""Independent rule reference and adversarial eligibility checks."""
import itertools
import math
import unittest

import numpy as np
import pandas as pd

from idgrm.classifier import classify_primary_pairs
from idgrm.config import IDGRMConfig


def pair_table():
    return pd.DataFrame([dict(pair_id="p", gene1="g1", gene2="g2", duplication_type="WGD",
                             role1="", role2="", ancestor_id="")])


def evidence(directions):
    return pd.DataFrame([dict(pair_id="p", tissue=f"t{i}", mean1=10., mean2=10.,
        active1=True, active2=True, method="DESeq2", pvalue=1e-5, qvalue=.01,
        log2fc_gene2_over_gene1={"gene1_high": -2., "gene2_high": 2., "no_difference": 0.}[d],
        direction=d) for i, d in enumerate(directions)])


def independent_reference(signs, fraction=1/3):
    """A plain-Python specification; no calls to iDGRM internals."""
    if len(signs) < 2:
        return "UNMAPPED"
    negative, positive = signs.count("gene1_high"), signs.count("gene2_high")
    if negative and positive:
        return "SUB_OR_NEOFUNCTIONALIZED"
    minimum = math.ceil(len(signs) * fraction)
    if max(negative, positive) >= minimum and min(negative, positive) == 0:
        return "ASYMMETRICALLY_EXPRESSED"
    return "NO_DIFFERENCE"


class ClassifierInvariantTests(unittest.TestCase):
    def setUp(self):
        self.config = IDGRMConfig()
        self.pairs = pair_table()

    def classify(self, frame):
        return classify_primary_pairs(frame, self.pairs, self.config).iloc[0]

    def test_reference_all_direction_combinations_through_four_tissues(self):
        for n in range(1, 5):
            for directions in itertools.product(("gene1_high", "gene2_high", "no_difference"), repeat=n):
                with self.subTest(directions=directions):
                    self.assertEqual(self.classify(evidence(directions)).primary_class,
                                     independent_reference(list(directions)))

    def test_copy_exchange_preserves_primary_and_exchanges_major(self):
        original = evidence(["gene1_high", "gene1_high", "no_difference"])
        swapped = original.copy()
        swapped["log2fc_gene2_over_gene1"] *= -1
        swapped["direction"] = swapped.direction.replace({"gene1_high": "gene2_high", "gene2_high": "gene1_high"})
        a, b = self.classify(original), self.classify(swapped)
        self.assertEqual(a.primary_class, b.primary_class)
        self.assertEqual((a.major_copy, b.major_copy), ("g1", "g2"))

    def test_evidence_order_invariance(self):
        frame = evidence(["gene1_high", "gene2_high", "no_difference"])
        pd.testing.assert_series_equal(self.classify(frame), self.classify(frame.iloc[::-1]))

    def test_missing_adjusted_p_excluded_even_if_raw_p_present(self):
        frame = evidence(["gene1_high", "gene2_high", "no_difference"])
        frame.loc[1, "qvalue"] = np.nan
        row = self.classify(frame)
        self.assertEqual(row.n_callable_tissues, 2)
        self.assertEqual(row.primary_class, "ASYMMETRICALLY_EXPRESSED")

    def test_nonfinite_effect_not_callable(self):
        frame = evidence(["gene1_high", "gene2_high"])
        frame.loc[1, "log2fc_gene2_over_gene1"] = np.inf
        self.assertEqual(self.classify(frame).primary_class, "UNMAPPED")

    def test_invalid_probability_not_callable(self):
        frame = evidence(["gene1_high", "gene2_high"])
        frame.loc[1, "qvalue"] = -0.1
        self.assertEqual(self.classify(frame).primary_class, "UNMAPPED")

    def test_local_difference_can_land_in_residual_class(self):
        row = self.classify(evidence(["gene1_high"] + ["no_difference"]*6))
        self.assertEqual(row.primary_class, "NO_DIFFERENCE")
        self.assertEqual(row.n_gene1_high, 1)

    def test_all_inactive_is_unmapped_not_conserved(self):
        frame = evidence(["no_difference", "no_difference"])
        frame[["active1", "active2"]] = False
        self.assertEqual(self.classify(frame).primary_class, "UNMAPPED")


if __name__ == "__main__":
    unittest.main()
