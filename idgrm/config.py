"""Configuration for the iDGRM inference engine."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from math import isfinite
from typing import Any


@dataclass(frozen=True)
class IDGRMConfig:
    """Thresholds used to construct tissue evidence and classify gene pairs.

    Defaults retain the main Science paper's twofold and 0.001 cutoffs while
    using Benjamini-Hochberg correction across gene pairs within each tissue.
    """

    log2fc_threshold: float = 1.0
    alpha: float = 0.001
    p_adjust: str = "bh"
    min_expression: float = 1.0
    ancestor_min_expression: float | None = None
    min_active_fraction: float = 0.5
    pseudocount: float = 0.1
    min_replicates: int = 2
    min_evaluable_tissues: int = 2
    aed_fraction: float = 1.0 / 3.0
    loss_fraction: float = 0.8
    detect_expression_loss: bool = True

    # Pair-only sub/neo proxy thresholds.
    sub_dominance_balance_min: float = 0.5
    sub_breadth_balance_min: float = 0.5
    sub_max_correlation: float = 0.0
    sub_max_active_overlap: float = 0.5
    neo_breadth_ratio_max: float = 0.5
    neo_retained_breadth_min: float = 2.0 / 3.0
    neo_novel_fraction_max: float = 1.0 / 3.0

    # Outgroup-supported sub/neo thresholds.
    ancestor_coverage_min: float = 0.8
    neo_parent_coverage_min: float = 0.6
    ancestor_similarity_delta: float = 0.10
    gain_fold: float = 4.0

    def validate(self) -> "IDGRMConfig":
        positive = {
            "log2fc_threshold": self.log2fc_threshold,
            "alpha": self.alpha,
            "min_expression": self.min_expression,
            "pseudocount": self.pseudocount,
            "gain_fold": self.gain_fold,
        }
        for name, value in positive.items():
            if not isfinite(value) or value <= 0:
                raise ValueError(f"{name} must be a finite value > 0; got {value!r}")

        if self.ancestor_min_expression is not None:
            if not isfinite(self.ancestor_min_expression) or self.ancestor_min_expression < 0:
                raise ValueError("ancestor_min_expression must be finite and >= 0")

        fractions = {
            "alpha": self.alpha,
            "min_active_fraction": self.min_active_fraction,
            "aed_fraction": self.aed_fraction,
            "loss_fraction": self.loss_fraction,
            "sub_dominance_balance_min": self.sub_dominance_balance_min,
            "sub_breadth_balance_min": self.sub_breadth_balance_min,
            "sub_max_active_overlap": self.sub_max_active_overlap,
            "neo_breadth_ratio_max": self.neo_breadth_ratio_max,
            "neo_retained_breadth_min": self.neo_retained_breadth_min,
            "neo_novel_fraction_max": self.neo_novel_fraction_max,
            "ancestor_coverage_min": self.ancestor_coverage_min,
            "neo_parent_coverage_min": self.neo_parent_coverage_min,
            "ancestor_similarity_delta": self.ancestor_similarity_delta,
        }
        for name, value in fractions.items():
            if not 0 <= value <= 1:
                raise ValueError(f"{name} must be between 0 and 1; got {value!r}")

        if self.p_adjust not in {"bh", "none"}:
            raise ValueError("p_adjust must be 'bh' or 'none'")
        if self.min_replicates < 2:
            raise ValueError("min_replicates must be at least 2")
        if self.min_evaluable_tissues < 1:
            raise ValueError("min_evaluable_tissues must be at least 1")
        return self

    @property
    def resolved_ancestor_min_expression(self) -> float:
        return (
            self.min_expression
            if self.ancestor_min_expression is None
            else self.ancestor_min_expression
        )

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)
