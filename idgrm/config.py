"""Configuration for the iDGRM inference engine."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from math import isfinite
from typing import Any


@dataclass(frozen=True)
class IDGRMConfig:
    """Thresholds used to construct tissue evidence and classify gene pairs.

    Defaults use a twofold effect threshold and a 0.05 false-discovery-rate threshold,
    with Benjamini-Hochberg correction across gene pairs within each tissue.
    """

    log2fc_threshold: float = 1.0
    alpha: float = 0.05
    p_adjust: str = "bh"
    min_expression: float = 1.0
    min_active_fraction: float = 0.5
    pseudocount: float = 0.1
    min_replicates: int = 2
    min_evaluable_tissues: int = 2
    asymmetry_tissue_fraction: float = 1.0 / 3.0

    def validate(self) -> "IDGRMConfig":
        positive = {
            "log2fc_threshold": self.log2fc_threshold,
            "alpha": self.alpha,
            "min_expression": self.min_expression,
            "pseudocount": self.pseudocount,
        }
        for name, value in positive.items():
            if not isfinite(value) or value <= 0:
                raise ValueError(f"{name} must be a finite value > 0; got {value!r}")

        fractions = {
            "alpha": self.alpha,
            "min_active_fraction": self.min_active_fraction,
            "asymmetry_tissue_fraction": self.asymmetry_tissue_fraction,
        }
        for name, value in fractions.items():
            if not 0 <= value <= 1:
                raise ValueError(f"{name} must be between 0 and 1; got {value!r}")

        if self.p_adjust not in {"bh", "none"}:
            raise ValueError("p_adjust must be 'bh' or 'none'")
        if isinstance(self.min_replicates, bool) or not isinstance(self.min_replicates, int) or self.min_replicates < 2:
            raise ValueError("min_replicates must be at least 2")
        if isinstance(self.min_evaluable_tissues, bool) or not isinstance(self.min_evaluable_tissues, int) or self.min_evaluable_tissues < 1:
            raise ValueError("min_evaluable_tissues must be at least 1")
        if self.asymmetry_tissue_fraction == 0 or self.min_active_fraction == 0:
            raise ValueError("Activity and asymmetry fractions must be > 0")
        return self

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)
