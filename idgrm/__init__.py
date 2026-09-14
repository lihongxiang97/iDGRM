"""iDGRM: inference of duplicated-gene retention mechanisms."""

from .config import IDGRMConfig
from .pipeline import (
    AnalysisResult,
    run_primary_classification,
    run_primary_classification_from_differential_expression,
    run_subtype_refinement,
)

__all__ = [
    "IDGRMConfig", "AnalysisResult", "run_primary_classification",
    "run_primary_classification_from_differential_expression", "run_subtype_refinement",
]
__version__ = "0.3.0"
