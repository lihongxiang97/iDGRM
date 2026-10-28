"""iDGRM: four-class duplicate-gene expression-pattern classification."""

from .config import IDGRMConfig
from ._version import __version__
from .pipeline import (
    AnalysisResult,
    run_primary_classification,
    run_primary_classification_from_differential_expression,
)

__all__ = [
    "IDGRMConfig", "AnalysisResult", "run_primary_classification",
    "run_primary_classification_from_differential_expression",
]
