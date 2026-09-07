"""iDGRM: inference of duplicated-gene retention mechanisms."""

from .config import IDGRMConfig
from .pipeline import AnalysisResult, run_analysis, run_legacy_analysis

__all__ = ["IDGRMConfig", "AnalysisResult", "run_analysis", "run_legacy_analysis"]
__version__ = "0.1.0"
