"""Four-class duplicate-gene expression-pattern workflows."""

from __future__ import annotations

import hashlib
import json
import platform
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

from .classifier import (
    PRIMARY_LABELS_CN,
    classify_primary_pairs,
)
from .config import IDGRMConfig
from ._version import __version__
from .io import normalize_expression, read_expression, read_pairs, read_samples
from .legacy import (
    DEFAULT_DIFFERENTIAL_EXPRESSION_FILENAME_REGEX,
    load_differential_expression_results,
)
from .statistics import compute_tissue_evidence


@dataclass
class AnalysisResult:
    classifications: pd.DataFrame
    evidence: pd.DataFrame
    summary: pd.DataFrame
    output_paths: dict[str, Path]
    warnings: list[str]


def _sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _input_record(path: str | Path | None) -> dict[str, str] | None:
    if path is None:
        return None
    resolved = Path(path).resolve()
    return {"path": str(resolved), "sha256": _sha256(resolved)}


def _metadata(config: IDGRMConfig, inputs: dict[str, Any], warnings: list[str], **extra: Any) -> dict[str, Any]:
    return {
        "software": "iDGRM",
        "version": __version__,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "python": platform.python_version(),
        "configuration": config.to_dict(),
        "classification_scope": "FOUR_CLASS_EXPRESSION_PATTERN",
        "method_profile": "ADAPTED_EXPRESSION_PATTERN_CLASSIFICATION",
        "reference_doi": "10.1126/science.aad8411",
        "statistical_scope": "Within-tissue significance; not final-category FDR or functional proof",
        "asymmetry_denominator": "Active tissues with finite effect and valid selected P statistic",
        "inputs": inputs,
        "warnings": warnings,
        **extra,
    }


def _summary(frame: pd.DataFrame, class_column: str, labels: dict[str, str], parent: str | None = None) -> pd.DataFrame:
    if parent is None:
        counts = frame[class_column].value_counts()
        return pd.DataFrame([
            {
                "class": code, "class_cn": labels[code], "pair_count": int(count),
                "proportion": float(count / len(frame)) if len(frame) else 0.0,
            }
            for code, count in counts.items()
        ])
    counts = frame.groupby([parent, class_column], dropna=False).size().rename("pair_count").reset_index()
    counts["class_cn"] = counts[class_column].map(labels)
    totals = frame.groupby(parent).size()
    counts["proportion_within_parent"] = counts.apply(
        lambda row: row.pair_count / totals[row[parent]], axis=1
    )
    return counts


def _write_primary(
    output_dir: str | Path,
    classifications: pd.DataFrame,
    evidence: pd.DataFrame,
    metadata: dict[str, Any],
) -> tuple[pd.DataFrame, dict[str, Path]]:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=False)
    summary = _summary(classifications, "primary_class", PRIMARY_LABELS_CN)
    paths = {
        "primary_classifications": output_dir / "primary_classifications.tsv",
        "tissue_evidence": output_dir / "tissue_evidence.tsv",
        "primary_summary": output_dir / "primary_summary.tsv",
        "run_metadata": output_dir / "primary_run_metadata.json",
    }
    if any(path.exists() for path in paths.values()):
        raise FileExistsError("Output files already exist; use a new output directory")
    classifications.to_csv(paths["primary_classifications"], sep="\t", index=False, na_rep="")
    evidence.to_csv(paths["tissue_evidence"], sep="\t", index=False, na_rep="")
    summary.to_csv(paths["primary_summary"], sep="\t", index=False, na_rep="")
    metadata["output_sha256"] = {path.name: _sha256(path) for name, path in paths.items() if name != "run_metadata"}
    paths["run_metadata"].write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8")
    return summary, paths


def run_primary_classification(
    expression_matrix_path: str | Path,
    duplicate_pairs_path: str | Path,
    output_dir: str | Path,
    sample_metadata_path: str | Path | None = None,
    config: IDGRMConfig | None = None,
    normalization: str = "none",
) -> AnalysisResult:
    """Classify all duplicate pairs into four primary expression-pattern classes."""

    config = (config or IDGRMConfig()).validate()
    warnings: list[str] = []
    expression = normalize_expression(read_expression(expression_matrix_path), normalization)
    pairs = read_pairs(duplicate_pairs_path)
    warnings.extend(pairs.attrs.get("warnings", []))
    samples = read_samples(sample_metadata_path, expression.columns) if sample_metadata_path else None
    if samples is None:
        raise ValueError("Biological replicate sample metadata is required; effect-only classification is not a supported public workflow")
    warnings.append("Log-expression paired t testing is a separate supported analysis, not a substitute for count-based DESeq2 or an exact reproduction of the original publication.")
    evidence = compute_tissue_evidence(expression, pairs, samples, config)
    classifications = classify_primary_pairs(evidence, pairs, config)
    metadata = _metadata(
        config,
        {
            "expression_matrix": _input_record(expression_matrix_path),
            "sample_metadata": _input_record(sample_metadata_path),
            "duplicate_pairs": _input_record(duplicate_pairs_path),
        },
        warnings,
        workflow="PRIMARY_EXPRESSION_PATTERN_CLASSIFICATION",
        classes=list(PRIMARY_LABELS_CN),
        normalization=normalization,
    )
    summary, paths = _write_primary(output_dir, classifications, evidence, metadata)
    return AnalysisResult(classifications, evidence, summary, paths, warnings)


def run_primary_classification_from_differential_expression(
    differential_expression_dir: str | Path,
    duplicate_pairs_path: str | Path,
    output_dir: str | Path,
    config: IDGRMConfig | None = None,
    file_pattern: str = "*.DESeq2.csv",
    tissue_name_regex: str = DEFAULT_DIFFERENTIAL_EXPRESSION_FILENAME_REGEX,
    activity_expression_matrix_path: str | Path | None = None,
    activity_sample_metadata_path: str | Path | None = None,
) -> AnalysisResult:
    """Classify from DE tests and measured expression for copy activity."""

    config = (config or IDGRMConfig()).validate()
    pairs = read_pairs(duplicate_pairs_path)
    if activity_expression_matrix_path is None or activity_sample_metadata_path is None:
        raise ValueError("Differential-expression input requires measured activity expression and sample metadata")
    activity_expression = read_expression(activity_expression_matrix_path)
    activity_samples = read_samples(activity_sample_metadata_path, activity_expression.columns)
    evidence, warnings = load_differential_expression_results(
        differential_expression_dir, pairs, config, pattern=file_pattern,
        filename_regex=tissue_name_regex,
        activity_expression=activity_expression,
        activity_samples=activity_samples,
    )
    warnings = list(pairs.attrs.get("warnings", [])) + warnings
    classifications = classify_primary_pairs(evidence, pairs, config)
    metadata = _metadata(
        config,
        {
            "differential_expression_directory": str(Path(differential_expression_dir).resolve()),
            "differential_expression_files": [
                _input_record(path)
                for path in sorted(Path(differential_expression_dir).glob(file_pattern))
                if path.is_file()
            ],
            "activity_expression_matrix": _input_record(activity_expression_matrix_path),
            "activity_sample_metadata": _input_record(activity_sample_metadata_path),
            "duplicate_pairs": _input_record(duplicate_pairs_path),
        },
        warnings,
        workflow="PRIMARY_CLASSIFICATION_FROM_DIFFERENTIAL_EXPRESSION",
        classes=list(PRIMARY_LABELS_CN),
        file_pattern=file_pattern,
        tissue_name_regex=tissue_name_regex,
    )
    summary, paths = _write_primary(output_dir, classifications, evidence, metadata)
    return AnalysisResult(classifications, evidence, summary, paths, warnings)
