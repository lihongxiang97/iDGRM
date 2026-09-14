"""Validated two-stage workflows for duplicate-gene expression classification."""

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
    REFINEMENT_LABELS_CN,
    classify_primary_pairs,
    refine_primary_candidates,
)
from .config import IDGRMConfig
from .io import aggregate_tissue_expression, normalize_expression, read_expression, read_pairs, read_samples
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
        "version": "0.3.0",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "python": platform.python_version(),
        "configuration": config.to_dict(),
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
    output_dir.mkdir(parents=True, exist_ok=True)
    summary = _summary(classifications, "primary_class", PRIMARY_LABELS_CN)
    paths = {
        "primary_classifications": output_dir / "primary_classifications.tsv",
        "tissue_evidence": output_dir / "tissue_evidence.tsv",
        "primary_summary": output_dir / "primary_summary.tsv",
        "run_metadata": output_dir / "primary_run_metadata.json",
    }
    classifications.to_csv(paths["primary_classifications"], sep="\t", index=False, na_rep="")
    evidence.to_csv(paths["tissue_evidence"], sep="\t", index=False, na_rep="")
    summary.to_csv(paths["primary_summary"], sep="\t", index=False, na_rep="")
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
        warnings.append("No sample metadata was supplied; columns are treated as tissue-level values.")
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
) -> AnalysisResult:
    """Classify from precomputed tissue-specific differential-expression results."""

    config = (config or IDGRMConfig()).validate()
    pairs = read_pairs(duplicate_pairs_path)
    evidence, warnings = load_differential_expression_results(
        differential_expression_dir, pairs, config, pattern=file_pattern,
        filename_regex=tissue_name_regex,
    )
    warnings = list(pairs.attrs.get("warnings", [])) + warnings
    classifications = classify_primary_pairs(evidence, pairs, config)
    metadata = _metadata(
        config,
        {
            "differential_expression_directory": str(Path(differential_expression_dir).resolve()),
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


def _read_table(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    return pd.read_csv(path, sep="," if path.suffix.lower() == ".csv" else "\t")


def _load_outgroup(
    expression_matrix_path: str | Path | None,
    sample_metadata_path: str | Path | None,
    normalization: str,
) -> pd.DataFrame | None:
    if expression_matrix_path is None:
        if sample_metadata_path is not None:
            raise ValueError("Outgroup sample metadata requires an outgroup expression matrix")
        return None
    expression = normalize_expression(read_expression(expression_matrix_path), normalization)
    samples = read_samples(sample_metadata_path, expression.columns) if sample_metadata_path else None
    return aggregate_tissue_expression(expression, samples)


def run_subtype_refinement(
    primary_classifications_path: str | Path,
    tissue_evidence_path: str | Path,
    duplicate_pairs_path: str | Path,
    output_dir: str | Path,
    config: IDGRMConfig | None = None,
    outgroup_expression_matrix_path: str | Path | None = None,
    outgroup_sample_metadata_path: str | Path | None = None,
    normalization: str = "none",
) -> AnalysisResult:
    """Refine eligible primary classes while retaining their primary assignments."""

    config = (config or IDGRMConfig()).validate()
    primary = _read_table(primary_classifications_path)
    evidence = _read_table(tissue_evidence_path)
    pairs = read_pairs(duplicate_pairs_path)
    outgroup = _load_outgroup(
        outgroup_expression_matrix_path, outgroup_sample_metadata_path, normalization
    )
    refined = refine_primary_candidates(primary, evidence, pairs, config, outgroup)
    summary = _summary(
        refined, "refined_class", REFINEMENT_LABELS_CN, parent="parent_primary_class"
    )
    warnings = list(pairs.attrs.get("warnings", []))
    if outgroup is None:
        warnings.append(
            "No outgroup expression was supplied; neofunctionalization is inferred from "
            "target-species expression patterns and is not evolutionarily polarized."
        )
    metadata = _metadata(
        config,
        {
            "primary_classifications": _input_record(primary_classifications_path),
            "tissue_evidence": _input_record(tissue_evidence_path),
            "duplicate_pairs": _input_record(duplicate_pairs_path),
            "outgroup_expression_matrix": _input_record(outgroup_expression_matrix_path),
            "outgroup_sample_metadata": _input_record(outgroup_sample_metadata_path),
        },
        warnings,
        workflow="HIERARCHICAL_SUBTYPE_REFINEMENT",
        eligible_primary_classes=["ASYMMETRICALLY_EXPRESSED", "SUB_OR_NEOFUNCTIONALIZED"],
    )
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "refined_classifications": output_dir / "refined_classifications.tsv",
        "refined_summary": output_dir / "refined_summary.tsv",
        "run_metadata": output_dir / "refinement_run_metadata.json",
    }
    refined.to_csv(paths["refined_classifications"], sep="\t", index=False, na_rep="")
    summary.to_csv(paths["refined_summary"], sep="\t", index=False, na_rep="")
    paths["run_metadata"].write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8")
    return AnalysisResult(refined, evidence, summary, paths, warnings)
