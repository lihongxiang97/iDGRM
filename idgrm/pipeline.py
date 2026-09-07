"""End-to-end iDGRM analysis workflows."""

from __future__ import annotations

import hashlib
import json
import platform
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

from .classifier import FATE_LABELS_CN, SCIENCE_LABELS_CN, classify_pairs
from .config import IDGRMConfig
from .io import (
    aggregate_tissue_expression,
    normalize_expression,
    read_expression,
    read_pairs,
    read_samples,
)
from .legacy import DEFAULT_FILENAME_REGEX, load_legacy_deseq2
from .report import write_html_report
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


def build_summary(classifications: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    total = len(classifications)
    for scheme, label_column, labels in (
        ("extended_fate", "extended_fate", FATE_LABELS_CN),
        ("science_class", "science_class", SCIENCE_LABELS_CN),
    ):
        counts = classifications[label_column].value_counts(dropna=False)
        for class_code, n_pairs in counts.items():
            rows.append(
                {
                    "scheme": scheme,
                    "class_code": class_code,
                    "label_cn": labels.get(str(class_code), str(class_code)),
                    "n_pairs": int(n_pairs),
                    "proportion": float(n_pairs / total) if total else 0.0,
                }
            )
    return pd.DataFrame(rows)


def _write_outputs(
    output_dir: str | Path,
    classifications: pd.DataFrame,
    evidence: pd.DataFrame,
    summary: pd.DataFrame,
    metadata: dict[str, Any],
) -> dict[str, Path]:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "classifications": output_dir / "classifications.tsv",
        "evidence": output_dir / "tissue_evidence.tsv",
        "summary": output_dir / "summary.tsv",
        "metadata": output_dir / "run_metadata.json",
        "report": output_dir / "report.html",
    }
    classifications.to_csv(paths["classifications"], sep="\t", index=False, na_rep="")
    evidence.to_csv(paths["evidence"], sep="\t", index=False, na_rep="")
    summary.to_csv(paths["summary"], sep="\t", index=False, na_rep="")
    paths["metadata"].write_text(
        json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    write_html_report(paths["report"], classifications, summary, evidence, metadata)
    return paths


def _metadata(config: IDGRMConfig, inputs: dict[str, Any], warnings: list[str], **extra: Any) -> dict[str, Any]:
    return {
        "software": "iDGRM",
        "version": "0.1.0",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "python": platform.python_version(),
        "config": config.to_dict(),
        "inputs": inputs,
        "warnings": warnings,
        **extra,
    }


def _load_ancestor(
    expression_path: str | Path | None,
    samples_path: str | Path | None,
    normalization: str,
) -> tuple[pd.DataFrame | None, pd.DataFrame | None]:
    if expression_path is None:
        if samples_path is not None:
            raise ValueError("--ancestor-samples requires --ancestor-expression")
        return None, None
    expression = normalize_expression(read_expression(expression_path), normalization)
    samples = read_samples(samples_path, expression.columns) if samples_path else None
    return expression, aggregate_tissue_expression(expression, samples)


def run_analysis(
    expression_path: str | Path,
    pairs_path: str | Path,
    output_dir: str | Path,
    samples_path: str | Path | None = None,
    ancestor_expression_path: str | Path | None = None,
    ancestor_samples_path: str | Path | None = None,
    config: IDGRMConfig | None = None,
    normalization: str = "none",
) -> AnalysisResult:
    """Run replicate-aware or tissue-level expression classification."""

    config = (config or IDGRMConfig()).validate()
    warnings: list[str] = []
    expression = normalize_expression(read_expression(expression_path), normalization)
    pairs = read_pairs(pairs_path)
    warnings.extend(pairs.attrs.get("warnings", []))
    samples = read_samples(samples_path, expression.columns) if samples_path else None
    if samples is None:
        warnings.append(
            "No sample metadata was supplied; columns were treated as tissue means and "
            "differential calls use the fold-change threshold without a p-value."
        )
    else:
        unused = sorted(set(expression.columns.astype(str)) - set(samples["sample_id"]))
        if unused:
            warnings.append(
                f"Ignored {len(unused)} expression column(s) not listed in sample metadata."
            )

    missing_gene1 = ~pairs["gene1"].isin(expression.index)
    missing_gene2 = ~pairs["gene2"].isin(expression.index)
    n_missing = int((missing_gene1 | missing_gene2).sum())
    if n_missing:
        warnings.append(f"{n_missing} pair(s) contain at least one gene absent from expression data.")

    _, ancestor_tissues = _load_ancestor(
        ancestor_expression_path, ancestor_samples_path, normalization
    )
    if ancestor_tissues is not None:
        requested = pairs["ancestor_id"].replace("", pd.NA).dropna()
        n_unmatched = int((~requested.isin(ancestor_tissues.index)).sum())
        if n_unmatched:
            warnings.append(f"{n_unmatched} non-empty ancestor ID(s) were absent from ancestor data.")

    evidence = compute_tissue_evidence(expression, pairs, samples, config)
    classifications = classify_pairs(evidence, pairs, config, ancestor_tissues)
    summary = build_summary(classifications)
    inputs = {
        "expression": _input_record(expression_path),
        "samples": _input_record(samples_path),
        "pairs": _input_record(pairs_path),
        "ancestor_expression": _input_record(ancestor_expression_path),
        "ancestor_samples": _input_record(ancestor_samples_path),
    }
    metadata = _metadata(
        config,
        inputs,
        warnings,
        analysis_mode="expression",
        normalization=normalization,
        n_pairs=int(len(pairs)),
        n_tissue_evidence_rows=int(len(evidence)),
    )
    paths = _write_outputs(output_dir, classifications, evidence, summary, metadata)
    return AnalysisResult(classifications, evidence, summary, paths, warnings)


def run_legacy_analysis(
    deseq_directory: str | Path,
    pairs_path: str | Path,
    output_dir: str | Path,
    config: IDGRMConfig | None = None,
    pattern: str = "*.DESeq2.csv",
    filename_regex: str = DEFAULT_FILENAME_REGEX,
    ancestor_expression_path: str | Path | None = None,
    ancestor_samples_path: str | Path | None = None,
    normalization: str = "none",
) -> AnalysisResult:
    """Classify original per-tissue DESeq2 outputs without rerunning alignment."""

    config = (config or IDGRMConfig()).validate()
    pairs = read_pairs(pairs_path)
    pair_warnings = list(pairs.attrs.get("warnings", []))
    evidence, warnings = load_legacy_deseq2(
        deseq_directory, pairs, config, pattern=pattern, filename_regex=filename_regex
    )
    warnings = pair_warnings + warnings
    _, ancestor_tissues = _load_ancestor(
        ancestor_expression_path, ancestor_samples_path, normalization
    )
    classifications = classify_pairs(evidence, pairs, config, ancestor_tissues)
    summary = build_summary(classifications)
    inputs = {
        "deseq_directory": str(Path(deseq_directory).resolve()),
        "pairs": _input_record(pairs_path),
        "ancestor_expression": _input_record(ancestor_expression_path),
        "ancestor_samples": _input_record(ancestor_samples_path),
    }
    metadata = _metadata(
        config,
        inputs,
        warnings,
        analysis_mode="legacy_deseq2",
        pattern=pattern,
        filename_regex=filename_regex,
        n_pairs=int(len(pairs)),
        n_tissue_evidence_rows=int(len(evidence)),
    )
    paths = _write_outputs(output_dir, classifications, evidence, summary, metadata)
    return AnalysisResult(classifications, evidence, summary, paths, warnings)
