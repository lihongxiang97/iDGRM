"""Classification rules for duplicate-gene expression fates."""

from __future__ import annotations

from math import ceil

import numpy as np
import pandas as pd

from .config import IDGRMConfig


PRIMARY_LABELS_CN = {
    "UNMAPPED": "未映射或不可定量",
    "NO_DIFFERENCE": "无显著差异",
    "ASYMMETRICALLY_EXPRESSED": "非对称表达",
    "SUB_OR_NEOFUNCTIONALIZED": "亚功能化或新功能化",
}


def _eligible_rows(rows: pd.DataFrame, config: IDGRMConfig) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Select active tissues with a valid effect and selected P statistic."""
    active = rows["active1"].fillna(False).astype(bool) | rows["active2"].fillna(False).astype(bool)
    effect = pd.to_numeric(rows.get("log2fc_gene2_over_gene1", pd.Series(index=rows.index, dtype=float)),
                           errors="coerce")
    active_rows = rows.loc[active & np.isfinite(effect)].copy()
    column = "qvalue" if config.p_adjust == "bh" else "pvalue"
    probabilities = pd.to_numeric(active_rows[column], errors="coerce")
    valid_statistic = np.isfinite(probabilities) & probabilities.between(0, 1)
    effect_only = active_rows["method"].astype(str).eq("effect_only")
    return active_rows, active_rows.loc[effect_only | valid_statistic].copy()


def _safe_ratio(numerator: float, denominator: float, default: float = 0.0) -> float:
    return float(numerator / denominator) if denominator else float(default)


def _safe_spearman(values1: np.ndarray, values2: np.ndarray) -> float:
    mask = np.isfinite(values1) & np.isfinite(values2)
    values1 = values1[mask]
    values2 = values2[mask]
    if len(values1) < 3:
        return float("nan")
    ranks1 = pd.Series(values1).rank(method="average").to_numpy(dtype=float)
    ranks2 = pd.Series(values2).rank(method="average").to_numpy(dtype=float)
    if np.std(ranks1) == 0 or np.std(ranks2) == 0:
        return float("nan")
    return float(np.corrcoef(ranks1, ranks2)[0, 1])


def _major_minor(
    gene1: str,
    gene2: str,
    n_gene1_high: int,
    n_gene2_high: int,
    rows: pd.DataFrame,
) -> tuple[str, str]:
    if n_gene1_high > n_gene2_high:
        return gene1, gene2
    if n_gene2_high > n_gene1_high:
        return gene2, gene1
    mean1 = float(rows["mean1"].mean(skipna=True))
    mean2 = float(rows["mean2"].mean(skipna=True))
    if np.isfinite(mean1) and np.isfinite(mean2):
        if mean1 > mean2:
            return gene1, gene2
        if mean2 > mean1:
            return gene2, gene1
    return "undetermined", "undetermined"


def classify_primary_pairs(
    evidence: pd.DataFrame,
    pairs: pd.DataFrame,
    config: IDGRMConfig,
) -> pd.DataFrame:
    """Classify duplicate pairs into the four primary expression-pattern classes."""

    config.validate()
    if evidence.duplicated(["pair_id", "tissue"]).any():
        raise ValueError("Duplicate pair/tissue evidence is not permitted")
    grouped = {pair_id: group.copy() for pair_id, group in evidence.groupby("pair_id")}
    records: list[dict[str, object]] = []

    for pair in pairs.itertuples(index=False):
        rows = grouped.get(pair.pair_id, pd.DataFrame()).copy()
        if rows.empty:
            rows = pd.DataFrame(
                columns=[
                    "tissue", "mean1", "mean2", "active1", "active2", "direction",
                    "pvalue", "qvalue", "method",
                ]
            )
        active_rows, analysis_rows = _eligible_rows(rows, config)
        n_active = int(len(active_rows))
        n_evaluable = int(len(analysis_rows))
        n_callable = n_evaluable
        n_gene1_high = int(analysis_rows["direction"].eq("gene1_high").sum())
        n_gene2_high = int(analysis_rows["direction"].eq("gene2_high").sum())
        n_no_difference = int(analysis_rows["direction"].eq("no_difference").sum())
        asymmetry_min_tissues = max(1, ceil(n_evaluable * config.asymmetry_tissue_fraction))
        reciprocal = n_gene1_high >= 1 and n_gene2_high >= 1

        breadth1 = _safe_ratio(int(analysis_rows.get("active1", False).sum()), n_evaluable)
        breadth2 = _safe_ratio(int(analysis_rows.get("active2", False).sum()), n_evaluable)
        max_breadth = max(breadth1, breadth2)
        breadth_balance = _safe_ratio(min(breadth1, breadth2), max_breadth)
        dominance_balance = _safe_ratio(
            min(n_gene1_high, n_gene2_high), max(n_gene1_high, n_gene2_high)
        )
        active_overlap = _safe_ratio(
            int((analysis_rows.get("active1", False) & analysis_rows.get("active2", False)).sum()),
            n_evaluable,
        )
        expression_correlation = _safe_spearman(
            np.log1p(pd.to_numeric(analysis_rows.get("mean1"), errors="coerce").to_numpy(float)),
            np.log1p(pd.to_numeric(analysis_rows.get("mean2"), errors="coerce").to_numpy(float)),
        )
        major, minor = _major_minor(
            pair.gene1, pair.gene2, n_gene1_high, n_gene2_high, analysis_rows
        )
        silencing1_fraction = _safe_ratio(
            int((~analysis_rows.get("active1", False) & analysis_rows.get("active2", False)).sum()),
            n_evaluable,
        )
        silencing2_fraction = _safe_ratio(
            int((~analysis_rows.get("active2", False) & analysis_rows.get("active1", False)).sum()),
            n_evaluable,
        )

        if n_evaluable < config.min_evaluable_tissues or n_callable < config.min_evaluable_tissues:
            primary_class = "UNMAPPED"
            primary_reason = (
                f"Only {n_callable} callable tissue(s) were available; at least "
                f"{config.min_evaluable_tissues} are required."
            )
        else:
            if reciprocal:
                primary_class = "SUB_OR_NEOFUNCTIONALIZED"
                primary_reason = (
                    "Both copies are significantly more highly expressed than their partner "
                    "in at least one tissue."
                )
            elif (
                n_gene1_high >= asymmetry_min_tissues and n_gene2_high == 0
            ) or (
                n_gene2_high >= asymmetry_min_tissues and n_gene1_high == 0
            ):
                primary_class = "ASYMMETRICALLY_EXPRESSED"
                primary_reason = (
                    f"One copy dominates {max(n_gene1_high, n_gene2_high)}/{n_evaluable} "
                    "evaluable tissues with no significant reversal."
                )
            else:
                primary_class = "NO_DIFFERENCE"
                primary_reason = "The pair does not meet either directional classification rule."

        record = {
            "pair_id": pair.pair_id,
            "gene1": pair.gene1,
            "gene2": pair.gene2,
            "duplication_type": pair.duplication_type,
            "role1": pair.role1,
            "role2": pair.role2,
            "primary_class": primary_class,
            "primary_class_cn": PRIMARY_LABELS_CN[primary_class],
            "primary_class_display": {
                "UNMAPPED": "Unmappable", "NO_DIFFERENCE": "No difference",
                "ASYMMETRICALLY_EXPRESSED": "Asymmetrically expressed",
                "SUB_OR_NEOFUNCTIONALIZED": "Sub-/neofunctionalized",
            }[primary_class],
            "assessment_status": "INSUFFICIENT_EXPRESSION_EVIDENCE" if primary_class == "UNMAPPED" else "ASSESSED",
            "mapping_status": "NOT_INFERRED_FROM_EXPRESSION",
            "active_tissue_test_coverage": _safe_ratio(n_evaluable, n_active),
            "major_copy": major,
            "minor_copy": minor,
            "n_evaluable_tissues": n_evaluable,
            "n_callable_tissues": n_callable,
            "n_active_tissues": n_active,
            "asymmetry_min_tissues": asymmetry_min_tissues,
            "n_gene1_high": n_gene1_high,
            "n_gene2_high": n_gene2_high,
            "n_no_difference": n_no_difference,
            "expression_breadth_gene1": round(breadth1, 6),
            "expression_breadth_gene2": round(breadth2, 6),
            "breadth_balance": round(breadth_balance, 6),
            "dominance_balance": round(dominance_balance, 6),
            "active_overlap": round(active_overlap, 6),
            "expression_correlation": expression_correlation,
            "silencing_fraction_gene1": round(silencing1_fraction, 6),
            "silencing_fraction_gene2": round(silencing2_fraction, 6),
            "primary_classification_reason": primary_reason,
        }
        records.append(record)

    return pd.DataFrame.from_records(records)
