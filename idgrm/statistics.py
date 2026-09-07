"""Replicate-aware tissue evidence construction."""

from __future__ import annotations

from math import isclose

import numpy as np
import pandas as pd

from .config import IDGRMConfig


def benjamini_hochberg(pvalues: pd.Series) -> pd.Series:
    """Benjamini-Hochberg adjusted p-values, preserving missing entries."""

    values = pd.to_numeric(pvalues, errors="coerce").to_numpy(dtype=float)
    result = np.full(values.shape, np.nan, dtype=float)
    valid = np.isfinite(values)
    if not valid.any():
        return pd.Series(result, index=pvalues.index, dtype=float)
    valid_values = values[valid]
    order = np.argsort(valid_values)
    ranked = valid_values[order]
    m = len(ranked)
    adjusted = ranked * m / np.arange(1, m + 1, dtype=float)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0.0, 1.0)
    restored = np.empty_like(adjusted)
    restored[order] = adjusted
    result[np.flatnonzero(valid)] = restored
    return pd.Series(result, index=pvalues.index, dtype=float)


def _paired_t_pvalue(values1: np.ndarray, values2: np.ndarray, pseudocount: float) -> float:
    mask = np.isfinite(values1) & np.isfinite(values2)
    values1 = values1[mask]
    values2 = values2[mask]
    if len(values1) < 2:
        return float("nan")
    differences = np.log2(values2 + pseudocount) - np.log2(values1 + pseudocount)
    if np.allclose(differences, differences[0], rtol=1e-12, atol=1e-12):
        return 1.0 if isclose(float(differences[0]), 0.0, abs_tol=1e-12) else 0.0
    try:
        from scipy.stats import ttest_rel
    except ImportError as exc:  # pragma: no cover - dependency error is explicit
        raise RuntimeError(
            "Replicate-aware analysis requires scipy. Install iDGRM with 'pip install -e .'"
        ) from exc
    result = ttest_rel(
        np.log2(values2 + pseudocount),
        np.log2(values1 + pseudocount),
        nan_policy="omit",
    )
    return float(result.pvalue) if np.isfinite(result.pvalue) else float("nan")


def compute_tissue_evidence(
    expression: pd.DataFrame,
    pairs: pd.DataFrame,
    samples: pd.DataFrame | None,
    config: IDGRMConfig,
) -> pd.DataFrame:
    """Build one auditable evidence row per duplicate pair and tissue."""

    config.validate()
    if samples is None:
        tissue_to_samples = {str(column): [str(column)] for column in expression.columns}
        method = "effect_only"
    else:
        tissue_to_samples = {
            str(tissue): subset["sample_id"].tolist()
            for tissue, subset in samples.groupby("tissue", sort=False)
        }
        method = "paired_t"

    records: list[dict[str, object]] = []
    for pair in pairs.itertuples(index=False):
        gene1_present = pair.gene1 in expression.index
        gene2_present = pair.gene2 in expression.index
        for tissue, sample_ids in tissue_to_samples.items():
            if gene1_present:
                values1 = expression.loc[pair.gene1, sample_ids].to_numpy(dtype=float)
            else:
                values1 = np.full(len(sample_ids), np.nan)
            if gene2_present:
                values2 = expression.loc[pair.gene2, sample_ids].to_numpy(dtype=float)
            else:
                values2 = np.full(len(sample_ids), np.nan)

            mean1 = float(np.nanmean(values1)) if np.isfinite(values1).any() else float("nan")
            mean2 = float(np.nanmean(values2)) if np.isfinite(values2).any() else float("nan")
            median1 = float(np.nanmedian(values1)) if np.isfinite(values1).any() else float("nan")
            median2 = float(np.nanmedian(values2)) if np.isfinite(values2).any() else float("nan")
            active_fraction1 = float(np.mean(values1 >= config.min_expression))
            active_fraction2 = float(np.mean(values2 >= config.min_expression))
            active1 = bool(
                np.isfinite(median1)
                and median1 >= config.min_expression
                and active_fraction1 >= config.min_active_fraction
            )
            active2 = bool(
                np.isfinite(median2)
                and median2 >= config.min_expression
                and active_fraction2 >= config.min_active_fraction
            )
            if np.isfinite(mean1) and np.isfinite(mean2):
                log2fc = float(
                    np.log2((mean2 + config.pseudocount) / (mean1 + config.pseudocount))
                )
            else:
                log2fc = float("nan")

            n_complete = int(np.sum(np.isfinite(values1) & np.isfinite(values2)))
            if not gene1_present or not gene2_present:
                row_method = "missing_gene"
                pvalue = float("nan")
            elif method == "effect_only":
                row_method = method
                pvalue = float("nan")
            elif n_complete < config.min_replicates:
                row_method = "insufficient_replicates"
                pvalue = float("nan")
            elif not (active1 or active2):
                row_method = method
                pvalue = 1.0
            else:
                row_method = method
                pvalue = _paired_t_pvalue(values1, values2, config.pseudocount)

            records.append(
                {
                    "pair_id": pair.pair_id,
                    "gene1": pair.gene1,
                    "gene2": pair.gene2,
                    "tissue": tissue,
                    "mean1": mean1,
                    "mean2": mean2,
                    "median1": median1,
                    "median2": median2,
                    "active1": active1,
                    "active2": active2,
                    "active_fraction1": active_fraction1,
                    "active_fraction2": active_fraction2,
                    "n_replicates": n_complete,
                    "log2fc_gene2_over_gene1": log2fc,
                    "pvalue": pvalue,
                    "method": row_method,
                }
            )

    evidence = pd.DataFrame.from_records(records)
    evidence["qvalue"] = np.nan
    for _, indices in evidence.groupby("tissue", sort=False).groups.items():
        evidence.loc[indices, "qvalue"] = benjamini_hochberg(
            evidence.loc[indices, "pvalue"]
        )

    active_union = evidence["active1"] | evidence["active2"]
    if method == "effect_only":
        evidence["significant"] = (
            active_union
            & evidence["log2fc_gene2_over_gene1"].abs().ge(config.log2fc_threshold)
        )
    else:
        test_values = evidence["qvalue"] if config.p_adjust == "bh" else evidence["pvalue"]
        evidence["significant"] = (
            active_union
            & test_values.le(config.alpha)
            & evidence["log2fc_gene2_over_gene1"].abs().ge(config.log2fc_threshold)
        )
    evidence["direction"] = "no_difference"
    evidence.loc[~active_union, "direction"] = "both_inactive"
    evidence.loc[
        evidence["significant"]
        & evidence["log2fc_gene2_over_gene1"].le(-config.log2fc_threshold),
        "direction",
    ] = "gene1_high"
    evidence.loc[
        evidence["significant"]
        & evidence["log2fc_gene2_over_gene1"].ge(config.log2fc_threshold),
        "direction",
    ] = "gene2_high"
    return evidence
