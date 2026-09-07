"""Classification rules for duplicate-gene expression fates."""

from __future__ import annotations

from math import ceil

import numpy as np
import pandas as pd

from .config import IDGRMConfig


FATE_LABELS_CN = {
    "SUBFUNCTIONALIZATION": "亚功能化",
    "NEOFUNCTIONALIZATION": "新功能化",
    "AMBIGUOUS_SUB_NEO": "亚/新功能化歧义",
    "AED": "非对称表达",
    "EXPRESSION_LOSS": "表达丢失/非功能化倾向",
    "NO_DIFFERENCE": "无显著差异",
    "INSUFFICIENT_DATA": "数据不足",
}

SCIENCE_LABELS_CN = {
    "SUB_OR_NEO": "亚/新功能化候选",
    "AED": "非对称表达",
    "NO_DIFFERENCE": "无显著差异",
    "INSUFFICIENT_DATA": "数据不足",
}


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


def _blank_ancestor_metrics() -> dict[str, object]:
    return {
        "n_ancestor_tissues": 0,
        "n_novel_gain_gene1": 0,
        "n_novel_gain_gene2": 0,
        "ancestor_coverage_gene1": float("nan"),
        "ancestor_coverage_gene2": float("nan"),
        "ancestor_combined_coverage": float("nan"),
        "ancestor_similarity_gene1": float("nan"),
        "ancestor_similarity_gene2": float("nan"),
    }


def _split_with_ancestor(
    pair: object,
    rows: pd.DataFrame,
    ancestor_profile: pd.Series,
    config: IDGRMConfig,
    dominance_balance: float,
) -> tuple[str, str, str, float, dict[str, object]]:
    profile = ancestor_profile.copy()
    profile.index = profile.index.astype(str)
    shared = rows[rows["tissue"].astype(str).isin(profile.index)].copy()
    shared["ancestor_mean"] = shared["tissue"].astype(str).map(profile)
    shared = shared[pd.to_numeric(shared["ancestor_mean"], errors="coerce").notna()]
    if len(shared) < config.min_evaluable_tissues:
        return "", "", "", float("nan"), _blank_ancestor_metrics()

    ancestor_values = pd.to_numeric(shared["ancestor_mean"], errors="coerce").to_numpy(float)
    ancestor_active = ancestor_values >= config.resolved_ancestor_min_expression
    direction1 = shared["direction"].eq("gene1_high").to_numpy(bool)
    direction2 = shared["direction"].eq("gene2_high").to_numpy(bool)
    active1 = shared["active1"].to_numpy(bool)
    active2 = shared["active2"].to_numpy(bool)
    means1 = pd.to_numeric(shared["mean1"], errors="coerce").to_numpy(float)
    means2 = pd.to_numeric(shared["mean2"], errors="coerce").to_numpy(float)
    fold1 = (means1 + config.pseudocount) / (ancestor_values + config.pseudocount)
    fold2 = (means2 + config.pseudocount) / (ancestor_values + config.pseudocount)

    gain1 = direction1 & active1 & ~ancestor_active & (fold1 >= config.gain_fold)
    gain2 = direction2 & active2 & ~ancestor_active & (fold2 >= config.gain_fold)
    n_gain1 = int(gain1.sum())
    n_gain2 = int(gain2.sum())
    n_ancestor_active = int(ancestor_active.sum())
    coverage1 = _safe_ratio(int((active1 & ancestor_active).sum()), n_ancestor_active, np.nan)
    coverage2 = _safe_ratio(int((active2 & ancestor_active).sum()), n_ancestor_active, np.nan)
    combined_coverage = _safe_ratio(
        int(((active1 | active2) & ancestor_active).sum()), n_ancestor_active, np.nan
    )
    ancestor_dom1 = int((direction1 & ancestor_active).sum())
    ancestor_dom2 = int((direction2 & ancestor_active).sum())
    similarity1 = _safe_spearman(np.log1p(means1), np.log1p(ancestor_values))
    similarity2 = _safe_spearman(np.log1p(means2), np.log1p(ancestor_values))
    metrics = {
        "n_ancestor_tissues": int(len(shared)),
        "n_novel_gain_gene1": n_gain1,
        "n_novel_gain_gene2": n_gain2,
        "ancestor_coverage_gene1": coverage1,
        "ancestor_coverage_gene2": coverage2,
        "ancestor_combined_coverage": combined_coverage,
        "ancestor_similarity_gene1": similarity1,
        "ancestor_similarity_gene2": similarity2,
    }

    if (n_gain1 > 0) ^ (n_gain2 > 0):
        if n_gain1:
            innovator = pair.gene1
            retained = pair.gene2
            retained_coverage = coverage2
            retained_dom = ancestor_dom2
            retained_similarity = similarity2
            innovator_similarity = similarity1
        else:
            innovator = pair.gene2
            retained = pair.gene1
            retained_coverage = coverage1
            retained_dom = ancestor_dom1
            retained_similarity = similarity1
            innovator_similarity = similarity2
        similarity_support = (
            not np.isfinite(retained_similarity)
            or not np.isfinite(innovator_similarity)
            or retained_similarity
            >= innovator_similarity + config.ancestor_similarity_delta
        )
        retention_support = (
            np.isfinite(retained_coverage)
            and retained_coverage >= config.neo_parent_coverage_min
            and retained_dom >= 1
        )
        if retention_support and (
            similarity_support or retained_coverage >= config.ancestor_coverage_min
        ):
            confidence = min(0.97, 0.82 + 0.04 * (n_gain1 + n_gain2))
            reason = (
                f"{innovator} has {n_gain1 + n_gain2} dominant expression gain(s) in "
                f"ancestor-inactive tissue(s), while {retained} retains "
                f"{retained_coverage:.0%} of the ancestral active-tissue breadth."
            )
            return "NEOFUNCTIONALIZATION", innovator, reason, confidence, metrics

    no_novel_gain = n_gain1 == 0 and n_gain2 == 0
    if (
        no_novel_gain
        and ancestor_dom1 >= 1
        and ancestor_dom2 >= 1
        and np.isfinite(combined_coverage)
        and combined_coverage >= config.ancestor_coverage_min
        and dominance_balance >= config.sub_dominance_balance_min
    ):
        confidence = min(0.95, 0.80 + 0.03 * min(ancestor_dom1, ancestor_dom2))
        reason = (
            "Both copies dominate different ancestor-expressed tissues, their union covers "
            f"{combined_coverage:.0%} of the ancestral expression domain, and no robust "
            "ancestor-inactive gain is detected."
        )
        return "SUBFUNCTIONALIZATION", "", reason, confidence, metrics

    reason = (
        "Reciprocal tissue dominance is present, but the outgroup pattern does not uniquely "
        "support complementary ancestral partition or a one-copy novel gain."
    )
    return "AMBIGUOUS_SUB_NEO", "", reason, 0.45, metrics


def _split_expression_only(
    pair: object,
    rows: pd.DataFrame,
    n_gene1_high: int,
    n_gene2_high: int,
    breadth1: float,
    breadth2: float,
    breadth_balance: float,
    dominance_balance: float,
    active_overlap: float,
    expression_correlation: float,
    config: IDGRMConfig,
) -> tuple[str, str, str, float]:
    onoff1 = int(
        (rows["direction"].eq("gene1_high") & rows["active1"] & ~rows["active2"]).sum()
    )
    onoff2 = int(
        (rows["direction"].eq("gene2_high") & rows["active2"] & ~rows["active1"]).sum()
    )
    n_tissues = max(1, len(rows))

    if breadth1 > breadth2:
        broad_gene, narrow_gene = pair.gene1, pair.gene2
        broad_breadth, narrow_breadth = breadth1, breadth2
        broad_dominance, narrow_dominance = n_gene1_high, n_gene2_high
        narrow_onoff = onoff2
    elif breadth2 > breadth1:
        broad_gene, narrow_gene = pair.gene2, pair.gene1
        broad_breadth, narrow_breadth = breadth2, breadth1
        broad_dominance, narrow_dominance = n_gene2_high, n_gene1_high
        narrow_onoff = onoff1
    else:
        broad_gene = narrow_gene = ""
        broad_breadth = narrow_breadth = breadth1
        broad_dominance = narrow_dominance = 0
        narrow_onoff = 0

    neo_proxy = (
        bool(narrow_gene)
        and broad_breadth >= config.neo_retained_breadth_min
        and breadth_balance < config.neo_breadth_ratio_max
        and narrow_onoff >= 1
        and broad_dominance >= 1
        and _safe_ratio(narrow_dominance, n_tissues) <= config.neo_novel_fraction_max
    )
    if neo_proxy:
        reason = (
            f"Expression-only proxy: {broad_gene} is broadly expressed ({broad_breadth:.0%} "
            f"of evaluable tissues), whereas {narrow_gene} has a restricted copy-specific "
            "dominant domain. An outgroup is required to prove that this domain is ancestral-new."
        )
        return "NEOFUNCTIONALIZATION", narrow_gene, reason, 0.65

    complementary_signal = (
        (onoff1 >= 1 and onoff2 >= 1)
        or (np.isfinite(expression_correlation) and expression_correlation <= config.sub_max_correlation)
        or active_overlap <= config.sub_max_active_overlap
    )
    sub_proxy = (
        dominance_balance >= config.sub_dominance_balance_min
        and breadth_balance >= config.sub_breadth_balance_min
        and complementary_signal
    )
    if sub_proxy:
        reason = (
            "Expression-only proxy: both copies dominate different tissues with balanced "
            f"contributions (dominance balance={dominance_balance:.2f}, breadth balance="
            f"{breadth_balance:.2f}). An outgroup is required to verify ancestral partition."
        )
        return "SUBFUNCTIONALIZATION", "", reason, 0.67

    reason = (
        "The Science reciprocal-dominance rule is met, but expression breadth and "
        "complementarity do not distinguish subfunctionalization from neofunctionalization."
    )
    return "AMBIGUOUS_SUB_NEO", "", reason, 0.40


def classify_pairs(
    evidence: pd.DataFrame,
    pairs: pd.DataFrame,
    config: IDGRMConfig,
    ancestor_tissue_expression: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Classify each pair under both the Science and extended iDGRM schemes."""

    config.validate()
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
        active_union = rows.get("active1", pd.Series(dtype=bool)).fillna(False).astype(bool) | rows.get(
            "active2", pd.Series(dtype=bool)
        ).fillna(False).astype(bool)
        evaluable = active_union & pd.to_numeric(
            rows.get("log2fc_gene2_over_gene1", pd.Series(index=rows.index, dtype=float)),
            errors="coerce",
        ).notna()
        active_rows = rows.loc[evaluable].copy()
        methods = active_rows.get("method", pd.Series(index=active_rows.index, dtype=str)).astype(str)
        has_stat = pd.to_numeric(
            active_rows.get("qvalue", pd.Series(index=active_rows.index, dtype=float)),
            errors="coerce",
        ).notna() | pd.to_numeric(
            active_rows.get("pvalue", pd.Series(index=active_rows.index, dtype=float)),
            errors="coerce",
        ).notna()
        callable_mask = methods.eq("effect_only") | has_stat
        analysis_rows = active_rows.loc[callable_mask].copy()
        n_active = int(len(active_rows))
        n_evaluable = int(len(analysis_rows))
        n_callable = n_evaluable
        n_gene1_high = int(analysis_rows["direction"].eq("gene1_high").sum())
        n_gene2_high = int(analysis_rows["direction"].eq("gene2_high").sum())
        n_no_difference = int(analysis_rows["direction"].eq("no_difference").sum())
        aed_min_tissues = max(1, ceil(n_evaluable * config.aed_fraction))
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
        loss1_fraction = _safe_ratio(
            int((~analysis_rows.get("active1", False) & analysis_rows.get("active2", False)).sum()),
            n_evaluable,
        )
        loss2_fraction = _safe_ratio(
            int((~analysis_rows.get("active2", False) & analysis_rows.get("active1", False)).sum()),
            n_evaluable,
        )

        ancestor_metrics = _blank_ancestor_metrics()
        innovating_copy = ""
        lost_copy = ""
        evidence_basis = "science_rule"

        if n_evaluable < config.min_evaluable_tissues or n_callable < config.min_evaluable_tissues:
            science_class = "INSUFFICIENT_DATA"
            fate = "INSUFFICIENT_DATA"
            confidence = 0.0
            reason = (
                f"Only {n_callable} callable tissue(s) were available; at least "
                f"{config.min_evaluable_tissues} are required."
            )
        else:
            if reciprocal:
                science_class = "SUB_OR_NEO"
            elif (
                n_gene1_high >= aed_min_tissues and n_gene2_high == 0
            ) or (
                n_gene2_high >= aed_min_tissues and n_gene1_high == 0
            ):
                science_class = "AED"
            else:
                science_class = "NO_DIFFERENCE"

            ancestor_available = (
                ancestor_tissue_expression is not None
                and bool(str(pair.ancestor_id).strip())
                and str(pair.ancestor_id) in ancestor_tissue_expression.index
            )
            if reciprocal:
                if ancestor_available:
                    outcome = _split_with_ancestor(
                        pair,
                        analysis_rows,
                        ancestor_tissue_expression.loc[str(pair.ancestor_id)],
                        config,
                        dominance_balance,
                    )
                    fate, innovating_copy, reason, confidence, ancestor_metrics = outcome
                    if fate:
                        evidence_basis = "outgroup_supported"
                    else:
                        fate, innovating_copy, reason, confidence = _split_expression_only(
                            pair,
                            analysis_rows,
                            n_gene1_high,
                            n_gene2_high,
                            breadth1,
                            breadth2,
                            breadth_balance,
                            dominance_balance,
                            active_overlap,
                            expression_correlation,
                            config,
                        )
                        evidence_basis = "expression_only_proxy"
                else:
                    fate, innovating_copy, reason, confidence = _split_expression_only(
                        pair,
                        analysis_rows,
                        n_gene1_high,
                        n_gene2_high,
                        breadth1,
                        breadth2,
                        breadth_balance,
                        dominance_balance,
                        active_overlap,
                        expression_correlation,
                        config,
                    )
                    evidence_basis = "expression_only_proxy"
            elif (
                config.detect_expression_loss
                and loss1_fraction >= config.loss_fraction
                and n_gene1_high == 0
            ):
                fate = "EXPRESSION_LOSS"
                lost_copy = pair.gene1
                confidence = min(0.95, 0.65 + 0.25 * loss1_fraction)
                reason = (
                    f"{pair.gene1} is inactive while {pair.gene2} is active in "
                    f"{loss1_fraction:.0%} of evaluable tissues and never significantly dominates."
                )
            elif (
                config.detect_expression_loss
                and loss2_fraction >= config.loss_fraction
                and n_gene2_high == 0
            ):
                fate = "EXPRESSION_LOSS"
                lost_copy = pair.gene2
                confidence = min(0.95, 0.65 + 0.25 * loss2_fraction)
                reason = (
                    f"{pair.gene2} is inactive while {pair.gene1} is active in "
                    f"{loss2_fraction:.0%} of evaluable tissues and never significantly dominates."
                )
            elif science_class == "AED":
                fate = "AED"
                confidence = min(
                    0.95,
                    0.65 + 0.25 * _safe_ratio(max(n_gene1_high, n_gene2_high), n_evaluable),
                )
                reason = (
                    f"One copy dominates {max(n_gene1_high, n_gene2_high)}/{n_evaluable} "
                    "evaluable tissues with no significant reversal."
                )
            else:
                fate = "NO_DIFFERENCE"
                confidence = 0.60
                reason = (
                    "The pair does not meet reciprocal-dominance, AED, or expression-loss criteria."
                )

        record = {
            "pair_id": pair.pair_id,
            "gene1": pair.gene1,
            "gene2": pair.gene2,
            "duplication_type": pair.duplication_type,
            "role1": pair.role1,
            "role2": pair.role2,
            "ancestor_id": pair.ancestor_id,
            "science_class": science_class,
            "science_class_cn": SCIENCE_LABELS_CN[science_class],
            "extended_fate": fate,
            "extended_fate_cn": FATE_LABELS_CN[fate],
            "evidence_basis": evidence_basis,
            "confidence": round(float(confidence), 3),
            "major_copy": major,
            "minor_copy": minor,
            "innovating_copy": innovating_copy,
            "lost_copy": lost_copy,
            "n_evaluable_tissues": n_evaluable,
            "n_callable_tissues": n_callable,
            "n_active_tissues": n_active,
            "aed_min_tissues": aed_min_tissues,
            "n_gene1_high": n_gene1_high,
            "n_gene2_high": n_gene2_high,
            "n_no_difference": n_no_difference,
            "expression_breadth_gene1": round(breadth1, 6),
            "expression_breadth_gene2": round(breadth2, 6),
            "breadth_balance": round(breadth_balance, 6),
            "dominance_balance": round(dominance_balance, 6),
            "active_overlap": round(active_overlap, 6),
            "expression_correlation": expression_correlation,
            "loss_fraction_gene1": round(loss1_fraction, 6),
            "loss_fraction_gene2": round(loss2_fraction, 6),
            "classification_reason": reason,
        }
        record.update(ancestor_metrics)
        records.append(record)

    return pd.DataFrame.from_records(records)
