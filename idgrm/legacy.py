"""Loader for tissue-specific differential-expression result tables."""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

from .config import IDGRMConfig
from .statistics import benjamini_hochberg


DEFAULT_DIFFERENTIAL_EXPRESSION_FILENAME_REGEX = r"^[^.]+\.[^.]+\.pairs\.(?P<tissue>.+)\.DESeq2\.csv$"


def _column_by_name(frame: pd.DataFrame, candidates: set[str]) -> str | None:
    mapping = {str(column).strip().lower(): str(column) for column in frame.columns}
    for candidate in candidates:
        if candidate.lower() in mapping:
            return mapping[candidate.lower()]
    return None


def load_differential_expression_results(
    directory: str | Path,
    pairs: pd.DataFrame,
    config: IDGRMConfig,
    pattern: str = "*.DESeq2.csv",
    filename_regex: str = DEFAULT_DIFFERENTIAL_EXPRESSION_FILENAME_REGEX,
    activity_expression: pd.DataFrame | None = None,
    activity_samples: pd.DataFrame | None = None,
) -> tuple[pd.DataFrame, list[str]]:
    """Convert tissue-specific DESeq2 tables to a normalized evidence table."""

    directory = Path(directory)
    files = sorted(directory.glob(pattern))
    if not files:
        raise FileNotFoundError(f"No differential-expression files matched {directory / pattern}")
    matcher = re.compile(filename_regex)
    config.validate()

    if activity_expression is None or activity_samples is None:
        raise ValueError(
            "Precomputed differential-expression results require a measured expression matrix "
            "and sample metadata for copy-activity calls; DESeq2 baseMean is not TPM"
        )
    tissue_activity: dict[str, pd.DataFrame] = {}
    genes = set(pairs["gene1"]) | set(pairs["gene2"])
    missing_genes = genes - set(activity_expression.index)
    if missing_genes:
        raise ValueError(f"Activity expression matrix lacks {len(missing_genes)} duplicate genes")
    for tissue, subset in activity_samples.groupby("tissue", sort=False):
        sample_ids = subset["sample_id"].tolist()
        values = activity_expression[sample_ids]
        median = values.median(axis=1)
        fraction = values.ge(config.min_expression).mean(axis=1)
        tissue_activity[str(tissue)] = pd.DataFrame({
            "mean": values.mean(axis=1),
            "median": median,
            "max": values.max(axis=1),
            "active_fraction": fraction,
            "active": median.ge(config.min_expression)
            & fraction.ge(config.min_active_fraction),
            "n_replicates": values.notna().sum(axis=1),
        })

    key_map: dict[str, tuple[object, bool]] = {}
    for pair in pairs.itertuples(index=False):
        candidate_keys = {
            pair.pair_id: (pair, False),
            f"{pair.gene1}-{pair.gene2}": (pair, False),
            f"{pair.gene2}-{pair.gene1}": (pair, True),
            f"{pair.gene1}__{pair.gene2}": (pair, False),
            f"{pair.gene2}__{pair.gene1}": (pair, True),
        }
        for key, value in candidate_keys.items():
            if key in key_map and key_map[key][0].pair_id != pair.pair_id:
                raise ValueError(f"Legacy pair key is ambiguous: {key}")
            key_map[key] = value

    records: list[dict[str, object]] = []
    warnings: list[str] = []
    tissues_without_adjusted_p: set[str] = set()
    for path in files:
        match = matcher.match(path.name)
        if not match or "tissue" not in match.groupdict():
            warnings.append(f"Skipped {path.name}: filename regex did not yield a tissue")
            continue
        tissue = match.group("tissue")
        if tissue_activity and tissue not in tissue_activity:
            raise ValueError(f"No activity-expression samples found for DE tissue {tissue!r}")
        frame = pd.read_csv(path)
        if frame.empty:
            warnings.append(f"Skipped empty file: {path.name}")
            continue
        lfc_column = _column_by_name(frame, {"log2FoldChange", "log2fc", "lfc"})
        p_column = _column_by_name(frame, {"pvalue", "p_value", "p"})
        q_column = _column_by_name(frame, {"padj", "qvalue", "fdr"})
        if q_column is None:
            tissues_without_adjusted_p.add(tissue)
        if lfc_column is None:
            raise ValueError(f"{path.name} has no log2FoldChange column")
        id_column = str(frame.columns[0])
        if id_column == lfc_column:
            raise ValueError(f"{path.name} needs a row-name/pair-key first column")

        for row in frame.itertuples(index=False, name=None):
            values = dict(zip(frame.columns, row))
            raw_key = str(values[id_column]).strip().strip('"')
            mapped = key_map.get(raw_key)
            if mapped is None:
                warnings.append(f"Unmatched pair key {raw_key!r} in {path.name}")
                continue
            pair, reversed_pair = mapped
            lfc = pd.to_numeric(pd.Series([values[lfc_column]]), errors="coerce").iloc[0]
            if not np.isfinite(lfc):
                lfc = float("nan")
            elif reversed_pair:
                lfc = -float(lfc)
            else:
                lfc = float(lfc)
            pvalue = (
                pd.to_numeric(pd.Series([values[p_column]]), errors="coerce").iloc[0]
                if p_column
                else float("nan")
            )
            qvalue = (
                pd.to_numeric(pd.Series([values[q_column]]), errors="coerce").iloc[0]
                if q_column
                else float("nan")
            )
            activity = tissue_activity[tissue]
            gene1, gene2 = str(pair.gene1), str(pair.gene2)
            mean1 = float(activity.at[gene1, "mean"])
            mean2 = float(activity.at[gene2, "mean"])
            median1 = float(activity.at[gene1, "median"])
            median2 = float(activity.at[gene2, "median"])
            active_fraction1 = float(activity.at[gene1, "active_fraction"])
            active_fraction2 = float(activity.at[gene2, "active_fraction"])
            active1 = bool(activity.at[gene1, "active"])
            active2 = bool(activity.at[gene2, "active"])
            n_replicates = int(min(activity.at[gene1, "n_replicates"], activity.at[gene2, "n_replicates"]))
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
                    "max1": float(activity.at[gene1, "max"]),
                    "max2": float(activity.at[gene2, "max"]),
                    "active1": active1,
                    "active2": active2,
                    "active_fraction1": active_fraction1,
                    "active_fraction2": active_fraction2,
                    "n_replicates": n_replicates,
                    "log2fc_gene2_over_gene1": lfc,
                    "pvalue": float(pvalue) if np.isfinite(pvalue) else float("nan"),
                    "qvalue": float(qvalue) if np.isfinite(qvalue) else float("nan"),
                    "method": "precomputed_differential_expression_with_tpm_activity",
                }
            )

    if not records:
        raise ValueError("Legacy files produced no matched evidence rows")
    evidence = pd.DataFrame.from_records(records)
    if evidence.duplicated(["pair_id", "tissue"]).any():
        raise ValueError("Duplicate pair/tissue evidence; merge biological libraries before differential testing")
    for column in ("pvalue", "qvalue"):
        values = evidence[column].dropna()
        if not values.between(0, 1).all():
            raise ValueError(f"{column} must contain probabilities between zero and one")
    # Insufficient measured biological replication must never become a callable test.
    evidence.loc[evidence.n_replicates < config.min_replicates, ["pvalue", "qvalue"]] = np.nan
    for tissue, indices in evidence.groupby("tissue", sort=False).groups.items():
        if tissue in tissues_without_adjusted_p:
            adjusted = benjamini_hochberg(evidence.loc[indices, "pvalue"])
            evidence.loc[indices, "qvalue"] = adjusted
    test_values = evidence["qvalue"] if config.p_adjust == "bh" else evidence["pvalue"]
    active_union = evidence["active1"] | evidence["active2"]
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
    warnings = list(dict.fromkeys(warnings))
    return evidence, warnings
