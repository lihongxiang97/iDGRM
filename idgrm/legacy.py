"""Compatibility loader for the original per-tissue DESeq2 CSV workflow."""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

from .config import IDGRMConfig
from .statistics import benjamini_hochberg


DEFAULT_FILENAME_REGEX = r"^[^.]+\.[^.]+\.pairs\.(?P<tissue>.+)\.DESeq2\.csv$"


def _column_by_name(frame: pd.DataFrame, candidates: set[str]) -> str | None:
    mapping = {str(column).strip().lower(): str(column) for column in frame.columns}
    for candidate in candidates:
        if candidate.lower() in mapping:
            return mapping[candidate.lower()]
    return None


def load_legacy_deseq2(
    directory: str | Path,
    pairs: pd.DataFrame,
    config: IDGRMConfig,
    pattern: str = "*.DESeq2.csv",
    filename_regex: str = DEFAULT_FILENAME_REGEX,
) -> tuple[pd.DataFrame, list[str]]:
    """Convert old DESeq2 CSV files to iDGRM's tidy evidence table."""

    directory = Path(directory)
    files = sorted(directory.glob(pattern))
    if not files:
        raise FileNotFoundError(f"No legacy DESeq2 files matched {directory / pattern}")
    matcher = re.compile(filename_regex)

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
    for path in files:
        match = matcher.match(path.name)
        if not match or "tissue" not in match.groupdict():
            warnings.append(f"Skipped {path.name}: filename regex did not yield a tissue")
            continue
        tissue = match.group("tissue")
        frame = pd.read_csv(path)
        if frame.empty:
            warnings.append(f"Skipped empty file: {path.name}")
            continue
        lfc_column = _column_by_name(frame, {"log2FoldChange", "log2fc", "lfc"})
        p_column = _column_by_name(frame, {"pvalue", "p_value", "p"})
        q_column = _column_by_name(frame, {"padj", "qvalue", "fdr"})
        base_column = _column_by_name(frame, {"baseMean", "basemean", "mean"})
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
            base_mean = (
                pd.to_numeric(pd.Series([values[base_column]]), errors="coerce").iloc[0]
                if base_column
                else float("nan")
            )
            if np.isfinite(base_mean) and np.isfinite(lfc):
                ratio = 2.0 ** float(lfc)
                mean1 = 2.0 * float(base_mean) / (1.0 + ratio)
                mean2 = mean1 * ratio
            else:
                mean1 = mean2 = float("nan")
            active1 = bool(np.isfinite(mean1) and mean1 >= config.min_expression)
            active2 = bool(np.isfinite(mean2) and mean2 >= config.min_expression)
            records.append(
                {
                    "pair_id": pair.pair_id,
                    "gene1": pair.gene1,
                    "gene2": pair.gene2,
                    "tissue": tissue,
                    "mean1": mean1,
                    "mean2": mean2,
                    "median1": mean1,
                    "median2": mean2,
                    "active1": active1,
                    "active2": active2,
                    "active_fraction1": float(active1),
                    "active_fraction2": float(active2),
                    "n_replicates": float("nan"),
                    "log2fc_gene2_over_gene1": lfc,
                    "pvalue": float(pvalue) if np.isfinite(pvalue) else float("nan"),
                    "qvalue": float(qvalue) if np.isfinite(qvalue) else float("nan"),
                    "method": "legacy_deseq2",
                }
            )

    if not records:
        raise ValueError("Legacy files produced no matched evidence rows")
    evidence = pd.DataFrame.from_records(records)
    for _, indices in evidence.groupby("tissue", sort=False).groups.items():
        missing_q = evidence.loc[indices, "qvalue"].isna()
        if missing_q.any():
            adjusted = benjamini_hochberg(evidence.loc[indices, "pvalue"])
            evidence.loc[indices, "qvalue"] = evidence.loc[indices, "qvalue"].fillna(adjusted)
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
