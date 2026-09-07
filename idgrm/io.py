"""Input parsing and validation for iDGRM."""

from __future__ import annotations

import csv
import re
from pathlib import Path
from typing import Iterable, Literal

import numpy as np
import pandas as pd


def _canonical(value: object) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(value).strip().lower())


def _unique_headers(headers: list[str]) -> list[str]:
    """Mangle duplicate headers in the same deterministic style as pandas."""

    counts: dict[str, int] = {}
    result: list[str] = []
    for raw in headers:
        name = str(raw).lstrip("\ufeff").strip()
        count = counts.get(name, 0)
        candidate = name if count == 0 else f"{name}.{count}"
        while candidate in counts:
            count += 1
            candidate = f"{name}.{count}"
        counts[name] = count + 1
        counts[candidate] = 1
        result.append(candidate)
    return result


def _read_with_separator(
    path: Path,
    sep: str,
    header: Literal["infer"] | None,
) -> pd.DataFrame:
    options = {
        "sep": sep,
        "dtype": str,
        "keep_default_na": True,
        "encoding": "utf-8-sig",
    }
    if header is None:
        return pd.read_csv(path, header=None, **options)

    # A few original parent/daughter tables contain a trailing role field that
    # has no header. Supplying complete names prevents pandas from silently
    # promoting the first gene ID to an index and shifting every named column.
    preview_rows: list[list[str]] = []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter=sep)
        for row in reader:
            if row:
                preview_rows.append(row)
            if len(preview_rows) >= 100:
                break
    if preview_rows:
        headers = _unique_headers(preview_rows[0])
        max_width = max(len(row) for row in preview_rows)
        if max_width > len(headers):
            headers.extend(
                f"__extra_{index + 1}" for index in range(max_width - len(headers))
            )
            return pd.read_csv(path, header=0, names=headers, **options)
    return pd.read_csv(path, header="infer", **options)


def read_delimited(
    path: str | Path,
    header: Literal["infer"] | None = "infer",
) -> pd.DataFrame:
    """Read CSV/TSV with a deterministic extension-first delimiter policy."""

    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)
    suffix = path.suffix.lower()
    sep = "," if suffix == ".csv" else "\t"
    frame = _read_with_separator(path, sep, header)
    if frame.shape[1] == 1:
        retry_sep = "\t" if sep == "," else ","
        retry = _read_with_separator(path, retry_sep, header)
        if retry.shape[1] > 1:
            frame = retry
    frame.columns = [str(column).lstrip("\ufeff").strip() for column in frame.columns]
    return frame


def _find_column(columns: Iterable[object], aliases: set[str]) -> str | None:
    for column in columns:
        if _canonical(column) in aliases:
            return str(column)
    return None


def _looks_like_gene_identifier(value: object) -> bool:
    text = str(value).strip()
    canonical = _canonical(text)
    if not text or any(character.isspace() for character in text):
        return False
    if re.fullmatch(r"(?:gene|copy|duplicate|paralog|id)\d*", canonical):
        return False
    return any(character.isdigit() for character in text)


def _infer_duplication_type(path: str | Path) -> str:
    name = Path(path).name.lower()
    for token, label in (
        ("tandem", "tandem"),
        ("proximal", "proximal"),
        ("transposed", "transposed"),
        ("segmental", "segmental/WGD"),
        ("wgd", "WGD"),
        ("dispersed", "dispersed"),
    ):
        if token in name:
            return label
    return "unspecified"


def read_expression(path: str | Path) -> pd.DataFrame:
    """Read a nonnegative genes-by-samples/tissues expression matrix."""

    frame = read_delimited(path)
    if frame.shape[1] < 2:
        raise ValueError("Expression matrix must have a gene column and at least one data column")

    gene_aliases = {"gene", "geneid", "id", "transcript", "transcriptid", "ancestorid"}
    gene_column = _find_column(frame.columns, gene_aliases) or str(frame.columns[0])
    frame = frame.rename(columns={gene_column: "gene_id"})
    if frame["gene_id"].isna().any():
        raise ValueError("Expression matrix contains a missing gene ID")
    frame["gene_id"] = frame["gene_id"].astype(str).str.strip()
    if frame["gene_id"].eq("").any():
        raise ValueError("Expression matrix contains an empty gene ID")
    duplicated = frame.loc[frame["gene_id"].duplicated(), "gene_id"].unique().tolist()
    if duplicated:
        preview = ", ".join(map(str, duplicated[:5]))
        raise ValueError(f"Expression matrix contains duplicate gene IDs: {preview}")

    value_columns = [column for column in frame.columns if column != "gene_id"]
    raw_values = frame[value_columns].copy()
    numeric_values = raw_values.apply(pd.to_numeric, errors="coerce")
    invalid = numeric_values.isna() & raw_values.notna()
    if invalid.any().any():
        row_index, column_index = np.argwhere(invalid.to_numpy())[0]
        bad_column = value_columns[int(column_index)]
        bad_gene = frame.iloc[int(row_index)]["gene_id"]
        bad_value = raw_values.iloc[int(row_index), int(column_index)]
        raise ValueError(
            f"Expression matrix contains a nonnumeric value at gene {bad_gene!r}, "
            f"column {bad_column!r}: {bad_value!r}"
        )
    frame[value_columns] = numeric_values
    values = frame[value_columns].to_numpy(dtype=float)
    if np.isinf(values).any():
        raise ValueError("Expression matrix contains infinite values")
    if not np.isfinite(values).any():
        raise ValueError("Expression matrix contains no finite expression values")
    if np.nanmin(values) < 0:
        raise ValueError("Expression values must be nonnegative")
    return frame.set_index("gene_id")


def normalize_expression(expression: pd.DataFrame, method: str = "none") -> pd.DataFrame:
    """Normalize columns using a lightweight, documented transformation."""

    method = method.lower().replace("_", "-")
    values = expression.astype(float).copy()
    if method == "none":
        return values
    if method == "cpm":
        totals = values.sum(axis=0, skipna=True)
        if (totals <= 0).any():
            bad = ", ".join(totals.index[totals <= 0].astype(str))
            raise ValueError(f"CPM normalization cannot use zero-sum columns: {bad}")
        return values.divide(totals, axis=1) * 1_000_000.0
    if method in {"median-ratio", "medianratio"}:
        matrix = values.to_numpy(dtype=float)
        with np.errstate(divide="ignore", invalid="ignore"):
            log_matrix = np.where(matrix > 0, np.log(matrix), np.nan)
            geometric_means = np.exp(np.nanmean(log_matrix, axis=1))
        valid_genes = np.isfinite(geometric_means) & (geometric_means > 0)
        if not valid_genes.any():
            raise ValueError("Median-ratio normalization found no genes with positive counts")
        ratios = matrix[valid_genes, :] / geometric_means[valid_genes, None]
        ratios[~np.isfinite(ratios) | (ratios <= 0)] = np.nan
        size_factors = np.nanmedian(ratios, axis=0)
        if np.any(~np.isfinite(size_factors)) or np.any(size_factors <= 0):
            raise ValueError("Could not estimate a positive median-ratio size factor for every sample")
        return values.divide(pd.Series(size_factors, index=values.columns), axis=1)
    raise ValueError("normalization must be one of: none, cpm, median-ratio")


def read_samples(path: str | Path, expression_columns: Iterable[str]) -> pd.DataFrame:
    """Read sample metadata and ensure every selected sample exists in the matrix."""

    frame = read_delimited(path)
    if frame.empty:
        raise ValueError("Sample metadata contains no samples")
    sample_column = _find_column(
        frame.columns,
        {"sample", "sampleid", "samplename", "library", "libraryid"},
    )
    tissue_column = _find_column(
        frame.columns,
        {"tissue", "organ", "group", "condition", "tissuename", "organname"},
    )
    if sample_column is None or tissue_column is None:
        raise ValueError("Sample metadata needs sample_id and tissue columns")
    frame = frame.rename(columns={sample_column: "sample_id", tissue_column: "tissue"})
    if frame["sample_id"].isna().any() or frame["tissue"].isna().any():
        raise ValueError("Sample metadata contains a missing sample ID or tissue name")
    frame["sample_id"] = frame["sample_id"].astype(str).str.strip()
    frame["tissue"] = frame["tissue"].astype(str).str.strip()
    if frame["sample_id"].eq("").any():
        raise ValueError("Sample metadata contains an empty sample ID")
    if frame["sample_id"].duplicated().any():
        duplicates = frame.loc[frame["sample_id"].duplicated(), "sample_id"].tolist()
        raise ValueError(f"Sample metadata contains duplicate sample IDs: {duplicates[:5]}")
    if frame["tissue"].eq("").any():
        raise ValueError("Sample metadata contains an empty tissue name")

    expression_set = set(map(str, expression_columns))
    missing = sorted(set(frame["sample_id"]) - expression_set)
    if missing:
        raise ValueError(
            "Samples listed in metadata are missing from the expression matrix: "
            + ", ".join(missing[:10])
        )
    return frame[["sample_id", "tissue"] + [
        column for column in frame.columns if column not in {"sample_id", "tissue"}
    ]]


def read_pairs(path: str | Path) -> pd.DataFrame:
    """Read modern or legacy duplicate-pair tables into a canonical schema."""

    frame = read_delimited(path)
    if frame.empty:
        raise ValueError("Pair table contains no duplicate-gene pairs")
    input_warnings: list[str] = []
    aliases1 = {
        "gene1", "copy1", "duplicate1", "paralog1", "parent", "parental",
        "parentcopy", "oldercopy",
    }
    aliases2 = {
        "gene2", "copy2", "duplicate2", "paralog2", "daughter", "daughtercopy",
        "newcopy", "transposed",
    }
    gene1_column = _find_column(frame.columns, aliases1)
    gene2_column = _find_column(frame.columns, aliases2)
    headerless = False
    if gene1_column is None or gene2_column is None:
        columns = list(frame.columns)
        if (
            len(columns) >= 3
            and _looks_like_gene_identifier(columns[0])
            and _looks_like_gene_identifier(columns[2])
        ):
            frame = read_delimited(path, header=None)
            gene1_column = "0"
            gene2_column = "2"
            headerless = True
            input_warnings.append(
                f"Detected a headerless pair table and interpreted columns 1 and 3 as genes: "
                f"{Path(path).name}"
            )
        else:
            raise ValueError(
                "Pair table needs gene1/gene2 columns (legacy 'Duplicate 1/2', "
                "'Parent copy/Daughter copy', 'Parental/Transposed', and headerless "
                "gene/type/gene/type layouts are supported)"
            )
    if frame[gene1_column].isna().any() or frame[gene2_column].isna().any():
        raise ValueError("Pair table contains a missing gene ID")

    pair_id_column = _find_column(frame.columns, {"pair", "pairid", "duplicatepairid"})
    type_column = _find_column(
        frame.columns,
        {"duplicationtype", "duplicatetype", "mechanism", "mode", "duplicationmode"},
    )
    ancestor_column = _find_column(
        frame.columns,
        {"ancestor", "ancestorid", "ancestralgene", "outgroup", "outgroupgene", "ortholog"},
    )
    role1_column = _find_column(frame.columns, {"role1", "gene1role", "copy1role"})
    role2_column = _find_column(frame.columns, {"role2", "gene2role", "copy2role"})

    canonical = pd.DataFrame()
    canonical["gene1"] = frame[gene1_column].astype(str).str.strip()
    canonical["gene2"] = frame[gene2_column].astype(str).str.strip()
    if canonical[["gene1", "gene2"]].eq("").any().any():
        raise ValueError("Pair table contains an empty gene ID")
    if (canonical["gene1"] == canonical["gene2"]).any():
        bad = canonical.loc[canonical["gene1"] == canonical["gene2"], "gene1"].iloc[0]
        raise ValueError(f"A duplicate pair contains the same gene twice: {bad}")

    if pair_id_column:
        if frame[pair_id_column].isna().any():
            raise ValueError("Pair table contains a missing pair ID")
        canonical["pair_id"] = frame[pair_id_column].astype(str).str.strip()
    else:
        canonical["pair_id"] = canonical["gene1"] + "__" + canonical["gene2"]
    if canonical["pair_id"].eq("").any():
        raise ValueError("Pair table contains an empty pair ID")
    canonical["duplication_type"] = (
        frame[type_column].fillna("unspecified").astype(str).str.strip()
        if type_column
        else _infer_duplication_type(path)
    )
    canonical.loc[canonical["duplication_type"].eq(""), "duplication_type"] = "unspecified"
    canonical["ancestor_id"] = (
        frame[ancestor_column].fillna("").astype(str).str.strip() if ancestor_column else ""
    )

    source1 = _canonical(gene1_column)
    source2 = _canonical(gene2_column)
    inferred_role1 = "parent" if "parent" in source1 or "older" in source1 else ""
    inferred_role2 = (
        "daughter"
        if "daughter" in source2 or "new" in source2 or "transposed" in source2
        else ""
    )
    canonical["role1"] = (
        frame[role1_column].fillna("").astype(str).str.strip()
        if role1_column
        else inferred_role1
    )
    canonical["role2"] = (
        frame[role2_column].fillna("").astype(str).str.strip()
        if role2_column
        else inferred_role2
    )

    duplicate_ids = canonical.loc[canonical["pair_id"].duplicated(keep=False), "pair_id"].unique()
    for pair_id in duplicate_ids:
        group = canonical[canonical["pair_id"].eq(pair_id)]
        if len(group[["gene1", "gene2"]].drop_duplicates()) > 1:
            raise ValueError(f"Pair table reuses pair ID {pair_id!r} for different gene pairs")
    n_duplicates = int(canonical["pair_id"].duplicated().sum())
    if n_duplicates:
        canonical = canonical.drop_duplicates("pair_id", keep="first")
        input_warnings.append(
            f"Removed {n_duplicates} exact duplicate pair row(s) from {Path(path).name}."
        )

    result = canonical[
        ["pair_id", "gene1", "gene2", "duplication_type", "role1", "role2", "ancestor_id"]
    ].reset_index(drop=True)
    result.attrs["warnings"] = input_warnings
    result.attrs["headerless"] = headerless
    return result


def aggregate_tissue_expression(
    expression: pd.DataFrame,
    samples: pd.DataFrame | None,
) -> pd.DataFrame:
    """Return gene-by-tissue median expression for sample- or tissue-level input."""

    if samples is None:
        result = expression.copy()
        result.columns = result.columns.astype(str)
        return result
    tissue_values: dict[str, pd.Series] = {}
    for tissue, subset in samples.groupby("tissue", sort=False):
        sample_ids = subset["sample_id"].tolist()
        tissue_values[str(tissue)] = expression[sample_ids].median(axis=1, skipna=True)
    return pd.DataFrame(tissue_values, index=expression.index)
