"""Optional, explicit outgroup-reference expression polarization (not ancestral reconstruction)."""
from __future__ import annotations

from dataclasses import asdict, dataclass
import json
from math import isfinite
from pathlib import Path

import numpy as np
import pandas as pd

from .io import read_delimited, read_expression, read_samples
from .pipeline import AnalysisResult, _input_record
from ._version import __version__


@dataclass(frozen=True)
class PolarizationConfig:
    on_threshold: float = 1.0
    off_threshold: float = 0.1
    minimum_replicates: int = 2
    minimum_matched_tissues: int = 2
    retained_coverage: float = 0.8
    allow_approximate_tissues: bool = False
    minimum_informative_fraction: float = 0.5

    def validate(self):
        if not all(isfinite(v) for v in (self.on_threshold, self.off_threshold, self.retained_coverage, self.minimum_informative_fraction)):
            raise ValueError("Polarization thresholds must be finite")
        if not 0 <= self.off_threshold < self.on_threshold:
            raise ValueError("Require 0 <= off_threshold < on_threshold")
        if not 0 < self.retained_coverage <= 1:
            raise ValueError("retained_coverage must be in (0,1]")
        if not 0 < self.minimum_informative_fraction <= 1:
            raise ValueError("minimum_informative_fraction must be in (0,1]")
        for value in (self.minimum_replicates, self.minimum_matched_tissues):
            if isinstance(value, bool) or not isinstance(value, int) or value < 2:
                raise ValueError("Minimum replicates and matched tissues must be integers >=2")
        return self


def _table(path, columns, unique):
    frame = read_delimited(path)
    if not set(columns).issubset(frame):
        raise ValueError(f"{Path(path).name} requires columns {columns}")
    if frame.empty or frame[list(columns)].isna().any().any():
        raise ValueError(f"{Path(path).name} is empty or contains missing required fields")
    for column in columns:
        frame[column] = frame[column].astype(str).str.strip()
        if frame[column].eq("").any():
            raise ValueError(f"Empty {column}")
    if frame.duplicated(unique).any():
        raise ValueError(f"Ambiguous or duplicate keys in {Path(path).name}: {unique}")
    return frame


def _state(values, config):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if len(values) < config.minimum_replicates:
        return "UNKNOWN"
    if np.median(values) >= config.on_threshold and np.mean(values >= config.on_threshold) >= .5:
        return "ON"
    # Absence requires all measured replicates low, not merely a low mean.
    if np.all(values <= config.off_threshold):
        return "OFF"
    return "UNKNOWN"


def polarize_patterns(states: pd.DataFrame, config: PolarizationConfig):
    """Pure rule function on qualified ON/OFF states; UNKNOWN never equals OFF."""
    config.validate()
    complete = states.loc[states[["copy1_state", "copy2_state", "outgroup_state"]].isin(["ON", "OFF"]).all(axis=1)]
    metrics = {"n_matched_tissues": len(states), "n_informative_tissues": len(complete)}
    metrics["informative_fraction"] = len(complete) / len(states) if len(states) else 0.
    if len(complete) < config.minimum_matched_tissues:
        return "UNRESOLVED", "Insufficient matched, qualified expression states", metrics
    if metrics["informative_fraction"] < config.minimum_informative_fraction:
        return "UNRESOLVED", "Insufficient informative coverage of mapped tissues", metrics
    a = complete.outgroup_state.eq("ON")
    x = complete.copy1_state.eq("ON")
    y = complete.copy2_state.eq("ON")
    n_ancestral = int(a.sum())
    # Expressions are compared as states, never as cross-species TPM fold changes.
    gain1 = int((~a & x & ~y & complete.direction.eq("gene1_high")).sum())
    gain2 = int((~a & y & ~x & complete.direction.eq("gene2_high")).sum())
    partition1 = int((a & x & ~y & complete.direction.eq("gene1_high")).sum())
    partition2 = int((a & y & ~x & complete.direction.eq("gene2_high")).sum())
    coverage1 = float((a & x).sum() / n_ancestral) if n_ancestral else 0.
    coverage2 = float((a & y).sum() / n_ancestral) if n_ancestral else 0.
    union = float((a & (x | y)).sum() / n_ancestral) if n_ancestral else 0.
    metrics.update(n_outgroup_on=n_ancestral, n_gain_copy1=gain1, n_gain_copy2=gain2,
        n_partition_copy1=partition1, n_partition_copy2=partition2,
        retention_copy1=coverage1, retention_copy2=coverage2, union_retention=union)
    if partition1 and partition2 and (gain1 or gain2):
        return "UNRESOLVED", "Mixed partition and gain evidence", metrics
    if partition1 and partition2 and not (gain1 or gain2) and union >= config.retained_coverage:
        return "OUTGROUP_SUPPORTED_SUB_CANDIDATE", "Complementary retention of outgroup-expressed tissues", metrics
    if bool(gain1) != bool(gain2):
        retained = coverage2 if gain1 else coverage1
        if n_ancestral and retained >= config.retained_coverage:
            metrics["candidate_gain_copy"] = "copy1" if gain1 else "copy2"
            return "OUTGROUP_SUPPORTED_NEO_CANDIDATE", "One-copy gain relative to qualified outgroup absence; partner retains outgroup expression", metrics
    return "UNRESOLVED", "No unique partition or one-copy gain interpretation", metrics


def run_outgroup_polarization(primary_results_dir, outgroup_expression_path,
        outgroup_samples_path, orthologs_path, tissue_map_path, detection_qc_path,
        output_dir, *, outgroup_species, reference_version, orthology_source,
        config=None):
    """Polarize only primary SUB/NEO pairs, with explicit mappings and detection QC.

    QC is a user declaration of gene/tissue detectability, not a validated probability.
    Primary tables and evidence are immutable inputs. No expression-only fallback.
    """
    config = (config or PolarizationConfig()).validate()
    for text in (outgroup_species, reference_version, orthology_source):
        if not str(text).strip():
            raise ValueError("Species, reference version and orthology source are required")
    root = Path(primary_results_dir)
    primary_path = root / "primary_classifications.tsv"
    evidence_path = root / "tissue_evidence.tsv"
    metadata_path = root / "primary_run_metadata.json"
    primary = _table(primary_path, ["pair_id", "gene1", "gene2", "primary_class"], ["pair_id"])
    required = ["pair_id", "tissue", "median1", "median2", "active1", "active2", "n_replicates", "direction"]
    evidence = _table(evidence_path, ["pair_id", "tissue", "direction"], ["pair_id", "tissue"])
    if not set(required).issubset(evidence):
        raise ValueError(f"Tissue evidence requires {required}")
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    if metadata.get("software") != "iDGRM":
        raise ValueError("Primary metadata must identify iDGRM")
    from .pipeline import _sha256
    recorded = metadata.get("output_sha256", {})
    for path in (primary_path, evidence_path):
        if recorded.get(path.name) != _sha256(path):
            raise ValueError("Primary outputs lack matching integrity hashes; regenerate primary results with the current iDGRM release")
    expected_on = metadata.get("configuration", {}).get("min_expression")
    if expected_on != config.on_threshold:
        raise ValueError("Polarization ON threshold must match primary expression threshold")
    if not set(evidence.pair_id).issubset(set(primary.pair_id)):
        raise ValueError("Evidence contains unknown pair IDs")
    if not set(primary.primary_class).issubset({"UNMAPPED", "NO_DIFFERENCE", "ASYMMETRICALLY_EXPRESSED", "SUB_OR_NEOFUNCTIONALIZED"}):
        raise ValueError("Unknown primary class")
    for column in ("active1", "active2"):
        converted = evidence[column].astype(str).str.lower().map({"true": True, "false": False})
        if converted.isna().any():
            raise ValueError(f"Invalid boolean {column}")
        evidence[column] = converted.astype(bool)
    for column in ("median1", "median2", "n_replicates"):
        evidence[column] = pd.to_numeric(evidence[column], errors="raise")
        values = evidence[column].dropna()
        if np.isinf(values).any() or values.lt(0).any():
            raise ValueError(f"Invalid numeric evidence: {column}")
    if not evidence.direction.isin(["gene1_high", "gene2_high", "no_difference", "both_inactive"]).all():
        raise ValueError("Unknown tissue direction")
    orthologs = _table(orthologs_path, ["pair_id", "outgroup_gene_id", "relationship"], ["pair_id"])
    if not orthologs.relationship.isin(["single_copy", "ambiguous"]).all():
        raise ValueError("Ortholog relationship must be single_copy or ambiguous")
    if not set(orthologs.pair_id).issubset(set(primary.pair_id)):
        raise ValueError("Ortholog table contains unknown pair IDs")
    mapping = _table(tissue_map_path, ["target_tissue", "outgroup_tissue", "match_quality"], ["target_tissue"])
    if mapping.outgroup_tissue.duplicated().any():
        raise ValueError("Tissue mapping must be one-to-one; do not reuse an outgroup tissue")
    if not mapping.match_quality.isin(["exact", "approximate"]).all():
        raise ValueError("match_quality must be exact or approximate")
    if not set(mapping.target_tissue).issubset(set(evidence.tissue)):
        raise ValueError("Tissue mapping contains unknown target tissues")
    expression = read_expression(outgroup_expression_path)
    samples = read_samples(outgroup_samples_path, expression.columns)
    if not set(mapping.outgroup_tissue).issubset(set(samples.tissue)):
        raise ValueError("Tissue mapping references absent outgroup samples")
    qc = _table(detection_qc_path, ["outgroup_gene_id", "outgroup_tissue", "assessment"], ["outgroup_gene_id", "outgroup_tissue"])
    if not qc.assessment.isin(["adequate", "inadequate"]).all():
        raise ValueError("Detection assessment must be adequate or inadequate")
    qc_lookup = qc.set_index(["outgroup_gene_id", "outgroup_tissue"]).assessment.to_dict()
    ortholog_lookup = orthologs.set_index("pair_id").to_dict("index")
    states_records, results = [], []
    lookup = evidence.set_index(["pair_id", "tissue"])
    for pair in primary.itertuples(index=False):
        result = dict(pair_id=pair.pair_id, gene1=pair.gene1, gene2=pair.gene2,
            parent_primary_class=pair.primary_class, outgroup_species=outgroup_species,
            polarization_class="NOT_ELIGIBLE", polarization_reason="Only primary SUB/NEO pairs are polarized")
        if pair.primary_class == "SUB_OR_NEOFUNCTIONALIZED":
            ortholog = ortholog_lookup.get(pair.pair_id)
            rows = []
            if ortholog and ortholog["relationship"] == "single_copy":
                gene = ortholog["outgroup_gene_id"]
                result["outgroup_gene_id"] = gene
                for match in mapping.itertuples(index=False):
                    if match.match_quality == "approximate" and not config.allow_approximate_tissues:
                        continue
                    if (pair.pair_id, match.target_tissue) not in lookup.index:
                        continue
                    row = lookup.loc[(pair.pair_id, match.target_tissue)]
                    states = []
                    for copy in (1, 2):
                        if row.n_replicates < config.minimum_replicates or not np.isfinite(row[f"median{copy}"]):
                            state = "UNKNOWN"
                        elif row[f"active{copy}"]:
                            state = "ON"
                        elif row[f"median{copy}"] <= config.off_threshold:
                            maximum = float(row.get(f"max{copy}", np.nan))
                            state = "OFF" if np.isfinite(maximum) and maximum <= config.off_threshold else "UNKNOWN"
                        else:
                            state = "UNKNOWN"
                        states.append(state)
                    out_state = "UNKNOWN"
                    adequate = qc_lookup.get((gene, match.outgroup_tissue)) == "adequate"
                    if gene in expression.index and adequate:
                        ids = samples.loc[samples.tissue.eq(match.outgroup_tissue), "sample_id"]
                        out_state = _state(expression.loc[gene, ids], config)
                    record = dict(pair_id=pair.pair_id, target_tissue=match.target_tissue,
                        outgroup_tissue=match.outgroup_tissue, match_quality=match.match_quality,
                        outgroup_gene_id=gene, detection_qc="adequate" if adequate else "not_adequate_or_missing",
                        copy1_state=states[0], copy2_state=states[1], outgroup_state=out_state,
                        direction=row.direction)
                    rows.append(record)
                    states_records.append(record)
            if rows:
                label, reason, metrics = polarize_patterns(pd.DataFrame(rows), config)
                result.update(polarization_class=label, polarization_reason=reason, **metrics)
            else:
                result.update(polarization_class="UNRESOLVED", polarization_reason="Missing unambiguous ortholog or usable tissue correspondence")
        results.append(result)
    classifications = pd.DataFrame(results)
    states = pd.DataFrame(states_records, columns=["pair_id", "target_tissue", "outgroup_tissue", "match_quality", "outgroup_gene_id", "detection_qc", "copy1_state", "copy2_state", "outgroup_state", "direction"])
    summary = classifications.polarization_class.value_counts().rename_axis("class").reset_index(name="pair_count")
    out = Path(output_dir)
    # Reserve a new directory, so concurrent or repeated runs cannot overwrite results.
    out.mkdir(parents=True, exist_ok=False)
    paths = {"polarization": out / "outgroup_polarization.tsv", "evidence": out / "outgroup_tissue_evidence.tsv",
        "summary": out / "outgroup_summary.tsv", "metadata": out / "polarization_run_metadata.json"}
    warnings = ["Outgroup expression is an ancestral-state proxy, not a reconstructed ancestor or proof of functional change.",
        "Orthology, tissue correspondence and detection adequacy are user supplied; no candidate-level FDR is asserted."]
    provenance = dict(software="iDGRM", version=__version__, module="OUTGROUP_EXPRESSION_POLARIZATION",
        configuration=asdict(config), outgroup_species=outgroup_species, reference_version=reference_version,
        orthology_source=orthology_source, warnings=warnings,
        inputs=[_input_record(p) for p in [primary_path, evidence_path, metadata_path, outgroup_expression_path,
            outgroup_samples_path, orthologs_path, tissue_map_path, detection_qc_path]])
    classifications.to_csv(paths["polarization"], sep="\t", index=False, na_rep="")
    states.to_csv(paths["evidence"], sep="\t", index=False, na_rep="")
    summary.to_csv(paths["summary"], sep="\t", index=False)
    paths["metadata"].write_text(json.dumps(provenance, ensure_ascii=False, indent=2), encoding="utf-8")
    return AnalysisResult(classifications, states, summary, paths, warnings)
