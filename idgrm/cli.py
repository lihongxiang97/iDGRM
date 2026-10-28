"""Command-line interface for the current iDGRM workflow."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import __version__
from .config import IDGRMConfig
from .legacy import DEFAULT_DIFFERENTIAL_EXPRESSION_FILENAME_REGEX
from .pipeline import (
    run_primary_classification,
    run_primary_classification_from_differential_expression,
)


def _add_configuration_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--absolute-log2-fold-change", type=float, default=1.0)
    parser.add_argument("--significance-threshold", type=float, default=0.05)
    parser.add_argument("--multiple-testing-correction", choices=["bh", "none"], default="bh")
    parser.add_argument("--expression-threshold", type=float, default=1.0)
    parser.add_argument("--minimum-active-replicate-fraction", type=float, default=0.5)
    parser.add_argument("--log-ratio-pseudocount", type=float, default=0.1)
    parser.add_argument("--minimum-biological-replicates", type=int, default=2)
    parser.add_argument("--minimum-evaluable-tissues", type=int, default=2)
    parser.add_argument("--asymmetry-tissue-fraction", type=float, default=1.0 / 3.0)


def _configuration(args: argparse.Namespace) -> IDGRMConfig:
    return IDGRMConfig(
        log2fc_threshold=args.absolute_log2_fold_change,
        alpha=args.significance_threshold,
        p_adjust=args.multiple_testing_correction,
        min_expression=args.expression_threshold,
        min_active_fraction=args.minimum_active_replicate_fraction,
        pseudocount=args.log_ratio_pseudocount,
        min_replicates=args.minimum_biological_replicates,
        min_evaluable_tissues=args.minimum_evaluable_tissues,
        asymmetry_tissue_fraction=args.asymmetry_tissue_fraction,
    ).validate()


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="idgrm",
        description="Four-class classification of duplicate-gene expression patterns.",
    )
    parser.add_argument("--version", action="version", version=f"iDGRM {__version__}")
    commands = parser.add_subparsers(dest="command", required=True)

    counts = commands.add_parser("fit-copy-counts", help="Fit raw integer copy counts with paired DESeq2 (requires R/DESeq2).")
    counts.add_argument("--counts",required=True,help="Whole-gene integer matrix used for size factors.")
    counts.add_argument("--duplicate-pairs",required=True)
    counts.add_argument("--sample-metadata",required=True)
    counts.add_argument("--copy1-counts")
    counts.add_argument("--copy2-counts")
    counts.add_argument("--observability-method",required=True,choices=["reciprocal-common-cds","equal-effective-length","uncorrected"])
    counts.add_argument("--test",choices=["wald","lrt"],default="wald")
    counts.add_argument("--rscript",default="Rscript")
    counts.add_argument("--significance-threshold",type=float,default=.05)
    counts.add_argument("--output-dir",required=True)

    primary = commands.add_parser(
        "classify-primary",
        help="Assign four primary expression-pattern classes to duplicate-gene pairs.",
    )
    source = primary.add_mutually_exclusive_group(required=True)
    source.add_argument("--expression-matrix")
    source.add_argument("--differential-expression-dir")
    primary.add_argument("--sample-metadata")
    primary.add_argument(
        "--activity-expression-matrix",
        help="Measured TPM matrix for copy activity with precomputed differential-expression tests.",
    )
    primary.add_argument(
        "--activity-sample-metadata", help="Sample-to-tissue metadata for the activity TPM matrix."
    )
    primary.add_argument("--duplicate-pairs", required=True)
    primary.add_argument("--output-dir", required=True)
    primary.add_argument("--normalization", choices=["none", "cpm", "median-ratio"], default="none")
    primary.add_argument("--differential-expression-file-pattern", default="*.DESeq2.csv")
    primary.add_argument(
        "--tissue-name-regex", default=DEFAULT_DIFFERENTIAL_EXPRESSION_FILENAME_REGEX
    )
    _add_configuration_arguments(primary)

    polarize = commands.add_parser("polarize-outgroup", help="Optional outgroup-supported expression candidates; never modifies primary classes.")
    for name in ("primary-results-dir", "outgroup-expression", "outgroup-samples", "orthologs",
                 "tissue-map", "detection-qc", "outgroup-species", "reference-version", "orthology-source", "output-dir"):
        polarize.add_argument("--" + name, required=True)
    polarize.add_argument("--on-threshold", type=float, default=1.)
    polarize.add_argument("--off-threshold", type=float, default=.1)
    polarize.add_argument("--minimum-replicates", type=int, default=2)
    polarize.add_argument("--minimum-matched-tissues", type=int, default=2)
    polarize.add_argument("--retained-coverage", type=float, default=.8)
    polarize.add_argument("--allow-approximate-tissues", action="store_true")
    polarize.add_argument("--minimum-informative-fraction", type=float, default=.5)

    return parser


def _print_result(result: object) -> None:
    print("iDGRM analysis completed.")
    for name, path in result.output_paths.items():
        print(f"  {name}: {Path(path).resolve()}")
    if result.warnings:
        print("Warnings:")
        for warning in result.warnings:
            print(f"  - {warning}")


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        if args.command == "polarize-outgroup":
            from .polarization import PolarizationConfig, run_outgroup_polarization
            result = run_outgroup_polarization(
                args.primary_results_dir, args.outgroup_expression, args.outgroup_samples,
                args.orthologs, args.tissue_map, args.detection_qc, args.output_dir,
                outgroup_species=args.outgroup_species, reference_version=args.reference_version,
                orthology_source=args.orthology_source,
                config=PolarizationConfig(args.on_threshold, args.off_threshold, args.minimum_replicates,
                    args.minimum_matched_tissues, args.retained_coverage, args.allow_approximate_tissues,
                    args.minimum_informative_fraction))
            _print_result(result)
            return 0
        if args.command=="fit-copy-counts":
            from .counts import fit_copy_counts
            path=fit_copy_counts(args.counts,args.duplicate_pairs,args.sample_metadata,args.output_dir,
                rscript=args.rscript,test=args.test,copy1_path=args.copy1_counts,copy2_path=args.copy2_counts,
                observability_method=args.observability_method,alpha=args.significance_threshold)
            print(f"Count fitting completed: {path}")
            print("Caution: high-variance small-sample error control and final-category FDR are not universally established.")
            return 0
        config = _configuration(args)
        if args.command == "classify-primary":
            if args.expression_matrix:
                if args.activity_expression_matrix or args.activity_sample_metadata:
                    raise ValueError(
                        "--activity-expression-matrix and --activity-sample-metadata "
                        "require --differential-expression-dir"
                    )
                if args.sample_metadata is None:
                    raise ValueError("--sample-metadata is required for replicate-aware inference")
                result = run_primary_classification(
                    expression_matrix_path=args.expression_matrix,
                    sample_metadata_path=args.sample_metadata,
                    duplicate_pairs_path=args.duplicate_pairs,
                    output_dir=args.output_dir,
                    config=config,
                    normalization=args.normalization,
                )
            else:
                if args.sample_metadata is not None:
                    raise ValueError("--sample-metadata is not used with --differential-expression-dir")
                if args.activity_expression_matrix is None or args.activity_sample_metadata is None:
                    raise ValueError(
                        "--differential-expression-dir requires --activity-expression-matrix "
                        "and --activity-sample-metadata"
                    )
                result = run_primary_classification_from_differential_expression(
                    differential_expression_dir=args.differential_expression_dir,
                    duplicate_pairs_path=args.duplicate_pairs,
                    output_dir=args.output_dir,
                    config=config,
                    file_pattern=args.differential_expression_file_pattern,
                    tissue_name_regex=args.tissue_name_regex,
                    activity_expression_matrix_path=args.activity_expression_matrix,
                    activity_sample_metadata_path=args.activity_sample_metadata,
                )
        _print_result(result)
        return 0
    except Exception as exc:
        print(f"iDGRM error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
