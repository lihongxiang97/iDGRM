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
    run_subtype_refinement,
)


def _add_configuration_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--absolute-log2-fold-change", type=float, default=1.0)
    parser.add_argument("--significance-threshold", type=float, default=0.05)
    parser.add_argument("--multiple-testing-correction", choices=["bh", "none"], default="bh")
    parser.add_argument("--expression-threshold", type=float, default=1.0)
    parser.add_argument("--outgroup-expression-threshold", type=float, default=None)
    parser.add_argument("--minimum-active-replicate-fraction", type=float, default=0.5)
    parser.add_argument("--log-ratio-pseudocount", type=float, default=0.1)
    parser.add_argument("--minimum-biological-replicates", type=int, default=2)
    parser.add_argument("--minimum-evaluable-tissues", type=int, default=2)
    parser.add_argument("--asymmetry-tissue-fraction", type=float, default=1.0 / 3.0)
    parser.add_argument("--silencing-tissue-fraction", type=float, default=0.8)
    parser.add_argument("--subfunctionalization-dominance-balance", type=float, default=0.5)
    parser.add_argument("--neofunctionalization-breadth-ratio", type=float, default=0.5)
    parser.add_argument("--outgroup-gain-fold-change", type=float, default=4.0)


def _configuration(args: argparse.Namespace) -> IDGRMConfig:
    return IDGRMConfig(
        log2fc_threshold=args.absolute_log2_fold_change,
        alpha=args.significance_threshold,
        p_adjust=args.multiple_testing_correction,
        min_expression=args.expression_threshold,
        ancestor_min_expression=args.outgroup_expression_threshold,
        min_active_fraction=args.minimum_active_replicate_fraction,
        pseudocount=args.log_ratio_pseudocount,
        min_replicates=args.minimum_biological_replicates,
        min_evaluable_tissues=args.minimum_evaluable_tissues,
        asymmetry_tissue_fraction=args.asymmetry_tissue_fraction,
        silencing_tissue_fraction=args.silencing_tissue_fraction,
        sub_dominance_balance_min=args.subfunctionalization_dominance_balance,
        neo_breadth_ratio_max=args.neofunctionalization_breadth_ratio,
        gain_fold=args.outgroup_gain_fold_change,
    ).validate()


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="idgrm",
        description="Hierarchical classification of duplicate-gene expression evolution.",
    )
    parser.add_argument("--version", action="version", version=f"iDGRM {__version__}")
    commands = parser.add_subparsers(dest="command", required=True)

    primary = commands.add_parser(
        "classify-primary",
        help="Assign four primary expression-pattern classes to duplicate-gene pairs.",
    )
    source = primary.add_mutually_exclusive_group(required=True)
    source.add_argument("--expression-matrix")
    source.add_argument("--differential-expression-dir")
    primary.add_argument("--sample-metadata")
    primary.add_argument("--duplicate-pairs", required=True)
    primary.add_argument("--output-dir", required=True)
    primary.add_argument("--normalization", choices=["none", "cpm", "median-ratio"], default="none")
    primary.add_argument("--differential-expression-file-pattern", default="*.DESeq2.csv")
    primary.add_argument(
        "--tissue-name-regex", default=DEFAULT_DIFFERENTIAL_EXPRESSION_FILENAME_REGEX
    )
    _add_configuration_arguments(primary)

    refine = commands.add_parser(
        "refine-subtypes",
        help="Refine eligible primary classes without altering the primary assignment.",
    )
    refine.add_argument("--primary-classifications", required=True)
    refine.add_argument("--tissue-evidence", required=True)
    refine.add_argument("--duplicate-pairs", required=True)
    refine.add_argument("--outgroup-expression-matrix")
    refine.add_argument("--outgroup-sample-metadata")
    refine.add_argument("--output-dir", required=True)
    refine.add_argument("--normalization", choices=["none", "cpm", "median-ratio"], default="none")
    _add_configuration_arguments(refine)
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
        config = _configuration(args)
        if args.command == "classify-primary":
            if args.expression_matrix:
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
                result = run_primary_classification_from_differential_expression(
                    differential_expression_dir=args.differential_expression_dir,
                    duplicate_pairs_path=args.duplicate_pairs,
                    output_dir=args.output_dir,
                    config=config,
                    file_pattern=args.differential_expression_file_pattern,
                    tissue_name_regex=args.tissue_name_regex,
                )
        else:
            result = run_subtype_refinement(
                primary_classifications_path=args.primary_classifications,
                tissue_evidence_path=args.tissue_evidence,
                duplicate_pairs_path=args.duplicate_pairs,
                outgroup_expression_matrix_path=args.outgroup_expression_matrix,
                outgroup_sample_metadata_path=args.outgroup_sample_metadata,
                output_dir=args.output_dir,
                config=config,
                normalization=args.normalization,
            )
        _print_result(result)
        return 0
    except Exception as exc:
        print(f"iDGRM error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
