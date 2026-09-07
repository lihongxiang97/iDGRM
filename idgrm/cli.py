"""Command-line interface for iDGRM."""

from __future__ import annotations

import argparse
import sys
from contextlib import ExitStack
from importlib.resources import as_file, files
from pathlib import Path

from . import __version__
from .config import IDGRMConfig
from .legacy import DEFAULT_FILENAME_REGEX
from .pipeline import run_analysis, run_legacy_analysis


def _add_config_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--log2fc", type=float, default=1.0, help="Absolute log2 fold-change cutoff (default: 1)")
    parser.add_argument("--alpha", type=float, default=0.001, help="P/q-value cutoff (default: 0.001)")
    parser.add_argument("--p-adjust", choices=["bh", "none"], default="bh", help="Multiple-testing policy")
    parser.add_argument("--min-expression", type=float, default=1.0, help="Expression on/off threshold")
    parser.add_argument("--ancestor-min-expression", type=float, default=None, help="Optional separate outgroup on/off threshold")
    parser.add_argument("--min-active-fraction", type=float, default=0.5, help="Minimum active-replicate fraction")
    parser.add_argument("--pseudocount", type=float, default=0.1, help="Pseudocount before log ratios")
    parser.add_argument("--min-replicates", type=int, default=2, help="Minimum paired replicates per tissue")
    parser.add_argument("--min-tissues", type=int, default=2, help="Minimum evaluable tissues per pair")
    parser.add_argument("--aed-fraction", type=float, default=1.0 / 3.0, help="AED dominance fraction")
    parser.add_argument("--loss-fraction", type=float, default=0.8, help="Expression-loss tissue fraction")
    parser.add_argument("--no-expression-loss", action="store_true", help="Do not emit EXPRESSION_LOSS")
    parser.add_argument("--sub-balance", type=float, default=0.5, help="Minimum reciprocal dominance balance")
    parser.add_argument("--neo-breadth-ratio", type=float, default=0.5, help="Maximum narrow/broad breadth ratio for neo proxy")
    parser.add_argument("--gain-fold", type=float, default=4.0, help="Minimum gain over an inactive ancestor")


def _config_from_args(args: argparse.Namespace) -> IDGRMConfig:
    return IDGRMConfig(
        log2fc_threshold=args.log2fc,
        alpha=args.alpha,
        p_adjust=args.p_adjust,
        min_expression=args.min_expression,
        ancestor_min_expression=args.ancestor_min_expression,
        min_active_fraction=args.min_active_fraction,
        pseudocount=args.pseudocount,
        min_replicates=args.min_replicates,
        min_evaluable_tissues=args.min_tissues,
        aed_fraction=args.aed_fraction,
        loss_fraction=args.loss_fraction,
        detect_expression_loss=not args.no_expression_loss,
        sub_dominance_balance_min=args.sub_balance,
        neo_breadth_ratio_max=args.neo_breadth_ratio,
        gain_fold=args.gain_fold,
    ).validate()


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="idgrm",
        description="Classify duplicated-gene expression fates across tissues.",
    )
    parser.add_argument("--version", action="version", version=f"iDGRM {__version__}")
    subparsers = parser.add_subparsers(dest="command", required=True)

    classify = subparsers.add_parser("classify", help="Analyze an expression matrix")
    classify.add_argument("--expression", required=True, help="Genes-by-samples/tissues CSV or TSV")
    classify.add_argument("--pairs", required=True, help="Duplicate-pair CSV or TSV")
    classify.add_argument("--samples", help="Sample metadata with sample_id and tissue")
    classify.add_argument("--ancestor-expression", help="Optional outgroup ortholog expression matrix")
    classify.add_argument("--ancestor-samples", help="Optional outgroup sample metadata")
    classify.add_argument("--output", required=True, help="Output directory")
    classify.add_argument(
        "--normalization",
        choices=["none", "cpm", "median-ratio"],
        default="none",
        help="Column normalization before testing",
    )
    _add_config_arguments(classify)

    legacy = subparsers.add_parser("legacy", help="Analyze original per-tissue DESeq2 CSV files")
    legacy.add_argument("--deseq-dir", required=True, help="Directory containing *.DESeq2.csv")
    legacy.add_argument("--pairs", required=True, help="Duplicate-pair CSV or TSV")
    legacy.add_argument("--output", required=True, help="Output directory")
    legacy.add_argument("--pattern", default="*.DESeq2.csv", help="File glob inside --deseq-dir")
    legacy.add_argument(
        "--filename-regex",
        default=DEFAULT_FILENAME_REGEX,
        help="Regex with a named 'tissue' capture group",
    )
    legacy.add_argument("--ancestor-expression", help="Optional outgroup ortholog expression matrix")
    legacy.add_argument("--ancestor-samples", help="Optional outgroup sample metadata")
    legacy.add_argument(
        "--normalization",
        choices=["none", "cpm", "median-ratio"],
        default="none",
    )
    _add_config_arguments(legacy)

    demo = subparsers.add_parser("demo", help="Run the bundled synthetic example")
    demo.add_argument("--output", default="results/demo", help="Output directory")
    demo.add_argument("--without-ancestor", action="store_true", help="Use pair-only proxy rules")
    _add_config_arguments(demo)
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
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        config = _config_from_args(args)
        if args.command == "classify":
            result = run_analysis(
                expression_path=args.expression,
                pairs_path=args.pairs,
                samples_path=args.samples,
                ancestor_expression_path=args.ancestor_expression,
                ancestor_samples_path=args.ancestor_samples,
                output_dir=args.output,
                config=config,
                normalization=args.normalization,
            )
        elif args.command == "legacy":
            result = run_legacy_analysis(
                deseq_directory=args.deseq_dir,
                pairs_path=args.pairs,
                output_dir=args.output,
                config=config,
                pattern=args.pattern,
                filename_regex=args.filename_regex,
                ancestor_expression_path=args.ancestor_expression,
                ancestor_samples_path=args.ancestor_samples,
                normalization=args.normalization,
            )
        else:
            data = files("idgrm.data")
            with ExitStack() as stack:
                example_paths = {
                    name: stack.enter_context(as_file(data / name))
                    for name in (
                        "expression.tsv",
                        "pairs.tsv",
                        "samples.tsv",
                        "ancestor_expression.tsv",
                        "ancestor_samples.tsv",
                    )
                }
                result = run_analysis(
                    expression_path=example_paths["expression.tsv"],
                    pairs_path=example_paths["pairs.tsv"],
                    samples_path=example_paths["samples.tsv"],
                    ancestor_expression_path=(
                        None
                        if args.without_ancestor
                        else example_paths["ancestor_expression.tsv"]
                    ),
                    ancestor_samples_path=(
                        None
                        if args.without_ancestor
                        else example_paths["ancestor_samples.tsv"]
                    ),
                    output_dir=args.output,
                    config=config,
                )
        _print_result(result)
        return 0
    except Exception as exc:
        print(f"iDGRM error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
