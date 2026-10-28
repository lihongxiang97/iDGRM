"""Release integration checks; synthetic fixtures are NOT accuracy benchmarks.

Run with an installed wheel from outside the source checkout.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    parser.add_argument("--examples", required=True)
    parser.add_argument("--rscript")
    args = parser.parse_args()
    output = Path(args.output).resolve()
    output.mkdir(parents=True, exist_ok=False)
    example = Path(args.examples).resolve()
    report = {"purpose": "Integration smoke test, not biological accuracy or error-control benchmark", "checks": []}

    def command(*arguments):
        result = subprocess.run([sys.executable, "-m", "idgrm", *map(str, arguments)],
            capture_output=True, text=True)
        report["checks"].append(dict(command=list(map(str, arguments)), returncode=result.returncode,
            stdout=result.stdout, stderr=result.stderr))
        (output / "integration_report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
        if result.returncode:
            raise RuntimeError(result.stderr)

    command("--version")
    command("classify-primary", "--expression-matrix", example / "expression.tsv",
        "--sample-metadata", example / "samples.tsv", "--duplicate-pairs", example / "pairs.tsv",
        "--output-dir", output / "primary")
    labels = pd.read_csv(output / "primary/primary_classifications.tsv", sep="\t").set_index("pair_id").primary_class
    assert labels["pair_sub"] == "SUB_OR_NEOFUNCTIONALIZED"
    assert labels["pair_neo"] == "SUB_OR_NEOFUNCTIONALIZED"
    if args.rscript:
        fixture = output / "counts_fixture"
        fixture.mkdir()
        rng = np.random.default_rng(20261028)
        n_pairs, n_samples = 300, 6
        ids = [f"p{i}" for i in range(n_pairs)]
        samples = [f"s{i}" for i in range(n_samples)]
        pairs = pd.DataFrame(dict(pair_id=ids, gene1=[f"g{i}a" for i in range(n_pairs)], gene2=[f"g{i}b" for i in range(n_pairs)]))
        means = np.exp(rng.uniform(np.log(20), np.log(1500), n_pairs))
        dispersions = rng.uniform(.03, .25, n_pairs)
        one = rng.negative_binomial(1 / dispersions[:, None], 1 / (1 + means[:, None] * dispersions[:, None]), size=(n_pairs, n_samples))
        two_means = means[:, None] * np.ones((1, n_samples))
        two_means[:60, :3] *= 4
        two_means[:30, 3:] *= .2
        two = rng.negative_binomial(1 / dispersions[:, None], 1 / (1 + two_means * dispersions[:, None]))
        pd.DataFrame(np.vstack([one, two]), index=list(pairs.gene1) + list(pairs.gene2), columns=samples).to_csv(fixture / "counts.tsv", sep="\t", index_label="gene_id")
        whole = pd.DataFrame(np.vstack([one, two]), index=list(pairs.gene1) + list(pairs.gene2), columns=samples)
        # All synthetic genes have the same assumed effective length.
        whole.divide(whole.sum(axis=0), axis=1).multiply(1e6).to_csv(fixture / "tpm.tsv", sep="\t", index_label="gene_id")
        pairs.to_csv(fixture / "pairs.tsv", sep="\t", index=False)
        pd.DataFrame(dict(sample_id=samples, tissue=["root"] * 3 + ["leaf"] * 3, replicate=[1, 2, 3] * 2)).to_csv(fixture / "samples.tsv", sep="\t", index=False)
        for label, values in (("one", one), ("two", two)):
            pd.DataFrame(values, index=ids, columns=samples).to_csv(fixture / f"{label}.tsv", sep="\t", index_label="pair_id")
        # Controlled equal observability is a fixture property, not a real-data claim.
        for test in ("wald", "lrt"):
            command("fit-copy-counts", "--counts", fixture / "counts.tsv", "--duplicate-pairs", fixture / "pairs.tsv",
                "--sample-metadata", fixture / "samples.tsv", "--copy1-counts", fixture / "one.tsv",
                "--copy2-counts", fixture / "two.tsv", "--observability-method", "equal-effective-length",
                "--test", test, "--rscript", args.rscript, "--output-dir", output / test)
            table = pd.read_csv(output / test / "idgrm.counts.pairs.root.DESeq2.csv")
            assert table.pvalue.notna().sum() > 200
            assert table.log2FoldChange[:60].median() > 1
            assert (output / test / "COMPLETE").exists()
            command("classify-primary", "--differential-expression-dir", output / test,
                "--activity-expression-matrix", fixture / "tpm.tsv", "--activity-sample-metadata", fixture / "samples.tsv",
                "--duplicate-pairs", fixture / "pairs.tsv", "--output-dir", output / f"{test}_primary")
            # These are simulated expression states, not biological observations.
    import idgrm
    report.update(status="PASS", software_version=idgrm.__version__, imported_from=idgrm.__file__,
        fixture_note="Synthetic data test executable plumbing and sign conventions only")
    report["output_sha256"] = {str(p.relative_to(output)): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in output.rglob("*") if p.is_file() and p.name != "integration_report.json"}
    (output / "integration_report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps({"status": "PASS", "checks": len(report["checks"]), "output": str(output)}))


if __name__ == "__main__":
    main()
