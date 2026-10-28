# iDGRM 1.0.0

iDGRM classifies duplicate-gene expression patterns across tissues. It provides two independent modules: primary four-class classification and optional outgroup-reference expression polarization.

## Scope

The primary module assigns four classes: **Unmappable**, **No difference**, **Asymmetrically expressed**, and **Sub-/neofunctionalized**.

The optional module analyzes only primary Sub-/neofunctionalized pairs and reports `OUTGROUP_SUPPORTED_SUB_CANDIDATE`, `OUTGROUP_SUPPORTED_NEO_CANDIDATE`, or `UNRESOLVED`. Other primary classes are marked `NOT_ELIGIBLE`. These outputs do not modify primary assignments. The module does not infer separate SUB/NEO labels without an outgroup, report probabilities of correctness, reconstruct a true ancestor, or establish changes in biochemical function.

The four-class framework follows Lan and Pritchard (2016), with explicitly documented statistical adaptations. Defaults are **BH-adjusted P ≤ 0.05** and **absolute log2FC ≥ 1**, rather than the publication's paired t-test with P < 0.001. Within-tissue BH adjustment does not guarantee false discovery rate control for the final expression classes.

## Dependencies and installation

Python **3.10 or later** is required. The Python package depends on NumPy, pandas, and SciPy.

**The raw-count workflow requires R and DESeq2.** Installing iDGRM with pip does not install R or DESeq2; install these separately and make `Rscript` available on your PATH, or provide its location with `--rscript`.

| Command or input route | R/DESeq2 required locally? |
|---|---|
| `fit-copy-counts`: differential testing of raw integer counts | **Yes**; calls the actual R/DESeq2 engine |
| `classify-primary --differential-expression-dir`: import completed differential tests | No; compatible test results must already exist |
| `classify-primary --expression-matrix`: log-expression paired t-test | No; this is a distinct statistical route, not a Python implementation of DESeq2 |
| `polarize-outgroup`: optional outgroup-reference analysis | No; uses existing primary results and expression evidence |

Download the wheel or source archive from the [1.0.0 release](https://github.com/lihongxiang97/iDGRM/releases/tag/v1.0.0). From an extracted source directory:

```bash
python -m pip install .
idgrm --version
idgrm --help
python -m unittest discover -s tests -v
```

## Module 1: primary four-class classification

Fit differential expression using raw integer counts from comparable, reliably assignable regions. Alternatively, import compatible differential-expression results generated beforehand.

```bash
idgrm fit-copy-counts \
  --counts whole_gene_counts.tsv \
  --copy1-counts comparable_copy1_counts.tsv \
  --copy2-counts comparable_copy2_counts.tsv \
  --duplicate-pairs pairs.tsv --sample-metadata samples.tsv \
  --observability-method reciprocal-common-cds \
  --output-dir results/count_tests

idgrm classify-primary \
  --differential-expression-dir results/count_tests \
  --activity-expression-matrix tpm.tsv \
  --activity-sample-metadata samples.tsv \
  --duplicate-pairs pairs.tsv --output-dir results/primary
```

The first column of the whole-gene count matrix contains gene IDs. The first column of each copy-specific matrix contains `pair_id`; the remaining columns contain sample counts in the same order as `samples.tsv`. Sample metadata must include `sample_id`, `tissue`, and `replicate`; the pair table must include `pair_id`, `gene1`, and `gene2`. Each tissue requires at least two independent biological libraries; three or more are recommended. Merge technical replicates before analysis. TPM is used for activity assessment and visualization, not as input to the DESeq2 model.

`--observability-method` records a user-supplied description of count provenance. It does not perform read swapping or verify mappability. The two copies must have comparable effective lengths or observable regions; uniquely mapping reads from one copy alone do not eliminate length and assignment biases. The `uncorrected` setting is intended for diagnostics and is not recommended for biological interpretation.

An independent log-expression paired t-test route is also available. It is not a reimplementation of DESeq2:

```bash
idgrm classify-primary --expression-matrix examples/expression.tsv \
  --sample-metadata examples/samples.tsv --duplicate-pairs examples/pairs.tsv \
  --output-dir results/example_primary
```

Outputs are `primary_classifications.tsv`, `tissue_evidence.tsv`, `primary_summary.tsv`, and `primary_run_metadata.json`. The machine-readable code `UNMAPPED` is retained for compatibility and displayed as Unmappable. This category also includes insufficient expression evidence: it does not, by itself, establish that a pair is genuinely unmappable. Consult `assessment_status`, `mapping_status`, and the classification reason. No difference means that the specified directional rules were not met, not that the copies have identical functions.

## Module 2: optional outgroup-reference expression polarization

```bash
idgrm polarize-outgroup \
  --primary-results-dir results/example_primary \
  --outgroup-expression examples/ancestor_expression.tsv \
  --outgroup-samples examples/ancestor_samples.tsv \
  --orthologs examples/orthologs.tsv --tissue-map examples/tissue_map.tsv \
  --detection-qc examples/detection_qc.tsv \
  --outgroup-species "Synthetic example" --reference-version fixture-v1 \
  --orthology-source "Synthetic fixture" --output-dir results/example_polarization
```

These synthetic examples verify program behavior; they are not accuracy benchmarks or biological findings. The historical filename `ancestor_expression.tsv` refers to an extant outgroup, not a reconstructed ancestor.

Users select the outgroup and supply orthology assignments, tissue correspondences, and detection-adequacy declarations. The module does not automatically validate these biological judgments. Select an outgroup with an appropriate phylogenetic position relative to the duplication event and a reliably identified single-copy ortholog; the closest species is not necessarily the most suitable reference. Match tissue identity, developmental stage, and experimental conditions where possible.

Required mapping tables are illustrated in `examples/`:

| File | Required columns | Accepted declarations |
|---|---|---|
| `orthologs.tsv` | `pair_id`, `outgroup_gene_id`, `relationship` | `single_copy` or `ambiguous` |
| `tissue_map.tsv` | `target_tissue`, `outgroup_tissue`, `match_quality` | `exact` or `approximate` |
| `detection_qc.tsv` | `outgroup_gene_id`, `outgroup_tissue`, `assessment` | `adequate` or `inadequate` |

Tissue correspondence must be one-to-one. Approximate matches require explicit opt-in with `--allow-approximate-tissues`. Missing genes, missing tissues, failed quantification, and insufficient detection evidence must not be encoded as absence of expression.

Defaults distinguish ON (TPM ≥ 1), OFF (all measured replicates ≤ 0.1 TPM), and UNKNOWN states. At least two valid replicates, two informative matched tissues, and a 50% informative fraction among included tissue mappings are required. The retention threshold is 0.8. These are configurable operational thresholds, not calibrated biological constants. See the [detailed outgroup guide (Chinese)](docs/OUTGROUP_GUIDE_CN.md) for the full input and interpretation requirements.

## Methods, limitations, and reproducibility

- [Methods and differences from the reference publication (Chinese)](docs/METHODOLOGY_CN.md)
- [Outgroup selection and input guide (Chinese)](docs/OUTGROUP_GUIDE_CN.md)
- [1.0.0 release notes and verification scope (Chinese)](docs/RELEASE_1.0.0_CN.md)

Run metadata records parameters and input SHA-256 checksums. Primary outputs include integrity hashes checked by the polarization module. Existing output files are not silently overwritten. Historical refinement code is retained internally for regression checks and is not a supported public interface; use `polarize-outgroup`. Historical results are not automatically relabeled.

In this release, stability refers to supported interfaces, rule execution, error handling, and reproducibility. It does not imply established biological accuracy or error control under every experimental condition. Small sample sizes, high dispersion, limited tissue coverage, and assignment bias between similar copies remain important limitations. Outgroup-supported candidates are expression-based interpretations, not direct proof of subfunctionalization or neofunctionalization.

## References

Lan X, Pritchard JK. Coregulation of tandem duplicate genes slows evolution of subfunctionalization in mammals. Science 352, 1009–1013 (2016). DOI: [10.1126/science.aad8411](https://doi.org/10.1126/science.aad8411).

For DESeq2 methods and installation guidance, see the [official DESeq2 documentation](https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html). Using iDGRM does not replace assessment of experimental design, data quality, and statistical assumptions.
