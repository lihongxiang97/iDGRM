# iDGRM 2.0.0

This release focuses iDGRM on four-class duplicate-gene expression-pattern classification: Unmappable, No difference, Asymmetrically expressed, and Sub-/neofunctionalized.

## Breaking changes and migration

- Removed the optional outgroup-reference expression polarization command, Python API, implementation, example inputs and tests.
- Removed the older internal subtype-refinement implementation and its unused configuration thresholds.
- The supported commands are `fit-copy-counts` and `classify-primary`. SUB and NEO remain combined. Existing four-class command lines remain supported.
- Ancestral-reference fields in pair tables are no longer interpreted or exported. Extra input columns can remain in user files but have no role in classification.
- Version 2.0.0 is a major release because it removes previously public interfaces. Existing 1.0.0 releases and tags remain available for reproducibility; old analysis files are not relabeled.

## Retained methods

The count workflow calls actual R/DESeq2 with within-library copy pairing; compatible completed differential tests can also be imported. A separate log-expression paired t-test route remains available. Default classification uses BH-adjusted P ≤0.05, |log2 fold change| ≥1, activity filtering and at least two callable tissues. These rules and the primary count model are unchanged from 1.0.0.

## Validation scope

Release checks cover the supported four-class interfaces, rejected removed commands, input validation, missing evidence, copy exchange and tissue order, clean wheel contents and installed-package examples. Optional synthetic count integration exercises actual R/DESeq2 Wald and LRT routes. These checks establish executable behavior rather than biological accuracy or universal error-rate control. See the validation guide and release verification records for completed checks.

Pre-release checks passed 39 core tests on Windows and Linux, installed-wheel CLI examples, and six integration commands including actual R/DESeq2 Wald/LRT fitting on Linux. Replaying the saved final tissue evidence for pear hapA (30,512 pairs), Arabidopsis (14,478) and pear hapB (31,743) preserved every primary label (0 changed among 76,733 dataset-specific pair entries). HapA and hapB share expression libraries; these are not independent biological replicates. The replay is a regression check, not a biological accuracy benchmark.
