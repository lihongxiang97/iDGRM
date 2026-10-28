# iDGRM 1.0.0

重复基因跨组织表达模式分类工具，包含两个独立模块：基础四分类与可选外群表达极化。Python ≥3.10；counts 拟合另需 R/DESeq2。

## 功能边界

基础四类：**Unmappable、No difference、Asymmetrically expressed、Sub-/neofunctionalized**。

可选模块只进一步分析最后一类，输出 **OUTGROUP_SUPPORTED_SUB_CANDIDATE、OUTGROUP_SUPPORTED_NEO_CANDIDATE、UNRESOLVED**；其他基础类别标为 NOT_ELIGIBLE。候选结果独立保存，不修改基础分类。无外群不做纯表达 SUB/NEO 猜测，不输出正确概率，不声称重建真实祖先或证明生化功能变化。

四分类参考 Lan & Pritchard (2016)，但默认采用 **BH-adjusted P ≤0.05、|log2FC| ≥1**。这是公开说明的统计调整，不是原论文配对 t 检验 / P <0.001 的原样复现。逐组织 BH 不能解释为最终类别的 FDR 保证。

## 安装与检查

```bash
python -m pip install .
idgrm --version
idgrm --help
python -m unittest discover -s tests -v
```

## 模块一：基础四分类

先从可比对区域的原始整数 counts 做差异检验（可使用其他软件产生的兼容结果）：

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

counts 矩阵首列为基因 ID；两张 copy 矩阵首列为 pair_id，其他列为 sample_id，顺序必须与 samples.tsv 一致。样本表包含 sample_id、tissue、replicate；pair 表包含 pair_id、gene1、gene2。每组织至少两份独立生物学文库，建议三份或以上；技术重复应先合并。TPM 仅用于活性判断与展示，不进入 DESeq2 模型。

`--observability-method` 为用户提供的计数来源声明，不会自动执行 read-swap 或验证可比对性。两 copy 必须具备可比的有效长度/可观测区域，不能仅保证单侧唯一比对就认为已经消除长度和归属偏倚。`uncorrected` 只供诊断，不建议用于正式生物学解释。

也可使用独立的 log-expression 配对 t 检验路径；它不是 DESeq2 的替代实现：

```bash
idgrm classify-primary --expression-matrix examples/expression.tsv \
  --sample-metadata examples/samples.tsv --duplicate-pairs examples/pairs.tsv \
  --output-dir results/example_primary
```

输出 primary_classifications.tsv、tissue_evidence.tsv、primary_summary.tsv、primary_run_metadata.json。保留机器代码 UNMAPPED 兼容既有结果，其显示名称为 Unmappable，但证据不足不等于证明不可比对，详见方法说明。

## 模块二：外群表达极化

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

示例仅验证程序行为，不是性能基准或真实生物学结论。兼容性历史名称 ancestor_expression.tsv 实际表示外群，不是祖先重建。

用户选择外群并提供正交关系、组织映射及检测条件声明。该模块不会自动认证这些生物学判断。请先阅读 [外群选择与输入规范](docs/OUTGROUP_GUIDE_CN.md)。默认只接受明确对应组织，近似对应须显式启用 `--allow-approximate-tissues`。

## 方法、版本与追溯

- [当前方法与原文差异](docs/METHODOLOGY_CN.md)
- [外群指南](docs/OUTGROUP_GUIDE_CN.md)
- [1.0.0 变更与发布验收](docs/RELEASE_1.0.0_CN.md)

输出带参数与输入 SHA-256；基础输出带完整性哈希供极化模块核对。重复运行不得覆盖同名结果文件。历史细分代码只为旧研究回归保留在内部模块，不是受支持接口；新的公开入口是 polarize-outgroup。旧结果不自动重标。

本版“稳定”指受支持接口、规则执行、异常处理与可重复性，不表示所有条件下都有已验证的生物学准确率或错误率。小样本、高离散度、组织覆盖不足和相似 copy 的归属偏倚仍需谨慎。

## 引用

Lan X, Pritchard JK. Coregulation of tandem duplicate genes slows evolution of subfunctionalization in mammals. Science 352, 1009–1013 (2016). DOI: [10.1126/science.aad8411](https://doi.org/10.1126/science.aad8411)。

DESeq2 方法参考其[官方说明](https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)。使用本工具不能替代对输入实验设计与假设的检查。
