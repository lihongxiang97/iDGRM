# iDGRM 与原脚本的一致性说明

## 结论

iDGRM 与原 Perl/R 流程在核心分类思想上兼容，但默认结果不承诺逐行完全一致。`science_class` 用于复现原方法的三类框架；`extended_fate` 是在该框架上的扩展结果。

## 保持一致的部分

- 以组织内两个 copy 的表达差异为基本证据。
- 默认使用两倍差异阈值，即 `|log2FC| >= 1`。
- 两个方向各至少出现一次时记为 Science 的 `SUB_OR_NEO` 候选。
- 一个方向跨至少约三分之一可评价组织持续占优、且无反向组织时记为 AED。
- 其余可检验基因对记为 `NO_DIFFERENCE`。
- 可直接导入旧流程产生的逐组织 `*.DESeq2.csv`，沿用其中的 `log2FoldChange`、`pvalue` 和 `padj`。

## 有意不同的部分

1. iDGRM 使用 `ceil(n/3)` 计算 AED 阈值；旧脚本的 `int(n/3)` 会向下取整，并可能在组织很少时得到零阈值。
2. AED 分母只包含任一 copy 有表达且真正可检验的组织；重复不足的组织保留在证据表中，但不稀释分母。
3. iDGRM 明确要求 reciprocal 两个方向均至少出现一次，避免旧版脚本在两个方向都为零时误判。
4. 原始表达入口默认在每个组织内做配对检验，并进行组织内跨基因对的 BH 校正；旧脚本主要直接使用 DESeq2 的 `padj`。
5. iDGRM 增加表达 on/off 过滤、缺失数据检查、完全重复 pair 去重和输入格式修复。
6. 原脚本只输出合并的 sub/neo；iDGRM 进一步输出 `SUBFUNCTIONALIZATION`、`NEOFUNCTIONALIZATION` 或 `AMBIGUOUS_SUB_NEO`。
7. iDGRM 还增加 `EXPRESSION_LOSS` 与 `INSUFFICIENT_DATA`，因此扩展类别不能与旧三分类直接一一对应。

## 当前可验证程度

原始目录中检测到 8 个最终 `*.pairs.RM.txt` 文件，但没有保留下来可供重算的逐组织 `*.DESeq2.csv`。因此目前已完成的是规则级、输入格式级和合成数据回归验证，而不是对原研究全部基因对的逐行一致率统计。

原目录中的 18 个 `.pairs` 文件均已通过 iDGRM 读取验证，共覆盖 37,953 条去重前后相关记录；标准表头、`Parent/Daughter`、`Parental/Transposed`、无表头和不齐整表头格式均已测试。

若要计算精确一致率，应提供当时的逐组织 DESeq2 输出，或用于生成这些输出的完整 count/TPM 表与样本—组织对应关系。建议同时报告：

- `science_class` 的总体一致率和 Cohen's kappa；
- 每类 precision、recall 和混淆矩阵；
- 因 AED 分母修正、BH 校正、表达过滤和 sub/neo 拆分造成的差异数量。
