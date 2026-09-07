# 原始脚本审计与 iDGRM 重构对应关系

## 审计范围

重点阅读了两个输入文件夹中的 `IDGRM-v2.pl`、`IDGRM-v3.pl`、`IDGRM-v3-new.pl`、拟南芥及 NG 版本、DESeq2 R 脚本、count/TPM 合并脚本、运行包装脚本和已有 `*.pairs.RM.txt` / `*.count.txt` 结果。

## 被保留的思想

- 以组织内 copy1/copy2 差异表达为最小证据单元。
- `|log2FC| > 1` 与严格显著性阈值。
- reciprocal 方向作为 sub/neo 候选。
- 单向、跨多个组织的稳定优势作为 AED。
- 显著占优组织数优先、总体表达用于平局时确定 major/minor。
- 保留 WGD、tandem、proximal、transposed、dispersed 等重复类型元数据。

## 旧实现中需要修正的地方

1. `int(tissues/3)` 向下取整；当组织数为 2 时阈值为 0，5 个组织时仍为 1。iDGRM 使用 `ceil(n_expressed/3)`。
2. 论文的 AED 分母是任一 copy 有表达的组织；旧脚本主要使用全部文件数。iDGRM 显式计算可评价且任一 copy 开启的组织。
3. v2 在显著方向次数相等时可把两个方向均为 0 的基因对误归入 sub/neo；v3-new 增加了非零条件。iDGRM 直接要求两个方向各至少 1 次。
4. 旧脚本把文件序号与 pair key 直接拼接，理论上可能碰撞，也丢失可读组织名。iDGRM 使用显式 `pair_id + tissue` 长表。
5. `split /-/` 不能安全处理本身含连字符的基因 ID。iDGRM 通过 pair 表精确映射，不拆分 ID 字符串。
6. 手工 `split /,/` 不能完整处理带引号/逗号的 CSV。iDGRM 使用标准 CSV 解析。
7. 旧流程固定支持 2 或 3 个重复、硬编码 TPMCalculator 列号和文件命名。iDGRM 支持任意重复数、显式样本元数据和可配置 legacy 文件正则。
8. `deseq2.R` 含固定工作目录，且设计式只有 copy condition，没有显式个体配对项。iDGRM 的归一化表达入口按同一样本对两个 copy 做 paired t-test；原 DESeq2 结果仍可通过 legacy 入口复用。
9. 旧输出没有保存实际组织证据、参数、软件版本和输入哈希。iDGRM 输出组织长表与 `run_metadata.json`。
10. 原方法无法区分 sub/neo。iDGRM 增加祖先参照规则和明确标注的 expression-only proxy，并保留歧义类别。
11. 原目录同时存在标准表头、`Transposed/Parental` 表头、无表头以及“数据列比表头多一列”的 pair 文件。iDGRM 对这些格式显式识别，避免基因 ID 被误当成行索引；完全重复的 pair 行会去重并写入运行警告。

## 兼容性说明

iDGRM 的 `science_class` 用于与旧结果比较，但默认又增加了三个提高可靠性的选择：BH 校正、表达 on/off 过滤、严格的可评价组织分母。因此结果不保证逐行等同于旧 Perl 输出。若需要最大程度接近论文原文，可使用 `--p-adjust none`；若导入旧 DESeq2 CSV，则直接使用其中的 `padj`。
