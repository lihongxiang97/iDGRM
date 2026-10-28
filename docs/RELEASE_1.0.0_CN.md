# iDGRM 1.0.0 首个正式版本

## 功能

- 基础四分类：Unmappable、No difference、Asymmetrically expressed、Sub-/neofunctionalized。
- 独立可选外群表达极化：SUB 候选、NEO 候选、无法极化，不覆盖基础分类。用户选择外群并提供正交关系、组织对应及检测条件声明。
- 原始整数 counts 的 R/DESeq2 拟合与逐组织证据导入；默认 BH-adjusted P ≤0.05、绝对 log2FC ≥1。这是对参考论文统计步骤的适配，不宣称原样复现。
- 输入检查、重复证据拒绝、缺失状态保留、输出覆盖保护和 SHA-256 追溯。

## 版本与边界

本版将本地 0.4.0 候选版正式编号为 1.0.0，不改变其分类规则。旧 refine-subtypes、网页界面和表达-only 的 SUB/NEO 猜测不属于受支持接口。当前公开命令为 fit-copy-counts、classify-primary、polarize-outgroup。

正式版本号表示公开接口的维护承诺，不意味着生物学准确率已完全验证。外群只是祖先状态参照；候选不是功能变化的实验证明。逐组织 BH 不保证最终类别 FDR；少重复、高变异及 copy 归属偏倚仍有局限。

## 验收

0.4.0 候选版已完成 Windows 75 项工作区测试、Linux 随包 57 项测试、安装包双模块示例及 R/DESeq2 Wald/LRT 集成检查。1.0.0 重新检查版本、随包测试、干净构建和安装包示例；示例与合成数据仅验证工程行为，不用于声称生物学准确率。

使用前阅读 README.md、METHODOLOGY_CN.md 和 OUTGROUP_GUIDE_CN.md。此次不重新计算或重标历史生物学结果。
