# iDGRM 2.0.0 方法及适用范围

## 基础分类与原文关系

参考 Lan & Pritchard (2016), DOI 10.1126/science.aad8411。保留四类表达模式框架，不把候选模式表述为已证实的生物学命运。

| 环节 | 参考论文正文 | iDGRM 2.0.0 |
|---|---|---|
| 差异检验 | 配对 t 检验 | 推荐整数 counts 的配对 DESeq2；另保留 log-expression 配对 t 检验 |
| 显著性 | P <0.001 | 默认逐组织 BH-adjusted P ≤0.05，可配置 |
| 效应值 | 至少两倍差异 | 绝对 log2FC ≥1 |
| SUB/NEO | 两 copy 各至少一个组织显著占优 | 相同方向反转框架，保持合并 |
| AED | 主要 copy 在至少三分之一任一 copy 表达的组织占优，无显著反向 | 可评价活跃组织中至少 ceil(N/3) 同向占优，无显著反向 |
| 不可评价 | 不可比对重复对 | 保留兼容代码 UNMAPPED，但另给 assessment_status、mapping_status 和原因，明确包含表达证据不足 |

这些调整具有明确统计和数据处理目的，但不是原文严格复现。修改显著性参数不能使所有其他步骤自动等同于原文。

## counts 模型

每组织分别拟合 `~biological_library + copy`。同一文库的两个 copy 共享由全基因 counts 的 poscounts 方法估计的 size factor。默认 Wald，也支持 LRT（约简模型 `~biological_library`）；总 counts <10 不检验，禁用独立过滤和自动 outlier replacement，保留 NA。R sessionInfo 记录真实运行版本。批次、处理或时间点不能不加区别地当作同一组织的生物学重复。

DESeq2 负二项模型及 BH 校正可以用于此适配分析，但其合理性依赖归属、有效长度、文库归一化和模型拟合假设。相互重叠的基因对及小样本高变异可能影响统计性质；不保证任何数据集都控制在 5%，更不保证最终类别的 FDR。BH 检验 H0: log2FC=0 后叠加效应筛选，不等同于直接检验 |log2FC|>1。

完全匹配且唯一归属是必要的输入质量要求之一，不自动保证两个 copy 的可观测区域对等。应使用共同同源且两侧可可靠归属的位置，或有依据的有效长度处理；软件声明字段不能替代实际定量处理。

## 可评价性及分类顺序

默认 TPM 中位数 ≥1 且至少一半重复达到阈值，判为表达活跃。任一 copy 活跃、效应有限、所选 P/q 为有效概率的组织才进入可评价分母。默认至少两个可评价组织。缺失 P/q 不填为 1，不因有原始 P 就替换已提供但缺失的 padj。

1. 证据不足：UNMAPPED；这是兼容性桶，不证明真实不可比对。mapping_status 明确标为 NOT_INFERRED_FROM_EXPRESSION。
2. 两方向各至少一次显著占优：SUB_OR_NEOFUNCTIONALIZED。
3. 一个方向在至少 ceil(N × asymmetry_tissue_fraction) 个组织占优且没有反向：ASYMMETRICALLY_EXPRESSED。
4. 其余：NO_DIFFERENCE。它表示未达到以上规则，不是证明表达完全相等或功能相同。

输出 active_tissue_test_coverage，提示有表达组织中多少能完成检验；低覆盖需谨慎，不能把缺少检验当作没有反向差异。

## 功能边界

本版只输出四类表达模式。SUB 与 NEO 保持合并，不提供外群极化、表达代理细分或祖先状态推断。2.0.0 删除了相关命令、Python 接口、内部历史细分实现和参数。四分类的 counts 模型、默认阈值和判定顺序不变。

## 科学依据

- 原论文：https://pmc.ncbi.nlm.nih.gov/articles/PMC5182070/
- DESeq2 官方手册（配对设计、NA、size factors、P 校正）：https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

本工具借鉴原文框架并公开统计适配，不声称新建生物学真值。工程验收与生物学验证分开报告。
