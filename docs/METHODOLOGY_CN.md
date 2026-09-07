# iDGRM 方法说明

## 1. 论文方法的可计算部分

来源：Xun Lan 与 Jonathan K. Pritchard，*Coregulation of tandem duplicate genes slows evolution of subfunctionalization in mammals*，Science 352:1009-1013 (2016)，DOI: [10.1126/science.aad8411](https://doi.org/10.1126/science.aad8411)。本地输入论文为 `lan2016.pdf`。

论文的分析链为：

1. 去除已注释假基因后，以 reciprocal best hit 构建 1,444 对高置信重复基因；要求编码序列可比对比例大于 80%，平均序列一致性大于 50%。
2. 以同义替换率 dS 近似重复时间，并用系统发育分布进行交叉校正。论文重点分析 dS < 0.7 的胎盘哺乳动物时期重复基因。
3. 为避免相似 copy 之间的 RNA-seq 错配，只使用两个 copy 都可唯一、无偏归属的同源位置估计表达比；代价是部分极年轻重复基因不可评价。
4. 使用 GTEx 的 46 个成人组织，每组织 10 个个体；并以 26 个小鼠组织复现主要结论。
5. 每对基因中，总体表达更高者为 major，另一者为 minor。
6. 若两个 copy 各自在至少一个组织显著高于对方，定义为 potential sub-/neofunctionalization：至少 2 倍差异、paired t-test P < 0.001。
7. 若 major 在“任一 copy 有表达”的组织中至少 1/3 显著更高，且没有任何组织发生显著反向，则定义为 AED（asymmetrically expressed duplicate）。
8. 其余为 no difference。论文另外通过外群单拷贝正交基因的总表达，检验年轻重复基因的 dosage sharing，而不是把 dosage sharing 直接作为上述三分类的一类。

## 2. 原脚本的实现

原始 Perl/R 流程大体包括：HISAT2 比对；TPMCalculator 计算 unique exon reads/TPM；按组织合并 2-3 个生物学重复；为每种重复模式提取基因对 count 矩阵；逐组织运行 DESeq2；最后由 `IDGRM-v3-new.pl` 汇总 `log2FoldChange` 与 `padj`。

脚本采用 `|log2FC| > 1`、`padj < 0.001` 识别组织内显著方向，并以显著方向出现次数确定 major/minor。一个 copy 至少一次显著占优且另一个 copy 也至少一次显著占优时，输出合并类别 `Sub-/neofunctionalized`。单向占优达到约 1/3 组织且无反向趋势时输出 AED，其余输出 No difference。

## 3. iDGRM 的组织级统计证据

对于带重复样本的归一化表达矩阵，iDGRM 在每个组织内对两个 copy 的同一样本表达做配对检验。默认在 `log2(expression + 0.1)` 上进行 paired t-test。每个组织内跨全部基因对进行 Benjamini-Hochberg 校正：

\[
\operatorname{log2FC}_{2/1,t}=\log_2\frac{\bar{x}_{2,t}+c}{\bar{x}_{1,t}+c}
\]

其中默认 `c=0.1`。显著方向要求：

\[
|\operatorname{log2FC}|\ge 1,\quad q\le 0.001.
\]

一个 copy 在组织内被判为 expressed，需要组织中位数达到 `min_expression`，且至少 `min_active_fraction` 的重复达到该阈值。AED 的分母严格使用“任一 copy 有表达、fold-change 有限且达到最少重复数，因而可实际检验”的组织数，阈值为 `ceil(n/3)`；重复不足的组织会保留在组织证据表中，但不会稀释分类分母。

如果输入列已经是组织均值，iDGRM 使用 effect-only 模式：满足表达 on/off 和 fold-change 即给出方向，不产生 P/q 值。该模式适合探索，不宜替代正式差异表达统计。

## 4. sub/neo 拆分规则

### 4.1 有外群单拷贝正交基因：推荐模式

首先必须满足 Science 的 reciprocal 条件，即 gene1-high 和 gene2-high 各至少出现一次。

亚功能化要求同时满足：

- 两个 copy 分别在至少一个“祖先有表达”组织占优；
- 两个 copy 的表达并集覆盖至少 80% 的祖先表达组织；
- 两个方向的组织数量较平衡，`min(n1,n2)/max(n1,n2) >= 0.5`；
- 没有 copy 在祖先未表达组织形成可靠新表达增益。

新功能化要求：

- 恰有一个 copy 在至少一个祖先未表达组织占优且表达开启；
- 该 copy 相对祖先表达至少增加 4 倍；
- 另一个 copy 在祖先表达组织中仍至少一次占优，并保留至少 60% 的祖先表达广度；
- 祖先相似度或覆盖度支持“一个 copy 保留、另一个 copy 创新”的方向。

若两个 copy 都出现祖先未表达组织增益，或祖先覆盖/相似度不足，则输出 `AMBIGUOUS_SUB_NEO`。

### 4.2 无外群：expression-only proxy

无外群时不能知道某个组织表达域是祖先保留还是后生获得，因此 iDGRM 明确使用代理规则：

- sub proxy：两个方向均存在，方向数量平衡度至少 0.5，表达广度平衡度至少 0.5，并且满足至少一种互补证据——两个 copy 各有 on/off 特异组织、跨组织表达相关不高于 0，或共同开启组织比例不高于 0.5。
- neo proxy：一个 copy 的表达广度至少为 2/3；另一个 copy 的广度不到前者的 1/2；窄谱 copy 至少有一个“自身表达而伙伴不表达”的占优组织；窄谱占优组织不超过全部可评价组织的 1/3；广谱 copy 在其他组织至少一次占优。
- 其他 reciprocal 模式保留为歧义型。

因此，expression-only 的 `NEOFUNCTIONALIZATION` 是“新表达域代理”，并不是功能获得的最终证明。

## 5. 其他扩展命运

- `EXPRESSION_LOSS`：某 copy 在至少 80% 可评价组织中关闭而伙伴开启，且从未显著占优。它表示表达层面的非功能化倾向，不等同于已形成假基因。
- `AED`：一侧显著占优达到 `ceil(n × 1/3)`，且无显著反向。
- `NO_DIFFERENCE`：未满足以上规则。
- `INSUFFICIENT_DATA`：可评价/可检验组织少于默认 2 个。

## 6. 解释与验证建议

1. 植物组织应尽量进行发育阶段和器官同源匹配；不对应的组织不能作为祖先 gain/loss 证据。
2. 跨物种 TPM 绝对量并不天然可比。祖先模式首先依赖表达/不表达状态，其次才使用 fold 和 profile similarity。
3. 年轻、高相似度重复基因必须处理多重比对；否则 minor copy 的低表达可能只是 read assignment 偏差。
4. neo 最好联合 parent/daughter 极性、多个外群、选择压力、结构域、互作和表型证据。
5. 阈值应做灵敏度分析；`run_metadata.json` 保存每次参数与输入哈希，便于完全复现。
