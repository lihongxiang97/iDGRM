# iDGRM 方法说明

## 1. 方法依据

初级分类规则依据 Lan 与 Pritchard（2016，Science 352:1009–1013，DOI: 10.1126/science.aad8411）描述的重复基因跨组织表达模式。期刊名称只用于文献溯源，不用于软件命令、参数、字段或文件名。

## 2. 严格定量要求

高度相似的重复基因容易发生read错误归属。正式分析应仅计入能够完全匹配且可唯一归属一个copy的reads，并保存比对规则、参考基因组版本、注释版本和计数方法。整数read counts应使用count-based差异表达模型；TPM只用于表达展示及表达状态辅助判定。

## 3. 组织级表达证据

每个组织记录两个copy的表达量、表达状态、log2 fold change、P值、校正后P值和显著方向。默认显著方向要求：

```text
absolute log2 fold change >= 1
Benjamini-Hochberg adjusted P value <= 0.05
```

若输入逐组织差异表达结果，iDGRM直接读取其效应值和显著性。若输入带生物学重复的表达矩阵，软件构建组织级证据；正式研究仍建议使用适用于整数计数的差异表达工具。

## 4. 初级表达模式分类

初级分类为四个互斥类别：

1. `UNMAPPED`：缺乏足够的可唯一归属表达证据，无法可靠分类；
2. `NO_DIFFERENCE`：未达到方向性表达差异规则；
3. `ASYMMETRICALLY_EXPRESSED`：同一copy在至少 `ceil(n_evaluable × asymmetry_tissue_fraction)` 个组织中显著占优，且不存在显著反向；
4. `SUB_OR_NEOFUNCTIONALIZED`：gene1-high和gene2-high各至少出现一次。

初级分类不包含任何iDGRM新增的一级命运类别。

## 5. 亚型细分

亚型细分读取 `primary_classifications.tsv` 和 `tissue_evidence.tsv`，只处理两个方向性类别，并保留 `parent_primary_class`。

### 5.1 非对称表达亚型

- `CONSISTENT_ASYMMETRIC_EXPRESSION`：一个copy呈稳定方向性优势，但伙伴copy不满足广泛表达沉默标准；
- `EXPRESSION_SILENCING_ASSOCIATED_ASYMMETRY`：一个copy在至少 `silencing_tissue_fraction` 的可评价组织中低于表达阈值，而伙伴copy处于表达状态，并且沉默copy从未显著占优。

后者描述表达模式，不证明基因功能丢失、假基因化或非功能化。

### 5.2 亚功能化与新功能化亚型

有外群单拷贝正交基因表达时：

- `SUBFUNCTIONALIZATION`：两个copy在不同祖先表达组织中分别占优，二者表达并集覆盖规定比例的祖先表达域，且没有可靠的新表达域增益；
- `NEOFUNCTIONALIZATION`：一个copy在祖先未表达组织获得显著表达域，而另一个copy保留足够的祖先表达域；
- `UNRESOLVED_SUB_OR_NEOFUNCTIONALIZATION`：现有证据不能唯一支持上述任一机制。

没有外群表达时，软件依据表达广度、组织互补性、方向平衡和表达相关性进行未极化推断，并将 `evidence_basis` 记录为 `TARGET_SPECIES_EXPRESSION_ONLY`。此时的新功能化判断不是祖先状态支持的最终结论。

## 6. 输出追溯性

初级分类和亚型细分分别保存输入文件SHA-256、参数、软件版本、警告及生成时间。亚型结果中的每一行都保留初级分类，因此不会因扩展分析改变参考分类结果。

## 7. 解释限制

新功能化的可靠证明需要结合外群表达、编码序列选择、蛋白结构或结构域变化、互作网络、表型和功能实验。表达沉默也需要结合更多组织、发育阶段、环境处理、基因组注释和可唯一比对长度进行验证。
