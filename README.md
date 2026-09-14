# iDGRM

iDGRM（Inference of Duplicated-Gene Retention Mechanisms）基于多组织表达证据，对重复基因对进行分层分类。软件将可复现的初级表达模式分类与研究性的进化命运亚型推断严格分离，避免将方法扩展误写为参考分类的一部分。

## 分类体系

### 初级表达模式分类

每个重复基因对只归入一个初级类别：

- `UNMAPPED`：缺乏足够的可唯一归属表达证据，无法可靠分类；
- `NO_DIFFERENCE`：未达到方向性表达差异规则；
- `ASYMMETRICALLY_EXPRESSED`：同一 copy 在规定比例的可评价组织中显著占优，且不存在显著反向；
- `SUB_OR_NEOFUNCTIONALIZED`：两个 copy 分别在至少一个组织中显著占优。

### 进化命运亚型细分

亚型细分只处理后两个方向性类别，并保留 `parent_primary_class`：

```text
ASYMMETRICALLY_EXPRESSED
├── CONSISTENT_ASYMMETRIC_EXPRESSION
└── EXPRESSION_SILENCING_ASSOCIATED_ASYMMETRY

SUB_OR_NEOFUNCTIONALIZED
├── SUBFUNCTIONALIZATION
├── NEOFUNCTIONALIZATION
└── UNRESOLVED_SUB_OR_NEOFUNCTIONALIZATION
```

`EXPRESSION_SILENCING_ASSOCIATED_ASYMMETRY` 是表达模式标签，不等同于假基因化或已证实的功能丢失。没有外群表达时，`NEOFUNCTIONALIZATION` 是未极化的表达模式推断；其 `evidence_basis` 会明确记录这一限制。

## 安装

```bash
python -m pip install -e .
```

浏览器界面：

```bash
python -m pip install -e ".[web]"
streamlit run app.py
```

## 初级分类

使用表达矩阵和生物学重复信息：

```bash
idgrm classify-primary \
  --expression-matrix expression.tsv \
  --sample-metadata samples.tsv \
  --duplicate-pairs pairs.tsv \
  --output-dir results/primary
```

使用已经完成的逐组织差异表达结果：

```bash
idgrm classify-primary \
  --differential-expression-dir tissue_differential_expression \
  --duplicate-pairs pairs.tsv \
  --output-dir results/primary
```

初级分类输出：

- `primary_classifications.tsv`；
- `primary_summary.tsv`；
- `tissue_evidence.tsv`；
- `primary_run_metadata.json`。

## 亚型细分

```bash
idgrm refine-subtypes \
  --primary-classifications results/primary/primary_classifications.tsv \
  --tissue-evidence results/primary/tissue_evidence.tsv \
  --duplicate-pairs pairs.tsv \
  --outgroup-expression-matrix outgroup_expression.tsv \
  --outgroup-sample-metadata outgroup_samples.tsv \
  --output-dir results/refined
```

亚型细分输出：

- `refined_classifications.tsv`；
- `refined_summary.tsv`；
- `refinement_run_metadata.json`。

## 输入要求

表达矩阵首列为基因ID，其余列为样本；样本信息至少包含 `sample_id` 和 `tissue`。重复基因对至少包含 `gene1` 和 `gene2`，建议提供稳定的 `pair_id`、重复类型、parent/daughter关系和外群单拷贝正交基因ID。

正式推断应使用能区分高度相似 copy 的定量结果。推荐仅计入完全匹配且可唯一归属一个 copy 的reads，并使用适用于整数计数的差异表达模型。TPM适合表达展示和表达阈值判断，不应直接替代count-based显著性模型。

默认方向性证据要求 `|log2FC| >= 1` 且Benjamini–Hochberg校正后P值不高于0.05。全部参数、输入文件哈希、软件版本和警告均写入运行元数据。

## 方法学边界

仅凭目标物种表达数据无法确定某个组织表达域是祖先保留还是后生获得。可靠区分亚功能化与新功能化应优先加入组织对应的外群单拷贝正交基因表达，并结合序列选择、蛋白结构、互作、表型或功能实验验证。

详细规则见 [方法说明](docs/METHODOLOGY_CN.md) 和 [验证说明](docs/VALIDATION_CN.md)。

## 测试

```bash
python -m unittest discover -s tests -v
```
