# iDGRM

iDGRM（Inference of Duplicated-Gene Retention Mechanisms）根据多组织转录组表达数据，对重复基因对的表达进化命运进行可复现分类。软件保留 Lan 与 Pritchard 2016 年 Science 论文的判定框架，同时将论文合并的 sub-/neofunctionalization 候选进一步拆分为亚功能化（sub）与新功能化（neo）。

## 核心输出

iDGRM 同时报告两套结果：

- `science_class`：`SUB_OR_NEO`、`AED`、`NO_DIFFERENCE`、`INSUFFICIENT_DATA`，用于与原论文及旧版 Perl 脚本对照。
- `extended_fate`：`SUBFUNCTIONALIZATION`、`NEOFUNCTIONALIZATION`、`AMBIGUOUS_SUB_NEO`、`AED`、`EXPRESSION_LOSS`、`NO_DIFFERENCE`、`INSUFFICIENT_DATA`。

每个结果都会给出组织级证据、major/minor copy、可能的 innovating/lost copy、置信度、判定依据和完整运行参数。

## 安装

建议使用 Python 3.10 或更高版本：

```bash
python -m pip install -e .
```

浏览器界面还需要 Streamlit：

```bash
python -m pip install -e ".[web]"
```

## 快速开始

运行附带的、含外群祖先表达的示例：

```bash
python -m idgrm demo --output results/demo
```

分析具有生物学重复的表达矩阵：

```bash
python -m idgrm classify \
  --expression expression.tsv \
  --samples samples.tsv \
  --pairs pairs.tsv \
  --output results/my_run
```

加入外群单拷贝正交基因，以获得祖先参照的 sub/neo 判定：

```bash
python -m idgrm classify \
  --expression expression.tsv \
  --samples samples.tsv \
  --pairs pairs.tsv \
  --ancestor-expression ancestor_expression.tsv \
  --ancestor-samples ancestor_samples.tsv \
  --output results/with_outgroup
```

若矩阵的每一列已经是组织均值，可不提供 `--samples`。此时只使用表达 on/off 与 fold-change，结果会标记为 effect-only，证据等级低于重复样本统计。

原始计数可选择 `--normalization cpm` 或 `--normalization median-ratio`。正式分析更推荐先使用 DESeq2/edgeR 获得组织内差异证据，或使用下面的旧流程兼容入口。

## 兼容原 Perl/DESeq2 流程

旧脚本产生的每组织 `*.DESeq2.csv` 可以直接读取：

```bash
python -m idgrm legacy \
  --deseq-dir path/to/deseq2_results \
  --pairs Ath.tandem.pairs \
  --output results/legacy_import
```

默认识别类似 `Ath.td.pairs.root.DESeq2.csv` 的文件名。自定义命名可通过 `--filename-regex` 提供含 `(?P<tissue>...)` 的正则表达式。

## 浏览器界面

Windows 可直接双击 `start_idgrm_web.bat`。也可以在终端运行：

```bash
streamlit run app.py
```

界面支持上传表达矩阵、样本信息、重复基因对和可选外群表达，运行后可下载全部 TSV、JSON 与 HTML 报告。

## 输入格式

### 1. 表达矩阵

首列是基因 ID，其余列是样本或组织，数值必须非负：

```text
gene_id    root_1    root_2    leaf_1    leaf_2
GeneA      12.1      11.8      0.2       0.1
GeneB      0.1       0.2       15.2      14.9
```

支持 CSV/TSV。首列名称也可为 `gene`、`id`、`transcript_id` 或 `ancestor_id`。

### 2. 样本信息

```text
sample_id    tissue    replicate
root_1       root      1
root_2       root      2
leaf_1       leaf      1
leaf_2       leaf      2
```

### 3. 重复基因对

```text
pair_id    gene1    gene2    duplication_type    ancestor_id
p001       GeneA    GeneB    tandem             OutgroupGeneA
```

只有 `gene1`、`gene2` 是必需的。旧版表头 `Duplicate 1/2`、`Parent copy/Daughter copy`、`Transposed/Parental` 会自动识别；后两种还会记录 parent/daughter 角色。原数据中的无表头 `gene/type/gene/type` 或 `gene/location/gene/location/...` 表也可直接读取，第一、三列解释为基因。完全重复的 pair 行只保留一次，并在运行警告中记录。

### 4. 外群表达

格式与表达矩阵相同，行 ID 应对应 pair 表的 `ancestor_id`。组织名称必须与目标物种匹配。跨物种表达量应在可比的归一化尺度上，至少应能可靠地区分表达/不表达。

## sub 与 neo 的规则

第一层沿用 Science 规则：两个 copy 分别在至少一个组织显著高于对方（默认 `|log2FC| >= 1` 且组织内 `q <= 0.001`），得到 reciprocal `SUB_OR_NEO` 候选。

随后拆分：

- 有外群时，`SUBFUNCTIONALIZATION` 要求两个 copy 在不同的祖先表达组织中分别占优势，二者并集覆盖大部分祖先表达域，且没有可靠的祖先未表达组织增益。
- 有外群时，`NEOFUNCTIONALIZATION` 要求恰有一个 copy 在祖先未表达组织形成明显、copy-specific 的表达增益，同时另一个 copy 保留足够多的祖先表达域。
- 无外群时，sub 是“平衡的组织间互补”代理；neo 是“一个广谱保留 copy + 一个窄谱、copy-specific 表达域 copy”代理。此时 `evidence_basis=expression_only_proxy`，不能当作对新生化功能的直接证明。
- reciprocal 候选不满足任一充分条件时保留为 `AMBIGUOUS_SUB_NEO`，不进行强制二分。

完整公式、阈值和方法边界见 [docs/METHODOLOGY_CN.md](docs/METHODOLOGY_CN.md)。

iDGRM 与原 Perl/R 脚本保持核心分类框架兼容，但默认不承诺逐行完全相同；具体一致项、修正项和验证限制见 [docs/VALIDATION_CN.md](docs/VALIDATION_CN.md)。

## 输出文件

- `classifications.tsv`：每个重复基因对一行的最终分类和证据指标。
- `tissue_evidence.tsv`：每个基因对 × 组织的均值、表达状态、log2FC、P/q 值与方向。
- `summary.tsv`：两套分类体系的数量和比例。
- `run_metadata.json`：输入文件 SHA-256、软件版本、参数、警告和时间戳。
- `report.html`：可搜索、可筛选的自包含浏览器报告。

## 重要解释边界

iDGRM 的 neo 是“表达调控层面的新功能化证据”。真正的蛋白质新功能仍需结合外群、编码序列选择、蛋白结构、互作、表型或实验功能验证。相似度很高的年轻重复基因还会受到多重比对偏差影响；应优先使用唯一可比对区域、可靠的定量方法或经模拟验证的 read assignment。

## 测试

```bash
python -m unittest discover -s tests -v
```
