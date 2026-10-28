# iDGRM 2.0.0 验证与质量控制

软件测试覆盖输入校验、重复对去重、旧格式差异结果导入、四类判定、copy 交换与组织重排不变性、缺失统计量、输出防覆盖以及原始整数 counts 接口。公开接口只包含 `fit-copy-counts` 和 `classify-primary`；发布检查确认已删除的细分命令和 Python 模块不可使用，安装包不包含相关示例。

随包示例仅验证程序行为，不是生物学准确率基准。counts 集成检查可使用 `scripts/verify_release.py --rscript <Rscript>` 运行真实 R/DESeq2 的 Wald/LRT 拟合。发布前还应验证安装后的 wheel、源代码归档、输入输出 SHA-256 与命令行帮助。

分类结果只能包含 `UNMAPPED`、`NO_DIFFERENCE`、`ASYMMETRICALLY_EXPRESSED` 和 `SUB_OR_NEOFUNCTIONALIZED`，每个输入重复对只有一个标签。SUB 与 NEO 保持合并。缺失检验不能作为无差异的证据，表达证据不足不等于独立证明不可比对。

实际数据应检查唯一归属与完全匹配要求、整数 counts、copy 可观测区域及有效长度、组织和生物学重复、离散度拟合、NA 比例、效应与 P 值分布。版本间标签一致只验证实现连续性；逐组织 BH 校正不保证最终四分类的错误发现率，也不证明功能命运。
