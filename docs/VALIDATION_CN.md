# iDGRM 验证与质量控制

## 软件级验证

自动测试覆盖输入格式、重复基因对去重、parent/daughter字段识别、组织级差异证据导入、四类初级分类、亚型细分和初级分类追溯。发布前要求全部测试通过，并检查Python字节码编译和Git差异格式。

## 分类不变量

1. 初级结果只能包含 `UNMAPPED`、`NO_DIFFERENCE`、`ASYMMETRICALLY_EXPRESSED` 和 `SUB_OR_NEOFUNCTIONALIZED`；
2. 每个重复基因对只能有一个初级类别；
3. 亚型细分只能接收 `ASYMMETRICALLY_EXPRESSED` 和 `SUB_OR_NEOFUNCTIONALIZED`；
4. 亚型结果必须保留 `parent_primary_class`；
5. 表达沉默相关模式不能成为独立一级类别；
6. 缺少外群时必须明确记录未极化表达推断的证据等级。

## 数据级验证

正式数据分析应检查：严格唯一比对比例、整数计数完整性、样本与组织对应、重复数量、低计数过滤、有效长度及可唯一比对长度、P值分布、离散度拟合、效应量分布和阈值敏感性。若大量基因对在极严格阈值下显著，应先排除模型设定或copy间可比性问题，再解释为生物学差异。

## 与参考实现的比较

推荐报告初级四分类的一致率、Cohen's kappa、各类precision/recall和混淆矩阵。iDGRM对亚型的扩展应单独报告，不能与初级分类混合计算一致率。
