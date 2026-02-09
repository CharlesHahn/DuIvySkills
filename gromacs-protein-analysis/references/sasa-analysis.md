# 溶剂可及表面积（SASA）分析

## 概述

溶剂可及表面积（SASA）测量蛋白质可被溶剂分子接触的表面积。SASA 分析提供了关于蛋白质-溶剂相互作用、疏水/亲水表面分布、结合位点可及性和蛋白质展开的见解。

## 何时使用 SASA

- 分析蛋白质表面的溶剂暴露
- 识别疏水和亲水区域
- 研究配体结合位点可及性
- 监测蛋白质展开或聚集
- 评估蛋白质溶解性
- 比较不同状态或突变体之间的表面性质
- 验证蛋白质结构稳定性

## 前提条件

- 轨迹文件（.xtc/.trr）- PBC 修正是可选的（如需要，请参阅 [PBC 修正指南](pbc-correction.md)）
- 拓扑文件（.tpr）
- 包含适当原子组的索引文件（.ndx）
- 支持 SASA 的 GROMACS 版本

## 工作流程

### 步骤 1：计算总 SASA

为不同原子组计算 SASA：

```bash
# 蛋白质主链的 SASA
echo -e "Backbone\n" | gmx sasa -s md.tpr -f md.xtc -o sas_backbone.xvg -tu ns

# C-alpha 原子的 SASA
echo -e "C-alpha\n" | gmx sasa -s md.tpr -f md.xtc -o sas_calpha.xvg -tu ns

# 所有蛋白质原子的 SASA
echo -e "Protein\n" | gmx sasa -s md.tpr -f md.xtc -o sas_protein.xvg -tu ns

# 蛋白质-配体复合物的 SASA
echo -e "Protein_Lig\n" | gmx sasa -s md.tpr -f md.xtc -o sas_complex.xvg -tu ns
```

**参数**：
- `-s md.tpr`：结构文件
- `-f md.xtc`：轨迹文件
- `-o sas_*.xvg`：SASA 数据的输出文件
- `-tu ns`：时间单位为纳秒

**输出文件**：
- `sas_backbone.xvg`：主链 SASA 随时间的变化
- `sas_calpha.xvg`：C-alpha SASA 随时间的变化
- `sas_protein.xvg`：完整蛋白质 SASA 随时间的变化
- `sas_complex.xvg`：复合物 SASA 随时间的变化

### 步骤 2：计算疏水和亲水 SASA

将 SASA 分为疏水和亲水成分：

```bash
# 带有疏水/亲水分解的总 SASA
echo -e "Protein\n" | gmx sasa -s md.tpr -f md.xtc -o sas_protein.xvg -output sas_protein_components.xvg -surface Protein -tu ns
```

**参数**：
- `-output`：带有成分分解的输出文件
- `-surface`：用于 SASA 计算的表面组

**输出文件**：
- `sas_protein_components.xvg`：包含多列：
  - 第 1 列：时间
  - 第 2 列：总 SASA
  - 第 3 列：疏水 SASA
  - 第 4 列：亲水 SASA

### 步骤 3：计算每个残基的 SASA（可选）

为每个残基计算 SASA：

```bash
# 每个残基的 SASA
echo -e "Protein\n" | gmx sasa -s md.tpr -f md.xtc -o sas_residue.xvg -or sas_per_residue.xvg -res -tu ns
```

**参数**：
- `-or`：输出每个残基的 SASA
- `-res`：计算每个残基的 SASA

**输出文件**：
- `sas_residue.xvg`：总 SASA 随时间的变化
- `sas_per_residue.xvg`：每个残基的 SASA（随时间平均）

### 步骤 4：计算特定结构域的 SASA（可选）

计算蛋白质结构域的 SASA：

```bash
# 创建包含结构域组的索引文件
echo -e "r 1-100\nname 10 Domain1\nr 101-200\nname 11 Domain2\nq\n" | gmx make_ndx -f md.tpr -o domains.ndx

# 计算 Domain1 的 SASA
echo -e "Domain1\n" | gmx sasa -s md.tpr -f md.xtc -n domains.ndx -o sas_domain1.xvg
```

### 步骤 5：调整探针半径（可选）

使用不同的探针半径进行 SASA 计算：

```bash
# 默认探针半径：0.14 nm（水）
# 使用更大的探针分析更大的溶剂分子
echo -e "Protein\n" | gmx sasa -s md.tpr -f md.xtc -o sas_large_probe.xvg -probe 0.2
```

**参数**：
- `-probe`：探针半径，单位为 nm（默认：0.14 nm）

### 步骤 6：可视化 SASA

生成 SASA 随时间变化的图：

```bash
# 基本 SASA 图
dit xvg_show -f sas_protein.xvg -x "Time (ns)" -y "SASA (nm²)" -t "Total SASA"

# 疏水/亲水分解
dit xvg_show -f sas_protein_components.xvg -c 1:4 -x "Time (ns)" -y "SASA (nm²)" -t "SASA Components" --multi

# 每个残基的 SASA
dit xvg_show -f sas_per_residue.xvg -x "Residue Number" -y "SASA (nm²)" -t "Per-Residue SASA"

# 比较多个 SASA 曲线
dit xvg_show -f sas_wt.xvg sas_mut1.xvg -x "Time (ns)" -y "SASA (nm²)" -t "SASA Comparison"
```

## 输出文件

- **sas_protein.xvg**：总 SASA 随时间的变化
- **sas_backbone.xvg**：主链 SASA 随时间的变化
- **sas_calpha.xvg**：C-alpha SASA 随时间的变化
- **sas_complex.xvg**：复合物 SASA 随时间的变化
- **sas_protein_components.xvg**：总、疏水和亲水 SASA
- **sas_per_residue.xvg**：每个残基的 SASA（平均）
- **sas_domain*.xvg**：特定结构域的 SASA

## 解释指南

### SASA 值

不同蛋白质大小的典型 SASA 范围：

- **< 50 nm²**：小蛋白质（< 50 个残基）
- **50-100 nm²**：小到中等蛋白质（50-100 个残基）
- **100-200 nm²**：中等蛋白质（100-200 个残基）
- **200-300 nm²**：中到大蛋白质（200-300 个残基）
- **300-500 nm²**：大蛋白质（300-500 个残基）
- **> 500 nm²**：非常大的蛋白质（> 500 个残基）

经验关系：SASA ≈ 6.3 × N^0.73（对于球形蛋白质）

### SASA 模式分析

- **稳定的 SASA**：围绕常数波动（~5-10% 的变化）
- **下降的 SASA**：蛋白质变得更加被掩埋（折叠或塌陷）
- **上升的 SASA**：蛋白质变得更加暴露（展开或变性）
- **突然变化**：大的构象变化或配体结合/解离
- **漂移**：连续增加表示展开

### 疏水与亲水 SASA

- **疏水 SASA**：非极性表面积（通常占总数的 40-60%）
- **亲水 SASA**：极性表面积（通常占总数的 40-60%）
- **疏水/亲水比率**：指示表面特征
  - > 1.0：主要为疏水（可能聚集）
  - ~1.0：平衡表面（典型于可溶性蛋白质）
  - < 1.0：主要为亲水（高度可溶）

### 每个残基的 SASA

- **低 SASA（< 0.1 nm²）**：被掩埋的残基（蛋白质核心）
- **中等 SASA（0.1-0.5 nm²）**：部分暴露的残基
- **高 SASA（> 0.5 nm²）**：表面暴露的残基

### 结合位点分析

- **活性位点 SASA**：通常中等（0.2-0.4 nm²）以保持可及性
- **配体结合**：配体结合时 SASA 减少
- **变构位点**：可能显示动态 SASA 变化

## 常见问题和解决方案

### 问题：SASA 显示突然跳跃

**可能原因**：
- PBC 伪影（分子跨越盒边界）
- 蛋白质居中问题
- 数值不稳定性

**解决方案**：
- 应用 PBC 修正：`gmx trjconv -pbc mol -center`
- 在 `trjconv` 中使用 `-pbc nojump` 标志
- 检查轨迹是否有伪影

### 问题：SASA 持续增加

**可能原因**：
- 蛋白质展开
- 模拟未收敛
- 力场不充分

**解决方案**：
- 检查蛋白质结构是否有展开
- 延长模拟时间
- 验证力场兼容性

### 问题：SASA 值异常低

**可能原因**：
- 原子选择不正确
- 蛋白质过度紧密
- 探针半径太小

**解决方案**：
- 验证原子选择组
- 检查蛋白质结构
- 如需要，调整探针半径

### 问题：疏水/亲水分解不工作

**可能原因**：
- GROMACS 版本不支持成分输出
- 缺少表面组规范

**解决方案**：
- 检查 GROMACS 版本（需要 2020 或更高版本）
- 使用 `-surface` 参数指定表面组

## 提示和最佳实践

- **原子选择**：全原子 SASA 是标准的。为提高效率使用主链。
- **探针半径**：默认 0.14 nm（水）对大多数分析是合适的。
- **时间选择**：仅使用生产阶段（排除平衡）。
- **多个选择**：为不同原子组计算 SASA。
- **成分分析**：分离疏水和亲水 SASA 以获得更深入的见解。
- **每个残基分析**：有助于识别结合位点或聚集倾向区域。
- **统计分析**：计算平均 SASA、标准差和置信区间。
- **与其他指标的相关性**：结合 Rg、RMSD 进行全面分析。

## 高级分析

### 时间依赖性 SASA 分布

为不同时间窗口计算 SASA 分布：

```bash
# 提取前半部分和后半部分的 SASA
# 使用直方图或 KDE 比较分布
# 识别溶剂暴露随时间的变化
```

### 不同模拟之间的 SASA 差异

比较不同条件之间的 SASA：

```bash
# 计算每个模拟的 SASA
# 并排绘制进行比较
# 使用统计检验评估显著性
```

### SASA 与 Rg 的相关性

分析 SASA 和 Rg 之间的关系：

```bash
# 绘制 SASA 与 Rg 的关系图
# 计算相关系数
- 识别展开中间体
```

### 结合位点可及性

监测结合位点残基的 SASA：

```bash
# 提取结合位点残基
# 计算这些残基随时间的 SASA
- 评估配体结合时的可及性变化
```

### 聚集倾向

识别聚集倾向区域：

```bash
# 计算每个残基的 SASA
# 识别具有高 SASA 的疏水残基
- 使用聚集预测工具（例如，TANGO、AGGRESCAN）
```

## 相关分析

- **RMSD**：提供结构稳定性，用 SASA 补充
- **Rg**：与紧致性和 SASA 相关
- **RMSF**：局部灵活性可能与 SASA 变化相关
- **氢键**：溶剂-蛋白质氢键与亲水 SASA 相关
- **FEL**：使用 SASA 作为构象景观的反应坐标

## 可视化增强

### 添加溶剂可及性阈值

在图上标记 SASA 阈值：

```bash
dit xvg_show -f sas_protein.xvg -x "Time (ns)" -y "SASA (nm²)" -t "Total SASA" --hline 150 --hlabel "Native State"
```

### 疏水/亲水比率

计算并绘制疏水/亲水比率：

```bash
# 从成分文件中提取疏水和亲水 SASA
# 计算比率：疏水 / 亲水
# 绘制比率随时间的变化
```

### 表面映射

将 SASA 映射到蛋白质结构上：

```bash
# 使用 PyMOL 或 VMD 按 SASA 对结构着色
# 红色 = 高 SASA（暴露）
# 蓝色 = 低 SASA（被掩埋）
```

### SASA 热图

创建多个模拟的每个残基 SASA 的热图：

```bash
# 组合来自多个模拟的每个残基 SASA 数据
# 生成显示表面暴露模式的热图
```

## 解释示例

### 示例 1：稳定的可溶性蛋白质

- SASA 稳定在 ~150 nm²
- 小波动（~10 nm²）
- 平衡的疏水/亲水比率（~1.0）
- **解释**：蛋白质折叠良好且可溶

### 示例 2：展开的蛋白质

- SASA 从 150 nm² 增加到 300 nm²
- 整个模拟过程中的连续漂移
- 疏水 SASA 显著增加
- **解释**：蛋白质正在展开，暴露疏水核心

### 示例 3：配体结合

- 配体结合时 SASA 减少约 20 nm²
- 减少发生在结合位点残基
- **解释**：配体结合到蛋白质，掩埋表面积

### 示例 4：聚集倾向蛋白质

- 高疏水 SASA（> 总数的 60%）
- 疏水表面残基的 SASA 大
- **解释**：蛋白质由于暴露的疏水表面而具有聚集倾向

## 与实验数据的比较

### 氢-氘交换（HDX）

将 SASA 与 HDX 保护因子进行比较：

- 低 SASA 区域应受到保护（慢交换）
- 高 SASA 区域应可接触（快交换）
- 验证模拟表面暴露

### 小角 X 射线散射（SAXS）

将 SASA 与 SAXS 散射相关联：

- SASA 影响散射剖面
- 将模拟 SASA 与 SAXS 导出的值进行比较
- 良好的一致性验证模拟

### NMR 化学位移

将 SASA 与 NMR 化学位移扰动进行比较：

- 表面暴露的残基显示较大的化学位移变化
- 被掩埋的残基显示最小变化
- 验证模拟溶剂可及性

## 理论背景

### SASA 计算

SASA 使用滚动探针方法计算：

1. 在每个原子周围放置探针球体（水半径为 0.14 nm）
2. 计算探针可接触的表面积
3. 汇总所有原子的贡献

### 与蛋白质大小的关系

对于球形蛋白质：

```
SASA ≈ A0 × N^α
```

其中：
- N = 残基数
- A0 ≈ 6.3 nm²
- α ≈ 0.73

### 疏水性标度

根据标度将残基分类为疏水或亲水：

- **疏水**：Ala、Val、Leu、Ile、Met、Phe、Trp、Pro、Tyr
- **亲水**：Asp、Glu、Asn、Gln、Lys、Arg、His、Ser、Thr、Cys

## 溶剂可及性和蛋白质功能

### 活性位点

- 通常具有中等 SASA（可及性和特异性的平衡）
- 可能在催化过程中显示动态 SASA 变化

### 膜蛋白质

- 具有独特的 SASA 模式（疏水跨膜区域）
- 可能需要使用膜模拟溶剂进行分析

### 蛋白质-蛋白质相互作用

- 复合物中的界面残基通常具有降低的 SASA
- 通过复合物形成时的 SASA 减少来识别结合位点

## 参考文献

有关理论背景，请参阅：
- Lee 和 Richards（1971）"Interpretation of protein structures: estimation of static accessibility"
- Shrake 和 Rupley（1973）"Environment and exposure to solvent of protein atoms"
- Richmond（1984）"Solvent accessible surface area and excluded volume in proteins"