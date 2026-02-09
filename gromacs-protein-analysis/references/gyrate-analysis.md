# 回转半径（Gyrate）分析

## 概述

回转半径（Rg）测量原子相对于质心的质量加权均方根距离。Gyrate 分析评估蛋白质紧致性，并提供关于折叠/展开转变、整体构象变化和整体蛋白质大小的见解。

## 何时使用 Gyrate

- 监测模拟期间的蛋白质紧致性
- 检测折叠/展开转变
- 评估整体构象变化
- 比较不同状态或突变体之间的紧致性
- 验证蛋白质稳定性
- 研究结构域运动或膨胀/收缩
- 监测蛋白质聚集倾向

## 前提条件

- 轨迹文件（.xtc/.trr）- PBC 修正是可选的（如需要，请参阅 [PBC 修正指南](pbc-correction.md)）
- 拓扑文件（.tpr）
- 包含适当原子组的索引文件（.ndx）

## 工作流程

### 步骤 1：计算回转半径

为不同原子组计算 Rg：

```bash
# 蛋白质主链的 Rg
echo -e "Backbone\n" | gmx gyrate -s md.tpr -f md.xtc -o gyrate_backbone.xvg

# C-alpha 原子的 Rg
echo -e "C-alpha\n" | gmx gyrate -s md.tpr -f md.xtc -o gyrate_calpha.xvg

# 所有蛋白质原子的 Rg
echo -e "Protein\n" | gmx gyrate -s md.tpr -f md.xtc -o gyrate_protein.xvg

# 蛋白质-配体复合物的 Rg
echo -e "Protein_Lig\n" | gmx gyrate -s md.tpr -f md.xtc -o gyrate_complex.xvg
```

**参数**：
- `-s md.tpr`：结构文件
- `-f md.xtc`：轨迹文件
- `-o gyrate_*.xvg`：Rg 数据的输出文件

**输出文件**：
- `gyrate_backbone.xvg`：主链 Rg 随时间的变化
- `gyrate_calpha.xvg`：C-alpha Rg 随时间的变化
- `gyrate_protein.xvg`：完整蛋白质 Rg 随时间的变化
- `gyrate_complex.xvg`：复合物 Rg 随时间的变化

### 步骤 2：计算每个轴的 Rg（可选）

分别计算每个主轴的 Rg：

```bash
echo -e "Protein\n" | gmx gyrate -s md.tpr -f md.xtc -o gyrate_axes.xvg -nm
```

**参数**：
- `-nm`：输出每个轴（x、y、z）的 Rg

**输出文件**：
- `gyrate_axes.xvg`：每个轴的 Rg 随时间的变化

**输出中的列**：
- 第 1 列：时间
- 第 2 列：总 Rg
- 第 3 列：沿 x 轴的 Rg
- 第 4 列：沿 y 轴的 Rg
- 第 5 列：沿 z 轴的 Rg

### 步骤 3：计算特定结构域的 Rg（可选）

计算蛋白质结构域的 Rg：

```bash
# 创建包含结构域组的索引文件
echo -e "r 1-100\nname 10 Domain1\nr 101-200\nname 11 Domain2\nq\n" | gmx make_ndx -f md.tpr -o domains.ndx

# 计算 Domain1 的 Rg
echo -e "Domain1\n" | gmx gyrate -s md.tpr -f md.xtc -n domains.ndx -o gyrate_domain1.xvg
```

### 步骤 4：可视化 Rg

生成 Rg 随时间变化的图：

```bash
# 基本 Rg 图
dit xvg_show -f gyrate_protein.xvg -x "Time (ns)" -y "Rg (nm)" -t "Radius of Gyration"

# 带平滑线的 Rg 图
dit xvg_show -f gyrate_protein.xvg -x "Time (ns)" -y "Rg (nm)" -t "Radius of Gyration" -m line --smooth

# 比较多个 Rg 曲线
dit xvg_show -f gyrate_backbone.xvg gyrate_calpha.xvg gyrate_protein.xvg -x "Time (ns)" -y "Rg (nm)" -t "Rg Comparison"

# 绘制每个轴的 Rg
dit xvg_show -f gyrate_axes.xvg -c 1:3 -x "Time (ns)" -y "Rg (nm)" -t "Rg per Axis" --multi
```

### 步骤 5：计算 Rg 分布（可选）

分析 Rg 分布以了解紧致性采样：

```bash
# 提取 Rg 值
# 计算直方图
# 拟合到分布模型
```

## 输出文件

- **gyrate_backbone.xvg**：主链 Rg 随时间的变化
- **gyrate_calpha.xvg**：C-alpha Rg 随时间的变化
- **gyrate_protein.xvg**：完整蛋白质 Rg 随时间的变化
- **gyrate_complex.xvg**：复合物 Rg 随时间的变化
- **gyrate_axes.xvg**：每个轴的 Rg 随时间的变化
- **gyrate_domain*.xvg**：特定结构域的 Rg

## 解释指南

### Rg 值

不同蛋白质大小的典型 Rg 范围：

- **< 1.0 nm**：小蛋白质（< 50 个残基）
- **1.0-1.5 nm**：小到中等蛋白质（50-100 个残基）
- **1.5-2.0 nm**：中等蛋白质（100-200 个残基）
- **2.0-2.5 nm**：中到大蛋白质（200-300 个残基）
- **2.5-3.0 nm**：大蛋白质（300-500 个残基）
- **> 3.0 nm**：非常大的蛋白质（> 500 个残基）

经验关系：Rg ≈ 0.395 × N^0.588（对于球形蛋白质）

### Rg 模式分析

- **稳定的 Rg**：围绕常数波动（~2-5% 的变化）
- **下降的 Rg**：蛋白质变得更加紧密（折叠或塌陷）
- **上升的 Rg**：蛋白质膨胀（展开或变性）
- **突然跳跃**：大的构象变化或 PBC 伪影
- **漂移**：连续增加/减少表示非平衡

### 原子选择的影响

- **主链**：反映整体蛋白质形状，常用
- **C-alpha**：类似于主链，计算效率高
- **所有原子**：由于侧链体积，值略大
- **重原子**：类似于所有原子，排除氢

### 每个轴分析

- **各向异性 Rg**：x、y、z 轴的 Rg 值不同
- **各向同性 Rg**：所有轴的 Rg 值相似（球形蛋白质）
- **形状信息**：轴的比率表示蛋白质形状（细长与球形）

### 紧致性评估

- **紧密球形**：相对于大小 Rg 低，各向同性
- **伸展构象**：相对于大小 Rg 高，各向异性
- **部分展开**：中等 Rg，可能显示漂移
- **完全展开**：高 Rg，各向同性（无规卷曲）

## 常见问题和解决方案

### 问题：Rg 显示突然跳跃

**可能原因**：
- PBC 伪影（分子跨越盒边界）
- 蛋白质居中问题
- 数值不稳定性

**解决方案**：
- 应用 PBC 修正：`gmx trjconv -pbc mol -center`
- 在 `trjconv` 中使用 `-pbc nojump` 标志
- 检查轨迹是否有伪影

### 问题：Rg 持续增加

**可能原因**：
- 蛋白质展开
- 模拟未收敛
- 力场不充分

**解决方案**：
- 检查蛋白质结构是否有展开
- 延长模拟时间
- 验证力场兼容性

### 问题：Rg 值异常低

**可能原因**：
- 原子选择不正确
- 蛋白质塌陷（过度紧密）
- 轨迹居中问题

**解决方案**：
- 验证原子选择组
- 检查蛋白质结构是否有不自然的塌陷
- 验证轨迹处理

### 问题：每个轴的 Rg 显示大的各向异性

**可能原因**：
- 蛋白质天然细长（例如，纤维状蛋白质）
- 轨迹未正确旋转
- 模拟盒形状效应

**解决方案**：
- 检查各向异性是否对蛋白质类型是预期的
- 验证轨迹拟合
- 检查模拟盒尺寸

## 提示和最佳实践

- **原子选择**：主链 Rg 是整体紧致性的标准。为提高效率使用 C-alpha。
- **轨迹拟合**：对于 Rg 并非严格必要，但建议保持一致性。
- **时间选择**：仅使用生产阶段（排除平衡）。
- **多个选择**：为不同原子组计算 Rg 以了解不同方面。
- **每个轴分析**：用于识别各向异性形状变化。
- **统计分析**：计算平均 Rg、标准差和置信区间。
- **与 RMSD 的相关性**：结合 Rg 和 RMSD 进行全面的稳定性评估。
- **目视检查**：检查轨迹以解释大的 Rg 变化。

## 高级分析

### 时间依赖性 Rg 分布

为不同时间窗口计算 Rg 分布：

```bash
# 提取前半部分和后半部分的 Rg
# 使用直方图或 KDE 比较分布
# 识别紧致性随时间的变化
```

### 不同模拟之间的 Rg 差异

比较不同条件之间的 Rg：

```bash
# 计算每个模拟的 Rg
# 并排绘制进行比较
# 使用统计检验评估显著性
```

### 与 RMSD 的相关性

分析 Rg 和 RMSD 之间的关系：

```bash
# 绘制 Rg 与 RMSD 的关系图
# 计算相关系数
# 识别构象转变
```

### 主轴分析

使用主轴分析蛋白质形状：

```bash
# 从协方差矩阵提取主轴
# 计算轴比率
- 与蛋白质形状和功能相关联
```

### 折叠/展开动力学

在折叠/展开期间监测 Rg：

```bash
# 使用 Rg 作为反应坐标
- 识别中间状态
- 计算转变速率
- 绘制折叠途径
```

## 相关分析

- **RMSD**：提供结构稳定性，用 Rg 补充
- **RMSF**：提供局部灵活性，Rg 提供整体紧致性
- **SASA**：与紧致性和溶剂暴露相关
- **PCA**：识别影响 Rg 的集体运动
- **FEL**：使用 Rg 作为构象景观的反应坐标

## 可视化增强

### 添加紧致性阈值

在 Rg 图上标记紧致性阈值：

```bash
dit xvg_show -f gyrate_protein.xvg -x "Time (ns)" -y "Rg (nm)" -t "Radius of Gyration" --hline 2.0 --hlabel "Compactness Threshold"
```

### 多个比较图

比较来自多个模拟的 Rg：

```bash
dit xvg_show -f gyrate_wt.xvg gyrate_mut1.xvg gyrate_mut2.xvg -x "Time (ns)" -y "Rg (nm)" -t "Rg Comparison"
```

### Rg 与 RMSD 散点图

分析 Rg 和 RMSD 之间的关系：

```bash
# 使用自定义脚本绘制 Rg 与 RMSD 的关系图
# 识别构象状态
```

### 形状分析

可视化蛋白质形状变化：

```bash
# 提取极端构象（最小/最大 Rg）
# 叠加并比较
- 识别负责膨胀/收缩的区域
```

## 解释示例

### 示例 1：稳定的球形蛋白质

- Rg 稳定在 ~1.8 nm
- 小波动（~0.05 nm）
- 每个轴 Rg 各向同性
- **解释**：蛋白质折叠良好且稳定

### 示例 2：展开的蛋白质

- Rg 从 1.8 nm 增加到 3.5 nm
- 整个模拟过程中的连续漂移
- 展开时变得更加各向同性
- **解释**：蛋白质正在展开，可能需要不同的条件

### 示例 3：构象变化

- Rg 在前 50 ns 内稳定在 1.8 nm
- 在 50 ns 时突然跳跃到 2.2 nm
- 在剩余的 50 ns 内稳定在新水平
- **解释**：蛋白质经历构象转变到更伸展的状态

### 示例 4：细长蛋白质

- Rg 稳定在 2.5 nm
- 强各向异性（x: 2.0, y: 2.8, z: 3.0 nm）
- **解释**：蛋白质天然细长（例如，纤维状蛋白质）

## 与实验数据的比较

### 小角 X 射线散射（SAXS）

将 Rg 与 SAXS 数据进行比较：

- 从 SAXS 散射曲线提取 Rg
- 与模拟 Rg 比较
- 良好的一致性验证模拟

### 流体动力学半径

将 Rg 与流体动力学半径（Rh）相关联：

```
Rh ≈ Rg / 0.775（对于球形蛋白质）
```

与来自 DLS 或粘度计的实验 Rh 进行比较。

### 晶体学堆积

将模拟 Rg 与晶体结构进行比较：

- 由于堆积，晶体结构可能更紧密
- 模拟在溶液中应显示略大的 Rg

## 理论背景

### Rg 计算

Rg 计算为：

```
Rg² = Σmi × ri² / Σmi
```

其中 mi 是原子 i 的质量，ri 是距离质心的距离。

### 与蛋白质大小的关系

对于球形蛋白质：

```
Rg ≈ R0 × N^ν
```

其中：
- N = 残基数
- R0 ≈ 0.395 nm
- ν ≈ 0.588（对于展开：ν ≈ 0.588，对于折叠：ν ≈ 0.33）

### 各向异性因子

计算各向异性因子：

```
A = 1 - 3 × (Rg_x × Rg_y × Rg_z) / (Rg_x² + Rg_y² + Rg_z²)^1.5
```

- A ≈ 0：各向同性（球形）
- A → 1：高度各向异性（细长）

## 参考文献

有关理论背景，请参阅：
- Flory（1953）"Principles of Polymer Chemistry"
- Miller 等人（2002）"Radius of gyration of proteins"
- Kohn 等人（2004）"Dynamics of protein folding"