# 蛋白质动力学的主成分分析（PCA）

## 概述

主成分分析（PCA）通过将蛋白质运动分解为主成分来识别蛋白质中的集体运动和主要构象变化。它在降低维度的同时保留最重要的动力学信息。

## 何时使用 PCA

- 识别集体运动和构象转变
- 提取蛋白质灵活性的主要模式
- 降低维度以进行进一步分析
- 分析蛋白质结构域运动
- 比较不同模拟之间的构象变化
- 理解功能运动

## 前提条件

- 轨迹文件（.xtc/.trr）- PBC 修正是可选的（如需要，请参阅 [PBC 修正指南](pbc-correction.md)）
- 拓扑文件（.tpr）
- 包含适当原子组的索引文件（.ndx）
- 轨迹应拟合以消除整体平移/旋转

## 工作流程

### 步骤 1：计算协方差矩阵

计算协方差矩阵并提取特征向量/特征值：

```bash
echo -e "C-alpha\nC-alpha\n" | gmx covar -s md.tpr -f md.xtc -o eigenvalues.xvg -v eigenvectors.trr -xpma covar.xpm
```

- 第一个输入：选择用于对齐的参考组（C-alpha）
- 第二个输入：选择用于协方差计算的组（C-alpha）

**输出文件**：
- `eigenvalues.xvg`：来自协方差矩阵对角化的特征值
- `eigenvectors.trr`：存储为轨迹的特征向量
- `covar.xpm`：协方差矩阵（可用于 DCCM）

### 步骤 2：分析特征值分布

检查特征值分布以了解集体运动：

```bash
# 绘制特征值
dit xvg_show -f eigenvalues.xvg -x "PC Number" -y "Eigenvalue" -t "Eigenvalue Spectrum"
```

计算贡献百分比：

```bash
# 读取 eigenvalues.xvg
# 第一列：PC 编号
# 第二列：特征值

# 总方差 = 所有特征值之和
# PC 贡献 = (单个特征值 / 总方差) × 100%

# 前 3 个 PC 的示例计算：
# PC1% = eigenvalue1 / total × 100%
# PC2% = eigenvalue2 / total × 100%
# PC3% = eigenvalue3 / total × 100%
# 累积% = (PC1 + PC2 + PC3) / total × 100%
```

**解释**：
- 前几个 PC 通常捕获 50-80% 的总运动
- 特征值的急剧下降表示主要的集体运动
- 均匀的特征值表明没有主要的集体运动

### 步骤 3：将轨迹投影到主成分上

将轨迹投影到前几个主成分上：

```bash
# 投影到 PC1 上
echo -e "C-alpha\nC-alpha\n" | gmx anaeig -s md.tpr -f md.xtc -v eigenvectors.trr -first 1 -last 1 -proj pc1.xvg

# 投影到 PC2 上
echo -e "C-alpha\nC-alpha\n" | gmx anaeig -s md.tpr -f md.xtc -v eigenvectors.trr -first 2 -last 2 -proj pc2.xvg

# 投影到 PC3 上
echo -e "C-alpha\nC-alpha\n" | gmx anaeig -s md.tpr -f md.xtc -v eigenvectors.trr -first 3 -last 3 -proj pc3.xvg
```

**输出文件**：
- `pc1.xvg`、`pc2.xvg`、`pc3.xvg`：投影到相应 PC 上的投影

### 步骤 4：可视化 PC 投影

生成 PC 投影随时间变化的图：

```bash
# PC1 投影
dit xvg_show -f pc1.xvg -x "Time (ps)" -y "PC1 Coordinate" -t "PC1 Projection"

# PC2 投影
dit xvg_show -f pc2.xvg -x "Time (ps)" -y "PC2 Coordinate" -t "PC2 Projection"

# PC1 vs PC2 散点图
dit xvg_show_scatter -f pc1.xvg pc2.xvg -c 1 1 -x "PC1" -y "PC2" -t "PC1 vs PC2"
```

### 步骤 5：提取极端构象

提取沿主成分的极值对应的结构：

```bash
# 提取沿 PC1 的极端构象
echo -e "C-alpha\nC-alpha\n" | gmx anaeig -s md.tpr -f md.xtc -v eigenvectors.trr -first 1 -last 1 -extreme pc1_extreme.pdb

# 提取平均结构
echo -e "C-alpha\nC-alpha\n" | gmx anaeig -s md.tpr -f md.xtc -v eigenvectors.trr -first 1 -last 1 -average average.pdb
```

**输出文件**：
- `pc1_extreme.pdb`：沿 PC1 具有最大偏差的结构
- `average.pdb`：轨迹上的平均结构

### 步骤 6：2D 投影（可选）

生成轨迹到 PC1-PC2 平面的 2D 投影：

```bash
echo -e "C-alpha\nC-alpha\n" | gmx anaeig -s md.tpr -f md.xtc -v eigenvectors.trr -first 1 -last 2 -2d 2dproj.xvg
```

**输出文件**：
- `2dproj.xvg`：随时间变化的 2D 投影（PC1 vs PC2）

## 输出文件

- **eigenvalues.xvg**：来自协方差矩阵的特征值
- **eigenvectors.trr**：存储为轨迹的特征向量
- **covar.xpm**：协方差矩阵
- **pc1.xvg、pc2.xvg、pc3.xvg**：投影到主成分上的投影
- **2dproj.xvg**：2D 投影（PC1 vs PC2）
- **pc1_extreme.pdb**：沿 PC1 的极端构象
- **average.pdb**：平均结构

## 解释指南

### 特征值分析

- **第一个 PC**：最大的特征值，捕获最主要的运动
- **累积方差**：前 N 个特征值之和 / 总方差
- **方差分数**：单个 PC 贡献（特征值 / 总方差）
- **维度**：解释 80-90% 方差所需的 PC 数量

### PC 投影

- **振幅**：PC 投影的范围表示运动的大小
- **时间尺度**：振荡频率与运动时间尺度相关
- **转变**：跳跃或位移表示构象变化
- **平稳性**：常数平均值表示平衡采样

### 极端构象

- **比较结构**：检查极端构象之间的差异
- **识别运动区域**：定位对 PC 运动贡献最大的原子
- **功能相关性**：将运动与蛋白质功能相关联

### 生物学解释

- **结构域运动**：结构域之间的大规模运动
- **环灵活性**：环和末端的局部运动
- **变构途径**：表示变构通讯的相关运动
- **功能运动**：与蛋白质活性相关的运动（例如，铰链弯曲、通道开放）

## 常见问题和解决方案

### 问题：特征值显示均匀分布

**可能原因**：
- 模拟时间不足
- 来自高频运动的过多噪声
- 没有主要的集体运动（蛋白质可能是刚性的）

**解决方案**：
- 延长模拟时间
- 应用时间平均或过滤
- 使用更大的原子组（主链而不是 C-alpha）

### 问题：PC 投影没有显示清晰模式

**解决方案**：
- 检查轨迹拟合（消除平移/旋转）
- 验证原子选择一致性
- 考虑使用不同的原子组

### 问题：极端构象看起来不切实际

**解决方案**：
- 验证模拟质量
- 检查数值不稳定性
- 考虑使用多个 PC 进行提取

### 问题：前几个 PC 没有捕获显著方差

**解决方案**：
- 可能表明蛋白质是刚性的或运动均匀分布
- 考虑使用更多 PC 进行分析
- 检查模拟时间和采样

## 提示和最佳实践

- **原子选择**：C-alpha 通常用于提高效率。使用主链以获得详细信息。
- **时间选择**：排除平衡。仅使用生产阶段。
- **轨迹拟合**：在 PCA 之前始终消除平移/旋转。
- **PC 数量**：分析前 3-5 个 PC，它们通常捕获大部分运动。
- **验证**：如有可用，与实验数据进行比较。
- **可视化**：使用多种可视化方法进行综合分析。

## 高级分析

### 时间依赖性 PCA

为不同时间窗口计算 PCA：

```bash
# 前半部分的 PCA
echo -e "C-alpha\nC-alpha\n" | gmx covar -s md.tpr -f md.xtc -b 0 -e 50000 -o eigenvalues_early.xvg -v eigenvectors_early.trr

# 后半部分的 PCA
echo -e "C-alpha\nC-alpha\n" | gmx covar -s md.tpr -f md.xtc -b 50000 -e 100000 -o eigenvalues_late.xvg -v eigenvectors_late.trr
```

比较特征值和特征向量以识别集体运动的变化。

### 重叠分析

计算来自不同模拟的特征向量之间的重叠：

```bash
echo -e "C-alpha\nC-alpha\n" | gmx anaeig -s md.tpr -f md.xtc -v eigenvectors1.trr -v2 eigenvectors2.trr -overlap overlap.xvg
```

量化集体运动之间的相似性。

### 准谐波分析

从 PCA 计算准谐波熵：

```bash
# 使用特征值计算熵
# S_qh = (k_B/2) * sum(ln(k_B*T/λ_i))
```

需要对特征值进行额外处理。

### 基本动力学

专注于基本子空间（前几个 PC）：

```bash
# 将轨迹投影到基本子空间
# 分析子空间内的采样
- 与实验数据比较
```

## 相关分析

- **DCCM**：分析相关运动，补充 PCA
- **FEL**：使用 PC 作为反应坐标映射构象景观
- **RMSF**：提供每个残基的灵活性信息
- **正模分析**：将理论模式与 PCA 导出的模式进行比较

## 可视化增强

### 3D 可视化

使用 PyMOL 或 VMD 将特征向量投影到蛋白质结构上：

```bash
# 使用 gmx anaeig 生成特征向量 PDB
echo -e "C-alpha\nC-alpha\n" | gmx anaeig -s md.tpr -f md.xtc -v eigenvectors.trr -first 1 -last 1 -filt eigenvector1.pdb
```

### 豪猪图

生成显示原子运动方向和幅度的豪猪图：

```bash
# 需要外部工具（例如，VMD 中的 NMWiz）
# 将特征向量可视化为蛋白质结构上的箭头
```

### PC 的交叉相关

分析 PC 之间的相关性：

```bash
# 计算 PC 投影的交叉相关矩阵
# 使用 DuIvyTools 或自定义脚本
```

## 与实验数据的比较

### 与 NMR 序列参数比较

将 NMR 序列参数投影到 PCA 模式上以验证运动。

### 与晶体结构比较

将 PC1 极端构象与不同的晶体结构（例如，开放/关闭状态）进行比较。

### 与 SAXS 比较

从 PC 导出的构象计算理论 SAXS 曲线并与实验数据进行比较。

## 参考文献

有关理论背景，请参阅：
- Amadei 等人（1993）"Essential dynamics of proteins"
- García（1992）"Large-amplitude collective motions in proteins"
- Berendsen 和 Hayward（2000）"Collective variables and molecular dynamics"