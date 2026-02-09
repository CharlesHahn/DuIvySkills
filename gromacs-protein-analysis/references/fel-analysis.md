# 自由能景观（FEL）分析

## 概述

自由能景观（FEL）绘制蛋白质构象的能量表面，识别稳定状态、能垒和转变途径。两种常用方法使用 RMSD + Gyrate 或主成分（PC）作为反应坐标。

## 何时使用 FEL

- 识别稳定的构象状态
- 量化状态之间的能垒
- 理解构象转变
- 分析蛋白质折叠/展开
- 研究结构域运动和变构
- 比较不同条件之间的自由能表面

## 两种方法

### 方法 1：RMSD + Gyrate

使用结构偏差和紧致性作为反应坐标。

**优点**：
- 直观的解释
- 直接与结构变化相关
- 易于计算

**缺点**：
- 可能无法捕获细微的构象变化
- 对参考结构选择敏感

### 方法 2：主成分

使用 PC 投影作为反应坐标。

**优点**：
- 捕获集体运动
- 基于动力学的自然坐标
- 揭示主要的构象变化

**缺点**：
- 需要先进行 PCA 分析
- 解释不够直观
- 可能遗漏非主导但重要的运动

## 前提条件

### 对于方法 1（RMSD + Gyrate）
- 轨迹文件（.xtc/.trr）- PBC 修正是可选的
- RMSD 数据（rmsd.xvg）
- Gyrate 数据（gyrate.xvg）
- 拓扑文件（.tpr）

### 对于方法 2（PCA）
- 轨迹文件（.xtc/.trr）- PBC 修正是可选的
- PCA 特征向量（eigenvectors.trr）
- PC 投影（pc1.xvg、pc2.xvg）
- 拓扑文件（.tpr）

## 工作流程：方法 1（RMSD + Gyrate）

### 步骤 1：计算 RMSD

计算相对于参考结构的 RMSD：

```bash
echo -e "Backbone\nProtein\n" | gmx rms -s md.tpr -f md.xtc -o rmsd.xvg
```

- 第一个输入：选择用于对齐的参考（Backbone）
- 第二个输入：选择用于 RMSD 计算的组（Protein）

**输出**：`rmsd.xvg`（时间、RMSD）

### 步骤 2：计算 Gyrate

计算回转半径：

```bash
echo -e "Protein\n" | gmx gyrate -s md.tpr -f md.xtc -o gyrate.xvg
```

- 输入：选择蛋白质组

**输出**：`gyrate.xvg`（时间、Rg、Rg_x、Rg_y、Rg_z）
- 第 1 列：时间
- 第 2 列：总 Rg
- 第 3-5 列：X、Y、Z 分量

### 步骤 3：组合 RMSD 和 Gyrate 数据

将 RMSD 和总 Rg 组合成单个文件以进行 sham 分析：

```bash
dit xvg_combine -f rmsd.xvg gyrate.xvg -c 0,1 1 -o sham.xvg
```

**参数**：
- `-c 0,1 1`：从第一个文件提取列（时间、RMSD），从第二个文件提取列 1（Rg）
- 创建 `sham.xvg`，包含列：时间、RMSD、Rg

### 步骤 4：生成自由能景观

使用 gmx sham 计算 FEL：

```bash
gmx sham -tsham 310 -nlevels 100 -f sham.xvg -ls gibbs.xpm -g gibbs.log -lsh enthalpy.xpm -lss entropy.xpm
```

**参数**：
- `-tsham 310`：开尔文温度
- `-nlevels 100`：能级数量（分辨率）
- `-f sham.xvg`：包含反应坐标的输入文件
- `-ls gibbs.xpm`：输出吉布斯自由能景观
- `-g gibbs.log`：包含能量极小值的日志文件
- `-lsh enthalpy.xpm`：焓景观
- `-lss entropy.xpm`：熵景观

**输出文件**：
- `gibbs.xpm`：吉布斯自由能景观（XPM 格式）
- `gibbs.log`：包含能量极小值和索引的日志文件
- `enthalpy.xpm`：焓景观
- `entropy.xpm`：熵景观
- `bindex.ndx`：带有帧分配的索引文件
- `ener.xvg`：能量值

### 步骤 5：可视化 FEL

生成自由能景观的可视化：

```bash
dit xpm_show -f gibbs.xpm -m contour -cmap jet -o fel.png
```

**参数**：
- `-m contour`：等高线图模式
- `-cmap jet`：Jet 颜色图（蓝色=低能量，红色=高能量）
- 设置轴标签：`-x "RMSD (nm)" -y "Rg (nm)"`

### 步骤 6：提取最低能量构象

识别并提取对应于能量极小值的结构：

```bash
# 查看能量极小值
cat gibbs.log

# 示例输出：
# Minimum 0 at index 26 energy 7.589
# Minimum 1 at index 46 energy 7.589
# Minimum 2 at index 50 energy 7.589
# Minimum 3 at index 56 energy 7.589
# Minimum 4 at index 141 energy 5.803

# 检查对应于极小值的帧
cat bindex.ndx

# 示例输出：
# [ 26 ]
# 1274
# [ 46 ]
# 2
# [ 141 ]
# 4
# 1282
```

在特定时间提取构象：

```bash
# 提取帧 1274
echo -e "Protein\n" | gmx trjconv -f md.xtc -s md.tpr -b <time_at_frame_1274> -e <time_at_frame_1274> -o min1.pdb
```

要查找特定帧的时间，请检查 sham.xvg：

```bash
# 帧编号对应于 sham.xvg 中的行号
# 帧 N 的时间 = sham.xvg 中第 N 行第 0 列的值
```

## 工作流程：方法 2（主成分）

### 步骤 1：执行 PCA 分析

完成 PCA 分析（参见 [PCA 分析指南](pca-analysis.md)）以获得：
- `eigenvectors.trr`：特征向量
- `pc1.xvg`：投影到 PC1 上
- `pc2.xvg`：投影到 PC2 上

### 步骤 2：组合 PC 投影

组合 PC1 和 PC2 投影：

```bash
dit xvg_combine -f pc1.xvg pc2.xvg -c 0,1 1 -o sham.xvg
```

**参数**：
- `-c 0,1 1`：从第一个文件提取列（时间、PC1），从第二个文件提取列 1（PC2）
- 创建 `sham.xvg`，包含列：时间、PC1、PC2

### 步骤 3：生成自由能景观

与方法 1 步骤 4 相同：

```bash
gmx sham -tsham 310 -nlevels 100 -f sham.xvg -ls gibbs.xpm -g gibbs.log -lsh enthalpy.xpm -lss entropy.xpm
```

### 步骤 4：可视化 FEL

生成可视化：

```bash
dit xpm_show -f gibbs.xpm -m contour -cmap jet -o fel.png -x "PC1" -y "PC2" -t "Free Energy Landscape"
```

### 步骤 5：提取最低能量构象

与方法 1 步骤 5-6 相同。

## 输出文件

- **sham.xvg**：组合的反应坐标数据
- **gibbs.xpm**：吉布斯自由能景观（XPM 格式）
- **gibbs.log**：能量极小值和索引
- **enthalpy.xpm**：焓景观
- **entropy.xpm**：熵景观
- **bindex.ndx**：能量极小值的帧分配
- **ener.xvg**：能量值
- **fel.png**：自由能景观的可视化
- **min*.pdb**：能量极小值处的结构

## 解释指南

### 能量极小值

- **深极小值**：稳定的构象状态
- **浅极小值**：亚稳态或瞬时构象
- **多个极小值**：存在多个稳定状态
- **单个极小值**：单个主导构象

### 能垒

- **高能垒**：缓慢转变，稀有事件
- **低能垒**：快速转变，频繁互变
- **能垒高度**：与转变速率相关（k ∝ exp(-ΔG/kT)）

### 景观拓扑

- **漏斗形状**：蛋白质折叠景观
- **多个漏斗**：多个折叠途径
- **粗糙景观**：许多局部极小值，玻璃态行为
- **平滑景观**：少数能垒，快速动力学

### 反应坐标

**RMSD + Gyrate**：
- **高 RMSD、高 Rg**：展开或伸展构象
- **低 RMSD、低 Rg**：紧密、类天然构象
- **低 RMSD、高 Rg**：展开但有序
- **高 RMSD、低 Rg**：错误折叠或塌陷结构

**主成分**：
- **PC1 轴**：主导运动
- **PC2 轴**：第二主导运动
- **区域**：沿集体运动的不同构象状态

## 常见问题和解决方案

### 问题：FEL 显示不切实际的能垒（> 20 kT）

**可能原因**：
- 采样不足
- 时间范围选择不正确
- 温度不匹配

**解决方案**：
- 延长模拟时间
- 使用一致的时间范围
- 验证 sham 命令中的温度设置

### 问题：FEL 显示噪声或碎片化

**解决方案**：
- 增加能级数量（`-nlevels`）
- 延长模拟时间以获得更好的采样
- 在分析之前平滑数据

### 问题：多个极小值合并为单个盆地

**解决方案**：
- 增加分辨率（更多能级）
- 使用不同的反应坐标
- 检查极小值是否在物理上是不同的

### 问题：FEL 没有显示清晰的极小值

**解决方案**：
- 蛋白质可能是刚性的（单个构象）
- 检查模拟质量
- 考虑使用不同的反应坐标

## 提示和最佳实践

- **时间选择**：仅使用生产阶段，排除平衡
- **温度**：在 sham 命令中匹配模拟温度
- **分辨率**：调整 `-nlevels` 以获得适当的细节
- **采样**：确保对构象空间的充分采样
- **比较**：比较来自不同条件或突变体的 FEL
- **验证**：将能量极小值与实验结构相关联

## 高级分析

### 时间依赖性 FEL

为不同时间窗口计算 FEL：

```bash
# 前半部分的 FEL
dit xvg_combine -f rmsd_early.xvg gyrate_early.xvg -c 0,1 1 -o sham_early.xvg
gmx sham -tsham 310 -nlevels 100 -f sham_early.xvg -ls gibbs_early.xpm

# 后半部分的 FEL
dit xvg_combine -f rmsd_late.xvg gyrate_late.xvg -c 0,1 1 -o sham_late.xvg
gmx sham -tsham 310 -nlevels 100 -f sham_late.xvg -ls gibbs_late.xpm
```

比较 FEL 以识别构象景观的变化。

### 差异 FEL

计算两个 FEL 之间的差异：

```bash
dit xpm_diff -f gibbs_early.xpm gibbs_late.xpm -o gibbs_diff.xpm
```

识别自由能发生显著变化的区域。

### 途径分析

识别状态之间的最小能量路径：

```bash
# 使用专用工具（例如，字符串方法、推动弹性带）
- 需要额外处理
```

### 多维 FEL

使用额外的反应坐标扩展到 3D：

```bash
# 组合三个坐标
dit xvg_combine -f pc1.xvg pc2.xvg pc3.xvg -c 0,1 1 1 -o sham_3d.xvg
gmx sham -tsham 310 -nlevels 100 -f sham_3d.xvg -ls gibbs_3d.xpm
```

使用 plotly 可视化 3D FEL：

```bash
dit xpm_show -f gibbs_3d.xpm -m 3d -eg plotly -o fel_3d.html
```

## 相关分析

- **PCA**：为 FEL 提供集体运动
- **DCCM**：识别能量极小值中的相关运动
- **RMSF**：分析能量极小值处的灵活性
- **接触分析**：识别不同状态中的接触

## 可视化增强

### 自定义颜色方案

使用不同的颜色图：

```bash
# Jet（默认）
dit xpm_show -f gibbs.xpm -m contour -cmap jet -o fel_jet.png

# Viridis（感知均匀）
dit xpm_show -f gibbs.xpm -m contour -cmap viridis -o fel_viridis.png

# Plasma（替代）
dit xpm_show -f gibbs.xpm -m contour -cmap plasma -o fel_plasma.png
```

### 3D 可视化

生成 3D 表面图：

```bash
dit xpm_show -f gibbs.xpm -m 3d -eg plotly -o fel_3d.html
```

### 叠加构象

使用 PyMOL 或 VMD 从不同极小值提取并叠加结构。

### 能量剖面

沿特定路径提取能量剖面：

```bash
# 沿特定 PC1 值提取能量值
- 需要对 gibbs.xpm 进行自定义处理
```

## 与实验数据的比较

### 与 NMR 集合比较

将 NMR 结构叠加在 FEL 上以验证构象采样。

### 与晶体结构比较

在 FEL 上定位晶体结构以了解其能量背景。

### 与 SAXS 比较

从 FEL 导出的构象计算理论 SAXS 曲线。

## 参考文献

有关理论背景，请参阅：
- 蛋白质折叠理论和能量景观
- 过渡态理论和能垒跨越
- 构象动力学和自由能计算