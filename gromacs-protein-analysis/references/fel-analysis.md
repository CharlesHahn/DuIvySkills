# FEL 分析

## 概述

自由能景观（FEL, Free Energy Landscape）绘制蛋白质构象的能量表面，识别稳定状态、能垒和转变途径。使用 RMSD + Gyrate 或主成分（PC）作为反应坐标。适用于识别稳定构象、量化能垒、理解构象转变、分析折叠/展开。

## 工作流程

### 方法 1：RMSD + Gyrate

```bash
# 步骤 1：计算 RMSD
gmx rms -s md.tpr -f md.xtc -o rmsd.xvg

# 步骤 2：计算 Gyrate
gmx gyrate -s md.tpr -f md.xtc -o gyrate.xvg

# 步骤 3：组合数据
dit xvg_combine -f rmsd.xvg gyrate.xvg -c 0,1 1 -o sham.xvg

# 步骤 4：生成 FEL
gmx sham -tsham 310 -nlevels 100 -f sham.xvg -ls gibbs.xpm -g gibbs.log
```

### 方法 2：主成分

```bash
# 步骤 1：执行 PCA 分析
gmx covar -s md.tpr -f md.xtc -o eigenvalues.xvg -v eigenvectors.trr
gmx anaeig -s md.tpr -f md.xtc -v eigenvectors.trr -first 1 -last 1 -proj pc1.xvg
gmx anaeig -s md.tpr -f md.xtc -v eigenvectors.trr -first 2 -last 2 -proj pc2.xvg

# 步骤 2：组合 PC 投影
dit xvg_combine -f pc1.xvg pc2.xvg -c 0,1 1 -o sham.xvg

# 步骤 3：生成 FEL
gmx sham -tsham 310 -nlevels 100 -f sham.xvg -ls gibbs.xpm -g gibbs.log
```

### 参数选择

**重要参数**：
- `-tsham`：温度（K），需匹配模拟温度
- `-nlevels`：能级数量（分辨率，推荐 100）
- `-f`：包含反应坐标的输入文件
- `-ls`：输出吉布斯自由能景观
- `-g`：输出能量极小值日志

### 输出
- **gibbs.xpm**：自由能景观（XPM 格式）
- **gibbs.log**：能量极小值和索引
- **bindex.ndx**：能量极小值的帧分配

### 可视化
```bash
dit xpm_show -f gibbs.xpm -m contour -cmap jet -o fel.png
```

## 结果解释

### 能量极小值
- **深极小值**：稳定的构象状态
- **浅极小值**：亚稳态或瞬时构象
- **多个极小值**：存在多个稳定状态
- **单个极小值**：单个主导构象

### 能垒
- **高能垒**：缓慢转变，稀有事件
- **低能垒**：快速转变，频繁互变
- **能垒高度**：与转变速率相关（k ∝ exp(-ΔG/kT)）

### 反应坐标（RMSD + Gyrate）
- **高 RMSD、高 Rg**：展开或伸展构象
- **低 RMSD、低 Rg**：紧密、类天然构象
- **低 RMSD、高 Rg**：展开但有序
- **高 RMSD、低 Rg**：错误折叠或塌陷结构

### 景观拓扑
- **漏斗形状**：蛋白质折叠景观
- **多个漏斗**：多个折叠途径
- **粗糙景观**：许多局部极小值
- **平滑景观**：少数能垒，快速动力学

## 常见问题

**Q: FEL 显示不切实际的能垒（> 20 kT）？**  
A: 可能采样不足或时间范围不正确，延长模拟时间或检查温度设置。

**Q: FEL 显示噪声或碎片化？**  
A: 增加能级数量（`-nlevels`）或延长模拟时间。

**Q: 没有显示清晰的极小值？**  
A: 蛋白质可能刚性（单个构象）或模拟质量有问题。