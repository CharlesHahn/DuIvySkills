# PCA 分析

## 概述

主成分分析（PCA, Principal Component Analysis）通过将蛋白质运动分解为主成分来识别集体运动和主要构象变化，降低维度同时保留重要动力学信息。适用于识别集体运动、提取主要模式、降低维度、分析结构域运动。

## 工作流程

### 输入文件
- `-s`：结构文件（.tpr）
- `-f`：轨迹文件（.xtc/.trr）

### 命令
```bash
# 步骤 1：计算协方差矩阵和特征向量
gmx covar -s md.tpr -f md.xtc -o eigenvalues.xvg -v eigenvectors.trr

# 步骤 2：投影到主成分
gmx anaeig -s md.tpr -f md.xtc -v eigenvectors.trr -first 1 -last 1 -proj pc1.xvg
gmx anaeig -s md.tpr -f md.xtc -v eigenvectors.trr -first 2 -last 2 -proj pc2.xvg

# 步骤 3：提取极端构象
gmx anaeig -s md.tpr -f md.xtc -v eigenvectors.trr -first 1 -last 1 -extreme pc1_extreme.pdb
```

### 参数选择

**原子组选择**（两次交互输入）：
- **C-alpha**：Cα 原子，最常用（推荐）
- **Backbone**：主链原子
- **Protein**：所有蛋白质原子

**重要参数**：
- `-o`：特征值输出文件
- `-v`：特征向量输出文件（.trr 格式）
- `-first/-last`：选择主成分编号
- `-proj`：投影到指定主成分
- `-extreme`：提取极端构象
- `-average`：提取平均结构

### 输出
- **eigenvalues.xvg**：特征值（PC 编号-特征值）
- **eigenvectors.trr**：特征向量
- **pc1.xvg, pc2.xvg**：投影到主成分
- **pc1_extreme.pdb**：极端构象

### 可视化
```bash
# 特征值谱
dit xvg_show -f eigenvalues.xvg -x "PC Number" -y "Eigenvalue"

# PC1 投影
dit xvg_show -f pc1.xvg -x "Time (ps)" -y "PC1 Coordinate"

# PC1 vs PC2 散点图
dit xvg_show_scatter -f pc1.xvg pc2.xvg -x "PC1" -y "PC2"
```

## 结果解释

### 特征值分析
- **第一个 PC**：最大特征值，最主要的运动
- **累积方差**：前 N 个特征值之和 / 总方差
- **贡献百分比**：单个 PC 贡献 = 特征值 / 总方差 × 100%
- **典型值**：前 3 个 PC 通常捕获 50-80% 的总运动

### PC 投影
- **振幅**：投影范围表示运动大小
- **时间尺度**：振荡频率与运动时间尺度相关
- **转变**：跳跃或位移表示构象变化

### 生物学解释
- **结构域运动**：大规模结构域间运动
- **环灵活性**：环和末端的局部运动
- **变构途径**：表示变构通讯的相关运动
- **功能运动**：与活性相关的运动（铰链弯曲、通道开放）

## 常见问题

**Q: 特征值显示均匀分布？**  
A: 可能模拟时间不足或蛋白质刚性，延长模拟时间。

**Q: PC 投影没有清晰模式？**  
A: 检查轨迹拟合和原子选择一致性。

**Q: 前几个 PC 没有捕获显著方差？**  
A: 可能蛋白质刚性或运动均匀分布，考虑使用更多 PC。