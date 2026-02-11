# RMSF 分析

## 概述

均方根涨落（RMSF, Root Mean Square Fluctuation）测量每个原子相对于其平均位置的平均偏差，提供残基水平的灵活性信息，识别刚性和柔性区域。适用于识别柔性环、评估局部稳定性、比较不同模拟、验证与实验 B 因子的一致性。

## 工作流程

### 输入文件
- `-s`：结构文件（.tpr）
- `-f`：轨迹文件（.xtc/.trr）

### 命令
```bash
gmx rmsf -s md.tpr -f md.xtc -o rmsf.xvg -res
```

### 参数选择

**原子组选择**：
- **C-alpha**：Cα 原子，最常用（推荐）
- **Backbone**：主链原子
- **Protein**：所有蛋白质原子

**重要参数**：
- `-o`：输出文件（.xvg 格式，残基-RMSF）
- `-res`：计算每个残基的 RMSF
- `-bf`：输出 B 因子 PDB 文件
- `-b`：起始时间
- `-e`：结束时间

### 输出
- **rmsf.xvg**：残基-RMSF 数据
  - 第 1 列：残基编号
  - 第 2 列：RMSF 值（单位：nm）
- **bfactor.pdb**：带有 B 因子的结构文件

### 可视化
```bash
# 残基-RMSF 曲线
dit xvg_show -f rmsf.xvg -x "Residue Number" -y "RMSF (nm)"

# 结构上 B 因子着色（PyMOL）
load bfactor.pdb; show cartoon; spectrum b, blue_white_red
```

## 结果解释

### RMSF 值范围

| RMSF (nm) | 含义 |
|-----------|------|
| < 0.05 | 非常刚性，结构化区域 |
| 0.05-0.1 | 刚性，稳定二级结构 |
| 0.1-0.2 | 中等灵活性，柔性环或末端 |
| 0.2-0.3 | 高灵活性，非结构化区域 |
| > 0.3 | 非常高灵活性，可能无序 |

### 二级结构模式
- **α-螺旋 (α-helix)**：低 RMSF（0.05-0.1 nm）
- **β-折叠 (β-sheet)**：低到中等 RMSF（0.05-0.15 nm）
- **环**：高 RMSF（0.15-0.3 nm）
- **N 端和 C 端 (N-terminus, C-terminus)**：极高 RMSF（0.2-0.5 nm）

### 功能区域
- **活性位点**：中等 RMSF（0.08-0.15 nm）
- **结合位点**：可变 RMSF，诱导契合
- **变构位点**：柔性 RMSF（0.15-0.25 nm）

## 常见问题

**Q: RMSF 值都很低？**  
A: 可能模拟时间不足或蛋白质刚性，延长模拟时间。

**Q: RMSF 值异常高？**  
A: 检查蛋白质是否展开或存在PBC周期性问题。

**Q: 与实验 B 因子相关性差？**  
A: 可能是晶体堆积效应或力场问题，模拟中的 B 因子通常略高于晶体。
