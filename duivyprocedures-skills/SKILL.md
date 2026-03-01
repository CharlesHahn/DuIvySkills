---
name: duivyprocedures-skills
description: "YAML-driven batch analysis tool for MD simulations. Use when Claude needs to: (1) Run multiple analyses in one command via YAML config, (2) Batch process multiple trajectories, (3) Auto-generate figures and data files. NOT for single-file visualization (use duivytools-skills instead)."
---

# DuIvyProcedures (DIP)

DuIvyProcedures 是一个 YAML 配置驱动的 MD 模拟批量分析工具。**核心价值：一条命令完成多种分析，自动出图出数据。**

## 快速判断

| 需求 | 使用 |
|------|------|
| 批量处理多条轨迹 | DIP ✓ |
| 一次运行多个分析 | DIP ✓ |
| 快速可视化单个 XVG/XPM | DIT (duivytools-skills) |
| 手动执行 GROMACS 命令 | gromacs-skills |
| 了解分析原理 | gromacs-protein-analysis |

## 安装

```bash
pip install DuIvyProcedures
```

**注意**：DIP 需要购买授权码使用，具体请联系微信公众号【杜艾维】。

## 命令行用法

### 生成配置文件

```bash
# 生成包含所有模块的配置模板
dip conf -o dip.yaml

# 指定模块和路径生成配置
dip conf -o dip.yaml -t gmx_RMSD DCCM -d MD0
```

### 运行分析

```bash
dip run -f dip.yaml
```

DIP 会按路径和模块顺序依次执行分析，结果保存在对应目录下。

## YAML 配置详解

### 完整配置结构

```yaml
Path:
  - MD0
  - MD1

Conf:
  gmx: gmx
  xtc: md.xtc
  tpr: md.tpr
  ndx: index.ndx
  fig: png        # 可选，输出图片格式：png/pdf/svg

Tasks:
  - gmx_RMSD:
      fit_group: Backbone
      calc_group: Protein
      rmsd_matrix: no
      gmx_parm:
        tu: ns

  - gmx_PCA:
      group: C-alpha

  - gmx_FEL:
      inputfile: ../gmx_PCA/pc12.xvg
      find_minimum: true
      minimum_num: 3
```

### 配置段说明

| 段 | 必需 | 说明 |
|------|------|------|
| `Path` | 是 | 分析路径列表，各路径下的文件名需一致 |
| `Conf` | 是 | 配置：gmx 路径、轨迹/拓扑/索引文件名、图片格式 |
| `Tasks` | 是 | 分析任务列表 |

### Conf 段参数

| 参数 | 说明 |
|------|------|
| `gmx` | GROMACS 可执行文件路径 |
| `xtc` | 轨迹文件名 |
| `tpr` | 拓扑文件名 |
| `ndx` | 索引文件名（可选） |
| `fig` | 图片格式：`png`/`pdf`/`svg`（v1.0.3+） |

### 模块通用结构

```yaml
- ModuleName:
    mkdir: OutputDir        # 输出目录（可选，默认为模块名）
    # ... 模块特定参数 ...
    gmx_parm:               # 仅 gmx_ 模块可用
      b: 0
      e: 10000
      tu: ns
    frame_start: 1000       # 帧选择（仅纯 Python 模块）
    frame_end: 5001
    frame_step: 10
```

### 模块间数据传递

后续模块可使用前序模块的输出：

```yaml
Tasks:
  - gmx_PCA:
      group: C-alpha

  - gmx_FEL:
      inputfile: ../gmx_PCA/pc12.xvg   # 相对路径
```

**路径规则**：每个模块的输出目录与配置文件同级，使用 `../ModuleName/` 访问。

### 原子选择语法

**GROMACS 模块（gmx_ 前缀）**：使用索引组名
- `Protein`、`C-alpha`、`Backbone` 等 GROMACS 预定义组
- 自定义组名（需在 .ndx 中定义）
- **组名必须以英文开头，不能以数字开头**（如 `6Lig` 会被识别为第 6 个组）

**纯 Python 模块**：使用 MDAnalysis 语法
- `protein` - 蛋白质
- `protein and name CA` - α-碳
- `resname LYS` - 赖氨酸残基
- `resname *ZIN` - 名称匹配 ZIN 的残基

语法参考：https://userguide.mdanalysis.org/2.7.0/selections.html

## 模块选择：gmx_ vs 纯Python

| 场景 | 选择 | 原因 |
|------|------|------|
| 已有 GROMACS 环境 | `gmx_` 模块 | 更快，与 GROMACS 生态兼容 |
| 无 GROMACS 或跨平台 | 纯 Python 模块 | 仅依赖 MDAnalysis，支持 AMBER 等格式 |
| 需要灵活原子选择 | 纯 Python 模块 | MDAnalysis 语法更灵活 |

## 帧选择参数（仅纯 Python 模块）

```yaml
frame_start: 1000    # 起始帧
frame_end: 5001      # 结束帧（不包含）
frame_step: 10       # 步长
```

## 前置处理

使用 DIP 前需自行完成：

1. **周期性校正**：DIP 不会自动校正轨迹，需保证分子完整性
2. **原子数对应**：拓扑文件和轨迹文件原子数必须一一对应
3. **索引文件**：gmx_ 模块需要正确的 .ndx 文件
4. **GROMACS 版本**：建议 2019-2023，2024 可能有兼容问题
5. **DSSP 安装**：GROMACS 2022 及以下运行 `gmx_DSSP` 需安装 DSSP 3.0
6. **RDKit 安装**：`PiStacking` 自动寻找芳香环需要 RDKit

## 重要注意事项

| 模块 | 注意事项 |
|------|----------|
| `RDCM` | `type_select: min` 计算非常慢，大体系需设置帧选择参数避免内存不足 |
| `MSM` | Demo 性质，需谨慎用于研究；参数不合适会报错，需反复调整 |
| `gmx_dPCA` | 使用了一些技巧，不能保证所有情况成功执行，建议检查结果 |
| `gmx_Hbond` | 组名不能以数字开头（如 `1ZIN` 应改为 `ZIN1`） |
| `SaltBridge` | 不同力场需修改原子命名（`NH3_atomnames`, `COO_atomnames`） |
| `gmx_cluster` | `dt` 设置过小会导致帧数多、计算量大 |
| `gmx_FEL` / `FEL` | 需要三列数据文件（时间, 数据1, 数据2），通常来自 PCA 或 dPCA |

## 模块速查

### 结构稳定性

| 模块 | 功能 | 关键参数 |
|------|------|----------|
| `gmx_RMSD` / `RMSD` | 均方根偏差 | `fit_group`, `calc_group` |
| `gmx_RMSF` / `RMSF` | 均方根涨落 | `group` / `calc_group` |
| `gmx_Gyrate` / `Gyrate` | 回转半径 | `calc_group` |
| `gmx_SASA` | 溶剂可及表面积 | `calc_group` |

### 构象分析

| 模块 | 功能 | 关键参数 |
|------|------|----------|
| `gmx_DCCM` / `DCCM` | 动态互相关矩阵 | `group` / `atom_selection` |
| `gmx_PCA` / `PCA` | 主成分分析 | `group` / `atom_selection` |
| `gmx_FEL` / `FEL` | 自由能景观 | `inputfile`（三列数据） |
| `gmx_dPCA` | 二面角 PCA | `group`（可能不稳定） |
| `RDCM` | 残基距离接触矩阵 | `atom_selection`, `type_select`（`min` 很慢） |

### 相互作用

| 模块 | 功能 | 关键参数 |
|------|------|----------|
| `gmx_Hbond` / `Hbond` | 氢键 | `group1`/`group2` 或 `donor_group`/`acceptor_group` |
| `SaltBridge` | 盐桥 | `byIndex`, `group` |
| `PiStacking` | Pi-堆叠 | `group1`, `group2` |
| `PiCation` | Pi-阳离子 | `group1`, `group2` |
| `Hydrophobic Contact` | 疏水接触 | `group1`, `group2` |
| `RDF` | 径向分布函数 | `center_group`, `calc_group` |

### 结构特征

| 模块 | 功能 | 关键参数 |
|------|------|----------|
| `gmx_DSSP` | 二级结构 | `group` |
| `gmx_cluster` | 聚类分析 | `fit_group`, `calc_group`, `method`, `cutoff` |
| `gmx_Mdmat` | 距离矩阵 | `group` |

### 密度与分布

| 模块 | 功能 | 关键参数 |
|------|------|----------|
| `gmx_Density` / `Density` | 一维密度分布 | `calc_group` / `atom_selection` |
| `DensityMap` | 三维密度映射 | `groups`, `byType`, `grid_bin` |
| `SPM` | 空间映射 | `atom_selection` |

### 降维与动力学

| 模块 | 功能 | 关键参数 |
|------|------|----------|
| `MSM` | 马尔可夫状态模型（Demo） | `atom_selection`, `lag4tica`, `lag4MSM` |
| `tICA` | 时间滞后独立成分分析 | `atom_selection`, `target`, `lag` |
| `tSNE` | t-SNE 降维 | `atom_selection`, `target` |
| `UMAP` | UMAP 降维 | `atom_selection`, `n_neighbors`, `min_dist` |

### 用户自定义

| 模块 | 功能 |
|------|------|
| `User_Mod` | 自定义分析模块 |

## 详细参考文档

| 文档 | 内容 |
|------|------|
| [gmx_RMSD.md](references/gmx_RMSD.md) | RMSD 分析 |
| [gmx_RMSF.md](references/gmx_RMSF.md) | RMSF 分析 |
| [gmx_Gyrate.md](references/gmx_Gyrate.md) | 回转半径 |
| [gmx_SASA.md](references/gmx_SASA.md) | 溶剂可及表面积 |
| [gmx_DCCM.md](references/gmx_DCCM.md) | 动态互相关矩阵 |
| [gmx_PCA.md](references/gmx_PCA.md) | 主成分分析 |
| [gmx_FEL.md](references/gmx_FEL.md) | 自由能景观 |
| [gmx_Hbond.md](references/gmx_Hbond.md) | 氢键分析 |
| [gmx_DSSP.md](references/gmx_DSSP.md) | 二级结构 |
| [gmx_cluster.md](references/gmx_cluster.md) | 聚类分析 |
| [gmx_Mdmat.md](references/gmx_Mdmat.md) | 距离矩阵 |
| [gmx_Density.md](references/gmx_Density.md) | 密度分析 |
| [gmx_dPCA.md](references/gmx_dPCA.md) | 二面角 PCA |
| [RMSD.md](references/RMSD.md) | RMSD（纯Python） |
| [RMSF.md](references/RMSF.md) | RMSF（纯Python） |
| [Gyrate.md](references/Gyrate.md) | 回转半径（纯Python） |
| [DCCM.md](references/DCCM.md) | DCCM（纯Python） |
| [PCA.md](references/PCA.md) | PCA（纯Python） |
| [FEL.md](references/FEL.md) | FEL（纯Python） |
| [Hbond.md](references/Hbond.md) | 氢键（纯Python） |
| [RDCM.md](references/RDCM.md) | 残基距离接触矩阵 |
| [RDF.md](references/RDF.md) | 径向分布函数 |
| [Density.md](references/Density.md) | 密度（纯Python） |
| [DensityMap.md](references/DensityMap.md) | 密度映射 |
| [MSM.md](references/MSM.md) | 马尔可夫状态模型 |
| [tICA.md](references/tICA.md) | tICA 降维 |
| [tSNE.md](references/tSNE.md) | t-SNE 降维 |
| [UMAP.md](references/UMAP.md) | UMAP 降维 |
| [SaltBridge.md](references/SaltBridge.md) | 盐桥分析 |
| [PiStacking.md](references/PiStacking.md) | Pi-堆叠 |
| [PiCation.md](references/PiCation.md) | Pi-阳离子 |
| [Hydrophobic_Contact.md](references/Hydrophobic_Contact.md) | 疏水接触 |
| [SPM.md](references/SPM.md) | 空间映射 |
| [User_Mod.md](references/User_Mod.md) | 用户自定义模块 |

## 官方文档

https://duivyprocedures-docs.readthedocs.io/