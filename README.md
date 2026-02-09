# DuIvySkills

GROMACS 分子动力学模拟的 iFlow CLI 技能集合，适用于 **Claude Code** 和 **iFlow CLI**。

## 安装技能

### 方式 1：复制到技能目录

将本仓库中的技能目录复制到 Claude Code 或 iFlow CLI 的技能目录：

```bash
# 克隆仓库
git clone https://github.com/CharlesHahn/DuIvySkills.git
cd DuIvySkills

# 复制技能到指定目录（根据您的安装位置调整）
# Claude Code 技能目录示例：
cp -r duivytools-skills ~/.claude-code/skills/
cp -r gromacs-skills ~/.claude-code/skills/
cp -r gromacs-protein-analysis ~/.claude-code/skills/

# 或 iFlow CLI 技能目录示例：
cp -r duivytools-skills ~/.iflow/skills/
cp -r gromacs-skills ~/.iflow/skills/
cp -r gromacs-protein-analysis ~/.iflow/skills/
```

### 方式 2：通过配置文件指定

在 Claude Code 或 iFlow CLI 的配置文件中添加技能路径：

```yaml
skills:
  - path: /path/to/DuIvySkills/duivytools-skills
  - path: /path/to/DuIvySkills/gromacs-skills
  - path: /path/to/DuIvySkills/gromacs-protein-analysis
```

### 验证安装

安装完成后，在 Claude Code 或 iFlow CLI 中测试：

```
"帮我使用 DuIvyTools 绘制 rmsd.xvg 文件"
```

如果系统能够识别并调用相关技能，说明安装成功。

## 项目简介

DuIvySkills 提供三个核心技能包，覆盖 GROMACS 分子动力学模拟的完整工作流：

- **duivytools-skills** - GROMACS 模拟结果的可视化工具
- **gromacs-skills** - GROMACS 软件的命令参考和工作流指南
- **gromacs-protein-analysis** - 蛋白质动力学分析的专门工作流

## 使用方法

### 安装 DuIvyTools

```bash
pip install DuIvyTools
```

或使用清华镜像（更快）：

```bash
pip install DuIvyTools -i https://pypi.tuna.tsinghua.edu.cn/simple
```

### 在 iFlow CLI 中使用

当您需要执行 GROMACS 相关任务时，iFlow CLI 会自动调用相应的技能：

1. **执行 GROMACS 命令**：调用 `gromacs-skills`
2. **分析蛋白质模拟结果**：调用 `gromacs-protein-analysis`
3. **可视化数据**：调用 `duivytools-skills`

### 示例

```bash
# 设置并运行 GROMACS 模拟
iFlow "帮我设置一个蛋白质的 GROMACS 模拟"

# 分析蛋白质动力学
iFlow "分析我的蛋白质 MD 模拟结果，包括 RMSD、PCA 和自由能景观"

# 可视化结果
iFlow "绘制我的 RMSD 和 DCCM 图"
```

## 技能包详情

### duivytools-skills

提供约 30 个命令用于处理 XVG、XPM、NDX 等 GROMACS 输出文件：

- 可视化 RMSD、RMSF、能量、氢键等数据
- 绘制 DCCM、FEL、DSSP 矩阵热图
- 操作 NDX 索引文件
- 统计分析和批量处理

### gromacs-skills

GROMACS 完整命令参考，涵盖：

- 拓扑和结构处理
- 模拟设置和运行
- 能量和轨迹分析
- 索引和原子选择

### gromacs-protein-analysis

蛋白质分析的五种核心工作流：

1. PBC 修正
2. DCCM 分析（动态互相关矩阵）
3. RDCM 分析（残基距离接触矩阵）
4. PCA 分析（主成分分析）
5. FEL 分析（自由能景观）

## 文档结构

```
DuIvySkills/
├── duivytools-skills/          # DuIvyTools 技能包
├── gromacs-skills/            # GROMACS 技能包
├── gromacs-protein-analysis/  # 蛋白质分析技能包
├── AGENTS.md                  # 项目上下文
└── README.md                  # 本文件
```

## 相关链接

- GitHub: https://github.com/CharlesHahn/DuIvySkills
- DuIvyTools: https://github.com/CharlesHahn/DuIvyTools
- GROMACS 文档: https://manual.gromacs.org/

## 许可证

GNU General Public License v3

---

# DuIvySkills

An iFlow CLI skill collection for GROMACS molecular dynamics simulations, compatible with **Claude Code** and **iFlow CLI**.

## Installing Skills

### Method 1: Copy to Skills Directory

Copy the skill directories from this repository to the Claude Code or iFlow CLI skills directory:

```bash
# Clone repository
git clone https://github.com/CharlesHahn/DuIvySkills.git
cd DuIvySkills

# Copy skills to specified directory (adjust according to your installation location)
# Claude Code skills directory example:
cp -r duivytools-skills ~/.claude-code/skills/
cp -r gromacs-skills ~/.claude-code/skills/
cp -r gromacs-protein-analysis ~/.claude-code/skills/

# Or iFlow CLI skills directory example:
cp -r duivytools-skills ~/.iflow/skills/
cp -r gromacs-skills ~/.iflow/skills/
cp -r gromacs-protein-analysis ~/.iflow/skills/
```

### Method 2: Specify via Configuration File

Add skill paths to the configuration file of Claude Code or iFlow CLI:

```yaml
skills:
  - path: /path/to/DuIvySkills/duivytools-skills
  - path: /path/to/DuIvySkills/gromacs-skills
  - path: /path/to/DuIvySkills/gromacs-protein-analysis
```

### Verify Installation

After installation, test in Claude Code or iFlow CLI:

```
"Help me use DuIvyTools to plot rmsd.xvg file"
```

If the system can recognize and invoke the relevant skills, the installation is successful.

## Project Overview

DuIvySkills provides three core skill packages covering the complete workflow of GROMACS molecular dynamics simulations:

- **duivytools-skills** - Visualization tool for GROMACS simulation results
- **gromacs-skills** - Command reference and workflow guide for GROMACS software
- **gromacs-protein-analysis** - Specialized workflow for protein dynamics analysis

## Usage

### Install DuIvyTools

```bash
pip install DuIvyTools
```

Or use Tsinghua mirror (faster):

```bash
pip install DuIvyTools -i https://pypi.tuna.tsinghua.edu.cn/simple
```

### Using in iFlow CLI

When you need to perform GROMACS-related tasks, iFlow CLI will automatically invoke the corresponding skills:

1. **Execute GROMACS commands**: Invoke `gromacs-skills`
2. **Analyze protein simulation results**: Invoke `gromacs-protein-analysis`
3. **Visualize data**: Invoke `duivytools-skills`

### Examples

```bash
# Set up and run GROMACS simulation
iFlow "帮我设置一个蛋白质的 GROMACS 模拟"

# Analyze protein dynamics
iFlow "分析我的蛋白质 MD 模拟结果，包括 RMSD、PCA 和自由能景观"

# Visualize results
iFlow "绘制我的 RMSD 和 DCCM 图"
```

## Skill Package Details

### duivytools-skills

Provides approximately 30 commands for processing GROMACS output files such as XVG, XPM, NDX:

- Visualize RMSD, RMSF, energy, hydrogen bonds and other data
- Plot DCCM, FEL, DSSP matrix heatmaps
- Manipulate NDX index files
- Statistical analysis and batch processing

### gromacs-skills

Complete GROMACS command reference, covering:

- Topology and structure processing
- Simulation setup and execution
- Energy and trajectory analysis
- Index and atom selection

### gromacs-protein-analysis

Five core workflows for protein analysis:

1. PBC correction
2. DCCM analysis (Dynamics Cross-Correlation Matrix)
3. RDCM analysis (Residue Distance Contact Matrix)
4. PCA analysis (Principal Component Analysis)
5. FEL analysis (Free Energy Landscape)

## Documentation Structure

```
DuIvySkills/
├── duivytools-skills/          # DuIvyTools skill package
├── gromacs-skills/            # GROMACS skill package
├── gromacs-protein-analysis/  # Protein analysis skill package
├── AGENTS.md                  # Project context
└── README.md                  # This file
```

## Related Links

- GitHub: https://github.com/CharlesHahn/DuIvySkills
- DuIvyTools: https://github.com/CharlesHahn/DuIvyTools
- GROMACS documentation: https://manual.gromacs.org/

## License

GNU General Public License v3

