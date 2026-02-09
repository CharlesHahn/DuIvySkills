# 周期性边界条件（PBC）修正

## 概述

周期性边界条件（PBC）修正消除了分子跨越模拟盒边界产生的伪影。当观察到轨迹质量问题（例如 RMSD 突然跳跃）时，可能需要进行此修正。许多分析在没有 PBC 修正的情况下也能正常工作。

## 何时应用 PBC 修正

- 分子（蛋白质或配体）跨越盒边界
- RMSD/RMSF 图显示突然跳跃或不连续性
- 可视化显示分子碎片化

注意：PBC 修正是可选的。仅在观察到上述问题时应用。分析问题可能有 PBC 伪影以外的其他原因。

## 前提条件

- 输入轨迹文件（.xtc）
- 拓扑文件（.tpr）
- 包含适当原子组的索引文件（.ndx）

## 工作流程

### 步骤 1：准备索引文件

创建或修改索引文件以包含必要的组：

```bash
# 创建基本索引文件
gmx make_ndx -f md.tpr -o index.ndx

# 如果需要，添加蛋白质-配体组
echo -e "Protein | Ligand\nq\n" | gmx make_ndx -f md.tpr -o index.ndx
```

### 步骤 2：找到中心原子

找到最接近几何中心的原子以进行适当的居中：

```bash
# 使用 DuIvyTools 查找中心原子
echo -e "Protein\n" | dit find_center -f npt.gro index.ndx
```

输出显示最接近中心的原子编号。记录此原子编号。

### 步骤 3：将中心组添加到索引文件

手动将中心组添加到索引文件中，包含识别的原子：

```bash
# 将中心组添加到 index.ndx
echo -e "\n[ center ]\n<atom_number>\n" >> index.ndx
```

将 `<atom_number>` 替换为步骤 2 中的实际原子编号。

### 步骤 4：居中轨迹

使用中心原子居中系统：

```bash
echo -e "center\nProtein_Lig\n" | gmx trjconv -s md.tpr -f md.xtc -o center.xtc -n index.ndx -pbc atom -center
```

- 第一个输入：选择中心组
- 第二个输入：选择输出组（Protein_Lig）

### 步骤 5：确保分子完整性

保持分子在周期性边界中的完整性：

```bash
echo -e "Protein_Lig\n" | gmx trjconv -s md.tpr -f center.xtc -o mol.xtc -n index.ndx -pbc mol -ur compact
```

- 输入：选择要保持完整的组（Protein_Lig）

### 步骤 6：消除整体平移和旋转

拟合轨迹以消除整体运动（可选但建议）：

```bash
echo -e "Backbone\nProtein_Lig\n" | gmx trjconv -s md.tpr -f mol.xtc -o fit.xtc -n index.ndx -fit rot+trans
```

- 第一个输入：选择用于拟合的参考结构（Backbone）
- 第二个输入：选择输出组（Protein_Lig）

### 步骤 7：更新拓扑文件

更新拓扑文件以匹配修正后的轨迹原子数：

```bash
echo -e "Protein_Lig\n" | gmx convert-tpr -s md.tpr -o fit.tpr -n index.ndx
```

- 输入：选择步骤 6 中使用的相同组

### 步骤 8：验证修正

目视检查修正后的轨迹：

```bash
# 转换为 PDB 以进行目视检查
echo -e "Protein_Lig\n" | gmx trjconv -s md.tpr -f fit.xtc -o fit.pdb -dt 1000 -n index.ndx
```

使用 PyMOL 或 VMD 检查 PDB 文件。检查以下内容：
- 分子完整且未碎片化
- 蛋白质/配体保持居中
- 没有分子跨越盒边界

## 输出文件

- **center.xtc**：居中的轨迹
- **mol.xtc**：具有完整分子的轨迹
- **fit.xtc**：最终修正的轨迹（无平移/旋转）
- **fit.tpr**：修正后的拓扑文件
- **fit.pdb**：用于目视检查的样本结构

## 常见问题和解决方案

### 问题：修正后分子仍然碎片化

**解决方案**：检查步骤 5 中是否选择了正确的原子组。确保该组包含所有应保持完整的原子。

### 问题：轨迹仍然显示 PBC 伪影

**解决方案**：验证步骤 2 中的中心原子选择。中心原子应接近系统的几何中心。

### 问题：tpr 和 xtc 原子数不匹配

**解决方案**：在任何组选择更改后，始终运行步骤 7 以更新拓扑文件。

### 问题：修正后 RMSD 仍然显示跳跃

**解决方案**：检查 RMSD 计算的参考结构是否合适。考虑使用平均结构作为参考。

## 提示和最佳实践

- **使用适当的时间间隔**：转换为 PDB 进行检查时，使用 `-dt` 参数限制输出大小
- **记录原子选择**：跟踪每个索引组中的原子
- **一致的组选择**：对所有相关分析使用相同的组
- **目视检查**：在下游分析之前始终目视检查修正后的轨迹
- **永远不要覆盖现有文件**：使用唯一的输出文件名以保留原始数据

## 相关分析

PBC 修正可能受益于 RMSD 显示突然跳跃的分析：
- RMSD/RMSF 分析
- DCCM 分析
- PCA 分析
- FEL 分析
- 任何需要对齐轨迹的结构分析

注意：PBC 修正是可选的，并非总是必要的。仅在观察到轨迹质量问题时应用。

## 替代方法

### 替代方法 1：单步修正

在一个步骤中合并居中和分子完整性：

```bash
echo -e "center\nProtein_Lig\n" | gmx trjconv -s md.tpr -f md.xtc -o fit.xtc -n index.ndx -pbc mol -center
```

### 替代方法 2：使用蛋白质中心

使用蛋白质质心而不是特定原子：

```bash
echo -e "Protein\nProtein_Lig\n" | gmx trjconv -s md.tpr -f md.xtc -o fit.xtc -n index.ndx -pbc mol -center
```

## 可视化

使用 DuIvyTools 技能可视化修正后的轨迹：

```bash
# 绘制 RMSD 以验证修正
dit xvg_show -f rmsd_corrected.xvg -x "Time (ns)" -y "RMSD (nm)"
```

修正后的 RMSD 应显示平滑进展，没有突然跳跃。