# 周期性校正

## 概述

周期性边界条件（PBC, Periodic Boundary Conditions）修正消除分子跨越模拟盒边界产生的问题。仅当轨迹质量问题（如 RMSD 突然跳跃）时应用，非必须步骤。适用于分子跨越盒边界、RMSD/RMSF 显示突变、可视化显示分子碎片化。

## 工作流程

### 输入文件
- `-s`：结构文件（.tpr）
- `-f`：轨迹文件（.xtc）

### 命令
```bash
# 步骤 1：居中系统
echo -e "Protein\nProtein\n" | gmx trjconv -s md.tpr -f md.xtc -o center.xtc -pbc atom -center

# 步骤 2：保持分子完整性
echo -e "Protein\n" | gmx trjconv -s md.tpr -f center.xtc -o mol.xtc -pbc mol -ur compact

# 步骤 3：消除整体平移和旋转（可选）
echo -e "Backbone\nProtein\n" | gmx trjconv -s md.tpr -f mol.xtc -o fit.xtc -fit rot+trans

# 步骤 4：更新拓扑文件
echo -e "Protein\n" | gmx convert-tpr -s md.tpr -o fit.tpr
```

### 参数选择

**原子组选择**：
- **Protein**：蛋白质原子（推荐）
- **Protein_Lig**：蛋白质-配体复合物
- **Backbone**：用于拟合

**重要参数**：
- `-pbc atom`：原子级别的 PBC 处理
- `-pbc mol`：保持分子完整性
- `-ur compact`：紧凑单位元
- `-center`：居中系统
- `-fit rot+trans`：消除平移和旋转

### 输出
- **fit.xtc**：修正后的轨迹文件
- **fit.tpr**：修正后的拓扑文件

### 验证
```bash
# 提取样本结构进行目视检查
echo -e "Protein\n" | gmx trjconv -s md.tpr -f fit.xtc -o fit.pdb -dt 1000
```

## 结果解释

### 修正效果
- **分子完整**：无碎片化
- **居中**：蛋白质/配体保持在盒中心
- **无跳跃**：RMSD 平滑无突变

### 何时需要
- **需要**：RMSD 显示突变、分子跨越盒边界、可视化显示碎片化
- **不需要**：RMSD 平滑、分子完整

## 常见问题

**Q: 修正后分子仍然碎片化？**  
A: 检查原子组选择，确保包含所有应保持完整的原子。

**Q: 轨迹仍然显示 PBC周期性问题？**  
A: 验证中心原子选择，中心原子应接近系统几何中心。

**Q: tpr 和 xtc 原子数不匹配？**  
A: 在任何组选择更改后运行 `gmx convert-tpr` 更新拓扑文件。

**Q: 修正后 RMSD 仍然显示跳跃？**  
A: 检查 RMSD 计算的参考结构，考虑使用平均结构。

