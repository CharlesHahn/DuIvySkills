# gmx_PCA

依赖 GROMACS 进行主成分分析。

## YAML 配置

```yaml
- gmx_PCA:
    group: C-alpha
    gmx_parm:
      tu: ns
```

## 参数

| 参数 | 必需 | 说明 |
|------|------|------|
| `group` | 是 | 分析组名，对于蛋白质一般选 `C-alpha` |
| `mkdir` | 否 | 输出目录名 |
| `gmx_parm` | 否 | 传递给 `gmx covar` 和 `gmx anaeig` 的共有参数（如 `-b`, `-e`） |

## 输出

- 前 3 主成分两两散点图（pc12.png, pc13.png, pc23.png）
- 前 10 主成分占比图
- 所有主成分占比图
- 主成分数据文件（pc12.xvg, pc13.xvg, pc23.xvg）- 可直接用于 `gmx_FEL`
- 余弦含量（检查采样收敛性）
- PC 极值投影 PDB 文件（pc1_proj.pdb 等）- 可用 PyMOL 可视化

## 余弦含量检查

前几个主成分余弦含量接近 1 时，说明该 PC 可能对应随机扩散，模拟未收敛。参考：Berk Hess. Phys. Rev. E 65, 031910 (2002).

## 官方文档

https://duivyprocedures-docs.readthedocs.io/en/latest/gmx_PCA.html
