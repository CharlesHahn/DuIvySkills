# gmx_RMSD

依赖 GROMACS 计算均方根偏差。

## YAML 配置

```yaml
- gmx_RMSD:
    calc_group: Protein
    fit_group: Backbone
    rmsd_matrix: no
```

## 参数

| 参数 | 必需 | 说明 |
|------|------|------|
| `calc_group` | 是 | 计算 RMSD 的组 |
| `fit_group` | 是 | 对齐组（RMSD 计算前先对齐） |
| `rmsd_matrix` | 否 | 是否输出帧间 RMSD 矩阵，`yes`/`no`，默认 `no`。计算矩阵较耗时 |
| `mkdir` | 否 | 输出目录名 |
| `gmx_parm` | 否 | 传递给 `gmx rms` 的参数 |

## gmx_parm 常用参数

```yaml
gmx_parm:
  tu: ns      # 时间单位
  dt: 0.1     # 帧间隔
  b: 1000     # 起始时间
  e: 50000    # 结束时间
```

## 输出

- RMSD 时间序列（.xvg + .png）
- RMSD 矩阵热图（如 `rmsd_matrix: yes`）

## 示例

```yaml
- gmx_RMSD:
    mkdir: RMSD_matrix
    calc_group: Protein
    fit_group: Backbone
    rmsd_matrix: yes
    gmx_parm:
      tu: ns
      dt: 0.1
```

## 官方文档

https://duivyprocedures-docs.readthedocs.io/en/latest/gmx_RMSD.html