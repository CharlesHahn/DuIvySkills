# gmx_cluster

依赖 GROMACS 进行聚类分析。

## YAML 配置

```yaml
- gmx_cluster:
    fit_group: Backbone
    calc_group: Protein
    gmx_parm:
      method: gromos
      cutoff: 0.1
      tu: ns
      dt: 0.1
```

## 参数

| 参数 | 必需 | 说明 |
|------|------|------|
| `fit_group` | 是 | 对齐组 |
| `calc_group` | 是 | 聚类计算组 |
| `mkdir` | 否 | 输出目录名 |
| `gmx_parm` | 否 | 传递给 `gmx cluster` 的参数 |

## gmx_parm 常用参数

| 参数 | 说明 |
|------|------|
| `method` | 聚类方法：`gromos`（默认）、`linkage`、`jarvis-patrick` 等 |
| `cutoff` | 聚类截断值（nm） |

## 输出

- RMSD 矩阵热图
- 聚类大小分布图
- 各聚类占比数据

## 注意

`dt` 设置过小会导致帧数多、计算量大。

## 官方文档

https://duivyprocedures-docs.readthedocs.io/en/latest/gmx_cluster.html