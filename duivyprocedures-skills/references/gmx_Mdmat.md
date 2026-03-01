# gmx_Mdmat

依赖 GROMACS 计算残基间最短距离矩阵（RDCM）。

## YAML 配置

```yaml
- gmx_Mdmat:
    group: Protein
    gmx_parm:
      t: 1.5
```

## 参数

| 参数 | 必需 | 说明 |
|------|------|------|
| `group` | 是 | 计算组名 |
| `mkdir` | 否 | 输出目录名 |
| `gmx_parm` | 否 | 传递给 `gmx mdmat` 的参数 |

## gmx_parm 常用参数

| 参数 | 说明 |
|------|------|
| `t` | 距离截断值（nm） |

## 输出

- 平均残基接触矩阵图
- 接触数量时间变化图

## 官方文档

https://duivyprocedures-docs.readthedocs.io/en/latest/gmx_Mdmat.html