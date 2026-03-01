# gmx_Density

依赖 GROMACS 计算沿轴向的密度分布。

## YAML 配置

```yaml
- gmx_Density:
    calc_group: Protein
    gmx_parm:
      d: Z
```

## 参数

| 参数 | 必需 | 说明 |
|------|------|------|
| `calc_group` | 是 | 计算组名 |
| `mkdir` | 否 | 输出目录名 |
| `gmx_parm` | 否 | 传递给 `gmx density` 的参数 |

## gmx_parm 常用参数

| 参数 | 说明 |
|------|------|
| `d` | 轴向方向：`X`、`Y`、`Z`（默认 Z） |

## 输出

- 密度分布图（.png + .xvg）

## 官方文档

https://duivyprocedures-docs.readthedocs.io/en/latest/gmx_Density.html