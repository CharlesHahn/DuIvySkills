# gmx_SASA

依赖 GROMACS 计算溶剂可及表面积。

## YAML 配置

```yaml
- gmx_SASA:
    calc_group: Protein
```

## 参数

| 参数 | 必需 | 说明 |
|------|------|------|
| `calc_group` | 是 | 计算组名 |
| `mkdir` | 否 | 输出目录名 |
| `gmx_parm` | 否 | 传递给 `gmx sasa` 的参数 |

## 输出

- SASA 时间序列（.xvg + .png）

## 示例

```yaml
- gmx_SASA:
    calc_group: Protein
    gmx_parm:
      tu: ns
```

## 官方文档

https://duivyprocedures-docs.readthedocs.io/en/latest/gmx_SASA.html