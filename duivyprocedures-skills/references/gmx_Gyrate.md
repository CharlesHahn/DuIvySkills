# gmx_Gyrate

依赖 GROMACS 计算回转半径。

## YAML 配置

```yaml
- gmx_Gyrate:
    calc_group: Protein
```

## 参数

| 参数 | 必需 | 说明 |
|------|------|------|
| `calc_group` | 是 | 计算组名 |
| `mkdir` | 否 | 输出目录名 |
| `gmx_parm` | 否 | 传递给 `gmx gyrate` 的参数 |

## 输出

- 回转半径时间序列（.xvg + .png）

## 示例

```yaml
- gmx_Gyrate:
    calc_group: Protein
    gmx_parm:
      tu: ns
```

## 官方文档

https://duivyprocedures-docs.readthedocs.io/en/latest/gmx_Gyrate.html