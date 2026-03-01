# gmx_DCCM

依赖 GROMACS 计算动态互相关矩阵。

## YAML 配置

```yaml
- gmx_DCCM:
    group: C-alpha
```

## 参数

| 参数 | 必需 | 说明 |
|------|------|------|
| `group` | 是 | 计算协方差矩阵的原子组 |
| `mkdir` | 否 | 输出目录名 |
| `gmx_parm` | 否 | 传递给 `gmx covar` 的参数 |

## 输出

- 协方差矩阵（.xpm）
- DCCM 矩阵（.xpm）
- DCCM 可视化图（.png）

## 示例

```yaml
- gmx_DCCM:
    group: C-alpha
    gmx_parm:
      tu: ns
```

## 官方文档

https://duivyprocedures-docs.readthedocs.io/en/latest/gmx_DCCM.html