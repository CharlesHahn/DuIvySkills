# gmx_RMSF

依赖 GROMACS 计算均方根涨落。

## YAML 配置

```yaml
- gmx_RMSF:
    group: C-alpha
```

## 参数

| 参数 | 必需 | 说明 |
|------|------|------|
| `group` | 是 | 计算组名 |
| `mkdir` | 否 | 输出目录名 |
| `gmx_parm` | 否 | 传递给 `gmx rmsf` 的参数 |

## 输出

- RMSF 数据文件（.xvg）
- RMSF 可视化图（.png）

## 示例

```yaml
- gmx_RMSF:
    group: C-alpha
    gmx_parm:
      tu: ns
```

## 官方文档

https://duivyprocedures-docs.readthedocs.io/en/latest/gmx_RMSF.html