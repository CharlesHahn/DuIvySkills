# DCCM

纯 Python 模块，计算动态互相关矩阵。使用 MDAnalysis 原子选择语法。

## YAML 配置

```yaml
- DCCM:
    atom_selection: protein and name CA
    byType: atom
    save_xpm: yes
```

## 参数

| 参数 | 必需 | 说明 |
|------|------|------|
| `atom_selection` | 是 | 原子选择，MDAnalysis 语法 |
| `byType` | 否 | 计算方式：`atom`（默认）、`res_com`、`res_cog`、`res_coc` |
| `save_xpm` | 否 | 是否保存 XPM 文件，默认 `yes` |
| `mkdir` | 否 | 输出目录名 |
| `frame_start` | 否 | 起始帧 |
| `frame_end` | 否 | 结束帧 |
| `frame_step` | 否 | 帧步长 |

## byType 说明

| 值 | 说明 |
|------|------|
| `atom` | 计算选中原子间的 DCCM |
| `res_com` | 残基质心间的 DCCM |
| `res_cog` | 残基几何中心间的 DCCM |
| `res_coc` | 残基电荷中心间的 DCCM |

使用 `res_*` 时，`atom_selection` 应选择残基的所有原子。

## 输出

- 协方差矩阵（.xpm + .csv）
- DCCM 矩阵（.xpm + .csv）
- 可视化图（.png）

## 官方文档

https://duivyprocedures-docs.readthedocs.io/en/latest/DCCM.html