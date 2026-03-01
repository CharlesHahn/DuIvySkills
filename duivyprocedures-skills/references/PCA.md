# PCA

纯 Python 模块，主成分分析。使用 MDAnalysis 原子选择语法。

## YAML 配置

```yaml
- PCA:
    atom_selection: protein and name CA
    byType: atom
```

## 参数

| 参数 | 必需 | 说明 |
|------|------|------|
| `atom_selection` | 是 | 原子选择，MDAnalysis 语法 |
| `byType` | 否 | 计算方式：`atom`、`res_com`、`res_cog`、`res_coc` |
| `mkdir` | 否 | 输出目录名 |
| `frame_start` | 否 | 起始帧 |
| `frame_end` | 否 | 结束帧 |
| `frame_step` | 否 | 帧步长 |

## 输出

- 主成分散点图
- 主成分占比图
- 主成分数据文件（.xvg）

## 官方文档

https://duivyprocedures-docs.readthedocs.io/en/latest/PCA.html