# UMAP

UMAP 降维分析。

## YAML 配置

```yaml
- UMAP:
    atom_selection: protein and name CA
    byType: atom
    target: coordinates
    n_neighbors: [50, 100, 200]
    min_dist: [1]
```

## 参数

| 参数 | 必需 | 说明 |
|------|------|------|
| `atom_selection` | 是 | 原子选择，MDAnalysis 语法 |
| `byType` | 否 | 计算方式：`atom`、`res_com`、`res_cog`、`res_coc` |
| `target` | 否 | 目标：`coordinates` 或 `dihedrals` |
| `n_neighbors` | 否 | 近邻数量列表，DIP 会遍历所有组合 |
| `min_dist` | 否 | 控制点堆积紧密程度的列表 |

## 输出

- 各参数组合的降维结果
- 散点图

## 官方文档

https://duivyprocedures-docs.readthedocs.io/en/latest/UMAP.html