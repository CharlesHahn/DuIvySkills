# tSNE

t-SNE 降维分析。

## YAML 配置

```yaml
- tSNE:
    atom_selection: protein and name CA
    byType: atom
    target: coordinates
```

## 参数

| 参数 | 必需 | 说明 |
|------|------|------|
| `atom_selection` | 是 | 原子选择，MDAnalysis 语法 |
| `byType` | 否 | 计算方式：`atom`、`res_com`、`res_cog`、`res_coc` |
| `target` | 否 | 目标：`coordinates` 或 `dihedrals` |

## 输出

- 降维结果数据文件
- 散点图

## 官方文档

https://duivyprocedures-docs.readthedocs.io/en/latest/tSNE.html