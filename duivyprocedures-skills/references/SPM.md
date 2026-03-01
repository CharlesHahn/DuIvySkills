# SPM

空间映射分析。

## YAML 配置

```yaml
- SPM:
    atom_selection: protein
    bin_size: 0.1
```

## 参数

| 参数 | 必需 | 说明 |
|------|------|------|
| `atom_selection` | 是 | 原子选择，MDAnalysis 语法 |
| `bin_size` | 否 | 格点大小（nm） |

## 输出

- 空间映射图
- 密度数据

## 官方文档

https://duivyprocedures-docs.readthedocs.io/en/latest/SPM.html