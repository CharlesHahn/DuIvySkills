# Hydrophobic Contact

疏水接触分析。

## YAML 配置

```yaml
- Hydrophobic Contact:
    distance_cutoff: 0.55
    group1: protein
    group2: resname *ZIN
    calc_lifetime: no
```

## 参数

| 参数 | 必需 | 说明 |
|------|------|------|
| `distance_cutoff` | 否 | 距离阈值（nm） |
| `group1` | 是 | 第一组，MDAnalysis 语法 |
| `group2` | 是 | 第二组，MDAnalysis 语法 |
| `calc_lifetime` | 否 | 计算生命周期 |

## 输出

- 接触时间序列
- 占有率图
- CSV 文件

## 官方文档

https://duivyprocedures-docs.readthedocs.io/en/latest/Hydrophobic%20Contact.html