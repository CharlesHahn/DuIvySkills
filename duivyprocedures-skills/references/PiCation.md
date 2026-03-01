# PiCation

Pi-阳离子相互作用分析。

## YAML 配置

```yaml
- PiCation:
    distance_cutoff: 0.60
    byIndex: no
    group1: protein
    group2: resname *ZIN
    only_aromatic_rings: yes
    calc_lifetime: no
```

## 参数

| 参数 | 必需 | 说明 |
|------|------|------|
| `distance_cutoff` | 否 | 距离阈值（nm） |
| `byIndex` | 否 | 是否通过索引定义 |
| `group1` / `group2` | 否 | 搜索组，MDAnalysis 语法 |
| `Pi_rings_Index` | 否 | 环索引列表 |
| `Cation_Index` | 否 | 阳离子索引列表 |
| `only_aromatic_rings` | 否 | 仅考虑芳香环 |
| `calc_lifetime` | 否 | 计算生命周期 |

## 输出

- 距离时间序列
- 占有率图
- CSV 文件：占有率、平均距离

## 官方文档

https://duivyprocedures-docs.readthedocs.io/en/latest/PiCation.html