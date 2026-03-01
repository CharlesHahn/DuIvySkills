# PiStacking

Pi-Pi 堆叠相互作用分析。

## YAML 配置

```yaml
- PiStacking:
    distance_max_cutoff: 0.55
    distance_min_cutoff: 0.05
    ring_center_offset: 0.20
    angle4T_stacking: [60, 90]
    angle4P_stacking: [0, 30]
    byIndex: no
    group1: protein
    group2: resname *ZIN
    only_aromatic_rings: yes
    calc_lifetime: no
```

## 参数

| 参数 | 必需 | 说明 |
|------|------|------|
| `distance_max_cutoff` | 否 | 最大距离（nm） |
| `distance_min_cutoff` | 否 | 最小距离（nm） |
| `ring_center_offset` | 否 | 环中心偏移阈值（nm） |
| `angle4T_stacking` | 否 | T-stacking 角度范围 |
| `angle4P_stacking` | 否 | P-stacking 角度范围 |
| `byIndex` | 否 | 是否通过索引定义环 |
| `group1` / `group2` | 否 | 搜索组（`byIndex: no`），MDAnalysis 语法 |
| `Pi_rings_Index` | 否 | 环索引列表（`byIndex: yes`） |
| `only_aromatic_rings` | 否 | 仅考虑芳香环 |
| `planarity_cutoff` | 否 | 平面度阈值（度） |
| `calc_lifetime` | 否 | 计算生命周期 |

## 输出

- 环结构 PDB 文件
- 距离/角度/偏移时间序列
- 占有率图
- CSV 文件：占有率、P-stacking/T-stacking 占有率

## 官方文档

https://duivyprocedures-docs.readthedocs.io/en/latest/PiStacking.html