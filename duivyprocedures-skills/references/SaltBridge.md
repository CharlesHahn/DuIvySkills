# SaltBridge

盐桥分析。

## YAML 配置

```yaml
- SaltBridge:
    dist_cutoff: 0.32
    byIndex: no
    group: protein
    ignore_chain_end: no
    calc_lifetime: no
```

## 参数

| 参数 | 必需 | 说明 |
|------|------|------|
| `dist_cutoff` | 否 | 距离阈值（nm），默认 0.32 |
| `byIndex` | 否 | 是否通过索引定义原子组 |
| `group` | 否 | 搜索组（`byIndex: no` 时使用），MDAnalysis 语法 |
| `positive_Index` | 否 | 正电基团索引列表（`byIndex: yes`） |
| `negative_Index` | 否 | 负电基团索引列表（`byIndex: yes`） |
| `NH3_atomnames` | 否 | NH3+ 原子名列表（不同力场需修改） |
| `COO_atomnames` | 否 | COO- 原子名列表（不同力场需修改） |
| `ignore_chain_end` | 否 | 忽略链端残基 |
| `calc_lifetime` | 否 | 计算生命周期 |

## 输出

- 带电基团 PDB 文件（确认选择）
- 盐桥距离时间序列
- 占有率矩阵图
- CSV 文件：占有率、平均距离

## 注意

不同力场的原子命名不同，需要根据体系修改 `NH3_atomnames` 和 `COO_atomnames`。

## 官方文档

https://duivyprocedures-docs.readthedocs.io/en/latest/SaltBridge.html