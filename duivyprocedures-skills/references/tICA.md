# tICA

时间滞后独立成分分析（Time-lagged Independent Component Analysis）。

## YAML 配置

```yaml
# 基于坐标
- tICA:
    atom_selection: protein and name CA
    byType: atom
    target: coordinates
    lag: 10

# 基于二面角
- tICA:
    atom_selection: protein
    byType: atom
    target: dihedrals
    lag: 10
```

## 参数

| 参数 | 必需 | 说明 |
|------|------|------|
| `atom_selection` | 是 | 原子选择，MDAnalysis 语法。二面角分析需包含骨架原子 |
| `byType` | 否 | 计算方式（仅 `coordinates` 有效）：`atom`、`res_com`、`res_cog`、`res_coc` |
| `target` | 否 | 目标：`coordinates` 或 `dihedrals` |
| `lag` | 是 | lag time 参数 |
| `mkdir` | 否 | 输出目录名 |

## byType 说明

| 值 | 说明 |
|------|------|
| `atom` | 选中原子坐标的 tICA |
| `res_com` | 残基质心的 tICA |
| `res_cog` | 残基几何中心的 tICA |
| `res_coc` | 残基电荷中心的 tICA |

使用 `res_*` 时，`atom_selection` 应包含残基所有原子。

## 二面角分析注意

角度会转换为 sin/cos 值再分析。请对照文献确认计算过程是否合适。

## 输出

- 前 3 成分数据文件（.xvg）
- 两两散点图

## 官方文档

https://duivyprocedures-docs.readthedocs.io/en/latest/tICA.html
