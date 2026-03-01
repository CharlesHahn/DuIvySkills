# DensityMap

三维密度映射，使用格点技术计算空间分布。

## YAML 配置

```yaml
- DensityMap:
    byType: Mass
    byIndex: no
    groups: [resname *ZIN, protein]
    grid_bin: 1.2
    doSplit_axis: "Y"
    doSplit_saveXPM: no
```

## 参数

| 参数 | 必需 | 说明 |
|------|------|------|
| `byType` | 否 | 计算类型：`Mass`、`Number`、`Charge` |
| `byIndex` | 否 | `no` 用 MDAnalysis 语法，`yes` 用 GROMACS 索引组 |
| `groups` | 是 | 原子组列表 |
| `grid_bin` | 否 | 格点大小（Å），默认 1.2 |
| `doSplit_axis` | 否 | 切片轴向：`X`、`Y`、`Z`、`XY` 等 |
| `doSplit_saveXPM` | 否 | 是否保存切片 XPM |

## 输出

- 各轴向上平均密度分布
- 各平面上密度热图
- 组间密度差值图
- 切片图（如设置）

## 官方文档

https://duivyprocedures-docs.readthedocs.io/en/latest/DensityMap.html