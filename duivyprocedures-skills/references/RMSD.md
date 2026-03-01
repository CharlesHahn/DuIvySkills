# RMSD

纯 Python 模块，计算均方根偏差。使用 MDAnalysis 原子选择语法。

## YAML 配置

```yaml
- RMSD:
    fit_group: backbone
    calc_group: resname *ZIN
    mass_weight: yes
```

## 参数

| 参数 | 必需 | 说明 |
|------|------|------|
| `fit_group` | 是 | 对齐组，MDAnalysis 语法 |
| `calc_group` | 是 | 计算 RMSD 的组，MDAnalysis 语法 |
| `mass_weight` | 否 | 是否质量加权，默认 `no` |
| `mkdir` | 否 | 输出目录名 |
| `frame_start` | 否 | 起始帧 |
| `frame_end` | 否 | 结束帧 |
| `frame_step` | 否 | 帧步长 |

## 输出

- RMSD 数据文件（.xvg）
- RMSD 可视化图（.png）

## 原子选择语法

参考：https://userguide.mdanalysis.org/2.7.0/selections.html

常用示例：
- `protein` - 蛋白质
- `protein and name CA` - α-碳
- `resname *ZIN` - 名称匹配 ZIN 的残基
- `backbone` - 主链

## 官方文档

https://duivyprocedures-docs.readthedocs.io/en/latest/RMSD.html