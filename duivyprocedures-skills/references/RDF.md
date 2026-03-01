# RDF

径向分布函数计算。

## YAML 配置

```yaml
- RDF:
    center_group: protein
    calc_group: resname *ZIN
    range: 4
```

## 参数

| 参数 | 必需 | 说明 |
|------|------|------|
| `center_group` | 是 | RDF 中心组，MDAnalysis 语法 |
| `calc_group` | 是 | 计算组，MDAnalysis 语法 |
| `range` | 否 | 计算半径范围（nm），默认 4 |
| `mkdir` | 否 | 输出目录名 |

## 输出

- RDF 数据文件（.xvg）
- RDF 可视化图（.png）

## 官方文档

https://duivyprocedures-docs.readthedocs.io/en/latest/RDF.html