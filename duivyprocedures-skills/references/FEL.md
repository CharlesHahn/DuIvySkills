# FEL

纯 Python 模块，从三列数据绘制自由能景观。

## YAML 配置

```yaml
- FEL:
    inputfile: ../PCA/pc12.xvg
    temperature: 300
    ngrid: 32
    find_minimum: true
    minimum_num: 5
```

## 参数

| 参数 | 必需 | 说明 |
|------|------|------|
| `inputfile` | 是 | 输入文件，三列数据（时间, 数据1, 数据2） |
| `temperature` | 是 | 体系温度（K），用于能量转换 |
| `ngrid` | 否 | 网格数，默认 32 |
| `find_minimum` | 否 | 是否寻找局部最小值 |
| `minimum_num` | 否 | 最小值数量 |
| `mkdir` | 否 | 输出目录名 |

## 输出

- Gibbs 自由能图（.xpm + .png）
- 最小值对应帧的 PDB 文件
- 最小值标记散点图

## 官方文档

https://duivyprocedures-docs.readthedocs.io/en/latest/FEL.html