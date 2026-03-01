# gmx_FEL

依赖 GROMACS 的 `gmx sham` 绘制自由能景观。

## YAML 配置

```yaml
- gmx_FEL:
    inputfile: ../gmx_PCA/pc12.xvg
    find_minimum: true
    minimum_num: 2
    gmx_parm:
      tsham: 310
```

## 参数

| 参数 | 必需 | 说明 |
|------|------|------|
| `inputfile` | 是 | 输入文件路径，三列数据（时间, 数据1, 数据2） |
| `find_minimum` | 否 | 是否寻找局部最小值，默认 `false` |
| `minimum_num` | 否 | 寻找的最小值数量 |
| `mkdir` | 否 | 输出目录名 |
| `gmx_parm` | 否 | 传递给 `gmx sham` 的参数 |

## gmx_parm 常用参数

| 参数 | 说明 |
|------|------|
| `tsham` | 体系温度（K） |

## 输出

- Gibbs 自由能图（gibbs.xpm + .png）
- 概率分布图（prob.xpm）
- 焓/熵图
- 最小值对应帧的 PDB 文件（如 `find_minimum: true`）

## 数据来源

通常使用 PCA 或 dPCA 的输出：
```yaml
inputfile: ../gmx_PCA/pc12.xvg
inputfile: ../gmx_dPCA/dpc12.xvg
```

## 官方文档

https://duivyprocedures-docs.readthedocs.io/en/latest/gmx_FEL.html